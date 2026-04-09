#' Get Folders with Results
#'
#' @description `r lifecycle::badge("stable")`
#'
#' @param results_dir A character string giving the directory in which results
#'  are stored (as produced by [regional_epinow()]).
#'
#' @return A named character vector containing the results to plot.
#' @keywords internal
get_regions <- function(results_dir) {
  # regions to include - based on folder names
  regions <- list.dirs(results_dir,
    recursive = FALSE,
    full.names = FALSE
  )

  # put into alphabetical order
  regions <- regions[regions != "runtimes.csv"]
  regions <- sort(regions)
  names(regions) <- regions
  regions
}

#' Get a Single Raw Result
#'
#' @description `r lifecycle::badge("stable")`
#'
#' @param file Character string giving the result files name.
#'
#' @param region Character string giving the region of interest.
#'
#' @param date Target date (in the format `"yyyy-mm-dd`).
#'
#' @param result_dir Character string giving the location of the target
#' directory.
#'
#' @return An R object read in from the targeted `.rds` file
#' @keywords internal
get_raw_result <- function(file, region, date,
                           result_dir) {
  file_path <- file.path(result_dir, region, date, file)
  readRDS(file_path)
}
#' Get Combined Regional Results
#'
#' @description `r lifecycle::badge("stable")`
#' Summarises results across regions either from input or from disk. See the
#' examples for details.
#'
#' @param regional_output A list of output as produced by [regional_epinow()]
#' and stored in the `regional` list.
#'
#' @param results_dir A character string indicating the folder containing the
#' `{EpiNow2}` results to extract.
#'
#' @param date A Character string (in the format "yyyy-mm-dd") indicating the
#' date to extract data for. Defaults to "latest" which finds the latest
#' results available.
#'
#' @param samples Logical, defaults to `TRUE`. Should samples be returned.
#'
#' @param forecast Logical, defaults to `FALSE`. Should forecast results be
#' returned.
#'
#' @return A list of estimates, forecasts and estimated cases by date of report.
#' @export
#' @importFrom purrr map safely
#' @importFrom data.table rbindlist
#' @examples
#' # get example multiregion estimates
#' regional_out <- readRDS(system.file(
#'   package = "EpiNow2", "extdata", "example_regional_epinow.rds"
#' ))
#'
#' # from output
#' results <- get_regional_results(regional_out$regional, samples = FALSE)
get_regional_results <- function(regional_output,
                                 results_dir, date,
                                 samples = TRUE,
                                 forecast = FALSE) {
  if (missing(regional_output)) {
    regional_output <- NULL
  }

  if (is.null(regional_output)) {
    # assign to latest likely date if not given
    if (missing(date)) {
      date <- "latest"
    }
    # find all regions
    regions <- get_regions(results_dir)

    load_data <- purrr::safely(get_raw_result) # nolint

    # get estimates
    get_estimates_file <- function(samples_path, summarised_path) {
      out <- list()

      if (samples) {
        samples <- purrr::map(
          regions, ~ load_data(
            samples_path, .,
            result_dir = results_dir,
            date = date
          )[[1]]
        )
        samples <- data.table::rbindlist(samples, idcol = "region", fill = TRUE)
        out$samples <- samples
      }
      # get incidence values and combine
      summarised <- purrr::map(
        regions, ~ load_data(
          summarised_path,
          .,
          result_dir = results_dir,
          date = date
        )[[1]]
      )
      summarised <- data.table::rbindlist(
        summarised,
        idcol = "region", fill = TRUE
      )
      out$summarised <- summarised
      out
    }
    out <- list()
    out$estimates <- get_estimates_file(
      samples_path = "estimate_samples.rds",
      summarised_path = "summarised_estimates.rds"
    )

    if (forecast) {
      out$estimated_reported_cases <- get_estimates_file(
        samples_path = "estimated_reported_cases_samples.rds",
        summarised_path = "summarised_estimated_reported_cases.rds"
      )
    }
  } else {
    out <- list()
    estimates_out <- list()

    if (samples) {
      samp <- purrr::map(regional_output, get_samples)
      samp <- data.table::rbindlist(samp, idcol = "region", fill = TRUE)
      estimates_out$samples <- samp
    }
    summarised <- purrr::map(
      regional_output, summary, type = "parameters"
    )
    summarised <- data.table::rbindlist(
      summarised,
      idcol = "region", fill = TRUE
    )
    estimates_out$summarised <- summarised
    out$estimates <- estimates_out

    if (forecast) {
      erc_out <- list()
      erc_data <- purrr::map(regional_output, estimates_by_report_date)
      if (samples) {
        samp <- purrr::map(erc_data, ~ .$samples)
        samp <- data.table::rbindlist(samp, idcol = "region", fill = TRUE)
        erc_out$samples <- samp
      }
      summarised <- purrr::map(erc_data, ~ .$summarised)
      summarised <- data.table::rbindlist(
        summarised,
        idcol = "region", fill = TRUE
      )
      erc_out$summarised <- summarised
      out$estimated_reported_cases <- erc_out
    }
  }
  out
}

#' Get Regions with Most Reported Cases
#'
#' @description `r lifecycle::badge("stable")`
#' Extract a vector of regions with the most reported cases in a set time
#' window.
#'
#' @param time_window Numeric, number of days to include from latest date in
#' data. Defaults to 7 days.
#'
#' @param no_regions Numeric, number of regions to return. Defaults to 6.
#'
#' @inheritParams regional_epinow
#'
#' @return A character vector of regions with the highest reported cases
#'
#' @importFrom data.table copy setorderv
#' @importFrom lubridate days
#' @keywords internal
get_regions_with_most_reports <- function(data,
                                          time_window = 7,
                                          no_regions = 6) {
  most_reports <- data.table::copy(data)
  most_reports <-
    most_reports[,
      .SD[date >= (max(date, na.rm = TRUE) - lubridate::days(time_window))],
      by = "region"
    ]
  most_reports <- most_reports[,
    .(confirm = sum(confirm, na.rm = TRUE)),
    by = "region"
  ]
  most_reports <- data.table::setorderv(
    most_reports,
    cols = "confirm", order = -1
  )
  most_reports[1:no_regions][!is.na(region)]$region
}

##' Estimate seeding time from delays and generation time
##'
##' The seeding time is set to the mean of the specified delays, constrained
##' to be at least the maximum generation time
##' @inheritParams estimate_infections
##' @return An integer seeding time
##' @keywords internal
get_seeding_time <- function(delays, generation_time, rt = rt_opts()) {
  # Estimate the mean delay -----------------------------------------------
  seeding_time <- sum(mean(delays, ignore_uncertainty = TRUE))
  if (!is.null(rt)) {
    ## make sure we have at least max(generation_time) seeding time
    seeding_time <- max(seeding_time, sum(max(generation_time)))
  }
  max(round(seeding_time), 1)
}

#' Get posterior samples from a fitted model
#'
#' @description `r lifecycle::badge("stable")`
#' Extracts posterior samples from a fitted model, combining all parameters
#' into a single data.table with dates and metadata.
#'
#' @param object A fitted model object (e.g., from `estimate_infections()`)
#' @param ... Additional arguments (currently unused)
#'
#' @return A `data.table` with columns: date, variable, strat, sample, time,
#'   value, type. Contains all posterior samples for all parameters.
#'
#' @export
#' @examples
#' \dontrun{
#' # After fitting a model
#' samples <- get_samples(fit)
#' # Filter to specific parameters
#' R_samples <- samples[variable == "R"]
#' }
get_samples <- function(object, ...) {
  UseMethod("get_samples")
}

#' @rdname get_samples
#' @export
get_samples.estimate_infections <- function(object, ...) {
  extract_enw_samples(object)
}

#' @rdname get_samples
#' @export
get_samples.epinow <- function(object, ...) {
  # If the epinow run failed (e.g., timeout), throw an informative error
  if (!is.null(object$error)) {
    cli_abort(c(
      "Cannot extract samples from a failed epinow run.",
      "i" = "The run failed with error: {object$error}"
    ))
  }
  # Otherwise delegate to the underlying estimate_infections method
  get_samples.estimate_infections(object, ...)
}

#' @rdname get_samples
#' @export
get_samples.estimate_truncation <- function(object, ...) {
  # Extract delay distribution parameters from epinowcast fit
  fit <- object$fit
  out <- list()

  # Reference delay mean and sd
  tryCatch({
    mean_draws <- fit$draws(
      variables = "refp_mean", format = "draws_matrix"
    )
    out$refp_mean <- enw_draws_to_dt(
      mean_draws, "refp_mean", seq_len(ncol(mean_draws))
    )
  }, error = function(e) NULL)

  tryCatch({
    sd_draws <- fit$draws(
      variables = "refp_sd", format = "draws_matrix"
    )
    out$refp_sd <- enw_draws_to_dt(
      sd_draws, "refp_sd", seq_len(ncol(sd_draws))
    )
  }, error = function(e) NULL)

  combined <- data.table::rbindlist(out, use.names = TRUE)
  if (nrow(combined) > 0) {
    combined[, strat := NA_character_]
    combined[, type := NA_character_]
  }
  combined[]
}

#' Format sample predictions
#'
#' Helper function to format posterior samples into the structure expected by
#' [scoringutils::as_forecast_sample()].
#'
#' @param samples Data.table with date, sample, and value columns
#' @param forecast_date Date when the forecast was made
#' @return Data.table with columns: forecast_date, date, horizon, sample,
#'   predicted
#' @keywords internal
format_sample_predictions <- function(samples, forecast_date) {
  predictions <- samples[, .(date, sample, predicted = value)]
  predictions[, forecast_date := forecast_date]
  predictions[, horizon := as.numeric(date - forecast_date)]
  data.table::setcolorder(
    predictions,
    c("forecast_date", "date", "horizon", "sample", "predicted")
  )
  predictions[]
}

#' Format quantile predictions
#'
#' Helper function to format posterior samples into quantiles in the structure
#' expected by [scoringutils::as_forecast_quantile()].
#'
#' @param samples Data.table with date and value columns
#' @param quantiles Numeric vector of quantile levels
#' @param forecast_date Date when the forecast was made
#' @return Data.table with columns: forecast_date, date, horizon,
#'   quantile_level, predicted
#' @keywords internal
format_quantile_predictions <- function(samples, quantiles, forecast_date) {
  predictions <- samples[
    ,
    .(predicted = quantile(value, probs = quantiles)),
    by = date
  ]
  predictions[, quantile_level := rep(quantiles, .N / length(quantiles))]
  predictions[, forecast_date := forecast_date]
  predictions[, horizon := as.numeric(date - forecast_date)]
  data.table::setcolorder(
    predictions,
    c("forecast_date", "date", "horizon", "quantile_level", "predicted")
  )
  predictions[]
}

#' Get predictions from a fitted model
#'
#' @description `r lifecycle::badge("stable")`
#' Extracts predictions from a fitted model. For `estimate_infections()` returns
#' predicted reported cases, for `estimate_secondary()` returns predicted
#' secondary observations. For `estimate_truncation()` returns reconstructed
#' observations adjusted for truncation.
#'
#' @param object A fitted model object (e.g., from `estimate_infections()`,
#'   `estimate_secondary()`, or `estimate_truncation()`)
#' @param format Character string specifying the output format:
#'   - `"summary"` (default): summary statistics (mean, sd, median, CrIs)
#'   - `"sample"`: raw posterior samples for
#'     [scoringutils::as_forecast_sample()]
#'   - `"quantile"`: quantile predictions for
#'     [scoringutils::as_forecast_quantile()]
#' @param CrIs Numeric vector of credible intervals to return. Defaults to
#'   c(0.2, 0.5, 0.9). Only used when `format = "summary"`.
#' @param quantiles Numeric vector of quantile levels to return. Defaults to
#'   c(0.05, 0.25, 0.5, 0.75, 0.95). Only used when `format = "quantile"`.
#' @param ... Additional arguments (currently unused)
#'
#' @return A `data.table` with columns depending on `format`:
#'   - `format = "summary"`: date, mean, sd, median, and credible intervals
#'   - `format = "sample"`: forecast_date, date, horizon, sample, predicted
#'   - `format = "quantile"`: forecast_date, date, horizon, quantile_level,
#'     predicted
#'
#' @export
#' @examples
#' \dontrun{
#' # After fitting a model
#' # Get summary predictions (default)
#' predictions <- get_predictions(fit)
#'
#' # Get sample-level predictions for scoringutils
#' samples <- get_predictions(fit, format = "sample")
#'
#' # Get quantile predictions for scoringutils
#' quantiles <- get_predictions(fit, format = "quantile")
#' }
get_predictions <- function(object, ...) {
  UseMethod("get_predictions")
}

#' @rdname get_predictions
#' @export
get_predictions.estimate_infections <- function(
    object,
    format = c("summary", "sample", "quantile"),
    CrIs = c(0.2, 0.5, 0.9),
    quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
    ...) {
  format <- rlang::arg_match(format)

  # Get samples for reported cases
  samples <- get_samples(object)
  reported_samples <- samples[variable == "reported_cases"]
  forecast_date <- max(object$observations$date, na.rm = TRUE)

  switch(format,
    summary = calc_summary_measures(
      reported_samples,
      summarise_by = "date",
      order_by = "date",
      CrIs = CrIs
    ),
    sample = format_sample_predictions(reported_samples, forecast_date),
    quantile = format_quantile_predictions(
      reported_samples, quantiles, forecast_date
    )
  )
}

#' @rdname get_predictions
#' @export
get_predictions.estimate_truncation <- function(
    object,
    format = c("summary", "sample", "quantile"),
    CrIs = c(0.2, 0.5, 0.9),
    quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
    ...) {
  format <- rlang::arg_match(format)

  # Extract nowcast posterior predictions
  nowcast <- summary(object$enw_fit, type = "nowcast")

  if (format == "summary") {
    return(nowcast)
  }

  # For sample/quantile, extract raw posterior draws
  fit <- object$fit
  pp_draws <- fit$draws(variables = "pp_inf_obs", format = "draws_matrix")
  n_samples <- nrow(pp_draws)
  n_cols <- ncol(pp_draws)

  # Build sample-level predictions
  latest <- object$enw_fit$latest[[1]]
  dates <- latest$reference_date[seq_len(min(n_cols, nrow(latest)))]
  forecast_date <- max(dates, na.rm = TRUE)

  samples_dt <- enw_draws_to_dt(
    pp_draws[, seq_len(length(dates)), drop = FALSE],
    "predicted", dates
  )
  data.table::setnames(samples_dt, "value", "predicted")
  samples_dt[, variable := NULL]
  samples_dt[, forecast_date := forecast_date]
  samples_dt[, horizon := as.numeric(date - forecast_date)]

  if (format == "sample") {
    data.table::setcolorder(
      samples_dt,
      c("forecast_date", "date", "horizon", "sample", "predicted")
    )
    return(samples_dt[])
  }

  # format == "quantile"
  quant_dt <- samples_dt[
    ,
    .(predicted = stats::quantile(predicted, probs = quantiles)),
    by = .(forecast_date, date, horizon)
  ]
  quant_dt[, quantile_level := rep(quantiles, .N / length(quantiles))]
  data.table::setcolorder(
    quant_dt,
    c("forecast_date", "date", "horizon", "quantile_level", "predicted")
  )
  quant_dt[]
}


#' @rdname get_parameters
#' @export
get_parameters.epinowfit <- function(x, ...) {
  result <- list()
  if (!is.null(x$args$generation_time)) {
    result$generation_time <- x$args$generation_time
  }
  if (!is.null(x$args$delays)) {
    result$reporting <- x$args$delays
  }
  if (!is.null(x$args$truncation)) {
    result$truncation <- x$args$truncation
  }
  result
}
