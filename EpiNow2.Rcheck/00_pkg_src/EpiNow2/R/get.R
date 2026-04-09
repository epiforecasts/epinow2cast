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
  raw_samples <- extract_samples(object$fit)
  # extract_delays returns data.table with variable column

  samples <- extract_delays(raw_samples, args = object$args)
  samples[]
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

  # Process input observations to get dates
  dirty_obs <- purrr::map(object$observations, data.table::as.data.table)
  earliest_date <- max(
    as.Date(
      purrr::map_chr(dirty_obs, function(x) x[, as.character(min(date))])
    )
  )
  dirty_obs <- purrr::map(dirty_obs, function(x) x[date >= earliest_date])
  nrow_obs <- order(purrr::map_dbl(dirty_obs, nrow))
  dirty_obs <- dirty_obs[nrow_obs]

  obs_sets <- object$args$obs_sets
  trunc_max <- object$args$delay_max[1]

  if (format == "summary") {
    # Extract reconstructed observations summary statistics
    recon_obs <- extract_stan_param(object$fit, "recon_obs",
      CrIs = CrIs,
      var_names = TRUE
    )
    recon_obs <- recon_obs[, id := variable][, variable := NULL]

    # Assign dataset index using modulo
    recon_obs <- recon_obs[, dataset := seq_len(.N)][
      ,
      dataset := dataset %% obs_sets
    ][
      dataset == 0, dataset := obs_sets
    ]

    # Link predictions to dates
    link_preds <- function(index) {
      target_obs <- dirty_obs[[index]][, idx := .N - 0:(.N - 1)]
      target_obs <- target_obs[idx < trunc_max]
      estimates <- recon_obs[dataset == index][, c("id", "dataset") := NULL]
      estimates <- estimates[, lapply(.SD, as.integer)]
      estimates <- estimates[, idx := .N - 0:(.N - 1)]
      if (!is.null(estimates$n_eff)) estimates[, "n_eff" := NULL]
      if (!is.null(estimates$Rhat)) estimates[, "Rhat" := NULL]

      result <- data.table::merge.data.table(
        target_obs[, .(date, idx)],
        estimates,
        by = "idx", all.x = TRUE
      )
      result[, report_date := max(target_obs$date)]
      result[order(date)][, idx := NULL]
    }

    predictions <- purrr::map(seq_len(obs_sets), link_preds)
    data.table::rbindlist(predictions)
  } else {
    # Both "sample" and "quantile" need raw samples first
    raw_samples <- extract_samples(object$fit, pars = "recon_obs")
    recon_samples <- data.table::as.data.table(raw_samples$recon_obs)
    recon_samples <- data.table::melt(recon_samples,
      measure.vars = seq_len(ncol(recon_samples)),
      variable.name = "obs_idx",
      value.name = "predicted"
    )
    recon_samples[, obs_idx := as.integer(obs_idx)]
    recon_samples[, sample := seq_len(.N), by = obs_idx]
    recon_samples[, dataset := ((obs_idx - 1) %% obs_sets) + 1]

    # Link samples to dates
    link_samples <- function(index) {
      target_obs <- dirty_obs[[index]][, idx := .N - 0:(.N - 1)]
      target_obs <- target_obs[idx < trunc_max]
      target_obs[, obs_idx := seq_len(.N)]

      samples_subset <- recon_samples[dataset == index]
      result <- data.table::merge.data.table(
        target_obs[, .(date, obs_idx)],
        samples_subset[, .(obs_idx, sample, predicted)],
        by = "obs_idx"
      )[, obs_idx := NULL]

      # Add forecast metadata
      forecast_date <- max(target_obs$date, na.rm = TRUE)
      result[, forecast_date := forecast_date]
      result[, horizon := as.numeric(date - forecast_date)]
      result[, dataset := index]

      result
    }

    predictions <- purrr::map(seq_len(obs_sets), link_samples)
    predictions <- data.table::rbindlist(predictions)

    if (format == "sample") {
      # Reorder columns for sample format
      data.table::setcolorder(
        predictions,
        c("dataset", "forecast_date", "date", "horizon", "sample", "predicted")
      )
    } else {
      # format == "quantile": aggregate to quantiles
      predictions <- predictions[
        ,
        .(predicted = quantile(predicted, probs = quantiles)),
        by = .(dataset, forecast_date, date, horizon)
      ]
      predictions[
        , quantile_level := rep(quantiles, .N / length(quantiles))
      ]
      data.table::setcolorder(
        predictions,
        c("dataset", "forecast_date", "date", "horizon",
          "quantile_level", "predicted")
      )
    }

    predictions[]
  }
}


#' @rdname get_parameters
#' @export
get_parameters.epinowfit <- function(x, ...) {
  result <- list()
  if (!is.null(x$args$generation_time)) {
    result$generation_time <- x$args$generation_time
  }
  result
}
