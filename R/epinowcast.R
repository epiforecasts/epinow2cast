# epinowcast backend
# epinow2cast translation layer
#' Convert EpiNow2 incidence data to epinowcast format
#'
#' @description Converts a simple incidence time series as used by
#'   [estimate_infections()] into the reporting triangle format required by
#'   [epinowcast::enw_preprocess_data()]. All observations are treated as
#'   fully reported on their reference date (delay 0). Actual reporting
#'   delays should be handled via `latent_reporting_delay` in
#'   [epinowcast::enw_expectation()].
#'
#' @param data A `data.frame` with columns `date` and `confirm` (daily
#'   incidence counts), as accepted by [estimate_infections()].
#'
#' @param horizon Integer; number of days to forecast beyond the last
#'   observation. Defaults to 0 (no forecast).
#'
#' @return An `enw_preprocess_data` object suitable for use with
#'   [epinowcast::epinowcast()].
#'
#' @importFrom data.table as.data.table copy
#' @keywords internal
incidence_to_enw <- function(data, horizon = 0L) {
  obs <- data.table::copy(data.table::as.data.table(data))

  # Convert to epinowcast's expected format:
  # reference_date, report_date, confirm (cumulative)
  # All cases reported immediately (delay 0), so max_delay = 1
  data.table::setnames(obs, "date", "reference_date")
  obs[, report_date := reference_date]
  obs[, confirm := cumsum(confirm)]

  # Extend for forecasting if needed
  if (horizon > 0L) {
    last_date <- max(obs$reference_date)
    last_confirm <- obs[reference_date == last_date, confirm]
    future_dates <- seq.Date(last_date + 1L, by = "day", length.out = horizon)
    future_obs <- data.table::data.table(
      reference_date = future_dates,
      report_date = future_dates,
      confirm = last_confirm,
      .observed = FALSE
    )
    obs[, .observed := TRUE]
    obs <- data.table::rbindlist(list(obs, future_obs), fill = TRUE)
  }

  epinowcast::enw_preprocess_data(obs, max_delay = 1)
}

#' Translate EpiNow2 generation time to epinowcast PMF
#'
#' @description Extracts a fixed probability mass function from a
#'   `<generation_time_opts>` object for use as the `generation_time`
#'   argument in [epinowcast::enw_expectation()].
#'
#' @param generation_time A `<generation_time_opts>` object as returned by
#'   [gt_opts()].
#'
#' @return A numeric vector summing to 1, representing the generation time
#'   PMF. The first element corresponds to delay 0.
#'
#' @keywords internal
gt_to_enw <- function(generation_time) {
  gt_fixed <- fix_parameters(generation_time)
  gt_disc <- discretise(gt_fixed)
  gt_collapsed <- collapse(gt_disc)
  pmf <- get_pmf(gt_collapsed)
  pmf / sum(pmf)
}

#' Translate EpiNow2 delays to epinowcast latent reporting delay
#'
#' @description Extracts a fixed PMF from a `<delay_opts>` object for use
#'   as the `latent_reporting_delay` argument in
#'   [epinowcast::enw_expectation()].
#'
#' @param delays A `<delay_opts>` object as returned by [delay_opts()].
#'
#' @return A numeric vector representing the reporting delay PMF.
#'
#' @keywords internal
delays_to_enw <- function(delays) {
  delay_fixed <- fix_parameters(delays)
  delay_disc <- discretise(delay_fixed)
  delay_collapsed <- collapse(delay_disc)
  pmf <- get_pmf(delay_collapsed)
  pmf / sum(pmf)
}

#' Translate EpiNow2 Rt options to epinowcast expectation formula
#'
#' @description Converts [rt_opts()] settings into an appropriate formula
#'   for the `r` argument of [epinowcast::enw_expectation()].
#'
#' @param rt An `<rt_opts>` object as returned by [rt_opts()], or `NULL`.
#'
#' @return A formula suitable for [epinowcast::enw_expectation()].
#'
#' @keywords internal
rt_to_enw_formula <- function(rt) {
  if (is.null(rt)) {
    cli::cli_abort(
      "Back-calculation mode (rt = NULL) is not supported with the
       epinowcast backend."
    )
  }

  if (rt$rw > 0) {
    # Random walk using epinowcast's rw() function
    rw_period <- rt$rw
    if (rw_period == 1) {
      # Daily random walk
      formula <- ~ 0 + rw(day, by = .group)
    } else {
      # Weekly or other period random walk
      formula <- ~ 0 + rw(week, by = .group)
    }
  } else {
    # Default: per-group daily random effects (epinowcast's own default)
    formula <- ~ 0 + (1 | day:.group)
  }
  formula
}

#' Translate EpiNow2 observation options to epinowcast
#'
#' @description Converts [obs_opts()] settings into arguments for
#'   [epinowcast::enw_obs()].
#'
#' @param obs An `<obs_opts>` object as returned by [obs_opts()].
#'
#' @return A list with `family` and `observation_indicator` suitable for
#'   [epinowcast::enw_obs()].
#'
#' @keywords internal
obs_to_enw <- function(obs) {
  # Observation formula for day-of-week effects
  if (isTRUE(obs$week_effect)) {
    observation_formula <- ~ 1 + day_of_week
  } else {
    observation_formula <- ~1
  }
  list(
    family = obs$family,
    observation_formula = observation_formula
  )
}

#' Run estimate_infections via epinowcast backend
#'
#' @description Internal function that translates [estimate_infections()]
#'   arguments into [epinowcast::epinowcast()] calls and converts the
#'   output back to EpiNow2's expected format.
#'
#' @param data A `data.frame` with columns `date` and `confirm`.
#' @param generation_time A `<generation_time_opts>` object.
#' @param delays A `<delay_opts>` object.
#' @param rt An `<rt_opts>` object or `NULL`.
#' @param obs An `<obs_opts>` object.
#' @param forecast A `<forecast_opts>` object.
#' @param stan A `<stan_opts>` object (used for sampling settings).
#' @param verbose Logical.
#'
#' @return An `estimate_infections` object matching EpiNow2's return
#'   structure.
#'
#' @importFrom data.table data.table
#' @keywords internal
run_epinowcast <- function(data, generation_time, delays, rt, obs,
                           forecast, stan, verbose = interactive()) {
  # Convert generation time to PMF

  gt_pmf <- gt_to_enw(generation_time)

  # Convert delays to latent reporting delay PMF
  delay_pmf <- delays_to_enw(delays)

  # Set up forecast horizon
  horizon <- if (!is.null(forecast)) forecast$horizon else 0L

  # Convert data to epinowcast format
  pobs <- incidence_to_enw(data, horizon = horizon)

  # Build epinowcast modules
  obs_args <- obs_to_enw(obs)
  expectation_module <- epinowcast::enw_expectation(
    r = rt_to_enw_formula(rt),
    generation_time = gt_pmf,
    latent_reporting_delay = delay_pmf,
    observation = obs_args$observation_formula,
    data = pobs
  )
  # Tighten the prior on the growth rate random effect SD to reduce
  # the hierarchical funnel
  expectation_module$priors[variable == "expr_beta_sd", sd := 0.2]

  # No reporting delay model — all delay handling is via
  # latent_reporting_delay in the expectation module
  reference_module <- epinowcast::enw_reference(
    parametric = ~0,
    non_parametric = ~0,
    data = pobs
  )

  report_module <- epinowcast::enw_report(~0, data = pobs)

  obs_module_args <- list(
    data = pobs,
    family = obs_args$family
  )
  if (horizon > 0L) {
    obs_module_args$observation_indicator <- ".observed"
  }
  obs_module <- do.call(epinowcast::enw_obs, obs_module_args)

  # Translate stan_opts to enw_fit_opts
  chains <- stan$chains %||% 4L
  samples <- stan$samples %||% 2000L
  warmup <- stan$warmup %||% 1000L
  fit_opts <- epinowcast::enw_fit_opts(
    sampler = epinowcast::enw_sample,
    chains = chains,
    iter_sampling = samples,
    iter_warmup = warmup,
    init_method = "prior",
    pp = TRUE,
    show_messages = verbose
  )

  model <- epinowcast::enw_model(threads = TRUE)

  # Run epinowcast
  enw_fit <- epinowcast::epinowcast(
    data = pobs,
    expectation = expectation_module,
    reference = reference_module,
    report = report_module,
    obs = obs_module,
    fit = fit_opts,
    model = model
  )

  # Convert output to EpiNow2 format
  enw_to_epinow2(enw_fit, data, generation_time, delays = delays)
}

#' Convert epinowcast output to EpiNow2 format
#'
#' @description Translates an epinowcast fit object into the
#'   `estimate_infections` S3 class structure expected by EpiNow2's
#'   downstream functions.
#'
#' @param enw_fit An epinowcast fit object.
#' @param original_data The original input data (with `date` and `confirm`).
#' @param generation_time The generation time options used.
#'
#' @return An `estimate_infections` object.
#'
#' @keywords internal
enw_to_epinow2 <- function(enw_fit, original_data, generation_time,
                           delays = NULL) {
  ret <- list(
    fit = enw_fit$fit[[1]],
    enw_fit = enw_fit,
    args = list(
      enw_data = enw_fit$data[[1]],
      generation_time = generation_time,
      delays = delays
    ),
    observations = original_data
  )
  class(ret) <- c("estimate_infections", "epinowfit", class(ret))
  ret
}

#' Extract posterior samples from an epinowcast-backed fit
#'
#' @description Extracts growth rates, infections, and reported cases from
#'   an epinowcast fit and formats them as a `data.table` matching
#'   EpiNow2's expected output structure.
#'
#' @param object An `estimate_infections` object produced by
#'   [run_epinowcast()].
#'
#' @return A `data.table` with columns: variable, time, date, sample,
#'   value, strat, type.
#'
#' @importFrom data.table data.table rbindlist
#' @keywords internal
extract_enw_samples <- function(object) {
  fit <- object$fit
  obs <- object$observations
  dates <- obs$date

  out <- list()

  # Growth rate (r) — epinowcast parameter name is "r"
  r_draws <- fit$draws(variables = "r", format = "draws_matrix")
  n_samples <- nrow(r_draws)
  n_r <- ncol(r_draws)
  out$growth_rate <- enw_draws_to_dt(
    r_draws, "growth_rate", dates[seq_len(n_r)]
  )

  # R = exp(r) when using generation time convolution
  out$R <- enw_draws_to_dt(
    exp(r_draws), "R", dates[seq_len(n_r)]
  )

  # Infections (exp_llatent = log expected latent observations)
  lat_draws <- fit$draws(variables = "exp_llatent", format = "draws_matrix")
  n_lat <- ncol(lat_draws)
  inf_dates <- dates[seq_len(min(n_lat, length(dates)))]
  out$infections <- enw_draws_to_dt(
    exp(lat_draws[, seq_len(length(inf_dates)), drop = FALSE]),
    "infections", inf_dates
  )

  # Reported cases (pp_inf_obs = posterior predictive for observations)
  tryCatch({
    pp_draws <- fit$draws(variables = "pp_inf_obs", format = "draws_matrix")
    n_pp <- ncol(pp_draws)
    pp_dates <- dates[seq_len(min(n_pp, length(dates)))]
    out$reported_cases <- enw_draws_to_dt(
      pp_draws[, seq_len(length(pp_dates)), drop = FALSE],
      "reported_cases", pp_dates
    )
  }, error = function(e) NULL)

  combined <- data.table::rbindlist(out, use.names = TRUE)
  combined[, strat := NA_character_]
  combined[, type := "estimate"]
  combined[]
}

#' Convert a draws matrix to a long data.table
#'
#' @param draws A matrix of posterior draws (samples x time points).
#' @param variable Character; the variable name.
#' @param dates Date vector for the time dimension.
#'
#' @return A `data.table` with columns: variable, time, date, sample, value.
#' @keywords internal
enw_draws_to_dt <- function(draws, variable, dates) {
  n_samples <- nrow(draws)
  n_times <- length(dates)
  data.table::data.table(
    variable = variable,
    time = rep(seq_len(n_times), each = n_samples),
    date = rep(dates, each = n_samples),
    sample = rep(seq_len(n_samples), times = n_times),
    value = as.vector(draws[, seq_len(n_times)])
  )
}

