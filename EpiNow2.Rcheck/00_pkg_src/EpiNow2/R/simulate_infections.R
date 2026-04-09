#' Simulate infections using the renewal equation
#'
#' Simulations are done from given initial infections and, potentially
#' time-varying, reproduction numbers. Delays and parameters of the observation
#' model can be specified using the same options as in [estimate_infections()].
#'
#' In order to simulate, all parameters that are specified such as the mean and
#' standard deviation of delays or observation scaling, must be fixed.
#' Uncertain parameters are not allowed.
#'
#' @param R a data frame of reproduction numbers (column `R`) by date (column
#'   `date`). Column `R` must be numeric and `date` must be in date format. If
#'   not all days between the first and last day in the `date` are present,
#'   it will be assumed that R stays the same until the next given date.
#' @param initial_infections numeric; the initial number of infections (i.e.
#'   before `R` applies). Note that results returned start the day after, i.e.
#'   the initial number of infections is not reported again. See also
#'   `seeding_time`
#' @param day_of_week_effect either `NULL` (no day of the week effect) or a
#'   numerical vector of length 7. Each element of the vector gives the weight
#'   given to reporting on this day (normalised to 1). The default is `NULL`.
#' @param seeding_time Integer; the number of days before the first time point
#'   of `R`; default is `NULL`, in which case it is set to the length of the
#'   generation time PMF.
#' @inheritParams estimate_infections
#' @inheritParams calc_CrIs
#' @importFrom checkmate assert_data_frame assert_date assert_numeric
#'   assert_subset assert_integerish
#' @importFrom data.table data.table merge.data.table nafill rbindlist
#' @importFrom cli cli_abort
#' @return A data.table of simulated infections (variable `infections`) and
#'   reported cases (variable `reported_cases`) by date.
#' @export
#' @examples
#' \donttest{
#' R <- data.frame(
#'   date = seq.Date(as.Date("2023-01-01"), length.out = 14, by = "day"),
#'   R = c(rep(1.2, 7), rep(0.8, 7))
#' )
#' sim <- simulate_infections(
#'   R = R,
#'   initial_infections = 100,
#'   generation_time = generation_time_opts(
#'     fix_parameters(example_generation_time)
#'   ),
#'   delays = delay_opts(fix_parameters(example_reporting_delay)),
#'   obs = obs_opts(family = "poisson")
#' )
#' }
simulate_infections <- function(R,
                                initial_infections,
                                day_of_week_effect = NULL,
                                generation_time = generation_time_opts(),
                                delays = delay_opts(),
                                obs = obs_opts(),
                                CrIs = c(0.2, 0.5, 0.9),
                                seeding_time = NULL) {
  assert_data_frame(R, any.missing = FALSE)
  assert_subset(c("date", "R"), colnames(R))
  assert_date(R$date)
  assert_numeric(R$R, lower = 0)
  assert_numeric(initial_infections, lower = 0)
  assert_numeric(day_of_week_effect, lower = 0, null.ok = TRUE)
  if (!is.null(seeding_time)) {
    assert_integerish(seeding_time, lower = 1)
  }
  assert_class(generation_time, "generation_time_opts")
  assert_class(delays, "delay_opts")
  assert_class(obs, "obs_opts")

  # Get generation time PMF
  gt_pmf <- gt_to_enw(generation_time)
  gt_len <- length(gt_pmf)

  # Get delay PMF
  delay_pmf <- delays_to_enw(delays)

  # Fill in R for all dates
  all_dates <- data.table::data.table(
    date = seq.Date(min(R$date), max(R$date), by = "day")
  )
  R <- data.table::merge.data.table(all_dates, R, by = "date", all.x = TRUE)
  R <- R[, R := data.table::nafill(R, type = "locf")]
  R <- R[!is.na(R)]

  if (is.null(seeding_time)) {
    seeding_time <- gt_len
  }

  # Seed infections using exponential growth implied by first R value
  r0 <- log(R$R[1]) / sum(seq_along(gt_pmf) * gt_pmf)
  seed_infections <- initial_infections * exp(r0 * seq(
    -(seeding_time - 1), 0
  ))

  # Simulate infections via renewal equation
  n_t <- nrow(R)
  infections <- numeric(seeding_time + n_t)
  infections[seq_len(seeding_time)] <- seed_infections

  for (t in seq_len(n_t)) {
    idx <- seeding_time + t
    past <- infections[max(1, idx - gt_len):(idx - 1)]
    gt_use <- rev(gt_pmf[seq_len(length(past))])
    infections[idx] <- R$R[t] * sum(past * gt_use)
  }

  # Convolve with reporting delay
  if (length(delay_pmf) > 1) {
    reported <- stats::convolve(infections, rev(delay_pmf), type = "open")
    reported <- reported[seq_len(length(infections))]
  } else {
    reported <- infections
  }

  # Apply day-of-week effect
  if (!is.null(day_of_week_effect)) {
    day_of_week_effect <- day_of_week_effect / sum(day_of_week_effect) *
      length(day_of_week_effect)
    all_dates_full <- seq.Date(
      min(R$date) - seeding_time, max(R$date), by = "day"
    )
    dow <- as.integer(format(all_dates_full, "%u"))
    reported <- reported * day_of_week_effect[dow]
  }

  # Apply observation model
  reported <- pmax(reported, 0)
  if (obs$family == "poisson") {
    reported_obs <- stats::rpois(length(reported), lambda = reported)
  } else {
    phi <- mean(fix_parameters(obs$dispersion))
    reported_obs <- stats::rnbinom(
      length(reported), mu = reported, size = 1 / phi^2
    )
  }

  # Build output
  dates <- c(
    seq(min(R$date) - seeding_time, min(R$date) - 1, by = "day"),
    R$date
  )
  reported_dates <- R$date

  inf_dt <- data.table::data.table(
    variable = "infections",
    date = dates,
    value = infections
  )
  rep_dt <- data.table::data.table(
    variable = "reported_cases",
    date = reported_dates,
    value = reported_obs[(seeding_time + 1):length(reported_obs)]
  )

  out <- data.table::rbindlist(list(inf_dt, rep_dt))
  out[]
}

#' Forecast infections
#'
#' @description `r lifecycle::badge("deprecated")`
#' @param estimates Unused.
#' @param ... Unused.
#' @export
forecast_infections <- function(estimates, ...) {
  lifecycle::deprecate_stop(
    "2.0.0", "forecast_infections()",
    details = paste(
      "This function is no longer supported since EpiNow2 v2.0.0.",
      "Use estimate_infections() with forecast = forecast_opts(horizon = N)",
      "to produce forecasts directly."
    )
  )
}
