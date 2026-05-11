skip_on_cran()
# Setup for testing -------------------------------------------------------

futile.logger::flog.threshold("FATAL")

reported_cases <- EpiNow2::example_confirmed[1:30]

default_estimate_infections <- function(..., gt = TRUE, delay = TRUE) {
  futile.logger::flog.threshold("FATAL")

  stan_args <- stan_opts(
    chains = 2, warmup = 100, samples = 100,
    backend = "cmdstanr"
  )

  suppressWarnings(estimate_infections(...,
    generation_time = if (gt) {
      gt_opts(fix_parameters(example_generation_time))
    } else {
      gt_opts()
    },
    delays = if (delay) {
      delay_opts(fix_parameters(example_reporting_delay))
    } else {
      delay_opts()
    },
    rt = rt_opts(rw = 7),
    stan = stan_args, verbose = FALSE
  ))
}

test_estimate_infections <- function(...) {
  out <- default_estimate_infections(...)
  expect_true(all(c("fit", "args", "observations") %in% names(out)))
  expect_true(nrow(get_samples(out)) > 0)
  expect_true(nrow(summary(out, type = "parameters")) > 0)
  expect_true(nrow(out$observations) > 0)

  invisible(out)
}

# Integration tests (MCMC-based) ------------------------------------------

# Run MCMC once and reuse across multiple tests to save time
get_default_fit <- local({
  cached <- NULL
  function() {
    testthat::skip_on_os("windows")
    if (is.null(cached)) {
      cached <<- default_estimate_infections(reported_cases)
    }
    cached
  }
})

# Core test: Core functionality with default settings (always runs)
test_that("estimate_infections successfully returns estimates using default settings", {
  default_fit <- get_default_fit()
  expect_true(all(c("fit", "args", "observations") %in% names(default_fit)))
  expect_true(nrow(get_samples(default_fit)) > 0)
  expect_true(nrow(summary(default_fit, type = "parameters")) > 0)
  expect_true(nrow(default_fit$observations) > 0)
})

# Variant tests: Only run in full test mode
test_that("estimate_infections works using no delays", {
  skip_integration()
  test_estimate_infections(reported_cases, delay = FALSE)
})

test_that("estimate_infections works using the poisson observation model", {
  skip_integration()
  test_estimate_infections(reported_cases, obs = obs_opts(family = "poisson"))
})

test_that("estimate_infections works using a random walk", {
  skip_integration()
  test_estimate_infections(reported_cases, rt = rt_opts(rw = 7))
})

test_that("estimate_infections works without setting a generation time", {
  skip_integration()
  df <- test_estimate_infections(reported_cases, gt = FALSE, delay = FALSE)
  ## check exp(r) == R
  samples <- get_samples(df)
  growth_rate <- samples[variable == "growth_rate"][
    ,
    list(date, sample, growth_rate = value)
  ]
  R <- samples[variable == "R"][
    ,
    list(date, sample, R = value)
  ]
  combined <- merge(growth_rate, R, by = c("date", "sample"), all = FALSE)
  expect_equal(exp(combined$growth_rate), combined$R)
})

test_that("estimate_infections produces no forecasts when forecast = NULL", {
  skip_integration()
  out <- test_estimate_infections(data = reported_cases, forecast = NULL)
  samples <- get_samples(out)
  expect_true(!"forecast" %in% unique(samples$type))
})

test_that("estimate_infections produces no forecasts when forecast_opts horizon is 0", {
  skip_integration()
  out <- test_estimate_infections(
    data = reported_cases, forecast = forecast_opts(horizon = 0)
  )
  samples <- get_samples(out)
  expect_true(!"forecast" %in% unique(samples$type))
})

# Non-integration tests (fast - use one MCMC fit for multiple checks) ----

test_that("summary with type='parameters' returns all dates by default", {
  out <- get_default_fit()

  summ <- summary(out, type = "parameters")
  summ_dates <- unique(summ$date)
  expect_gt(length(summ_dates), 1)

  expect_true("infections" %in% summ$variable)
  expect_true("R" %in% summ$variable)

  # When target_date is explicitly provided, should filter to that date
  target <- summ_dates[length(summ_dates) %/% 2]
  summ_filtered <- summary(out, type = "parameters", target_date = target)
  expect_equal(unique(summ_filtered$date), target)
})

test_that("summary with type='parameters' has variable with semantic names", {
  out <- get_default_fit()

  summ <- summary(out, type = "parameters")
  expect_true("variable" %in% names(summ))
  expect_false("parameter" %in% names(summ))
  expect_true("R" %in% summ$variable)
})

test_that("get_predictions works with format='summary'", {
  out <- get_default_fit()

  preds <- get_predictions(out, format = "summary")

  expect_s3_class(preds, "data.table")
  expect_true("date" %in% names(preds))
  expect_true("mean" %in% names(preds))
  expect_true("median" %in% names(preds))
  expect_false("confirm" %in% names(preds))
})

test_that("get_predictions works with format='sample'", {
  out <- get_default_fit()

  preds <- get_predictions(out, format = "sample")

  expect_s3_class(preds, "data.table")
  expect_true(all(c(
    "forecast_date", "date", "horizon", "sample", "predicted"
  ) %in% names(preds)))
  expect_false("observed" %in% names(preds))
  expect_true(nrow(preds) > 0)
  expect_true(is.numeric(preds$predicted))
  expect_true(is.integer(preds$sample))
  expect_true(is.numeric(preds$horizon))
})

test_that("get_predictions works with format='quantile'", {
  out <- get_default_fit()

  preds <- get_predictions(out, format = "quantile")

  expect_s3_class(preds, "data.table")
  expect_true(all(c(
    "forecast_date", "date", "horizon", "quantile_level", "predicted"
  ) %in% names(preds)))
  expect_true(all(c(0.05, 0.25, 0.5, 0.75, 0.95) %in% preds$quantile_level))
})

test_that("get_predictions default format is 'summary'", {
  out <- get_default_fit()

  preds_default <- get_predictions(out)
  preds_explicit <- get_predictions(out, format = "summary")

  expect_equal(names(preds_default), names(preds_explicit))
  expect_equal(preds_default, preds_explicit)
})

test_that("get_predictions forecast_date equals last observation date", {
  out <- get_default_fit()

  preds <- get_predictions(out, format = "sample")

  expected_forecast_date <- max(out$observations$date, na.rm = TRUE)
  expect_equal(unique(preds$forecast_date), expected_forecast_date)
})

test_that("get_predictions horizon is correctly calculated", {
  out <- get_default_fit()

  preds <- get_predictions(out, format = "sample")

  expected_horizon <- as.numeric(preds$date - preds$forecast_date)
  expect_equal(preds$horizon, expected_horizon)
})

test_that("get_predictions format='sample' compatible with scoringutils", {
  skip_integration()
  skip_if_not_installed("scoringutils")

  fit <- estimate_infections(
    reported_cases,
    generation_time = gt_opts(fix_parameters(example_generation_time)),
    delays = delay_opts(fix_parameters(example_reporting_delay)),
    rt = rt_opts(rw = 7),
    stan = stan_opts(
      samples = 100, warmup = 100,
      chains = 2, backend = "cmdstanr"
    ),
    forecast = forecast_opts(horizon = 7),
    verbose = FALSE
  )

  preds <- get_predictions(fit, format = "sample")
  forecasts <- preds[horizon > 0]
  forecasts <- merge(
    forecasts, data.table::as.data.table(example_confirmed), by = "date"
  )

  forecast_obj <- scoringutils::as_forecast_sample(
    forecasts,
    forecast_unit = "horizon",
    observed = "confirm",
    sample_id = "sample"
  )
  expect_s3_class(forecast_obj, "forecast_sample")

  scores <- scoringutils::score(forecast_obj)
  expect_s3_class(scores, "data.table")
  expect_true("crps" %in% names(scores))
})

test_that("get_predictions format='quantile' compatible with scoringutils", {
  skip_integration()
  skip_if_not_installed("scoringutils")

  fit <- estimate_infections(
    reported_cases,
    generation_time = gt_opts(fix_parameters(example_generation_time)),
    delays = delay_opts(fix_parameters(example_reporting_delay)),
    rt = rt_opts(rw = 7),
    stan = stan_opts(
      samples = 100, warmup = 100,
      chains = 2, backend = "cmdstanr"
    ),
    forecast = forecast_opts(horizon = 7),
    verbose = FALSE
  )

  preds <- get_predictions(fit, format = "quantile")
  forecasts <- preds[horizon > 0]
  forecasts <- merge(
    forecasts, data.table::as.data.table(example_confirmed), by = "date"
  )

  forecast_obj <- scoringutils::as_forecast_quantile(
    forecasts,
    forecast_unit = "horizon",
    observed = "confirm",
    quantile_level = "quantile_level"
  )
  expect_s3_class(forecast_obj, "forecast_quantile")

  scores <- scoringutils::score(forecast_obj)
  expect_s3_class(scores, "data.table")
  expect_true("wis" %in% names(scores))
})

# Deprecation tests -------------------------------------------------------

test_that("summary.estimate_infections with type = 'samples' errors", {
  out <- get_default_fit()
  expect_error(summary(out, type = "samples"), "get_samples")
})

test_that("$samples accessor errors", {
  out <- get_default_fit()
  expect_error(out$samples, "get_samples")
})

test_that("$summarised accessor errors", {
  out <- get_default_fit()
  expect_error(out$summarised, "summary")
})

test_that("[[ accessor handles deprecated elements", {
  out <- get_default_fit()
  expect_error(out[["samples"]], "get_samples")
  expect_error(out[["summarised"]], "summary")

  expect_no_error(out[["fit"]])
  expect_no_error(out[["args"]])
  expect_no_error(out[["observations"]])
})

test_that("deprecated gp argument errors", {
  expect_error(
    estimate_infections(reported_cases, gp = NULL),
    "deprecated"
  )
})

test_that("deprecated backcalc argument errors", {
  expect_error(
    estimate_infections(reported_cases, backcalc = NULL),
    "deprecated"
  )
})

test_that("deprecated truncation argument errors", {
  expect_error(
    estimate_infections(reported_cases, truncation = NULL),
    "deprecated"
  )
})
