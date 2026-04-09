# Support `max_delay = 1` in epinowcast

## Problem

`enw_preprocess_data(obs, max_delay = 1)` fails when `reference_date == report_date` for all rows. This blocks using epinowcast for Rt estimation from simple incidence data without reporting delays.

## What to change

### 1. `R/check.R` — `check_timestep_by_date()`

The check at line 695 requires ≥2 rows per reference_date to infer timestep from report_date spacing. With `max_delay = 1` there's only one row per reference_date (delay 0). Skip this check when `max_delay = 1` — add a `max_delay` argument and pass it from `enw_preprocess_data()`.

### 2. `R/model-modules.R` — `enw_reference()`

Lines 70-76 error when both `parametric = ~0` and `non_parametric = ~0`. Remove this error. The Stan model already handles `model_refp = 0` and `model_refnp = 0` — both branches are gated by `if (model_refp)` / `if (model_refnp)`.

### 3. `R/preprocess.R` — `enw_preprocess_data()`

The call to `check_timestep_by_date()` at line 1207 happens before `max_delay` is resolved. Either move it after, or resolve `max_delay` earlier, so it can be passed to the check.

### 4. Verify downstream

Check that `enw_reporting_triangle()`, `enw_latest_data()`, and `enw_filter_delay()` handle a single-column reporting triangle (delay 0 only) without error.

## Test

```r
obs <- data.frame(
  reference_date = as.Date("2021-01-01") + 0:9,
  report_date = as.Date("2021-01-01") + 0:9,
  confirm = cumsum(rpois(10, 100))
)
# Should succeed:
pobs <- enw_preprocess_data(obs, max_delay = 1)
# Should succeed with no delay model:
ref <- enw_reference(parametric = ~0, non_parametric = ~0, data = pobs)
```
