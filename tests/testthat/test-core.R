test_that("registry: register, retrieve, and list models", {
  register_cr_model(
    key = "test_trivial",
    fit = function(x, y_time, y_event, args = list()) {
      list(intercept = mean(y_time), causes = sort(unique(y_event[y_event != 0])))
    },
    predict_cif = function(fit_obj, x_new, times) {
      n  <- nrow(x_new)
      K  <- length(fit_obj$causes)
      Tm <- length(times)
      array(0.1, dim = c(n, K, Tm))
    },
    info = function() list(
      name         = "Trivial test model",
      supports     = "CIF",
      needs_tuning = FALSE,
      default_grid = function() tibble::tibble()
    )
  )

  expect_true("test_trivial" %in% list_cr_models())

  mdl <- get_cr_model("test_trivial")
  expect_named(mdl, c("fit", "predict_cif", "info"))

  info <- mdl$info()
  expect_equal(info$name, "Trivial test model")
  expect_false(info$needs_tuning)
})

test_that("registry: error on unknown key", {
  expect_error(get_cr_model("no_such_model"), "not registered")
})

test_that("prepare_data renames columns and sorts by time", {
  df <- data.frame(t = c(3, 1, 2), e = c(1, 0, 2), x = 1:3)
  out <- prepare_data(df, time_col = "t", event_col = "e")
  expect_named(out, c("time", "event", "x"))
  expect_equal(out$time, c(1, 2, 3))
})

test_that("prepare_data offsets zero minimum time", {
  df <- data.frame(time = c(0, 1, 2), event = c(0, 1, 0))
  out <- prepare_data(df, "time", "event")
  expect_true(min(out$time) > 0)
})

test_that("pool_cifs_mean returns element-wise mean", {
  a <- array(1, dim = c(5, 2, 3))
  b <- array(3, dim = c(5, 2, 3))
  res <- pool_cifs_mean(list(a, b))
  expect_equal(res[1, 1, 1], 2)
})

test_that("trapezoidal.integration is correct", {
  # integral of f(x)=x from 0 to 1 = 0.5
  x <- seq(0, 1, by = 0.01)
  y <- x
  expect_equal(trapezoidal.integration(x, y), 0.5, tolerance = 1e-4)
})

test_that("WeightedBrierScore returns values in [0, 1]", {
  skip_if_not_installed("prodlim")
  set.seed(1)
  n      <- 100
  time   <- rexp(n)
  status <- sample(c(0L, 1L, 2L), n, replace = TRUE)
  preds  <- matrix(runif(n), ncol = 1)
  res    <- WeightedBrierScore(preds, tau = 0.5, time = time,
                                status = status, cause = 1L,
                                cens.code = 0L, cmprsk = TRUE)
  expect_true(all(res$weighted.brier.score >= 0))
})

test_that("sim_cmprks returns a data frame with correct columns", {
  df <- sim_cmprks(n = 200, seed = 42)
  expect_s3_class(df, "data.frame")
  expect_true("Tobs"  %in% names(df))
  expect_true("cause" %in% names(df))
  expect_true(all(df$cause %in% c(0L, 1L, 2L)))
})

test_that("built-in models are registered on package load", {
  expect_true("FGR"   %in% list_cr_models())
  expect_true("FGRP"  %in% list_cr_models())
  expect_true("csCPH" %in% list_cr_models())
  expect_true("RSF"   %in% list_cr_models())
})
