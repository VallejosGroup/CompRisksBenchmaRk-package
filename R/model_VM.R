################################################################################
# Vertical Model for Competing Risks
# Nicolaie, van Houwelingen & Putter (2010), Stat Med 29:1190-1205
#
# Registered keys:
#   "VM_spline"    - natural cubic spline time basis (recommended)
#   "VM_piecewise" - piecewise-constant time basis
#
# NOTE: This model requires `vertical_modeling.R` to be sourced (defines
# vertical_model(), predict_vertical(), vm_to_df()).
# Place vertical_modeling.R in the package as R/vertical_modeling.R and
# remove the source() call, or source() it manually before using these keys.
################################################################################

# source("../VerticalModelling_Trials/vertical_modeling.R")
# Uncomment the line above OR ensure vertical_model() / predict_vertical()
# are available in the search path before calling these registered models.

register_cr_model(
  key = "VM_spline",

  fit = function(x, y_time, y_event, args = list()) {
    x       <- as.data.frame(x)
    y_time  <- as.numeric(y_time)
    y_event <- as.integer(as.character(y_event))
    causes  <- sort(unique(y_event[y_event != 0L]))
    K       <- length(causes)

    if (K < 2L)
      stop("VM_spline: need at least 2 competing causes; found ", K, ".")

    args <- utils::modifyList(
      list(basis              = "spline",
           knots              = NULL,
           interact_with_time = FALSE,
           scale_continuous   = TRUE,
           cause_labels       = NULL),
      args)

    if (!is.null(args$cause_labels) && length(args$cause_labels) != K) {
      warning("VM_spline: args$cause_labels length != number of causes; ignored.")
      args$cause_labels <- NULL
    }

    has_cov <- ncol(x) > 0L
    if (has_cov) {
      is_const <- vapply(x, function(v) {
        uv <- unique(suppressWarnings(as.numeric(v)))
        length(uv[is.finite(uv)]) <= 1L
      }, logical(1L))
      if (any(is_const)) {
        warning("VM_spline: dropping constant column(s): ",
                paste(names(x)[is_const], collapse = ", "))
        x <- x[, !is_const, drop = FALSE]
      }
      has_cov <- ncol(x) > 0L
    }

    dat        <- cbind(data.frame(.time = y_time, .event = y_event,
                                   stringsAsFactors = FALSE), x)
    vm_formula <- prodlim::Hist(.time, .event) ~ 1

    vm <- vertical_model(
      formula            = vm_formula,
      data               = dat,
      multinom_formula   = if (has_cov) stats::reformulate(names(x)) else NULL,
      cause_labels       = args$cause_labels,
      basis              = args$basis,
      knots              = args$knots,
      interact_with_time = args$interact_with_time,
      scale_continuous   = args$scale_continuous
    )

    structure(
      list(vm           = vm,
           causes       = causes,
           cause_labels = vm$cause_labels,
           has_cov      = has_cov,
           covar_names  = if (has_cov) names(x) else character(0L)),
      class = c("cr_model_vm", "cr_model")
    )
  },

  predict_cif = function(fit_obj, x_new, times) {
    x_new <- as.data.frame(x_new)
    times <- as.numeric(times)
    n     <- nrow(x_new)
    K     <- length(fit_obj$causes)
    Tm    <- length(times)
    vm    <- fit_obj$vm

    out <- array(0, dim     = c(n, K, Tm),
                 dimnames = list(NULL, fit_obj$cause_labels, as.character(times)))

    if (!fit_obj$has_cov) {
      for (j in seq_len(K)) {
        cif_j <- stats::approx(vm$event_times, vm$F_j[, j],
                               xout = times, method = "constant",
                               f = 0, rule = 2)$y
        out[, j, ] <- matrix(cif_j, nrow = n, ncol = Tm, byrow = TRUE)
      }
    } else {
      for (i in seq_len(n)) {
        nd <- as.list(x_new[i, fit_obj$covar_names, drop = FALSE])
        pr <- predict_vertical(vm, newdata = nd)
        for (j in seq_len(K)) {
          out[i, j, ] <- stats::approx(pr$event_times, pr$F_j[, j],
                                       xout = times, method = "constant",
                                       f = 0, rule = 2)$y
        }
      }
    }
    out
  },

  info = function() list(
    name         = "Vertical Model \u2014 spline basis (Nicolaie et al. 2010)",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble(
      basis              = "spline",
      knots              = list(NULL),
      interact_with_time = FALSE,
      scale_continuous   = TRUE
    )
  )
)


register_cr_model(
  key = "VM_piecewise",

  fit = function(x, y_time, y_event, args = list()) {
    x       <- as.data.frame(x)
    y_time  <- as.numeric(y_time)
    y_event <- as.integer(as.character(y_event))
    causes  <- sort(unique(y_event[y_event != 0L]))
    K       <- length(causes)

    if (K < 2L)
      stop("VM_piecewise: need at least 2 competing causes; found ", K, ".")

    args <- utils::modifyList(
      list(basis              = "piecewise",
           knots              = NULL,
           interact_with_time = FALSE,
           scale_continuous   = TRUE,
           cause_labels       = NULL),
      args)

    if (!is.null(args$cause_labels) && length(args$cause_labels) != K) {
      warning("VM_piecewise: args$cause_labels length != number of causes; ignored.")
      args$cause_labels <- NULL
    }

    has_cov <- ncol(x) > 0L
    if (has_cov) {
      is_const <- vapply(x, function(v) {
        uv <- unique(suppressWarnings(as.numeric(v)))
        length(uv[is.finite(uv)]) <= 1L
      }, logical(1L))
      if (any(is_const)) {
        warning("VM_piecewise: dropping constant column(s): ",
                paste(names(x)[is_const], collapse = ", "))
        x <- x[, !is_const, drop = FALSE]
      }
      has_cov <- ncol(x) > 0L
    }

    dat        <- cbind(data.frame(.time = y_time, .event = y_event,
                                   stringsAsFactors = FALSE), x)
    vm_formula <- prodlim::Hist(.time, .event) ~ 1

    vm <- vertical_model(
      formula            = vm_formula,
      data               = dat,
      multinom_formula   = if (has_cov) stats::reformulate(names(x)) else NULL,
      cause_labels       = args$cause_labels,
      basis              = args$basis,
      knots              = args$knots,
      interact_with_time = args$interact_with_time,
      scale_continuous   = args$scale_continuous
    )

    structure(
      list(vm           = vm,
           causes       = causes,
           cause_labels = vm$cause_labels,
           has_cov      = has_cov,
           covar_names  = if (has_cov) names(x) else character(0L)),
      class = c("cr_model_vm", "cr_model")
    )
  },

  predict_cif = function(fit_obj, x_new, times) {
    x_new <- as.data.frame(x_new)
    times <- as.numeric(times)
    n     <- nrow(x_new)
    K     <- length(fit_obj$causes)
    Tm    <- length(times)
    vm    <- fit_obj$vm

    out <- array(0, dim     = c(n, K, Tm),
                 dimnames = list(NULL, fit_obj$cause_labels, as.character(times)))

    if (!fit_obj$has_cov) {
      for (j in seq_len(K)) {
        cif_j <- stats::approx(vm$event_times, vm$F_j[, j],
                               xout = times, method = "constant",
                               f = 0, rule = 2)$y
        out[, j, ] <- matrix(cif_j, nrow = n, ncol = Tm, byrow = TRUE)
      }
    } else {
      for (i in seq_len(n)) {
        nd <- as.list(x_new[i, fit_obj$covar_names, drop = FALSE])
        pr <- predict_vertical(vm, newdata = nd)
        for (j in seq_len(K)) {
          out[i, j, ] <- stats::approx(pr$event_times, pr$F_j[, j],
                                       xout = times, method = "constant",
                                       f = 0, rule = 2)$y
        }
      }
    }
    out
  },

  info = function() list(
    name         = "Vertical Model \u2014 piecewise basis (Nicolaie et al. 2010)",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble(
      basis              = "piecewise",
      knots              = list(NULL),
      interact_with_time = FALSE,
      scale_continuous   = TRUE
    )
  )
)

message("Vertical Model registered: keys = 'VM_spline', 'VM_piecewise'")
