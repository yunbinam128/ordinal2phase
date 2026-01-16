#' @title Sieve Maximum Likelihood Estimation for Ordinal Outcomes
#'
#' @description Implements the SMLE approach for two-phase sampling with ordinal outcomes.
#'
#' @param formula Model formula.
#' @param data_R1 A data frame for subjects selected into Phase 2.
#' @param data_R0 A data frame for subjects not selected into Phase 2.
#' @param Bspline_R1 Spline basis for subjects selected into Phase 2.
#' @param Bspline_R0 Spline basis for subjects not selected into Phase 2.
#' @param Xterm Name of the expensive covariate.
#' @param se_calc Logical, calculate standard errors?
#' @param theta_init Initial values for parameters.
#' @param hnscale Scale for the perturbation bandwidth.
#' @param tol Convergence tolerance.
#' @param maxiter Maximum EM iterations.
#' @param print_nsteps Frequency of printing progress.
#' @param verbose Logical, print detailed output?
#'
#' @importFrom MASS polr
#' @export
# ---------- Main fitter (EM + optional SE) ----------
smle_po <- function (formula, data_R1, data_R0, Bspline_R1, Bspline_R0, Xterm,
                     se_calc = TRUE, theta_init = NULL, hnscale = 1,
                     tol = 1e-4, maxiter = 1000, print_nsteps = 10, verbose = FALSE) {

  ws <- smle_build_workspace(formula, data_R1, data_R0, Bspline_R1, Bspline_R0, Xterm)

  # theta init
  nbetas <- ncol(ws$Xmat_R1)
  theta_hat <- if (is.null(theta_init)) c(seq_len(ws$nYcats - 1), rep(0, nbetas)) else theta_init

  # p_kj init
  p_kj <- ws$p_kj_init

  # indices to update weights
  r0_idx <- (ws$n_R1 + 1L):nrow(ws$design_df)

  # -- EM loop ----
  CONV_T <- CONV_P <- FALSE; nstep <- 1L
  while (nstep <= maxiter && !(CONV_T && CONV_P)) {
    ## -- E-step: calculate psi and q ----
    prob_po <- exp(po_loglik_vec_cpp(theta_hat, ws$X_terms, ws$cumY_mat))

    ws$p_mat[] <- as.matrix(p_kj[rep(seq_len(ws$m), times = ws$n_R0), ws$p_cols, drop = FALSE])
    resE <- estep_qvec_cpp(ws$B_mat, ws$p_mat, prob_po, ws$id_vec, ws$k_vec, ws$n_R0, ws$m)
    ws$design_df$w[r0_idx] <- resE$q_vec

    # -- M-step: update theta and p ----
    if (ws$nYcats < 2) stop("Need at least 2 outcome categories.")
    mod <- suppressWarnings(MASS::polr(formula, data = ws$design_df, weights = w))
    theta_new <- c(mod$zeta, -mod$coefficients)

    p_kj_new_num <- as.matrix(ws$p_kj_num_R1[, ws$p_cols, drop = FALSE]) + resE$p_kj_num_R0
    p_kj_new_denom <- ws$p_kj_denom_R1 + colSums(resE$p_kj_num_R0)
    p_kj_new_denom[p_kj_new_denom == 0] <- 1
    p_new_mat <- sweep(p_kj_new_num, 2, p_kj_new_denom, "/")
    p_new <- cbind.data.frame(X = ws$p_kj_num_R1[[ws$Xterm]], p_new_mat)
    colnames(p_new) <- c(ws$Xterm, ws$p_cols)

    # convergence checks
    CONV_T <- (max(abs(theta_hat - theta_new)) < tol)
    CONV_P <- (max(abs(p_kj[, ws$p_cols, drop = FALSE] - p_new[, ws$p_cols, drop = FALSE])) < tol)
    if (CONV_T && CONV_P) cat("Converged at nstep=", nstep, "\n")

    theta_hat <- theta_new
    p_kj <- p_new
    if (verbose) print(theta_hat)

    if (print_nsteps > 0 && nstep %% print_nsteps == 0 && verbose) cat(nstep, "th step done!\n")
    nstep <- nstep + 1L
  }
  if (!(CONV_T && CONV_P) && verbose) cat("Maximum number of iteration reached without convergence\n")

  if (!se_calc) {
    return(list(coef = theta_hat, nstep = nstep, se = rep(NA_real_, length(theta_hat))))
  }

  ws$p_mat[] <- as.matrix(p_kj[rep(seq_len(ws$m), times = ws$n_R0), ws$p_cols, drop = FALSE])

  # -- SE estimation via profile likelihood ----
  pl_0d <- smle_obs_loglik(theta_hat, p_kj, ws)
  nparams <- length(theta_hat)
  n_tot <- ws$n_R1 + ws$n_R0
  hn <- hnscale * (n_tot^(-1/2))
  e_mat <- diag(rep(hn, nparams))
  Hess_mat <- matrix(NA_real_, nrow = nparams, ncol = nparams)
  pl_2d <- matrix(NA_real_, nrow = nparams, ncol = nparams)
  pl_1d <- numeric(nparams)

  for (i in seq_len(nparams)) {
    for (j in i:nparams) {
      theta_pert <- as.numeric(theta_hat + e_mat[i, ] + e_mat[j, ])
      pl_2d[i, j] <- pl_2d[j, i] <- smle_profile_loglik(
        theta_hat_fixed = theta_pert, p_kj_init = p_kj, ws = ws,
        tol = tol, maxiter = maxiter, verbose = FALSE)
    }
  }
  for (i in seq_len(nparams)) {
    theta_pert <- as.numeric(theta_hat + e_mat[i, ])
    pl_1d[i] <- smle_profile_loglik(
      theta_hat_fixed = theta_pert, p_kj_init = p_kj, ws = ws,
      tol = tol, maxiter = maxiter, verbose = FALSE)
  }

  for (i in seq_len(nparams)) {
    for (j in i:nparams) {
      Hess_mat[i, j] <- Hess_mat[j, i] <- (pl_2d[i, j] - pl_1d[i] - pl_1d[j] + pl_0d)
    }
  }

  Cov_mat <- solve(-Hess_mat / (hn^2))
  se <- sqrt(diag(Cov_mat))
  list(coef = theta_hat, nstep = nstep, se = se, vcov = Cov_mat)
}
# ---- end smle_po.R ------------------------------------------------------
