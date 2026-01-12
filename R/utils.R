##### log-likelihood #####
# log of subject-specific AC part
#' @keywords internal
logac_i <- function(ods_pi, gammaj, gammajp1){
  # Calculate the probability mass function
  pmf <- gammajp1 - gammaj
  pmf <- c(gammaj[1], pmf)
  # Return AC part in individual log-likelihood
  log(sum(ods_pi * pmf))
}


# Subject-specific (ascertainment corrected) log-likelihood
#' @keywords internal
po_loglik_i <- function(theta_hat, Xmati, cumYi, ac, ods_pi, Yi=NULL, k=NULL){
  nbetas <- length(Xmati)
  if(is.null(cumYi) & !is.null(k)) {
    cumYi <- sapply(1:k, function(y) as.integer(Yi <= y))
  } else {
    k <- length(cumYi)
  }
  # If sampling strategy depends on a covariate term
  if (ac == TRUE) {
    if (is.list(ods_pi)) {
      conds <- strsplit(names(ods_pi), "=")
      for (c in 1:length(conds)) {
        cond <- gsub(" ", "", conds[[c]])
        if (Xmati[cond[1]] == cond[2]) {
          ods_pi <- ods_pi[[c]]
          break
        }}}}

  # Format theta_hat matrix
  a_ix <- seq(1, k-1)
  b_ix <- seq(k, k+nbetas-1)
  alphaj <- theta_hat[a_ix]
  betas <- theta_hat[b_ix]
  theta_mat <- cbind(alphaj,  # (k-1, 1 + n_betas)
                     matrix(rep(betas, k-1), nrow=k-1, byrow=T))
  Xidesign <- c(1, Xmati)

  # Define each component
  etaj <- c(theta_mat %*% Xidesign)
  gammaj <- 1 / (1+exp(-etaj))
  gammajp1 <- c(gammaj[-1], 1)
  phj <- log(gammaj / (gammajp1-gammaj))
  gphj <- log(1 + exp(phj))
  cumYj <- cumYi[-k]
  cumYjp1 <- c(cumYj[-1], 1)

  loglik <- sum(cumYj*phj - cumYjp1*gphj)

  if (ac == TRUE) {
    # Return log-likelihood for a single person with AC correction
    loglik - logac_i(ods_pi, gammaj, gammajp1) + log(ods_pi[min(which(cumYi==1))])
  } else {
    # Return log-likelihood for a single person without AC correction
    loglik
  }
}


# (Ascertainment corrected) Log-likelihood
#' @keywords internal
po_loglik <- function(theta_hat, Xmat, cumY, ac, ods_pi=NULL, weights=NULL, Y=NULL, k=NULL){
  n <- nrow(Xmat)
  if (!is.null(weights)) {
    if (is.matrix(ods_pi)) {
      sum(sapply(1:n, function (i) {
        weights[i] * po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi[i,], Yi=Y[i], k=k)
      }))
    } else {
      sum(sapply(1:n, function (i) {
        weights[i] * po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi, Yi=Y[i], k=k)
      }))
    }
  } else {
    if (is.matrix(ods_pi)) {
      sum(sapply(1:n, function (i) {
        po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi[i,], Yi=Y[i], k=k)
      }))
    } else {
      sum(sapply(1:n, function (i) {
        po_loglik_i(theta_hat, Xmati=as.numeric(Xmat[i,]), cumYi=cumY[i,], ac, ods_pi, Yi=as.numeric(Y[i]), k=k)
      }))
    }
  }
}


# Negative (ascertainment corrected) log-likelihood
#' @keywords internal
neg_po_loglik <- function (theta_hat, Xmat, cumY, ac, ods_pi, weights) {
  -po_loglik(theta_hat, Xmat, cumY, ac, ods_pi, weights)
}


##### Score #####
# Gradient of log of subject-specific AC part
#' @keywords internal
gr_logac_i <- function(ods_pi, gammaj, gammajp1, dejdt, dejp1dt){
  nparams <- ncol(dejdt)
  # Calculate denominator
  pmf <- gammajp1 - gammaj
  pmf <- c(gammaj[1], pmf)
  denom <- sum(ods_pi * pmf)

  # Calculate numerator (nparams)
  gammam <- c(gammaj, 1)
  gammamm1 <- c(0, gammaj)
  demdt <- rbind(dejdt, matrix(0, ncol=nparams))
  demm1dt <- rbind(matrix(0, ncol=nparams), dejdt)
  gammam_rep <- matrix(rep(gammam, nparams), ncol=nparams)
  gammamm1_rep <- matrix(rep(gammamm1, nparams), ncol=nparams)
  ods_pi_rep <- matrix(rep(ods_pi, nparams), ncol=nparams)
  num <- colSums(
    ods_pi_rep *
      (demdt*gammam_rep*(1-gammam_rep) - demm1dt*gammamm1_rep*(1-gammamm1_rep))
  )

  return(num/c(denom))
}


# Subject-specific gradient of (ascertainment corrected) log-likelihood
#' @keywords internal
gr_po_loglik_i <- function(theta_hat, Xmati, cumYi, ac, ods_pi){
  nbetas <- length(Xmati)
  k <- length(cumYi)
  nparams <- (k-1) + nbetas
  # If sampling strategy depends on a covariate term
  if (ac == TRUE) {
    if (is.list(ods_pi)) {
      conds <- strsplit(names(ods_pi), "=")
      for (c in 1:length(conds)) {
        cond <- gsub(" ", "", conds[[c]])
        if (Xmati[cond[1]] == cond[2]) {
          ods_pi <- ods_pi[[c]]
          break
        }}}}

  # Format theta_hat matrix
  a_ix <- seq(1, k-1)
  b_ix <- seq(k, nparams)
  alphaj <- theta_hat[a_ix]
  betas <- theta_hat[b_ix]
  theta_mat <- cbind(alphaj,  # (k-1, 1+n_betas)
                     matrix(rep(betas, k-1), nrow=k-1, byrow=T))
  Xidesign <- c(1, Xmati)

  # Define each component
  etaj <- c(theta_mat %*% Xidesign)
  gammaj <- 1 / (1+exp(-etaj))
  gammajp1 <- c(gammaj[-1], 1)
  phj <- log(gammaj / (gammajp1-gammaj))
  gphj <- log(1+exp(phj))
  egphj <- exp(gphj)
  cumYj <- cumYi[-k]
  cumYjp1 <- c(cumYj[-1], 1)

  dejdt <- cbind(diag(k-1),  # for parameter alphaa: dejdt=1 for j=a and 0 for j/=a)
                 matrix(t(Xmati)[rep(1, k-1),], nrow=k-1))  # for parmeter betas: dejdt=x
  dejp1dt <- rbind(dejdt[-1,],
                   matrix(0, ncol=ncol(dejdt)))

  egphj_rep <- matrix(rep(egphj, nparams), ncol=nparams)
  gammaj_rep <- matrix(rep(gammaj, nparams), ncol=nparams)
  gammajp1_rep <- matrix(rep(gammajp1, nparams), ncol=nparams)
  cumYj_rep <- matrix(rep(cumYj, nparams), ncol=nparams)
  cumYjp1_rep <- matrix(rep(cumYjp1, nparams), ncol=nparams)

  score_i <- colSums(egphj_rep *
                       (dejdt*(1-gammaj_rep) - dejp1dt*(1-gammajp1_rep)) *
                       (cumYj_rep - cumYjp1_rep*gammaj_rep/gammajp1_rep))

  if (length(b_ix) > 0) {
    names(score_i) <- c(paste0("alpha", 1:length(a_ix)), paste0("beta", names(Xmati)))
  } else {
    names(score_i) <- c(paste0("alpha", 1:length(a_ix)))
  }

  if (ac == TRUE) {
    # Return score for a single person with AC correction
    score_i - gr_logac_i(ods_pi, gammaj, gammajp1, dejdt, dejp1dt)
  } else {
    # Return score for a single person without AC correction
    score_i
  }
}


# Gradient of (ascertainment corrected) log-likelihood
#' @keywords internal
gr_po_loglik <- function (theta_hat, Xmat, cumY, ac, ods_pi, weights) {
  n <- nrow(Xmat)
  if (!is.null(weights)) {
    if (is.matrix(ods_pi)) {
      rowSums(sapply(1:n, function (i) {
        weights[i] * gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi[i,])
      }))
    } else {
      rowSums(sapply(1:n, function (i) {
        weights[i] * gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi)
      }))
    }
  } else {
    if (is.matrix(ods_pi)) {
      rowSums(sapply(1:n, function (i) {
        gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi[i,])
      }))
    } else {
      rowSums(sapply(1:n, function (i) {
        gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi)
      }))
    }
  }
}


#' @keywords internal
cov_gr_po_loglik <- function (theta_hat, Xmat, cumY, ac, ods_pi, weights) {
  n <- nrow(Xmat)
  if (!is.null(weights)) {
    Reduce("+", lapply(1:n, function (i) {
      gi <- weights[i] * gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi)
      gi %*% t(gi)
    }))
  } else {
    Reduce("+", lapply(1:n, function (i) {
      gi <- gr_po_loglik_i(theta_hat, Xmati=Xmat[i,], cumYi=cumY[i,], ac, ods_pi)
      gi %*% t(gi)
    }))
  }
}


# Negative gradient of (ascertainment corrected) log-likelihood
#' @keywords internal
neg_gr_po_loglik <- function (theta_hat, Xmat, cumY, ac, ods_pi, weights) {
  -gr_po_loglik(theta_hat, Xmat, cumY, ac, ods_pi, weights)
}


##### Model fitting #####
# Negative PO log-likelihood with attribute of gradient
#' @keywords internal
neg_po_loglik_nlm <- function (theta_hat, Xmat, cumY, ac, ods_pi, weights) {
  res <- neg_po_loglik(theta_hat, Xmat, cumY, ac, ods_pi, weights)
  attr(res, "gradient") <- neg_gr_po_loglik(theta_hat, Xmat, cumY, ac, ods_pi, weights)
  res
}


# Calculate bread by using numerical hessian
#' @keywords internal
vcov_po_ac <- function (theta_hat, formula=NULL, data=NULL, Xmat, cumY,
                        ac=FALSE, ods_pi, s=1e-6, robust_se=FALSE,
                        N_subsamples=NULL, n_subsamples=NULL) {
  # Extract variables from formula
  if (!is.null(formula)) {
    response <- all.vars(formula)[1]
    Y <- as.numeric(data[[response]])  # (n_subjects)
    covariates <- all.vars(formula)[-1]
    Xmat <- stats::model.matrix(formula, data)[,-1]  # (n_subjects, n_covariates)
  }
  k <- length(unique(Y))  # number of ordinal outcome categories
  # Create ordinal outcome cumulative indicator
  cumY <- sapply(1:k, function(j) {  # (n_subjects, k)
    indicator_cumYj <- as.integer(Y <= j)
    return(indicator_cumYj)
  }); colnames(cumY) <- paste0("cumY", 1:k)

  score_theta_hat <- gr_po_loglik(theta_hat, Xmat, cumY, ac, ods_pi, weights)

  # Numerical hessian
  hess_mat <- matrix(0, length(theta_hat), length(theta_hat))
  for (i in 1:length(theta_hat)) {
    theta_hat_s <- theta_hat
    theta_hat_s[i] <- theta_hat_s[i] + s
    score_theta_hat_s <- gr_po_loglik(theta_hat_s, Xmat, cumY, ac, ods_pi, weights)
    hess_mat[, i] <- (score_theta_hat - score_theta_hat_s)/s
  }

  if (is.null(n_subsamples)) {
    solve(hess_mat)
  }
}


#' @keywords internal
mvrnorm_restricted <- function(mu, Sigma, k) {  # mvrnorm function where drawn intercepts are monotonically increasing
  n_repeat <- 0
  repeat {
    coef <- MASS::mvrnorm(1, mu = mu, Sigma = Sigma)
    n_repeat <- n_repeat + 1
    if (all(diff(coef[1:(k-1)]) > 0)) {return(list(coef = coef, n_repeat = n_repeat))}
  }
}


#' @keywords internal
impute_Xobs <- function (Xmati_ipw, Xmati_acml, cumYi, coef_acml_m, coef_ipw_m, sigma_hat_ipw, Xcandidates, Xcandidates_mat) {

  # Calculate P(X=Xcandidatesi|Z, W)
  mui <- c(coef_ipw_m %*% Xmati_ipw)
  prob_ipw <- stats::dnorm(Xcandidates, mean = mui, sd = sigma_hat_ipw)

  # Calculate P(Y|X=Xval_R1i, Z, W)
  k <- length(cumYi)

  n_candidates <- length(Xcandidates)
  if (is.null(Xcandidates_mat)) {
    Xmati_acml_all <- cbind(Xcandidates, matrix(rep(Xmati_acml[names(Xmati_ipw)[-1]], n_candidates), nrow = n_candidates, byrow = TRUE))
  } else {
    Xmati_acml_all <- cbind(Xcandidates_mat, matrix(rep(Xmati_acml[-c(1:3)], n_candidates), nrow = n_candidates, byrow = TRUE))
  }
  prob_acml <- exp(apply(Xmati_acml_all, 1, function (row) {
    po_loglik_i(coef_acml_m, row, cumYi, ac=FALSE)
  }))

  # Calculate the distribution of Xcandidates and the expected X
  num <- prob_acml * prob_ipw
  denom <- sum(num)
  distn <- num / denom

  sample(Xcandidates, size = 1, prob = distn)
}


# -- Workspace builder (allocates big objects once) ----
#' @keywords internal
smle_build_workspace <- function(formula, data_R1, data_R0,
                                 Bspline_R1, Bspline_R0, Xterm) {
  covariates <- setdiff(all.vars(formula)[-1], "parms_list")
  response <- all.vars(formula)[1]

  Xmat_R1 <- stats::model.matrix(formula, data_R1)[, -1, drop = FALSE]  # (n_R1, n_covariates)
  if (nrow(Xmat_R1) < nrow(data_R1)) stop("NA detected in Xmat_R1")
  terms <- colnames(Xmat_R1)

  Y_R1 <- as.integer(data_R1[[response]])  # (n_R1)
  Y_R0 <- as.integer(data_R0[[response]])  # (n_R0)
  Y <- c(Y_R1, Y_R0)
  Y_uniq <- sort(unique(Y))
  nYcats <- length(Y_uniq)  # number of categories of Y
  if (length(Y_R1) < nrow(data_R1)) stop("NA detected in Y_R1")

  # keep only needed cols
  data_R1_small <- data_R1[, c(response, covariates), drop = FALSE]
  data_R0_small <- data_R0[, c(response, covariates[covariates != Xterm]), drop = FALSE]

  # splines
  n_sieve <- ncol(Bspline_R1)
  B_cols <- paste0("B_", seq_len(n_sieve))
  colnames(Bspline_R1) <- colnames(Bspline_R0) <- B_cols
  data_R1_small <- cbind(data_R1_small, Bspline_R1)  # (n_R1, + n_sieve)
  data_R0_small <- cbind(data_R0_small, Bspline_R0)  # (n_R0, + n_sieve)

  # cumY for R0 (`nYcats` columns)
  cumY_R0 <- outer(Y_R0, Y_uniq, "<=") * 1L
  cumY_cols <- paste0("cumY", Y_uniq)
  colnames(cumY_R0) <- cumY_cols
  data_R0_small <- cbind(data_R0_small, cumY_R0)

  # R1 info
  B_R1_mat <- as.matrix(data_R1_small[, B_cols, drop = FALSE])
  X_vec_R1_num <- data_R1_small[[Xterm]]
  Y_R1_int <- Y_R1

  # Initial p grid levels (k)
  pj_num_mat <- rowsum(as.matrix(
    data_R1_small[, B_cols, drop = FALSE]), group = X_vec_R1_num, reorder = FALSE)
  X_levels <- as.numeric(rownames(pj_num_mat))  # character keeps safe type
  p_cols <- paste0("p_k", seq_len(n_sieve))
  p_kj_num_R1 <- cbind.data.frame(X = X_levels, pj_num_mat, row.names = NULL)
  colnames(p_kj_num_R1) <- c(Xterm, p_cols)
  p_kj_denom_R1 <- colSums(p_kj_num_R1[, -1, drop = FALSE])
  p_kj_denom_R1[p_kj_denom_R1 == 0] <- 1
  p_init_mat <- sweep(p_kj_num_R1[, -1, drop = FALSE], 2, p_kj_denom_R1, "/")
  p_kj_init <- cbind.data.frame(p_kj_num_R1[, Xterm, drop = FALSE], p_init_mat)  # (m, n_sieve), Initial p_kj values

  # Replicate R0 over k
  m <- length(X_levels)
  n_R0 <- nrow(data_R0_small)
  id_rep <- rep(seq_len(n_R0), each = m)
  k_rep <- rep(seq_len(m), times = n_R0)

  data_R0_rep <- cbind(  # replicate data_R0 rows by k
    id = id_rep, k = k_rep,
    data_R0_small[id_rep, , drop = FALSE],
    p_kj_init[k_rep, , drop = FALSE]
  )
  Xmat_R0_rep <- cbind(  # model matrix for the replicated R0
    data_R0_rep[, !names(data_R0_rep) %in% covariates, drop = FALSE],
    stats::model.matrix(formula, data_R0_rep)[, -1, drop = FALSE]
  )

  # slices for kernels
  idx_terms <- match(terms, colnames(Xmat_R0_rep))
  idx_cumY <- match(cumY_cols, colnames(Xmat_R0_rep))
  X_terms <- as.matrix(Xmat_R0_rep[, idx_terms, drop = FALSE])
  cumY_mat <- as.matrix(Xmat_R0_rep[, idx_cumY, drop = FALSE])
  B_mat <- as.matrix(data_R0_rep[, B_cols, drop = FALSE])
  p_mat <- matrix(0.0, nrow = nrow(B_mat), ncol = ncol(B_mat))
  colnames(p_mat) <- colnames(B_mat)

  # modeling df + weights (built once)
  design_df <- rbind(
    data_R1_small[, c(response, covariates), drop = FALSE],
    data_R0_rep[, c(response, covariates), drop = FALSE]
  )
  design_df[[response]] <- factor(design_df[[response]])
  # add weight column (R0 rows initialized to 0, R1 rows to 1)
  design_df$w <- 0
  n_R1 <- nrow(data_R1_small)
  design_df$w[seq_len(n_R1)] <- 1

  structure(list(
    response = response, covariates = covariates, Xterm = Xterm,
    Y_uniq = Y_uniq, nYcats = nYcats,
    B_cols = B_cols, p_cols = p_cols, cumY_cols = cumY_cols,
    n_R0 = n_R0, n_R1 = n_R1, m = m, n_sieve = n_sieve,
    # R1 caches
    Xmat_R1 = Xmat_R1, B_R1_mat = B_R1_mat, X_vec_R1_chr = as.character(X_vec_R1_num), Y_R1_int = Y_R1_int,
    # replicated R0 slices
    data_R0_rep = data_R0_rep,
    X_terms = X_terms, cumY_mat = cumY_mat,
    B_mat = B_mat, p_mat = p_mat,
    id_vec = as.integer(data_R0_rep$id),
    k_vec = as.integer(data_R0_rep$k),
    # p normalizers (R1 part)
    p_kj_init = p_kj_init, p_kj_num_R1 = p_kj_num_R1, p_kj_denom_R1 = p_kj_denom_R1,
    # modeling data frame
    design_df = design_df
  ), class = "smle_ws")
}


# -- Observed log-likelihood using workspace ----
#' @keywords internal
smle_obs_loglik <- function(theta_hat, p_kj, ws) {
  # R1 contribution
  idx <- match(ws$X_vec_R1_chr, as.character(p_kj[[ws$Xterm]]))
  p_R1 <- as.matrix(p_kj[idx, ws$p_cols, drop = FALSE])
  p_R1[p_R1 == 0] <- 1
  cumY_R1 <- outer(ws$Y_R1_int, ws$Y_uniq, "<=") * 1L
  loglik_R1_rows <- po_loglik_vec_cpp(theta_hat, ws$Xmat_R1, cumY_R1)
  loglik_R1 <- sum(loglik_R1_rows) + sum(ws$B_R1_mat * log(p_R1))

  # R0 contribution
  ws$p_mat[] <- as.matrix(p_kj[rep(seq_len(ws$m), times = ws$n_R0), ws$p_cols, drop = FALSE])
  prob_po <- exp(po_loglik_vec_cpp(theta_hat, ws$X_terms, ws$cumY_mat))
  loglik_R0 <- loglik_R0_cpp(prob_po, ws$B_mat, ws$p_mat, ws$id_vec, ws$n_R0)

  loglik_R0 + loglik_R1
}


# -- Profile log-likelihood (theta fixed) using workspace ----
#' @keywords internal
smle_profile_loglik <- function(theta_hat_fixed, p_kj_init, ws,
                                tol = 1e-4, maxiter = 1000, verbose = FALSE) {
  p_kj <- p_kj_init
  ws$p_mat[] <- as.matrix(p_kj[rep(seq_len(ws$m), times = ws$n_R0), ws$p_cols, drop = FALSE])
  CONVERGED <- FALSE; nstep <- 1L

  while (nstep <= maxiter && !CONVERGED) {
    prob_po <- exp(po_loglik_vec_cpp(theta_hat_fixed, ws$X_terms, ws$cumY_mat))
    resE <- estep_qvec_cpp(ws$B_mat, ws$p_mat, prob_po, ws$id_vec, ws$k_vec, ws$n_R0, ws$m)
    q_vec <- resE$q_vec
    pkj_R0 <- resE$p_kj_num_R0

    # update p
    p_kj_new_num <- as.matrix(ws$p_kj_num_R1[, ws$p_cols, drop = FALSE]) + pkj_R0
    p_kj_new_denom <- ws$p_kj_denom_R1 + colSums(pkj_R0)
    p_kj_new_denom[p_kj_new_denom == 0] <- 1
    p_new_mat <- sweep(p_kj_new_num, 2, p_kj_new_denom, "/")
    p_kj_new  <- cbind.data.frame(X = ws$p_kj_num_R1[[ws$Xterm]], p_new_mat)
    colnames(p_kj_new) <- c(ws$Xterm, ws$p_cols)

    CONVERGED <- (max(abs(p_kj[, ws$p_cols, drop = FALSE] - p_kj_new[, ws$p_cols, drop = FALSE])) < tol)
    if (verbose && nstep %% 25L == 0L) cat(nstep, " steps...\n")

    p_kj <- p_kj_new
    ws$p_mat[] <- as.matrix(p_kj[rep(seq_len(ws$m), times = ws$n_R0), ws$p_cols, drop = FALSE])
    nstep <- nstep + 1L
  }
  if (!CONVERGED && verbose) cat("profile: reached maxiter without convergence\n")

  smle_obs_loglik(theta_hat_fixed, p_kj, ws)
}
