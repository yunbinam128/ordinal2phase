#' @title Multiple Imputation for Two-Phase Ordinal Data
#'
#' @description Implements the MI approach for ordinal outcomes in two-phase designs.
#'
#' @param formula The main model formula for the outcome.
#' @param formula_exposure The imputation model formula for the expensive exposure.
#' @param data_R1 A data frame for subjects selected into Phase 2 (complete-data).
#' @param data_R0 A data frame for subjects not selected into Phase 2 (missing expensive covariate).
#' @param ods_pi Sampling probabilities used for ascertainment correction.
#' @param n_candidates Number of grid points for expensive covariate in imputation.
#' @param n_impute Number of imputation datasets (default is 50).
#' @param verbose_TF Logical; if TRUE, prints progress every 25 imputations.
#' @param srs_TF Logical; set to TRUE if the design is Simple Random Sampling.
#'
#' @return A list containing the pooled coefficients and standard errors calculated via Rubin's Rules.
#'
#' @importFrom MASS mvrnorm
#' @importFrom sandwich vcovHC
#' @importFrom stats lm coef vcov as.formula
#' @importFrom rms orm
#' @export
mi_po <- function (formula, formula_exposure, data_R1, data_R0, ods_pi, n_candidates, n_impute=50, verbose_TF=FALSE, srs_TF=FALSE) {

  if (srs_TF == FALSE) {
    # 1. Fit ACML [Y|X,Z,W]
    mdl_acml <- suppressWarnings(
      acml_po(formula, data=data_R1, ac=TRUE, ods_pi=ods_pi))
    coef_acml <- mdl_acml$coef
    vcov_acml <- mdl_acml$vcov

    # 2. Fit IPW [X|Z,W]
    mdl_ipw <- lm(formula_exposure, data=data_R1, weights=ipw_weights)
    coef_ipw <- coef(mdl_ipw)
    vcov_ipw <- sandwich::vcovHC(mdl_ipw)
    sigma_hat_ipw <- sqrt(sum(mdl_ipw$residuals^2) / mdl_ipw$df.residual)
  } else {
    # 1. Fit ML [Y|X,Z,W]
    mdl_acml <- orm(formula, data_R1)
    coef_acml <- -mdl_acml$coef
    vcov_acml <- vcov(mdl_acml, intercepts="all")

    # 2. Fit Linear Regression [X|Z,W]
    mdl_ipw <- lm(formula_exposure, data=data_R1)
    coef_ipw <- coef(mdl_ipw)
    vcov_ipw <- vcov(mdl_ipw)
    sigma_hat_ipw <- sqrt(sum(mdl_ipw$residuals^2) / mdl_ipw$df.residual)
  }

  # Observed X values
  expensive_var <- all.vars(formula_exposure)[1]
  Xval_R1 <- unique(data_R1[[expensive_var]])
  Xcandidates <- seq(min(Xval_R1), max(Xval_R1), length.out=n_candidates)

  response <- all.vars(formula)[1]
  Y_R1 <- as.integer(data_R1[[response]]); Y_R0 <- as.integer(data_R0[[response]]); Y <- c(Y_R1, Y_R0)
  if (length(unique(Y)) > length(unique(Y_R0))) {
    k <- length(unique(Y))
    Y_uniq_R0 <- unique(Y_R0); Y_uniq_R0 <- Y_uniq_R0[order(Y_uniq_R0)]
    Y_uniq_R1 <- unique(Y_R1); Y_uniq_R1 <- Y_uniq_R1[order(Y_uniq_R1)]
    Y_uniq <- unique(Y); Y_uniq <- Y_uniq[order(Y_uniq)]
    cumY_R0 <- sapply(Y_uniq_R1, function(y) as.integer(as.integer(data_R0[[response]]) <= y))
  } else {
    Y_uniq <- unique(Y); Y_uniq <- Y_uniq[order(Y_uniq)]
    k <- length(Y_uniq)
    cumY_R0 <- sapply(Y_uniq, function(y) as.integer(as.integer(data_R0[[response]]) <= y))
  }

  Xmat_R0_acml <- stats::model.matrix(formula, data_R0)
  expensive_terms <- grep(expensive_var, colnames(Xmat_R0_acml), value = TRUE)
  if (length(expensive_terms) > 1) {
    pos <- regexpr(")", expensive_terms, fixed = TRUE)[1]
    term_tmp <- substr(expensive_terms, 1, pos)[1]
    formula_tmp <- as.formula(paste(
      response, "~",  term_tmp
    ))
    data_tmp <- data_R1[rep(1, n_candidates),]
    data_tmp[[expensive_var]] <- Xcandidates
    Xcandidates_mat <- stats::model.matrix(formula_tmp, data_tmp)[,-1]
  }
  Xmat_R0_ipw <- stats::model.matrix(formula_exposure, data_R0)

  # Create empty matrix and array to save results
  p <- length(coef_acml)
  coef_mat <- matrix(NA, nrow=n_impute, ncol=p)
  vcov_arr <- array(NA, dim=c(p, p, n_impute))
  coef_acml_m_repeats <- 0
  # 3. For m=1,...,M
  for (m in 1:n_impute) {
    # a. Draw coefficients
    coef_acml_m_sampled <- mvrnorm_restricted(mu=coef_acml, Sigma=vcov_acml, k)
    coef_acml_m <- coef_acml_m_sampled$coef
    coef_acml_m_repeats <- coef_acml_m_repeats + coef_acml_m_sampled$n_repeat

    coef_ipw_m <- MASS::mvrnorm(1, mu=coef_ipw, Sigma=vcov_ipw)

    # b. Perform imputation
    if (length(expensive_terms) > 1) {
      data_R0[[expensive_var]] <- sapply(seq_len(nrow(data_R0)), function(rowidx) {
        impute_Xobs(Xmati_ipw = Xmat_R0_ipw[rowidx,], Xmati_acml = Xmat_R0_acml[rowidx,], cumYi = cumY_R0[rowidx,],
                    coef_acml_m, coef_ipw_m, sigma_hat_ipw, Xcandidates, Xcandidates_mat)
      })
    } else {
      data_R0[[expensive_var]] <- sapply(seq_len(nrow(data_R0)), function(rowidx) {
        impute_Xobs(Xmati_ipw = Xmat_R0_ipw[rowidx,], Xmati_acml = Xmat_R0_acml[rowidx,], cumYi = cumY_R0[rowidx,],
                    coef_acml_m, coef_ipw_m, sigma_hat_ipw, Xcandidates, Xcandidates_mat = NULL)
      })
    }
    data_imputed <- rbind(data_R1[,colnames(data_R0)], data_R0)
    mdl_m <- orm(formula, data_imputed)
    coef_mat[m,] <- -coef(mdl_m)
    vcov_arr[,,m] <- vcov(mdl_m, intercepts="all")

    if (verbose_TF == TRUE & m %% 25 == 0) {cat(m, "th imputation done! \n")}
  }

  # 4. Aggregate results from `n_impute` sets of estimations using Rubin's rules
  coef_pool <- colMeans(coef_mat)
  varw <- apply(vcov_arr, 1:2, mean)
  e <- coef_mat - matrix(coef_pool, nrow=n_impute, ncol=p, byrow=TRUE)
  varb <- (t(e) %*% e)/(n_impute - 1)
  vcov_pool <- varw + (1 + 1/n_impute) * varb
  se_pool <- sqrt(diag(vcov_pool))
  # names(coef_pool) <- names(se_pool) <- c(paste0("y<=", Y_uniq[-k]), covariates)

  list(coef = coef_pool, se = se_pool, vcov = vcov_pool, coef_acml_repeats = coef_acml_m_repeats)
}

