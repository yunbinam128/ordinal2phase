#' @title Ascertainment Corrected Maximum Likelihood for Ordinal Outcomes
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param data A data frame for subjects selected into Phase 2 (complete-data).
#' @param Y Numeric vector of the ordinal outcome.
#' @param Xmat Numeric matrix of predictors.
#' @param ac Logical; if TRUE, performs ascertainment correction.
#' @param ods_pi Probability of being sampled into Phase 2 based on outcome/covariates.
#' @param minimizer String specifying the optimization method: "BFGS" (using optim) or "nlm".
#' @param theta_init Numeric vector of initial parameter values.
#' @param weights Numeric vector of sampling weights.
#'
#' @return A list containing:
#' \item{coef}{Estimated parameters (cutpoints and coefficients)}
#' \item{se}{Standard errors}
#' \item{vcov}{Variance-covariance matrix}
#' \item{loglik}{Log-likelihood at convergence}
#'
#' @importFrom stats optim nlm model.matrix complete.cases
#' @export
acml_po <- function (formula=NULL, data=NULL, Y, Xmat,
                     ac=FALSE, ods_pi=NULL, minimizer="BFGS", theta_init=NULL,
                     weights=NULL) {
  # Extract variables from formula
  if (!is.null(formula)) {
    data <- data[complete.cases(data),]
    covariates <- all.vars(formula)[-1]
    Xmat <- stats::model.matrix(formula, data)[,-1]  # (n_subjects, n_covariates)
    response <- all.vars(formula)[1]
    Y <- as.numeric(data[[response]])  # (n_subjects)
    if (length(Y) < nrow(data)) {stop("NA detected in Y")}
  }
  # If length(unique(Y_R1)) < length(unique(Y))
  Y_uniq <- unique(Y); Y_uniq <- Y_uniq[order(Y_uniq)]
  k <- length(unique(Y_uniq))
  if (!identical(seq_len(k), Y_uniq)) {
    Y_raw <- Y
    if(is.matrix(ods_pi)) {
      ods_pi <- ods_pi[,Y_uniq]
    } else {ods_pi <- ods_pi[Y_uniq]}
    Y <- as.numeric(factor(Y_raw, levels=Y_uniq, labels=1:k))
  }

  # Create ordinal outcome cumulative indicator
  cumY <- sapply(seq_len(k), function(j) {  # (n_subjects, k)
    indicator_cumYj <- as.integer(Y <= j)
    return(indicator_cumYj)
  }); colnames(cumY) <- paste0("cumY", 1:k)

  # Initialize theta
  nbetas <- ncol(Xmat)
  if (is.null(theta_init)) {
    theta_hat <- c(seq_len(k-1), rep(0, nbetas))
  } else {theta_hat <- theta_init}
  # Find MLE using minimizer
  if (minimizer == "BFGS") {
    fit_bfgs <- stats::optim(par=theta_hat, fn=neg_po_loglik, gr=neg_gr_po_loglik, method="BFGS",
                             Xmat=Xmat, cumY=cumY, ac=ac, ods_pi=ods_pi, weights=weights, hessian=TRUE)
    coef <- fit_bfgs$par
    loglik <- -fit_bfgs$value
    vcov <- solve(fit_bfgs$hessian)
  } else if(minimizer == "nlm"){
    fit_nlm <- stats::nlm(f=neg_po_loglik_nlm, p=theta_hat,
                          Xmat=Xmat, cumY=cumY, ac=ac, ods_pi=ods_pi, weights=weights, hessian=TRUE)
    coef <- fit_nlm$estimate
    loglik <- -fit_nlm$minimum
    vcov <- solve(fit_nlm$hessian)
  }

  # Calculate SE
  se <- sqrt(diag(vcov))
  list(
    coef = coef, se = se, vcov = vcov, loglik = loglik
  )
}
