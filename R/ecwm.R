#' A logistic regression model for the ECWM
#'
#' @description Fits a logistic regression model to data collected with the
#'   extended crosswise model, and predicts the probability of the sensitive
#'   characteristic conditional on the predictors.
#'
#' @param formula a "\code{\link{formula}}" in the form
#' \code{response ~ terms}, where \code{response} is the response
#' variable and \code{terms} describes the linear predictor. The
#' response variable should be numeric with the values 1 for
#' "One Yes answer" and 2 for "Two No or two Yes answers".
#' @param data data frame with individual records containing the response
#' and \code{prob_var} variables, and the predictors.
#' @param prob_var name of the variable with the individual probabilities
#' of answering "Yes" to the unrelated question.
#'
#' @details Maximization of the log-likelihood may fail when the population
#' prevalence estimate is close to the boundary, resulting in a
#' non-invertible Hessian matrix. In that case the standard errors are set
#' to \code{NA}s, and a warning is issued that the solution is not reliable.
#'
#' @return A list with the following objects:
#' \describe{
#' \item{coefs}{data frame with the parameter estimates.}
#' \item{gof}{fit statistics including the log-likelihood, AIC and a likelihood-ratio
#' test against the intercept-only model.}
#' \item{fitted}{vector with the estimated probabilities.}
#' \item{cormat}{correlation matrix of the parameter estimates.}
#' \item{D}{design matrix of the fitted model.}
#' }
#'
#' @examples
#' logreg(formula = Q1 ~ Time + Difficult, data = doping, prob_var = P1)
#'
#' # Separate prevalence estimates for each sub-sample, which yields the
#' # same LR statistic as the function 'prev'
#'
#' logreg(formula = Q1 ~ factor(P1), data = doping, prob_var = P1)
#'
#' @importFrom stats model.matrix optim pchisq pt complete.cases
#' @importFrom dplyr mutate case_when %>%
#'
#' @export

logreg <- function(formula, data, prob_var){



  y <- data[[formula[[2]]]]

  if(!is.numeric(y))         stop("response variable is not of class numeric.")

  if(!all(y %in% 1:2))       stop("response variable has values other than 1 and 2.")

  if(sum(is.na(y)) > 0)      stop("response variable has missing values.")

  y <- 2 - y

  p <- data[[substitute(prob_var)]]

  if(!is.numeric(p))         stop("prob_var is not of class numeric.")

  if(sum(is.na(p)) > 0)      stop("prob_var variable has missing values.")

  if(length(unique(p)) > 2)  warning("prob_var variable has more than two unique values")

  p2 <- mutate(data.frame(y, p),
               prc = case_when(p == min(p) & y == 0 ~ min(p),
                               p == min(p) & y == 1 ~ max(p),
                               p == max(p) & y == 0 ~ max(p),
                               p == max(p) & y == 1 ~ min(p)),
               .keep = "unused") %>%
    as.matrix()

  D     <- model.matrix(formula, data)

  reg0  <- optim(0, logl_prev, p2 = p2, method = "Brent", lower = 0, upper = 1)

  reg   <- tryCatch(optim(rep(0, ncol(D)), logl_reg, p2 = p2, D = D, hessian = TRUE),
                    error = function(e){0})

  if(!is.list(reg)) stop("Maximization failed, the fitted model may be too complex for the data")

  se    <- tryCatch(sqrt(diag(solve(reg$hessian))), error = function(e) NA)

  coefs <- data.frame(est       = reg$par,
                      se        = se,
                      t.value   = reg$par / se,
                      p.value   = 2 * pt(-abs(reg$par / se), df = nrow(p2) - ncol(D)),
                      row.names = colnames(D))

  LR  <-  2 * (reg0$value - reg$value)

  gof <-  data.frame(loglike = -reg$value %>%  round(1),
                     AIC     = 2 * (reg$value + ncol(D)) %>%  round(1),
                     LR      = LR %>%  round(2),
                     df      = ncol(D) - 1,
                     p       = pchisq(LR, ncol(D) - 1, lower.tail = FALSE) %>%  round(4),
                     row.names = "")

  if(anyNA(se)){

    warning("Unreliable solution due a non-invertible Hessian")

    cormat <- matrix(nrow = ncol(D), ncol = ncol(D))

  } else {

    cormat <- solve(reg$hessian)/(se %*% t(se))

  }

  dimnames(cormat) <- list(colnames(D), colnames(D))


  cat("\n")
  cat("Parameter estimates")
  cat("\n")
  print(coefs %>% round(3))
  cat("\n")
  cat("\n")
  cat("Fit statistics")
  cat("\n")
  print(gof)
  cat("\n")
  if(reg$convergence != 0) warning("The model did not successfully converge")

  invisible(list(coefs  = coefs,
                 gof    = gof,
                 fitted = c(logistic(D %*% reg$par)),
                 cormat = round(cormat, 3),
                 D      = D
  )
  )
}


#' Prevalence estimation for ECWM data
#'
#' @description Point and interval estimates of the
#' Method of Moments Estimator (MME) and the
#' Maximum Likelihood Estimator (MLE) for data collected
#' with the extended crosswise model.
#'
#' @param data data frame with individual records containing the
#' numeric \code{response} the \code{prob_var} variables.
#' @param response name of the response variable,
#' with value 1 denoting "One Yes answer" and value
#' 2 denoting "Two No or two Yes answers".
#' @param prob_var name of the variable with the probabilities
#' of answering "Yes" to the unrelated question.
#'
#' @details When the MME is outside the parameter space, a
#' non-parametric bootstrap percentile confidence interval is
#' reported for the MLE.
#'
#'
#' @return A list with the elements:
#' \describe{
#' \item{est}{data frame with point and interval estimates for the MME and MLE.}
#' \item{LR}{likelihood-ratio test against the model with separate
#' prevalence estimates for each sub-sample.}
#' \item{freqs}{data frame with the observed and estimated response frequencies.}
#' }
#'
#' @examples
#' prev(data = doping, response = Q1, prob_var = P1)
#'
#' @importFrom stats model.matrix optim pchisq pt xtabs quantile
#' @importFrom dplyr mutate case_when slice_sample %>%
#' @importFrom MASS ginv
#'
#' @export


prev <- function(data, response, prob_var){


  y <- data[[substitute(response)]]

  if(!is.numeric(y))         stop("response variable is not of class numeric.")

  if(!all(y %in% 1:2))       stop("response variable has values other than 1 and 2.")

  if(sum(is.na(y)) > 0)      stop("response variable has missing values.")

  y <- 2 - y

  p <- data[[substitute(prob_var)]]

  if(!is.numeric(p))         stop("prob_var is not of class numeric.")

  if(sum(is.na(p)) > 0)      stop("prob_var has missing values.")

  if(length(unique(p)) != 2) stop("prob_var should have two unique values")

  ### MME

  d <- data.frame(xtabs(~ y + p))

  p1 <- min(p)

  pistar <- d$Freq/sum(d$Freq)

  P <- matrix(c(p1,    1 - p1,
                1 - p1,    p1,
                1 - p1,    p1,
                p1,    1 - p1), 4, 2)

  mme  <- 2 * ginv(P) %*% pistar

  semm <- sqrt((4 * (ginv(P) %*% (diag(c(pistar)) - pistar %*% t(pistar)) %*% t(ginv(P))) / length(y))[1])



  ### MLE

  p2 <- mutate(data.frame(y, p),
               prc = case_when(p == min(p) & y == 0 ~ min(p),
                               p == min(p) & y == 1 ~ max(p),
                               p == max(p) & y == 0 ~ max(p),
                               p == max(p) & y == 1 ~ min(p)),
               .keep = "unused")


  sat <- optim(c(0, 0), logl_reg, p2 = p2, D = model.matrix(~ p))

  mle <- optim(0.5, logl_prev, p2 = p2, method = "Brent", lower = 0, upper = 1, hessian = TRUE)

  seml <- sqrt(solve(mle$hessian))

  bootci <- rep(0, 2)

  if(mme[1] <= 0 || mme[1] >= 1){

    bootci <- bootstrap(y = y, p = p)

  }

  est  <- data.frame(est       = c(mme[1], mle$par),
                     se        = c(semm, ifelse(mme[1] <= 0 || mme[1] >= 1, NA, seml)),
                     min95     = c(mme[1] - 1.96 * semm, ifelse(mme[1] <= 0 || mme[1] >= 1, bootci[1],
                                                                ifelse(mle$par - 1.96 * seml < 0, 0, mle$par - 1.96 * seml))),
                     max95     = c(mme[1] + 1.96 * semm, ifelse(mme[1] <= 0 || mme[1] >= 1, bootci[2],
                                                                ifelse(mle$par + 1.96 * seml > 1, 1, mle$par + 1.96 * seml))),
                     row.names = c("mme", "mle"))

  LR <- data.frame(G  = -2 * (sat$value - mle$value),
                   df = 1,
                   p  = pchisq(-2 * (sat$value - mle$value), df = 1, lower.tail = FALSE),
                   row.names = "")
  colnames(LR)[1] <- "X-squared"

  est_mme = c(rep(sum(d$Freq[1:2]), 2), rep(sum(d$Freq[3:4]), 2)) * P %*% mme
  est_mle = c(rep(sum(d$Freq[1:2]), 2), rep(sum(d$Freq[3:4]), 2)) * P %*% c(mle$par, 1 - mle$par)


  freqs <- data.frame(d, est_mme, est_mle)

  freqs <- mutate(freqs, y = 3 - as.numeric(y))

  colnames(freqs)[1:3] = c(substitute(response), substitute(prob_var), "obs")

  cat("prevalence estimates \n \n")
  print(est %>% round(5))
  cat("\n")
  cat("LR test \n \n")
  print(LR %>% round(4))
  cat("\n")
  cat("observed and estimated frequencies \n \n")
  print(freqs %>% mutate(est_mme = round(est_mme, 2),
                         est_mle = round(est_mle, 2)), row.names = FALSE)

  invisible(list(est   = est,
                 LR    = LR,
                 freqs = freqs))

}
