logistic <- function(x) exp(x) / ( 1 + exp(x))


# log-likelihood for prevalence estimation (on probability scale)

logl_prev <- function(b, p2){

  lik <- p2 * b + (1 - p2) * (1 - b)

  -sum(log(lik))
}


# log-likelihood for regression model (on logistic scale)

logl_reg <- function(b, p2, D){

  pr   <- logistic(D %*% b)

  lik <- p2 * pr + (1 - p2) * (1 - pr)

  -sum(log(lik))
}


bootstrap <- function(y, p){

  set.seed(1)

  d1 <- data.frame(y, p)[p == min(p), ]
  d2 <- data.frame(y, p)[p == max(p), ]

  piboot <- NULL

  for(i in 1:500){

    bootsample <- rbind(slice_sample(d1, n = nrow(d1), replace = TRUE),
                        slice_sample(d2, n = nrow(d2), replace = TRUE))

    p2 <- mutate(bootsample,
                 prc = case_when(p == min(p) & y == 0 ~ min(p),
                                 p == min(p) & y == 1 ~ max(p),
                                 p == max(p) & y == 0 ~ max(p),
                                 p == max(p) & y == 1 ~ min(p)),
                 .keep = "unused")

    piboot <- c(piboot,
                    optim(0.5, logl_prev, p2 = p2, method = "Brent", lower = 0, upper = 1)$par)
  }
  quantile(piboot, probs=c(.025, .975))
}
