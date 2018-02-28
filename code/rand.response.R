# Perturb network
make.private.rr <- function (formula, dp.epsilon) {
  # get original network
  y.orig <- ergm.getnetwork(formula)
  
  p <- 1 - exp(dp.epsilon)/(1+exp(dp.epsilon))
  y <- rand.response.nw(y.orig, p)
  formula <- nonsimp.update.formula(formula, y ~ ., from.new=TRUE)
  
  noise <- list("dp.epsilon" = dp.epsilon, "p" = p, "method"="rr")
  return(list("formula" = formula, "noise" = noise))
}

rand.response.nw <- function(nw, p) {
  # flip each dyad with probabilitiy p_ij
  n <- length(nw$val)
  nw.mat <- as.matrix(nw)
  to.flip <- forceSymmetric(matrix(rbinom(n*n, 1, p), n, n))
  nw.mat <- as.matrix(abs(nw.mat - to.flip))
  diag(nw.mat) <- 0
  nw.new <- network(nw.mat, directed = FALSE)
  return(nw.new)
} 
