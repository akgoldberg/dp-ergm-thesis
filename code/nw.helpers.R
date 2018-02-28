################################################################################
##                                  Helpers                                   ##
################################################################################

# compute the max degree node of a network
network.maxdegree <- function(network) {
  nw.mat <- as.matrix.csr(as.matrix.network(network))
  nw.degs <- diag((nw.mat %*% nw.mat))
  return(max(nw.degs))
}

# compute max shared partner of network
network.maxsharedpartners <- function(network) {
    nw.mat <- as.matrix.csr(as.matrix.network(network))
    nw.prod <- as.matrix(nw.mat %*% nw.mat)
    diag(nw.prod) <- 0
    return(max(nw.prod))
}

network.numedges <- function(network) {
  return(length(network$mel))
}

network.altkstar <- function(network, lambda=0.5) {
  return(summary(network ~ altkstar(lambda, fixed=TRUE))[[1]])
}

network.altktri <- function(network, lambda=0.5) {
    return(summary(network ~ gwesp(lambda, fixed=TRUE))[[1]])
}

network.altktwopath <- function(network, lambda=0.5) {
    return(summary(network ~ gwdsp(lambda, fixed=TRUE))[[1]])
}