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

network.triangles <- function(network) {
    return(summary(network ~ triangles)[[1]])
}

get.statname <- function(ergm.name) {
  if (ergm.name == 'edges') {
    return('Edges')
  }
  if (str_count(ergm.name, 'altkstar') == 1) {
    return('Alt k-Star')
  }
  if (str_count(ergm.name, 'gwesp') == 1) {
    return('Alt k-Triangle')
  }
  if (str_count(ergm.name, 'gwdsp') == 1) {
    return('Alt k-Two-Path')
  }
  if (str_count(ergm.name, 'nodefactor') == 1) {
    return(sprintf('P-%s', tail(str_split(ergm.name, "[.]")[[1]],1)))
  }
  if (str_count(ergm.name, 'nodematch') == 1) {
    return(sprintf('H-%s', tail(str_split(ergm.name, "[.]")[[1]],1)))
  }
}

copy.vertex.attrs <- function(nw.old, nw.new) {
    for (attr in list.vertex.attributes(nw.old)) {
    set.vertex.attribute(nw.new, attr, get.vertex.attribute(nw.old, attr)) 
  }
  delete.vertex.attribute(nw.new, 'vertex.names')
  return(nw.new)
}
