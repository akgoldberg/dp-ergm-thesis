################################################################################
##                                Restricted Sensitivity                      ##
################################################################################

# draw Laplace noise under restricted sensitivty
draw.lap.noise.restricted <- function(terms, dp.epsilonTot, dp.k, privacy.type) {
  noise.level <- rep(NA, length(terms))
  # split epsilon evenly between terms
  dp.epsilon <- dp.epsilonTot/length(terms)
  # calculate amount of noise to add for each term
  for (t in 1L:length(terms)) {
    name <- terms[[t]]$name
    param <- tail(terms[[t]]$inputs, 1)
    # edge-level privacy
    if (privacy.type == 'edge') {
      if (name == 'edges') {
        noise.level[[t]] <- 1/dp.epsilon
      }
      if (name == 'altkstar') {
        noise.level[[t]] <- (2*(1./param))/dp.epsilon
      }
      if (name == 'gwesp') {
        noise.level[[t]] <- 3*(2*(dp.k-1) + 1./param)/dp.epsilon
      }
      if (name == 'gwdsp') {
        noise.level[[t]] <- 3*(2*(dp.k-1))/dp.epsilon
      }
    }

    # node-level privacy
    # TODO #
    if (privacy.type == 'node') {
      noise.level[[t]] <- 10.
    }
  }
  # draw Laplace noise to use
  noise.draw <- rlaplace(n = length(terms), s = noise.level)
  return(list("level" = noise.level, "draw" = noise.draw))
}

# projection into H_k for edge-level privacy
projection.edge <- function(nw, dp.k) {
  edge.list <- as.edgelist(nw)
  edge.ordering <- sample(1L:length(edge.list[,1]))
  # fix ordering of edges and iterate over it
  for (i in edge.ordering) {
    u <- edge.list[i, 1]
    v <- edge.list[i, 2]
    u.deg <- length(get.edgeIDs(nw, u))
    v.deg <- length(get.edgeIDs(nw, v))
    # if too many edges
    if (max(u.deg, v.deg) > dp.k) {
      # get edge id associated w/ edge between u and v
      e.id <- get.edgeIDs(nw, u, v)
      # delete that edge
      delete.edges(nw, e.id)
    }
  }
  return(nw)
}