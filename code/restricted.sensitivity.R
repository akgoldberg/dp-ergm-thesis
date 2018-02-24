#####################################################################
##                        Privatize a Network                      ##
#####################################################################

# Project the network to H_k and add noise to suff stats
make.private <- function (formula, dp.epsilon, dp.k, privacy.type="edge", dp.delta=NULL) {
  # project network to H_k
  # get original network
  y.orig <- ergm.getnetwork(formula)
  
  d.hat <- NULL
  # project network to H_k
  if (privacy.type=="edge") y <- projection.edge(y.orig, dp.k)
  if (privacy.type=="node") {
    out <- projection.node(y.orig, dp.k)
    y <- out$y
    d.hat <- out$d.hat
  }
  
  # update formula with projection
  formula <- nonsimp.update.formula(formula, y ~ ., from.new=TRUE)
  model <- ergm.getmodel(formula, y)
  # draw laplace noise
  noise <- draw.lap.noise.restricted(model$terms, dp.epsilon, dp.k, privacy.type, dp.delta, d.hat)
  noise$dp.epsilon <- dp.epsilon
  if (privacy.type == "node") noise$dp.delta <- dp.delta
  noise$dp.k <- dp.k
  return(list("formula" = formula, "noise" = noise))
}

################################################################################
##                                Restricted Sensitivity                      ##
################################################################################

# Draw Laplace noise under restricted sensitivty
draw.lap.noise.restricted <- function(terms, dp.epsilonTot, dp.k, privacy.type, dp.deltaTot=NULL, d.hat=NULL) {
  
  # output vector of level of laplace noise to add
  noise.level <- rep(NA, length(terms))
  # split epsilon evenly between terms if it is given as scalar
  if (length(dp.epsilonTot == 1)) {
    dp.epsilon <- rep(dp.epsilonTot/length(terms), length(terms))
  } else {
    if (length(dp.epsilonTot) != length(terms)) {
      print("Epsilon must be either scalar or specified for all terms.")
      return()
    }
    dp.epsilon <- dp.epsilonTot
  }

  # split delta evenly between terms if it is given as scalar
  if (privacy.type == "node") {
      # split delta evenly between terms if it is given as scalar
      if (length(dp.deltaTot == 1)) {
        dp.delta <- rep(dp.deltaTot/length(terms), length(terms))
      } else {
        if (length(dp.deltaTot) != length(terms)) {
          print("Delta must be either scalar or specified for all terms.")
          return()
        }
        dp.delta <- dp.deltaTot
      }
  }


  # calculate amount of noise to add for each term
  for (t in 1L:length(terms)) {
    name <- terms[[t]]$name
    param <- tail(terms[[t]]$inputs, 1)
    # edge-level privacy
    if (privacy.type == 'edge') {
      if (name == 'edges') {
        noise.level[[t]] <- 1/dp.epsilon[t]
      }
      if (name == 'altkstar') {
        noise.level[[t]] <- (2*(1./param))/dp.epsilon[t]
      }
      if (name == 'gwesp') {
        noise.level[[t]] <- 3*(2*(dp.k-1) + 1./param)/dp.epsilon[t]
      }
      if (name == 'gwdsp') {
        noise.level[[t]] <- 3*(2*(dp.k-1))/dp.epsilon[t]
      }
    }

    # node-level privacy
    if (privacy.type == 'node') {
      # restricted sensitivity is computed over H_2k
      dp.k <- 2*dp.k

      # set up multiplicative factor for smooth sensitivity
      beta <- dp.epsilon[t]/(2*log(1/dp.delta))
      if ((beta/4.) <= 2./5.) {
        x <- beta/4
        g <- (2./x) * exp((5/2)*x - 1)
      } else {
        g <- 5
      }
      C <- g*exp(beta*d.hat/4.)
      if (name == 'edges') {
        noise.level[[t]] <- 2*C*dp.k/dp.epsilon[t]
      }
      if (name == 'altkstar') {
        noise.level[[t]] <- 2*C*(3*dp.k*(1./param))/dp.epsilon[t]
      }
      if (name == 'gwesp') {
        noise.level[[t]] <- 2*C*((dp.k**2) + ((1./param) - 1))/dp.epsilon[t]
      }
      if (name == 'gwdsp') {
        noise.level[[t]] <-  2*C*(dp.k**2)
      }
    }
  }
  # draw Laplace noise to use
  noise.draw <- rlaplace(n = length(terms), s = noise.level)
  return(list("level" = noise.level, "draw" = noise.draw))
}

########################################################################
##                            Projections to H_k                      ##
########################################################################

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

# projection into H_k for node-level privacy
projection.node <- function(nw, dp.k) {
  tic("LP setup")
  
  nw.mat <- as.matrix(nw)
  n <- dim(nw.mat)[1]
  num.vars <- n + choose(n,2)
  # list of edges
  a <- nw.mat[upper.tri(nw.mat, diag=FALSE)]
  # setup linear program
  # objective function
  objective.in <- c(rep(1, n), rep(0, choose(n,2)))

  # non-negativity constraints
  const.mat.nonneg <- diag(num.vars)
  const.dir.nonneg <- rep(">=", num.vars)
  const.rhs.nonneg <- rep(0,num.vars)

  # projection constraints
  const.mat.proj1 <- cbind(matrix(0, nrow=choose(n,2), ncol=n), diag(choose(n,2)))
  const.dir.proj1 <- rep("<=", choose(n,2))
  const.rhs.proj1 <- a
  
  # set up matrix of u's and v's for w_uv + x_u + x_v >= a_uv
  x.uv.mat <- matrix(0, nrow=choose(n,2), ncol=n)
  row <- 1
  for (v in 2:n) {
    for (u in 1:(v-1)) {
       x.uv.mat[row, u] <- 1
       x.uv.mat[row, v] <- 1
       row <- row + 1
    }
  }
  const.mat.proj2 <- cbind(x.uv.mat, diag(choose(n,2)))
  const.dir.proj2 <- rep(">=", choose(n,2))
  const.rhs.proj2 <- a

  # restricted degree constraint
  # set up matrix of all edges incident to vertex u
   u.mat <- matrix(0, nrow=n, ncol=choose(n,2))
   x <- 1
   for (v in 2:n) {
    for (u in 1:(v-1)) {
       u.mat[u, x] = 1
       u.mat[v, x] = 1
       x <- x + 1
    }
  }
  const.mat.rs <- cbind(matrix(0, nrow=n, ncol=n), u.mat)
  const.dir.rs <- rep("<=", n)
  const.rhs.rs <- rep(dp.k, n)

  const.mat <- rbind(const.mat.nonneg, const.mat.proj1, const.mat.proj2, const.mat.rs)
  const.dir <- c(const.dir.nonneg, const.dir.proj1, const.dir.proj2, const.dir.rs)
  const.rhs <- c(const.rhs.nonneg, const.rhs.proj1, const.rhs.proj2, const.rhs.rs)
  
  toc()
  
  tic("LP solution")
  LP <- lp(direction = "min", objective.in, const.mat, const.dir, const.rhs)
  toc()
  
  d.hat <- 4*LP$objval
  to.remove <- LP$solution[1:n] >= 0.25
  nw.mat[to.remove,] <- 0
  nw.mat[,to.remove] <- 0
  y <- network(nw.mat, directed=FALSE)
  
  return(list("y" = y, "d.hat" = d.hat, "lp" = LP))
}

