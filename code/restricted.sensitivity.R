#####################################################################
##                        Privatize a Network                      ##
#####################################################################

# Project the network to H_k and add noise to suff stats
make.private.restr <- function (formula, dp.epsilon, dp.k, privacy.type="edge", dp.delta=NULL, proj=TRUE) {
  # project network to H_k
  # get original network
  y.orig <- ergm.getnetwork(formula)
  
  if (proj) {
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
  } else {
    y <- y.orig
  }
 
  model <- ergm.getmodel(formula, y)
  # draw laplace noise
  noise <- draw.lap.noise.restricted(model$terms, dp.epsilon, dp.k, privacy.type, dp.delta, d.hat)
  noise$dp.epsilon <- dp.epsilon
  if (privacy.type == "node") noise$dp.delta <- dp.delta
  noise$dp.k <- dp.k
  noise$method <- "restr"
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
  noise.draw <- rlaplace(n = length(terms), scale = noise.level)
  return(list("level" = noise.level, "draw" = noise.draw))
}

########################################################################
##                            Projections to H_k                      ##
########################################################################

projection.edge <- function(nw, dp.k) {
  edge.list <- as.edgelist(nw)
  edge.ordering <- sample(1L:length(edge.list[,1]))
  nw.mat <- as.matrix(nw)
  
  # fix ordering of edges and iterate over it
  for (i in edge.ordering) {
    u <- edge.list[i, 1]
    v <- edge.list[i, 2]
    u.deg <- sum(nw.mat[u,])
    v.deg <- sum(nw.mat[v,])
    # if too many edges
    if (max(u.deg, v.deg) > dp.k) {
      nw.mat[u,v] <- 0
      nw.mat[v,u] <- 0
    }
    if(i%%1000==0) {
      nw.mat.sparse <- as.matrix.csr(nw.mat)
      if (max(diag((nw.mat.sparse %*% nw.mat.sparse))) <= dp.k) return(network(nw.mat, directed=FALSE))
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

  # set up matrix of all edges incident to vertex u/v (load from file if possible)
  fname <- sprintf("obj/lp_setup/const.%d", n)
  if (!file.exists(fname)) {
    adj.mat <- Matrix(0, nrow=n, ncol=choose(n,2), sparse=TRUE)
    x <- 1
    for (v in 2:n) {
      for (u in 1:(v-1)) {
        adj.mat[c(u,v), x] = 1
        x <- x + 1
      }
    }
    # save matrix to file
    save(adj.mat, file = fname)
  } else {
    load(fname)
  }
  
  
  # projection constraints
  const.mat.proj1 <- cbind(Matrix(0, nrow=choose(n,2), ncol=n, sparse=TRUE), .sparseDiagonal(choose(n,2)))
  const.dir.proj1 <- rep("<=", choose(n,2))
  const.rhs.proj1 <- a
  
  const.mat.proj2 <- cbind(t(adj.mat), .sparseDiagonal(choose(n,2)))
  const.dir.proj2 <- rep(">=", choose(n,2))
  const.rhs.proj2 <- a

  const.mat.rs <- cbind(matrix(0, nrow=n, ncol=n), adj.mat)
  const.dir.rs <- rep("<=", n)
  const.rhs.rs <- rep(dp.k, n)

  # setup of constraints
  const.mat <- rbind(const.mat.proj1, const.mat.proj2, const.mat.rs)
  const.dir <- c(const.dir.proj1, const.dir.proj2, const.dir.rs)
  const.rhs <- c(const.rhs.proj1, const.rhs.proj2, const.rhs.rs)
  
  toc()
  
  tic("LP solution")
   # all vars are non-negative + real by default 
  LP <- Rglpk_solve_LP(objective.in, const.mat, const.dir, const.rhs, max=FALSE)
  toc()
  
  d.hat <- 4*LP$optimum
  to.remove <- LP$solution[1:n] >= 0.25
  nw.mat[to.remove,] <- 0
  nw.mat[,to.remove] <- 0
  y <- network(nw.mat, directed=FALSE)
  
  return(list("y" = y, "d.hat" = d.hat, "lp" = LP))
}


