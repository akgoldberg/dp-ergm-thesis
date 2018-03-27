#####################################################################
##                        Privatize a Network                      ##
#####################################################################

# Project the network to H_k and add noise to suff stats
make.private.restr <- function (formula, dp.epsilon, dp.k, privacy.type="edge",
                                 dp.delta=NULL, proj=TRUE, proj.type=NULL, noise.type='lap', labels.priv=FALSE, attrs=NULL, 
                                 eps.labels=1.) {
  # project network to H_k
  # get original network
  y.orig <- ergm.getnetwork(formula)
  
  # project network to H_k
  if (proj) {
    proj.aux.info <- NULL
    
    # edge projection
    if (privacy.type=="edge") y <- projection.edge(y.orig, dp.k)
    # node projection
    if (privacy.type=="node" && proj.type =='LP') {
      out <- projection.node.LP(y.orig, dp.k)
      y <- out$y
      proj.aux.info <- out$d.hat
    }
    if (privacy.type=="node" && proj.type =='trunc') {
      out <- projection.node.trunc(y.orig, dp.k)
      y <- out$y
      proj.aux.info <- out$Cs
    }
  } else {
    y <- y.orig
  }
  
  # privatize labels
  if (labels.priv) {
    y <- private.labels(y, eps.labels, attrs=attrs)
  }
  
  # update formula with projection
  formula <- nonsimp.update.formula(formula, y ~ ., from.new=TRUE)
 
  model <- ergm.getmodel(formula, y)
  # draw laplace noise
  noise <- draw.lap.noise.restricted(model$terms, dp.epsilon, dp.k, privacy.type, 
                                     dp.delta, proj.aux.info=proj.aux.info, proj.type=proj.type, noise.type=noise.type, labels.priv=labels.priv)
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
draw.lap.noise.restricted <- function(terms, dp.epsilonTot,
                                      dp.k, privacy.type, labels.priv=FALSE,
                                      dp.deltaTot=NULL, proj.type=NULL, noise.type='lap', proj.aux.info=NULL) {
  
  # output vector of level of laplace noise to add
  noise.level <- rep(NA,  length(unlist(lapply(terms, function(x) {x$coef.names}))))
  # split epsilon evenly between terms if it is given as scalar
  if (length(dp.epsilonTot) == 1) {
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
      if (noise.type == 'lap') {
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
      } else {
        dp.delta <- rep(NULL, length(terms))
      }
      

      # SMOOTH SENSITIVITY SETUP
      # (assuming even split of noise)
      beta <- get.beta(noise.type, dp.epsilon[1], dp.delta[1])
      if (proj.type == "LP") {
        # restricted sensitivity is computed over H_2k
        dp.k <- 2*dp.k
        C <- smooth.sens.LP(proj.aux.info, beta)
      }
      if (proj.type == "trunc") {
        C <- smooth.sens.trunc(proj.aux.info, beta)
      }
      if (noise.type == "cauchy") {
        C <- C*sqrt(2)
      }
      if (noise.type == "lap") {
        C <- C*2
      }
  } else {
    C <- NULL
  }

  # calculate amount of noise to add for each term
  i <- 1
  for (t in 1L:length(terms)) {
    name <- terms[[t]]$name
    param <- tail(terms[[t]]$inputs, 1)
    num.terms <- length(terms[[t]]$coef.names)
    
    # edge-level privacy 
    if (privacy.type == 'edge') {
      if (name == 'edges') {
        noise.level[i] <- 1/dp.epsilon[t]
      }
      if (name == 'altkstar') {
        noise.level[i] <- (2*(1./param))/dp.epsilon[t]
      }
      if (name == 'gwesp') {
        noise.level[i] <- 3*(2*(dp.k-1) + 1./param)/dp.epsilon[t]
      }
      if (name == 'gwdsp') {
        noise.level[i] <- 3*(2*(dp.k-1))/dp.epsilon[t]
      }
      # label-dependant terms
      if (name == 'nodematch') {
        if (!labels.priv) {
          noise.level[i:(i+num.terms)] <- 1/dp.epsilon[t]
        } else {
          noise.level[i:(i+num.terms)] <- 3*dp.k/dp.epsilon[t]
        }
      } 
      if (name == 'nodefactor') {
        if (!labels.priv) {
          noise.level[i:(i+num.terms)] <- 2/dp.epsilon[t]
        } else {
          noise.level[i:(i+num.terms)] <- 3*1*dp.k/dp.epsilon[t]
        }
      }
    }

    # node-level privacy structural terms
    if (privacy.type == 'node') {

      if (name == 'edges') {
        noise.level[i] <- C*dp.k/dp.epsilon[t]
      }
      if (name == 'altkstar') {
        noise.level[i] <- C*(dp.k*(1./param))/dp.epsilon[t]
      }
      if (name == 'gwesp') {
        noise.level[i] <- C*((dp.k**2) + ((1./param) - 1))/dp.epsilon[t]
      }
      if (name == 'gwdsp') {
        noise.level[i] <-  C*(dp.k**2)/dp.epsilon[t]
      }
      # label-dependant terms
      if (name == 'nodematch') {
          noise.level[i:(i+num.terms)] <- C*dp.k/dp.epsilon[t]
      } 
      if (name == 'nodefactor') {
          noise.level[i:(i+num.terms)] <- C*dp.k/dp.epsilon[t]
      }
    }
    
    i <- i+num.terms
  }
  # draw Laplace noise to use
  if (noise.type == 'lap')  noise.draw <- rlaplace(n = length(noise.level), scale = noise.level)
  if (noise.type == 'cauchy') noise.draw <- rcauchy(n = length(noise.level), scale=noise.level)
  return(list("level" = noise.level, "draw" = noise.draw, "C" = C))
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
      if (max(diag((nw.mat.sparse %*% nw.mat.sparse))) <= dp.k) break
    } 
  }
  nw.out <- network(nw.mat, directed=FALSE)
  nw.out <- copy.vertex.attrs(nw, nw.out)

  return(nw.out)
}

projection.node.trunc <- function(nw, dp.k) {
   nw.mat <- as.matrix(nw)

   n <- dim(nw.mat)[1]
   to.del <- c()
   for (i in 1:n) {
     if (sum(nw.mat[i,]) > dp.k) to.del <- append(to.del, i)
   }
   # truncate vertices
   nw.mat[to.del,] <- 0
   nw.mat[,to.del] <- 0

    y <- network(nw.mat, directed=FALSE)
    y <- copy.vertex.attrs(nw, y)

    # find Ct for smooth sens
    nw.matorig <- as.matrix(nw)
    deg.dist <- data.table(table(apply(nw.matorig, 1, sum)))
    colnames(deg.dist) <- c("deg", "num")
    deg.dist$deg <- as.numeric(deg.dist$deg)
    C_t <- function(t) {1 + t + deg.dist[(deg >= (dp.k - t)) & (deg <= (dp.k + t + 1)), sum(num)]}
    Cs <- sapply(seq(0,1+round(2*dp.k)), function(t) C_t(t))
    return(list("y" = y, "Cs"=Cs))
}

get.beta <- function(noise.type, dp.epsilon, dp.delta=NULL) {
  if (noise.type == 'lap') {
    return(dp.epsilon/(2*log(1/dp.delta)))
  } 
  if (noise.type == 'cauchy') {
    return(dp.epsilon/sqrt(2))
  }
}

smooth.sens.trunc <- function(Cs, dp.beta) {
   vals <- unlist(sapply(seq(0,length(Cs)), function(t) exp(-dp.beta*t)*Cs[t]))
   beta.smooth <- max(vals)
   return(beta.smooth)
}

smooth.sens.LP <- function(d.hat, dp.beta) {
   if ((dp.beta/4.) <= 2./5.) {
        x <- dp.beta/4
        g <- (2./x) * exp((5/2)*x - 1)
      } else {
        g <- 5
      }
      C <- g*exp(dp.beta*d.hat/4.)
}

# projection into H_k for node-level privacy
projection.node.LP <- function(nw, dp.k) {

  #tic("LP setup")
  
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
  
  #toc()
  
  #tic("LP solution")
   # all vars are non-negative + real by default 
  LP <- Rglpk_solve_LP(objective.in, const.mat, const.dir, const.rhs, max=FALSE)
  #toc()
  
  d.hat <- 4*LP$optimum
  to.remove <- LP$solution[1:n] >= 0.25
  nw.mat[to.remove,] <- 0
  nw.mat[,to.remove] <- 0
  y <- network(nw.mat, directed=FALSE)
  y <- copy.vertex.attrs(nw, y)
  
  return(list("y" = y, "d.hat" = d.hat))
}

# privatize labels in a network
private.labels <- function(nw, eps, attrs=NULL) {
  if (is.null(attrs)) {
    label.df <- subset(rbindlist(lapply(nw$val, data.frame)), select=-na)
  } else {
    label.df <- subset(rbindlist(lapply(nw$val, data.frame)), select=attrs)
  }
  attrs <- colnames(label.df)
  n <- dim(label.df)[1]
  # generate histogram of labels
  label.df <- data.frame(table(label.df))
  colnames(label.df) <- c("attr", "Freq")
  
  # add noise
  noisyFreq <- sapply(label.df$Freq + rlaplace(n=dim(label.df)[1], scale=1./eps), max, 0)
  postproc <-round(noisyFreq)
  label.df <- cbind(label.df, noisyFreq, postproc)
  sort.order <- order(abs(label.df$noisyFreq - label.df$postproc), decreasing=TRUE)
  label.df <- label.df[sort.order, ]
  n.diff <- n-sum(label.df$postproc)
  # postprocess
  i <- 1
  num.fixed <- 0
  while (num.fixed < abs(n.diff)) {
    i.touse <- i
    if (i > dim(label.df)[1]) i.touse <- 1 + (i %% dim(label.df)[1])
    # do not get negative counts
    if (label.df$postproc[i.touse] + sign(n.diff) > 0) {
      label.df$postproc[i.touse] <- label.df$postproc[i.touse] + sign(n.diff)
      num.fixed <- num.fixed + 1
    }
    i <- i + 1
  }
  # make nw
  nw.out <- copy(nw)
  #delete.vertex.attribute(nw.out, list.vertex.attributes(nw.out))
  for (attr in attrs) {
    set.vertex.attribute(nw.out, attr, as.vector(rep(label.df$attr, label.df$postproc)))
  }
  return(nw.out)
}


