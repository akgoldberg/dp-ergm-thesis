################################################################################
##                          Modified Private Inference                        ##
################################################################################

# use modified bergm procedure
bergm.modified.private <- function (formula, 
                        noise,
                        burn.in=100,
                        main.iters=1000,
                        aux.iters=1000, 
                        m.prior = NULL, 
                        sigma.prior = NULL, 
                        nchains = 1, 
                        gamma = 0.5, 
                        sigma.epsilon = NULL,
                        print.out = 1000,
                        ...){ 	
  
  # setup inputs to ergm
  y <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, y)
  Clist <- ergm.Cprepare(y, model)

  if (noise$method == "rr") y.mat <- as.matrix(y)
  
  # add noise to sufficient statistics
  stats0 <- summary(formula)
  if (noise$method == "smooth" || noise$method == "restr") {
    stats0 <- stats0 + noise$draw
  } 

  print(sprintf("Noisy Stats: %s", sstats.to.str(matrix(stats0))))
  control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                      MCMC.interval = 0)
  control$MCMC.samplesize <- 1
  
  MHproposal <- MHproposal.ergm(object= model, 
                                constraints = ~., arguments = control$MCMC.prop.args, 
                                nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                reference = ~Bernoulli, response = NULL)     
  
  snooker <- 0
  # initialize priors
  if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
  if (is.null(sigma.prior)) sigma.prior <- diag(5, Clist$nstats)
  if (is.null(nchains)) nchains <- 2 * Clist$nstats
  if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
  if (Clist$nstats == 1) {
    nchains <- 1
    sigma.epsilon <- diag(gamma, Clist$nstats)
  }
  Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
  theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
  acc.counts <- rep(0L, nchains)
  theta1 <- rep(NA, Clist$nstats)
  tot.iters <- burn.in + main.iters
  
  # initialize best guess of model to the noisy suff stats
  suffstats.old <-  matrix(rep(stats0, nchains), Clist$nstats, nchains)
  # initialize best guess of how far off probability is 
  if (noise$method == "rr") edge.diff.old <- rep(0, nchains)
  
  for (k in 1L:tot.iters) {
    if (k%%print.out == 0) {
      print(sprintf("Completed %d iterations. Thetas = %s", k, theta.to.str(theta)))
      if (noise$method != "rr") {
        print(sprintf("Best guess of suff stats is now: %s", sstats.to.str(suffstats.old)))
      }
    }  
    for (h in 1L:nchains) {
      if (Clist$nstats > 1 && nchains > 1) {
        snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
      }
      
      theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
      
      pr <- dmvnorm(rbind(theta1, theta[, h]), 
                    mean = m.prior, 
                    sigma = sigma.prior,
                    log=TRUE)
      
      # delta = f(x) - f(y)
      z <-ergm.mcmcslave(Clist, 
                              MHproposal, 
                              eta0 = theta1, 
                              control, 
                              verbose = FALSE)

      delta <- z$s
      if (noise$method != "rr") delta <- z$s - noise$draw
      
      suffstats <- stats0 + delta

      # change in suff stats from x to x' -- f(x') - f(x)
      delta.x <- suffstats - suffstats.old[, h]
      # change in unnormalized likelihood from x to x'
      beta.lik <- as.vector(theta[, h] - theta1) %*% as.vector(delta.x)
      
      # theta acceptance probability
      beta.theta <- beta.lik + pr[1] - pr[2]
      
      if (beta.theta >= log(runif(1))) {
        theta[, h] <- theta1
        beta.lik <- 0 # set to 0 for appropriate gibbs update
        if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1
      }

      # x acceptance probability
      if (noise$method == "rr") {
        nw.sampled <- newnw.extract(y, z)
        edge.diff <- sum(abs(as.matrix(nw.sampled) - y.mat))/2.
        beta.x <- beta.lik + (edge.diff - edge.diff.old[h])*log(noise$p) 
      } 
      else {
        beta.x <- beta.lik + (1./noise$level) %*% (abs(suffstats.old[, h] - as.vector(stats0)) - as.vector(abs(suffstats - stats0)))
      }
      
      # replace current best guess of x
      if (beta.x >=  log(runif(1))) {
            suffstats.old[, h] <- suffstats
            if (noise$method == "rr") edge.diff.old[h] <- edge.diff
      }

    }
    if (k > burn.in) Theta[k - burn.in, , ] <- theta
  }
  if (nchains == 1) Theta <- as.matrix(Theta[, , 1])
  
  
  out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
             formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
             dim = Clist$nstats, nchains = nchains, stats = stats0, 
             Theta = Theta, nchains = nchains, AR = acc.counts / main.iters, 
             m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters, noise = noise)
  out
}


#####################################################################
##                        Non-Private Inference                    ##
#####################################################################
# basic inference using bergm (modified to include printing)
bergm.orig <- function (formula, 
                   burn.in=100,
                   main.iters=1000,
                   aux.iters=1000, 
                   m.prior = NULL, 
                   sigma.prior = NULL, 
                   nchains = 1, 
                   gamma = 0.5, 
                   sigma.epsilon = NULL,
                   print.out = 1000,
                   ...){ 	
    
    # setup inputs to ergm
    y <- ergm.getnetwork(formula)
    model <- ergm.getmodel(formula, y)
    Clist <- ergm.Cprepare(y, model)
        
    stats0 <- summary(formula)
    control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                        MCMC.interval = 0)
    control$MCMC.samplesize <- 1
    
    MHproposal <- MHproposal.ergm(object= model, 
                                  constraints = ~., arguments = control$MCMC.prop.args, 
                                  nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                  reference = ~Bernoulli, response = NULL)     
    
    snooker <- 0
    # initialize priors
    if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
    if (is.null(sigma.prior)) sigma.prior <- diag(5, Clist$nstats)
    if (is.null(nchains)) nchains <- 2 * Clist$nstats
    if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
    if (Clist$nstats == 1) {
    	 nchains <- 1
        sigma.epsilon <- diag(gamma, Clist$nstats)
    }
    Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
    theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
    acc.counts <- rep(0L, nchains)
    theta1 <- rep(NA, Clist$nstats)
    tot.iters <- burn.in + main.iters

    for (k in 1L:tot.iters) {
        if (k%%print.out == 0)  print(sprintf("Completed %d iterations. Thetas = %s", k, theta.to.str(theta)))
        for (h in 1L:nchains) {
            if (Clist$nstats > 1 && nchains > 1) {
                snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
            }
            
            theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
            
            pr <- dmvnorm(rbind(theta1, theta[, h]), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)

            delta <- ergm.mcmcslave(Clist, 
                                    MHproposal, 
                                    eta0 = theta1, 
                                    control, 
                                    verbose = FALSE)$s

            beta <- as.vector(theta[, h] - theta1) %*% as.vector(delta) + pr[1] - pr[2]
            
            if (beta >= log(runif(1))) {
                
                theta[, h] <- theta1
                if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1
            
            }
        
        }
        if (k > burn.in) Theta[k - burn.in, , ] <- theta
    }
    if (nchains == 1) Theta <- as.matrix(Theta[, , 1])


    out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
        formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
        dim = Clist$nstats, nchains = nchains, stats = stats0, 
        Theta = Theta, nchains = nchains, AR = acc.counts / main.iters, 
        m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters)
    out

}


#####################################################################
##                        Naive Private Inference                  ##
#####################################################################

# use private statistics directly in original bergm
bergm.orig.private <- function (formula, 
                        noise,
                        burn.in=100,
                        main.iters=1000,
                        aux.iters=1000, 
                        m.prior = NULL, 
                        sigma.prior = NULL, 
                        nchains = 1, 
                        gamma = 0.5, 
                        sigma.epsilon = NULL,
                        print.out = 1000,
                        ...){ 	
  
  # setup inputs to ergm
  y <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, y)
  Clist <- ergm.Cprepare(y, model)
  
  stats0 <- summary(formula)
  control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                      MCMC.interval = 0)
  control$MCMC.samplesize <- 1
  
  MHproposal <- MHproposal.ergm(object= model, 
                                constraints = ~., arguments = control$MCMC.prop.args, 
                                nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                reference = ~Bernoulli, response = NULL)     
  
  snooker <- 0
  # initialize priors
  if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
  if (is.null(sigma.prior)) sigma.prior <- diag(5, Clist$nstats)
  if (is.null(nchains)) nchains <- 2 * Clist$nstats
  if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
  if (Clist$nstats == 1) {
    nchains <- 1
    sigma.epsilon <- diag(gamma, Clist$nstats)
  }
  Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
  theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
  acc.counts <- rep(0L, nchains)
  theta1 <- rep(NA, Clist$nstats)
  tot.iters <- burn.in + main.iters
  
  for (k in 1L:tot.iters) {
    if (k%%print.out == 0)  print(sprintf("Completed %d iterations. Thetas = %s", k, theta.to.str(theta)))
    for (h in 1L:nchains) {
      if (Clist$nstats > 1 && nchains > 1) {
        snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
      }
      
      theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
      
      pr <- dmvnorm(rbind(theta1, theta[, h]), 
                    mean = m.prior, 
                    sigma = sigma.prior,
                    log=TRUE)
      
      delta <- ergm.mcmcslave(Clist, 
                              MHproposal, 
                              eta0 = theta1, 
                              control, 
                              verbose = FALSE)$s
      
      # add noise to delta
      if (noise$method == "smooth" || noise$method == "restr") delta <- delta - noise$draw
      
      beta <- as.vector(theta[, h] - theta1) %*% as.vector(delta) + pr[1] - pr[2]
      
      if (beta >= log(runif(1))) {
        
        theta[, h] <- theta1
        if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1
        
      }
      
    }
    if (k > burn.in) Theta[k - burn.in, , ] <- theta
  }
  if (nchains == 1) Theta <- as.matrix(Theta[, , 1])
  
  
  out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
             formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
             dim = Clist$nstats, nchains = nchains, stats = stats0, 
             Theta = Theta, nchains = nchains, AR = acc.counts / main.iters, 
             m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters)
  out
  
}


#####################################################################
##                            Helpers                              ##
#####################################################################

# turn theta matrix into string for printing
theta.to.str <- function(theta) {
  sprintf("[%s]", paste(apply(round(theta,3), 2, paste, collapse=","), collapse="] ["))
}

sstats.to.str <- function(sstats) {
   sprintf("[%s]", paste(apply(round(sstats,3), 2, paste, collapse=","), collapse="] ["))
}