source("./restricted.sensitivity.R", chdir=TRUE)

################################################################################
##                           Bayesian Inference with Noise                    ##
################################################################################


# apply Bayesian inference algorithm to training
dp.bergm <- function (formula, 
                   dp.epsilon,
                   dp.k,
                   privacy.type='edge',
                   dp.delta = NULL,
                   burn.in=100,
                   main.iters=1000,
                   aux.iters=1000, 
                   m.prior = NULL, 
                   sigma.prior = NULL, 
                   nchains = 1, 
                   gamma = 0.5, 
                   sigma.epsilon = NULL,
                   ...){ 	
  
  # project network to H_k
  # get original network
  y.orig <- ergm.getnetwork(formula)
  # project network to H_k
  if(privacy.type == 'edge') y <- projection.edge(y.orig, dp.k)
  # update formula with projection
  formula <- nonsimp.update.formula(formula, y ~ ., from.new=TRUE)

  model <- ergm.getmodel(formula, y)
  Clist <- ergm.Cprepare(y, model)
  
  # get statistics
  stats0 <- summary(formula)
  print("Stats: ")
  print(stats0)

  # draw noise for edge-level privacy under k-restricted sensitivity
  noise <- draw.lap.noise(model$terms, dp.epsilon, dp.k, privacy.type)
  print("Noise Level: ")
  print(noise)

  # set up MCMC control
  control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                      MCMC.interval = 0)
  control$MCMC.samplesize <- 1
  control$parallel <- 4

  MHproposal <- MHproposal.ergm(object= model, 
                                constraints = ~., arguments = control$MCMC.prop.args, 
                                nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                reference = ~Bernoulli, response = NULL)     
  
  snooker <- 0

  # initialize priors
  if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
  if (is.null(sigma.prior)) sigma.prior <- diag(100, Clist$nstats)
  if (is.null(nchains)) nchains <- 2 * Clist$nstats
  if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
  if (Clist$nstats == 1) {
    nchains <- 1
    sigma.epsilon <- diag(gamma, Clist$nstats)
  }
  Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
  # initialize thetas uniformly at random
  theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
  acc.counts <- rep(0L, nchains)
  theta1 <- rep(NA, Clist$nstats)
  tot.iters <- burn.in + main.iters
  
  # draw initial x as very far off from correct stats
  delta.y.old <- matrix(runif(Clist$nstats * nchains, min = -10.*noise$level, max = -5.*noise$level), Clist$nstats, nchains)

  print("Beginning MCMC...")
  # iterations of exchange algorithm
  for (k in 1L:tot.iters) {
    # number of chains to run for population MCMC
    for (h in 1L:nchains) {
      # set up 
      if (Clist$nstats > 1 && nchains > 1) {
        snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
      }
      
      # new proposal of theta
      theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
      
      # gets probabilities for new theta and old theta under prior on params
      pr <- dmvnorm(rbind(theta1, theta[, h]), 
                    mean = m.prior, 
                    sigma = sigma.prior,
                    log=TRUE)
      
      # sample a new network x
      delta <- ergm.mcmcslave(Clist, 
                              MHproposal, 
                              eta0 = theta1, 
                              control, 
                              verbose = FALSE)$s
        
      
      #### ADJUST delta to reflect correct delta ### 
      # delta.y = (suff stats new x) - (suff stats y noisy) 
      #         = (suff stats new x - suff stats y) - noise
      delta.y <- as.vector(delta - noise$draw)

      # delta.x     = (suff stats new x) - (suff stats old)
      #             = ((suff stats new x) - (suff stats y noisy)) - ((suff stats old x) - (suff stats y noisy))
      delta.x <- delta.y - delta.y.old[, h]
      
      # MH acceptance ratio 
      beta <- (theta[, h] - theta1) %*% delta.x
      # MH acceptance prob for theta
      beta.theta <- beta + pr[1] - pr[2] 
      # MH acceptance prob for x (#### BUG HERE ??? ###)
      beta.x <- beta + ((1./noise$level) %*% (abs(delta.y.old[, h]) - abs(delta.y)))

      # decide whether to accept new theta
      if (beta.theta >= log(runif(1))) {
        theta[, h] <- theta1
        if (k > burn.in) {
          acc.counts[h] <- acc.counts[h] + 1
          # make proposed theta steps smaller after burn in
          # sigma.epsilon <- diag(c(0.0025, 0.0025,0.0025, 0.0025), Clist$nstats)
        }
      }

      # decide whether to accept new x
      if (beta.x >= log(runif(1))) {
        # update deltas to reflect new x being used
        delta.y.old[, h] <- delta.y
      }
      
    
    if (k > burn.in) Theta[k - burn.in, , ] <- theta
    
    }
    if (k %% 1000 == 0) {
      print(sprintf("MCMC iter %d", k))
      print(delta.y.old)
      print(theta)
    } 
  }
  if (nchains == 1) Theta <- as.matrix(Theta[, , 1])
  
  out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
             formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
             dim = Clist$nstats, nchains = nchains, stats = stats0, 
             Theta = Theta, nchains = nchains, AR = acc.counts / main.iters, 
             m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters, noise = noise)
  out
  
}





################################################################################
##                                      Standard Inference                    ##
################################################################################
standard.bergm <- function (formula, 
                   dp.epsilon,
                   dp.k,
                   privacy.type='edge',
                   dp.delta = NULL, 
                   burn.in=100,
                   main.iters=1000,
                   aux.iters=1000, 
                   m.prior = NULL, 
                   sigma.prior = NULL, 
                   nchains = 1, 
                   gamma = 0.5, 
                   sigma.epsilon = NULL,
                   ...){ 	
    
  # project network to H_k
  # get original network
  y.orig <- ergm.getnetwork(formula)
  # project network to H_k
  if(privacy.type == 'edge') y <- projection.edge(y.orig, dp.k)
  # update formula with projection
  formula <- nonsimp.update.formula(formula, y ~ ., from.new=TRUE)

  model <- ergm.getmodel(formula, y)
  Clist <- ergm.Cprepare(y, model)

  noise <- draw.lap.noise(model$terms, dp.epsilon, dp.k, privacy.type)
        
    stats0 <- summary(formula)
    control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                        MCMC.interval = 0)
    control$MCMC.samplesize <- 1
    
    # constrain samples to have bounded degree
    MHproposal <- MHproposal.ergm(object= model, 
                                  constraints = ~bd(maxout=dp.k), arguments = control$MCMC.prop.args, 
                                  nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                  reference = ~Bernoulli, response = NULL)     
    
    snooker <- 0
    # initialize priors
    if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
    if (is.null(sigma.prior)) sigma.prior <- diag(100, Clist$nstats)
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
        if (k%%100 == 1) {
          print(sprintf("Completed %d iterations. Theta = ", k, theta))
          
        } 
        for (h in 1L:nchains) {
            if (Clist$nstats > 1 && nchains > 1) {
                snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
            }

            theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
            
            # gets probs for old theta and new theta both
            pr <- dmvnorm(rbind(theta1, theta[, h]), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)

            ### CAN INCORPORATE NOISE TO DELTA BY SUBTRACTING OFF NOISE ARGUMENT ###            
            delta <- ergm.mcmcslave(Clist, 
                                    MHproposal, 
                                    eta0 = theta1, 
                                    control, 
                                    verbose = FALSE)$s
            
            #delta <- as.vector(delta - noise$draw)

            # MH Acceptance probability <- CAN CHANGE THIS ACCEPTANCE PROBABILITY                       
            beta <- (theta[, h] - theta1) %*% delta + pr[1] - pr[2]
            
            if (beta >= log(runif(1))) {
                
                theta[, h] <- theta1
                if (k > burn.in) {
                  acc.counts[h] <- acc.counts[h] + 1
                }
            
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