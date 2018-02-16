
################################################################################
##                                  Sample G(n,p)                             ## 
################################################################################
sample_bernoulli <- function(n, p, t, nsim=20) {
  MCMC_control <- control.simulate.formula(MCMC.burnin=1e6,
                                           MCMC.interval=1e4,
                                           MCMC.prop.weights="default",
                                           MCMC.prop.args=list(),
                                           MCMC.packagenames=c(),
                                           MCMC.runtime.traceplot=FALSE,
                                           network.output="network",
                                           
                                           parallel=0,
                                           parallel.type=NULL,
                                           parallel.version.check=TRUE)
  # generate samples
  # gwesp - same as alt-k-triangle, gwdsp - same as alt-k-two-path
  theta <- log(p/(1-p))
  sim <- simulate(network(n, directed=FALSE) 
                  ~ edges,
                  coef = theta,
                  sequential = TRUE,
                  nsim = nsim,
                  control = MCMC_control)
  # calculate degree of simulated graphs
  sim.deg <- unlist(lapply(sim, network.maxdegree))
  sim.sp <- unlist(lapply(sim, network.maxsharedpartners))
  #hist(sim.deg)
  return(sim)
}

################################################################################
##                                  Sample ERGM                               ##
################################################################################
sample_ergm <- function(n, theta, nsim=50) {
    MCMC_control <- control.simulate.formula(MCMC.burnin=4*choose(n,2),
                                         MCMC.interval=0.5*choose(n,2),
                                         MCMC.prop.weights="default",
                                         MCMC.prop.args=list(),
                                         MCMC.packagenames=c(),
                                         MCMC.runtime.traceplot=TRUE,
                                         network.output="network",
                                         
                                         parallel=4,
                                         parallel.type=NULL,
                                         parallel.version.check=TRUE)

  # generate samples
  # gwesp - same as alt-k-triangle, gwdsp - same as alt-k-two-path
  sim <- simulate(network(n, directed=FALSE) 
                  ~ edges + altkstar(0.5, fixed=TRUE) 
                  + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE),
           coef = theta,
           sequential = TRUE,
           nsim = nsim*MCMC_control$parallel,
           control = MCMC_control,
           verbose=TRUE)
  
  sample <- filter.null(sim)
  
  # compute degrees and shared partners of sampled networks
  deg <- dist.samples(sample, "deg")
  sp <- dist.samples(sample, "sp")
  
  return (list("sample" = sample, "deg" = deg, "sp" = sp, "n" = n, "theta" = theta))
}

################################################################################
##                                  Helpers                                   ##
################################################################################

# compute the max degree node of a network
network.maxdegree <- function(network) {
  dd <- degreedist(network, print=FALSE)
  max_deg <- names(tail(dd, 1))
  return(as.numeric(str_extract_all(max_deg, "[0-9]+")[[1]]))
}

# compute max shared partner of network
network.maxsharedpartners <- function(network) {
    nw.mat <- as.matrix.network(network)
    n.nodes <- length(nw.mat[1,])
    max.val <- 0
    # iterate over all dyads
    for (u in 1L:(n.nodes-1)) {
        for (v in (u+1):n.nodes) {
            # shared partners of u and v
            C <- nw.mat[u, ] %*% nw.mat[, v]
            max.val <- max(C, max.val)
        }
    }
    return(max.val)
}

dist.samples <- function(nws, type) {
    if (type == "deg") return(unlist(lapply(nws, network.maxdegree)))
    if (type == "sp") return(unlist(lapply(nws, network.maxsharedpartners)))
}

filter.null <- function(samp) {
  return(samp[!sapply(samp, is.null)]) 
}

get.levels <- function(noise) {
  sapply(noise, function (o) {o$level})
}

################################################################################
##             Compute Noise Added on Sample for Smooth vs. Restricted        ##
################################################################################
test.noise.sample <- function(nws.sample, dp.epsilon, dp.delta, dp.k) {
    # alt k triangle noise
    ktri.noise.smooth <- lapply(nws.sample$deg, draw.lap.noise.smooth.term,
                                 term = "gwesp", param = 0.5, dp.epsilon = dp.epsilon, dp.delta = dp.delta, max.deg=NULL)
    
    terms <- list(c(list("name" = c("gwesp"), "inputs" = c(0.5))))
    ktri.noise.restr <- draw.lap.noise.restricted(terms, dp.epsilon, dp.k, "edge")
    ktri.noise.level <- list("smooth" = mean(get.levels(ktri.noise.smooth)), "restricted" = ktri.noise.restr$level)

    # gwdsp noise
    ktwop.noise.smooth <- lapply(nws.sample$sp, draw.lap.noise.smooth.term, 
                                    term = "gwdsp", param = 0.5, dp.epsilon = dp.epsilon, dp.delta = dp.delta, max.sp=NULL)
    terms <-list(c(list("name" = c("gwdsp"), "inputs" = c(0.5))))
    ktwop.noise.restr <- draw.lap.noise.restricted(terms, dp.epsilon, dp.k, "edge")
    ktwop.noise.level <- list("smooth" = mean(get.levels(ktwop.noise.smooth)), "restricted" = ktwop.noise.restr$level)

    return(list("ktri" = ktri.noise.level, "ktwop" = ktwop.noise.level, 
                "max.deg" = mean(nws.maxdegrees), "max.sp" = mean(nws.maxsharedpartners)))
}