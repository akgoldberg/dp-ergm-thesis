################################################################################
##                                  Sample ERGM                               ##
################################################################################
sample_ergm <- function(n, theta, nsim=50, verbose=TRUE) {
    MCMC_control <- control.simulate.formula(MCMC.burnin=2*choose(n,2),
                                         MCMC.interval=0.1*choose(n,2),
                                         MCMC.prop.weights="TNT",
                                         MCMC.prop.args=list(),
                                         MCMC.packagenames=c(),
                                         MCMC.runtime.traceplot=TRUE,
                                         network.output="network")

  # generate samples
  # gwesp - same as alt-k-triangle, gwdsp - same as alt-k-two-path
  sim <- simulate(network(n, directed=FALSE) 
                  ~ edges + altkstar(0.5, fixed=TRUE) 
                  + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE),
           coef = theta,
           sequential = TRUE,
           nsim = nsim,
           control = MCMC_control,
           verbose=verbose)
  
  sample <- filter.null(sim)
  
  # compute degrees and shared partners of sampled networks
  deg <- dist.samples(sample, "deg")
  sp <- dist.samples(sample, "sp")
  altktri <- dist.samples(sample, "altktri")
  altktwopath <- dist.samples(sample, "altktwopath")
  altkstar <- dist.samples(sample, "altkstar")
  edge <- dist.samples(sample, "edge")
  
  return (list("sample" = sample, "deg" = deg, "sp" = sp,
               "n" = n, "theta" = theta, "altktri" = altktri,
               "altktwopath" = altktwopath, "edge" = edge))
}

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

# compute distribution of graph statistics
dist.samples <- function(nws, type) {
    if (type == "deg") return(unlist(lapply(nws, network.maxdegree)))
    if (type == "sp") return(unlist(lapply(nws, network.maxsharedpartners)))
    if (type == "edge") return(unlist(lapply(nws, network.numedges)))
    if (type == "altkstar") return(unlist(lapply(nws, network.altkstar)))
    if (type == "altktri") return(unlist(lapply(nws, network.altktri)))
    if (type == "altktwopath") return(unlist(lapply(nws, network.altktwopath)))
}

# filter nulls from outputted samples
filter.null <- function(samp) {
  return(samp[!sapply(samp, is.null)]) 
}

# extract level field of list of noise
get.levels <- function(noise) {
  sapply(noise, function (o) {o$level})
}

# extract draw field of list of noise
get.draws  <- function(noise) {
  sapply(noise, function (o) {o$draw})
}

################################################################################
##             Compute Noise Added on Sample for Smooth vs. Restricted        ##
################################################################################
test.noise.sample <- function(nws.sample, dp.epsilon, dp.delta, dp.k) {
    # alt k triangle noise
    ktri.noise.smooth <- lapply(nws.sample$sp, draw.lap.noise.smooth.term,
                                 term = "gwesp", param = 0.5, dp.epsilon = dp.epsilon, dp.delta = dp.delta, max.deg=NULL)
    
    terms <- list(c(list("name" = c("gwesp"), "inputs" = c(0.5))))
    ktri.noise.restr <- draw.lap.noise.restricted(terms, dp.epsilon, dp.k, "edge")
    ktri.noise.level <- list("smooth" = mean(get.levels(ktri.noise.smooth)), "restricted" = ktri.noise.restr$level)

    # gwdsp noise
    ktwop.noise.smooth <- lapply(nws.sample$deg, draw.lap.noise.smooth.term, 
                                    term = "gwdsp", param = 0.5, dp.epsilon = dp.epsilon, dp.delta = dp.delta, max.sp=NULL)
    terms <-list(c(list("name" = c("gwdsp"), "inputs" = c(0.5))))
    ktwop.noise.restr <- draw.lap.noise.restricted(terms, dp.epsilon, dp.k, "edge")
    ktwop.noise.level <- list("smooth" = mean(get.levels(ktwop.noise.smooth)), "restricted" = ktwop.noise.restr$level)
    
    return(list("ktri" = ktri.noise.level, "ktwop" = ktwop.noise.level))
}

dist.noise.samples <- function(nws.samples, dp.epsilon, dp.delta) {
  out <- data.frame()
  for (nws.sample in nws.samples) {
    dp.k <- max(nws.sample$deg) + 1
    new <- test.noise.sample(nws.sample, dp.epsilon, dp.delta, dp.k)
    new$n <- nws.sample$n
    new$altktri <- mean(nws.sample$altktri)
    new$altktwopath <- mean(nws.sample$altktwopath)
    new$dp.k <- dp.k
    out <- rbind(out, data.frame(rbind(unlist(new))))
  }
  return(out)
}

plot.noise.samples <- function(out, dp.epsilon) {
  xrange <- range(out$n)
  yrange <- range(out$ktri.restr, out$ktri.smooth)
  yrange[2] <- yrange[2] + 10
  yrange[1] <- yrange[1] - 10
  plot(xrange, yrange,  type="n", xlab="n", ylab="Noise Level", main=)
  pchs <- c(16,17)
  colors = c("dark green", "dark red")
  lines(out$n, out$ktri.restr, type='b',  col=colors[1], pch=pchs[1])
  lines(out$n, out$ktri.smooth, type='b',  col=colors[2], pch=pchs[2])
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  legend("topright", inset=c(-0.3, 0), legend = c("Restricted", "Smooth"), 
         col=colors, pch=pchs, cex=0.8)
  title.text <- bquote(paste(epsilon, "=",.(dp.epsilon)))
  title(main=title.text, font.main=4)
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}



################################################################################
##                                Generate/Load Samples                       ##
################################################################################

# generate samples from n = start to 1000
generate.samples <- function(samp.id, start=100) {
  for (n in seq(start,1000, 100)) {
    tic(sprintf("Sample n=%d", n))
    tic(sprintf("n=%d - sampling", n))
    print(sprintf("Starting sample with %d nodes...", n))
    p <- log(n)/(2*n)
    theta1 <- log(p/(1-p))
    theta3 <- 0.
    theta <- c(theta1, 0, theta3, 0)
    samp <- sample_ergm(n, theta, 50, verbose=TRUE)
    print(sprintf("Samples have avg. edges: %g, max degree: %d, max shared partners: %d", mean(samp$edge), max(samp$deg), max(samp$sp)))
    toc()
    tic(sprintf("n=%d - saving", n))
    print(sprintf("Saving sample with %d nodes...", n))
    samp_name <- sprintf("obj/samples/sample%d-%d", samp.id, n)
    save(samp, file=samp_name)
    remove(samp)
    toc()
    toc()
  }
}

# load generated samples
load.samples <- function(samp.id, start=100, stop=1000) {
  i <- 1
  samples <- list()
  for (n in seq(start,1000, 100)) {
    samp_name <- sprintf("obj/samples/sample%d-%d", samp.id, n)
    load(samp_name)
    samples[[i]] <- samp
    remove(samp)
    i <- i+1
  }
  return(samples)
}
