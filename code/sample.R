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
##                                Generate/Load Samples                       ##
################################################################################

# generate samples from n = start to 1000
generate.samples <- function(samp.id, start=100, stop=1000, theta=NULL, verbose=FALSE) {
  for (n in seq(start,stop, 100)) {
    tic(sprintf("Sample n=%d", n))
    tic(sprintf("n=%d - sampling", n))
    print(sprintf("Starting sample with %d nodes...", n))
    p <- log(300)/(2*300)
    theta[1] <- log(p/(1-p))
    samp <- sample_ergm(n, theta, 50, verbose=verbose)
    print(sprintf("Samples have avg. edges: %g, max degree: %d, max shared partners: %d", mean(samp$edge), max(samp$deg), max(samp$sp)))
    toc()
    tic(sprintf("n=%d - saving", n))
    print(sprintf("Saving sample with %d nodes...", n))
    samp_name <- sprintf("obj/samples.final/sample%d-%d", samp.id, n)
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
  for (n in seq(start, stop, 100)) {
    samp_name <- sprintf("obj/samples.final/sample%d-%d", samp.id, n)
    load(samp_name)
    samples[[i]] <- samp
    remove(samp)
    i <- i+1
  }
  return(samples)
}



# #### GENERATE SAMPLES ####
# theta1.factor <- 2.
# theta1 <- 0
# theta2 <- 0.
# theta3 <- 1.
# theta4 <- 0.
# theta <- c(theta1, theta2, theta3, theta4)
# generate.samples(1, theta=theta)

# theta1.factor <- 2.
# theta1 <- 0
# theta2 <- 0.
# theta3 <- -1.
# theta4 <- 0.
# theta <- c(theta1, theta2, theta3, theta4)
# generate.samples(2, theta=theta)

# theta1.factor <- 2.
# theta1 <- 0.
# theta2 <- 0.
# theta3 <- 2.0
# theta4 <- -0.1
# theta <- c(theta1, theta2, theta3, theta4)
# generate.samples(3, theta=theta)

# theta1.factor <- 2.
# theta1 <- 0.
# theta2 <- 0.
# theta3 <- -2.0
# theta4 <- 0.1
# theta <- c(theta1, theta2, theta3, theta4)
# generate.samples(4, theta=theta)

# theta1.factor <- 2.
# theta1 <- 0.
# theta2 <- 2.0
# theta3 <- 2.0
# theta4 <- -0.5
# theta <- c(theta1, theta2, theta3, theta4)
# generate.samples(5, theta=theta)


################################################################################
##             Compute Noise Added on Sample for Smooth vs. Restricted        ##
################################################################################
test.noise.sample <- function(nws.sample, dp.epsilon, dp.delta, dp.ks, dp.knames) {
    out <- vector("list", length(dp.ks))

    for (i in 1:length(dp.ks)) {
      dp.k <- dp.ks[i]
      dp.kname <- dp.knames[i]

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

      noise.out <- list("ktri" = ktri.noise.level, "ktwop" = ktwop.noise.level, "ktype" = dp.kname)
      out[[i]] <- noise.out
    }
    
    return(out)
}

dist.noise.samples <- function(nws.samples, dp.epsilon, dp.delta) {
  out <- data.frame()
  for (nws.sample in nws.samples) {
    ######## TEST FOR DIFFERENT k's####
    dp.k.min <- min(nws.sample$deg)
    dp.k.med <- median(nws.sample$deg)
    dp.k.max <- max(nws.sample$deg)
    dp.k.cons <- 1.5*max(nws.sample$deg)
    dp.ks <- c(dp.k.min, dp.k.med, dp.k.max, dp.k.cons)
    dp.knames <- c("min", "med", "max", "conservative")
    
    new <- test.noise.sample(nws.sample, dp.epsilon, dp.delta, dp.ks, dp.knames)
    new$n <- nws.sample$n
    new$altktri <- mean(nws.sample$altktri)
    new$altktwopath <- mean(nws.sample$altktwopath)
    new$dp.k <- dp.k
    out <- rbind(out, data.frame(rbind(unlist(new))))
  }
  return(out)
}

plot.noise.sample <- function(nws.samples, dp.epsilon, stat) {
  df <- dist.noise.samples(nws.samples, dp.epsilon, 1e-6)
  if (stat == 'ktri') {
    restr.noise <- df$ktri.restr
    smooth.noise <- df$ktri.smooth
  }
  if (stat == 'ktwopath') {
    restr.noise <- df$ktwop.restr
    smooth.noise <- df$ktwop.smooth
  }
  title.text <- bquote(paste(epsilon, "=",.(dp.epsilon)))
  plt <- ggplot(df, aes(x=n)) + 
    geom_point(aes(y = restr.noise, color='restricted'), size=2) +
    geom_point(aes(y = smooth.noise, color='smooth'), size=2) +
    geom_line(aes(y = restr.noise, color='restricted')) +
    geom_line(aes(y = smooth.noise, color='smooth')) +
    ggtitle(title.text) +
    ylab('Noise Level') +
    theme(axis.text.x = element_text(size=6, angle=45),
          axis.text.y = element_text(size=6),
          axis.title.x = element_text(size=8, margin=unit(c(-5, 0, 0, 0), "pt")),
          axis.title.y = element_text(size=8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 9)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) +
    labs(color='Sensitivity:') +
    theme(plot.margin = unit(c(-15,1,5,5), "pt")) +
    theme(legend.position="top", legend.title = element_text(size=10),
          legend.text = element_text(size=10), 
          legend.key.size = unit(0.25, "line"),
          legend.background = element_rect(colour = 'black', size = 0.1, linetype='solid'),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'))
  return(plt)
}

# Plot the average noise level per sample
plot.noise.samples <- function(nws.samples, stat, title.string, dp.epsilons = c(0.1, 0.5, 1.)) {
    plts <- list()
    i <- 1
    for (dp.epsilon in dp.epsilons) {
      plts[[i]] <- plot.noise.sample(nws.samples, dp.epsilon, stat)
      i <- i + 1
    }
    df.summary <- dist.noise.samples(nws.samples, 1.0, 1e-6)
    
    # generate plot with statistic value for sense of scale
    if (stat == 'ktri') {
      plt.summary <- ggplot(df.summary, aes(x=n, y=altktri))
    }
    else if (stat == 'ktwopath') {
      plt.summary <- ggplot(df.summary, aes(x=n, y=altktwopath))
    } 
    else {
      print("Valid stats: ktri or ktwopath")
      return()
    }
    plts[[i+1]] <- 
      plt.summary + 
      geom_point(size=2) +
      geom_line() +
      ggtitle("Computed Statisic Value") +
      ylab('Value') +
      theme(axis.text.x = element_text(size=6, angle=45),
            axis.text.y = element_text(size=6),
            axis.title.x = element_text(size=8, margin=unit(c(-5, 0, 0, 0), "pt")),
            axis.title.y = element_text(size=8),
            plot.title = element_text(size=12, face='plain')) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
    
    legend <- get_legend(plts[[1]])
    plts <- lapply(plts, function (plt)  plt + theme(legend.position="none"))
    
    prow <- plot_grid(plotlist=plts) 
    prowlegend <- plot_grid(legend, prow, ncol=1, rel_heights = c(.2, 1))
    title <- ggdraw() + draw_label(title.string, fontface='bold', size = 12)
    plot_grid(title, prowlegend + theme(plot.margin = unit(c(-15,1,5,5), "pt")), ncol=1, rel_heights=c(0.1, 1))
}


# Summarize samples
summarize.samples <- function(samp.ids = c(0,1,3,5)) {
    df.rows <- data.frame()
    for (id in samp.ids) {
      samples <- load.samples(id)
      print(sprintf("Loaded sample %d", id))
      for (i in seq(1, 10)) {
        row = list("n" = i*100, "samp.id" = sample.map.id(id))
        row$density.mean <- mean(samples[[i]]$edge/choose(i*100,2))
        row$triangles.mean <- mean(as.numeric(lapply(samples[[i]]$sample, network.triangles)))
        row$deg.min <- min(samples[[i]]$deg)
        row$deg.median <- median(samples[[i]]$deg)
        row$deg.max <- max(samples[[i]]$deg)
        df.row <- data.frame(row)
        df.rows <- rbindlist(list(df.row, df.rows))
      }
    }
    #write.table(df.summary.samples, file = sprintf("obj/df/samples.summarydf.txt",i), sep = ",", col.names = colnames(df.summary.samples))
    return(df.rows)
}

plot.summary.samples <- function() {
  df.summary.samples <- read.table(file = sprintf("obj/df/samples.summarydf.txt",i), sep = ",", header=TRUE)
  plt.degs <- ggplot(data=df.summary.samples, mapping=aes(fill = factor(samp.id,levels = rev(unique(samp.id))))) +
    geom_line(aes(x=n, y=deg.median), color="grey", linetype=2) +
    geom_ribbon(aes(x=n, ymax=deg.max, ymin=deg.min), alpha=0.4) +
    scale_x_continuous(name='n', breaks=seq(200, 1000, 200)) +
    scale_y_continuous(name='Degree') +
    scale_linetype_manual('', values = ltps) +
    ggtitle("Degree of Simulated Networks") +
    scale_fill_discrete(name = "Model") +
    theme_gray() + 
    theme(plot.title = element_text(size=12, face='plain', hjust = 0.5))
  plt.tris <- ggplot(data=df.summary.samples, mapping=aes(color = factor(samp.id,levels = rev(unique(samp.id))))) +
    geom_line(aes(x=n, y=triangles.mean)) +
    scale_x_continuous(name='n', breaks=seq(200, 1000, 200)) +
    scale_y_continuous(name='Number of Triangles') +
    scale_linetype_manual('', values = ltps) +
    ggtitle("Triangle Count of Simulated Networks") +
    scale_color_discrete(name = "Model") +
    theme_gray() + 
    theme(plot.title = element_text(size=12, face='plain', hjust = 0.5))
  plts <- plot_grid(plt.degs, plt.tris)
  ggsave(filename = "plots/samples.summary.png", plot = plts, width=10, height=5)
  return(plts)
}

sample.map.id <- function(sample.id) {
  switch (as.character(sample.id),
    '0' = 0,
    '1' = 1,
    '3' = 2,
    '5' = 3,
    '6' = 4,
    '7' = 5,
    '8' = 6
  )
}
