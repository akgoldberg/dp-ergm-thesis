################################################################################
##                                 Testing Noise                              ##
################################################################################

# run tests comparing noise of smooth and restricted
run.noise.tests <- function(samples, dp.epsilons, dp.delta=1e-6, stop=1000) {
  # dataframe to output
  df <- data.frame("n"=numeric(), "type"=character(), "stat.name"=character(),
                   "dp.k"=numeric(), "eps"=numeric(), "delta"=numeric(),
                   "scale"=numeric(), "stat.value"=numeric(), "bias"=numeric(),
                   "rmse"=numeric(), "dp.klevel"=character())
  
  # iterate over all samples (over synthetic networks of size n=100,200,...)
  for (sample in samples) {
    if (sample$n > stop) {
      break
    }
    
    print(sprintf("Testing samples of size %d", sample$n))
    for (dp.epsilon in dp.epsilons) {
          print(sprintf("Epsilon=%g", dp.epsilon))
          tic("Restricted")
          df <- rbindlist(list(df, run.noise.tests.restr(sample, dp.epsilon)),
                          use.names=TRUE, fill=TRUE)
          toc()

          tic("Smooth")
          df <- rbindlist(list(df, run.noise.tests.smooth(sample, dp.epsilon, dp.delta)),
                          use.names=TRUE, fill=TRUE)
          toc()
      }
    }
  return(df)
}

# Test restricted sensitivty for model with different values of k
run.noise.tests.restr <- function(sample, dp.epsilon) {
  # test (min, median, max, 1.5*max)
  dp.ks <- c(0.75*min(sample$deg), min(sample$deg), median(sample$deg), max(sample$deg), 1.5*max(sample$deg))
  dp.klevels <- c("aggressive", "min", "median", "max", "conservative")
 
  # collect data frames for different k
  df.rows <- data.frame()

  for (iter in 1:length(dp.ks)) {
    dp.k <- dp.ks[iter]
    
    # set up row to add to dataframe
    row <- list()
    row$n <- sample$n
    row$type <- "restricted"
    row$dp.k <- dp.k
    row$dp.klevel <-dp.klevels[iter] 
    row$eps <- dp.epsilon
    row$delta <- NULL

    # scale of laplace noise for value of k
    scale <- make.private.restr(sample$sample[[1]] ~ edges + altkstar(0.5, fixed=TRUE) 
                                      + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE),
                                      rep(dp.epsilon, 4), dp.k, proj=FALSE)$noise$level
    
    row$scale <- scale
  
    # averages over all datapoints
    true.stats.tot <- 0
    rmse.tot <- 0
    bias.tot <- 0

    N <- length(sample)
    # iterate over nws in the sample
    for (i in 1:N) {
      nw <- sample$sample[[i]]
      true.stats <- summary(nw  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
      
      ### restricted sensitivity test ###
      # get scale of laplace noise
        
        var <- 2*(scale)**2
        # if over-estimated degree, then don't need to get expected bias
        if (dp.k >= sample$deg[i]) {
          bias.sq.mean <- 0
          bias.mean <- 0
        } else {
            # get average bias^2 from projections
            n.proj <- 10
            bias.sq.mean <- 0
            bias.mean <- 0
            for (i in 1:n.proj) {
              nw.proj <- projection.edge(nw, dp.k)
              proj.stats <- summary(nw.proj  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
              bias.sq.mean <- bias.sq.mean + (true.stats - proj.stats)**2
              bias.mean <- bias.mean + (proj.stats - true.stats)
            }
            bias.sq.mean <- bias.sq.mean/n.proj
            bias.mean <- bias.mean/n.proj
        }
      
        # calculate rmse and relative rmse (bias^2 + variance)
        rmse <- sqrt(bias.sq.mean + var)

        # update aggregate stats
        true.stats.tot <- true.stats.tot + true.stats
        rmse.tot <- rmse.tot + rmse
        bias.tot <- bias.tot + bias.mean
    }

    row$stat.value <- true.stats.tot/N
    row$bias <- bias.tot/N
    row$rmse <- rmse.tot/N

    df.row <- data.frame(row)
    setDT(df.row, keep.rownames = TRUE)
    setnames(df.row, 1, "stat.name") 
    df.rows <- rbindlist(list(df.row, df.rows))
  }

  return(df.rows)
}

# Test smooth sensitivity draws
run.noise.tests.smooth <- function(sample, dp.epsilon, dp.delta) {
    row <- list()
    row$n <- sample$n
    row$type <- "smooth"
    row$dp.k <- NULL
    row$eps <- dp.epsilon
    row$delta <- dp.delta
    row$bias <- 0
    
  
    # averages over all datapoints
    scale.tot <- 0
    true.stats.tot <- 0
    rmse.tot <- 0

    N <- length(sample)
    # iterate over nws in the sample
    for (i in 1:N) {
      nw <- sample$sample[[i]]
      sp <- sample$sp[[i]]
      deg <- sample$deg[[i]]
      true.stats <- summary(nw  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))

      n.draws <- 50
      rmse <- 0
      scale.mean <- 0
      for (i in 1:n.draws) {
        # scale of laplace noise for value of k
        noise <- make.private.smooth(nw ~ edges + altkstar(0.5, fixed=TRUE) 
                                          + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE),
                                          rep(dp.epsilon, 4), rep(dp.delta, 4),
                                          max.deg=deg, max.sp=sp)$noise
        scale.mean <- scale.mean + noise$level
        rmse <- rmse + (noise$draw)**2
      }
      scale.mean <- scale.mean/n.draws
      rmse <- sqrt(rmse/n.draws)
      
      true.stats.tot <- true.stats.tot + true.stats
      rmse.tot <- rmse.tot + rmse
      scale.tot <- scale.tot + scale.mean
  }

  row$stat.value <- true.stats.tot/N
  row$scale <- scale.tot/N
  row$rmse <- rmse.tot/N

  df.row <- data.frame(row)
  setDT(df.row, keep.rownames = TRUE)
  setnames(df.row, 1, "stat.name") 
  return(df.row)
}

# visualize results of test
visualize.noise.tests <- function(tests.id,
                                  dp.klevels=c('min','median','max'),
                                  stat.names=c('edges', 'altkstar.0.5', 'gwesp.fixed.0.5', 'gwdsp.fixed.0.5')) {
  df.samples <- data.table(read.table(file=sprintf("obj/df/df.samples%d.txt",tests.id), sep = ",", header=TRUE))
  
  # refer to smooth as 'private local'
  df.samples$type = revalue(df.samples$type, c('smooth' = 'private local'))
  
  # get epsilon values in dataframe
  eps.values <- unique(df.samples$eps)
  
  # get line type to use (dotted for smooth)
  ltps <- rep(1, length(dp.klevels)+1)
  ltps[1] <- 3
  
  # iterate over stats to generate plots
  for (this.stat.name in stat.names) {
    # store plots for each epsilon value
    plts <- vector(mode='list', length(eps.values)+1)
    # generate plots for each value of eps
    for (i in 1:length(eps.values)) {
        this.eps <- eps.values[i]
        df.data <- df.samples[(eps==this.eps) & (stat.name==this.stat.name) 
                                & (dp.klevel %in% dp.klevels | is.na(dp.klevel))]
        plts[[i]] <- 
          ggplot(data=df.data, mapping=aes(x =n, y=rmse, color=gsub(" NA", "", paste(type, dp.klevel)), linetype=gsub(" NA", "", paste(type, dp.klevel)))) +
           geom_point() +
           geom_line() +
           guides(shape=FALSE) +
           scale_x_continuous(name='n', breaks=seq(min(df.data$n),max(df.data$n), 200)) +
           scale_y_continuous(name='RMSE') +
           scale_color_brewer('', type='seq', palette = 'YlGnBu', direction=0.5) + 
           scale_linetype_manual('', values = ltps) +
           ggtitle(bquote(paste(epsilon, "=",.(this.eps)))) +
           theme_gray() + 
           theme(plot.title = element_text(size=12, face='plain', hjust = 0.5))
    }
    df.statvalues <- unique(df.samples[stat.name==this.stat.name, n, stat.value])
    # make plot of stat value
    plts[[length(eps.values)+1]] <- 
      ggplot(data=df.data, mapping=aes(x =n, y=stat.value)) +
        geom_point() +
        geom_line() +
        scale_x_continuous(name='n', breaks=seq(min(df.data$n),max(df.data$n), 200)) +
        scale_y_continuous(name='Value') +
        ggtitle("True Statistic") +
        theme_gray() +
        theme(plot.title = element_text(size=12, face='plain', hjust = 0.5))
    
    # get legend from plot
    legend <- get_legend(plts[[1]] +
                         theme(legend.position="bottom",
                                legend.justification="center",
                                legend.background = element_rect(colour = 'black', size = 0.1, linetype='solid'),
                                legend.box.margin=margin(-15, -25, 0, 0, "pt")) +
                         guides(color=guide_legend(ncol=2))
    )
    # hide legends on plots
    plts <- lapply(plts, function (plt)  plt + theme(legend.position="none"))
    # put everything together
    plot.name <- sprintf("Model %d - %s", tests.id, get.statname(this.stat.name))
    title <- ggdraw() + draw_label(plot.name, fontface='bold', size = 12)
    prow <- plot_grid(plotlist=plts)
    final.plt <- plot_grid(title, prow, legend, ncol=1, rel_heights = c(.1, 1., .15))
    ggsave(filename = sprintf("plots/noise.final/%s.png", plot.name), plot = final.plt)
  }
}

################################################################################
################################################################################
##                                 Testing Inference                          ##
################################################################################
################################################################################

run.inference.tests <- function(samp.id, n, formula.rhs, dp.epsilon = 1.0, num.tests = 10,
                               burn.in=10000, main.iters=10000, sigma.epsilon=NULL, parallel=TRUE, non.private=TRUE,
                               method='restr') {
  
  
  samples <- load.samples(samp.id, start=n, stop=n)[[1]]
  nw <- samples$sample[[25]]
  dp.k <- median(samples$deg)
  true.theta <- samples$theta

  formula <- nonsimp.update.formula(formula.rhs, nw ~ ., from.new=TRUE)
  print(sprintf("Running %d tests with sample=%d, n=%d, k=%d, eps=%g, method=%s", num.tests, samp.id, n, dp.k, dp.epsilon, method))

  tic("Runtime")
  
  if (non.private) {
    print("Non-Private")
    nonprivate.out <- bergm.orig(formula,
                        burn.in=burn.in, main.iters = main.iters, aux.iters = 0.1*choose(n,2),
                        sigma.epsilon = sigma.epsilon,
                        print.out=2500, nchains = 3)
    
  } else {
    nonprivate.out <- NULL
  }
  
  # setup up closure for parallel processing
  test.run <- function(t) {
    run.one.test(t, formula, n, dp.epsilon, dp.k, burn.in, main.iters, 
                  sigma.epsilon, parallel, method)
  }
  
  if (parallel) {
    num.cores <- min(num.tests+1, detectCores()-1)
    print(sprintf("Can use %d cores.", num.cores))
    outs <- mclapply(as.list(1:num.tests), test.run, mc.preschedule=FALSE, mc.cores	= num.cores)
  } else {
    # run tests sequentially
    for (t in 1:num.tests) {
      outs <- list()
      outs[[t]] <- test.run(t)
    }
  }
  
  toc()
    
  return(list("private"=outs, "nonprivate"=nonprivate.out, "true"=true.theta))
}

run.one.test <- function(t,
                         formula, 
                         n,
                         dp.epsilon, 
                         dp.k, 
                         burn.in, 
                         main.iters, 
                         sigma.epsilon,
                         parallel,
                         method,
                         print.out = 2500,
                         nchains = 3) {
  
  print(sprintf("Test: %d", t))
  if (!parallel) {
    tic()
  }
  if (method == "smooth") {
    dp.delta <- 1e-6
    nw.private <- make.private.smooth(formula, dp.epsilon, dp.delta)
  } 
  if (method == "rr") {
     nw.private <- make.private.rr(formula, dp.epsilon)
  }
  if (method == "restr") {
     nw.private <- make.private.restr(formula, dp.epsilon, dp.k, privacy.type = "edge")
  }
  
  private.out <- bergm.modified.private(nw.private$formula, nw.private$noise,
                                        burn.in=burn.in, main.iters = main.iters, aux.iters = 0.2*choose(n,2),
                                        sigma.epsilon = sigma.epsilon,
                                        print.out=1000, nchains = 3)
  
  if (!parallel) toc()

  return(private.out)
}

load.inference.tests <- function(i, method, dp.epsilon) {
  fname = sprintf("obj/inference.tests/inference.tests%d%s-eps%g", i, method, dp.epsilon)
  load(fname)
  return(inference.tests)
}

get.summary.tests <- function(private.tests) {
  num.tests <- length(private.tests)
  test.summary <- rep(NA, num.tests)
  for (i in 1:num.tests) {
    x <- private.tests[[i]]
    
    # get acceptance rates
    rates <- matrix(x$AR,x$nchains,1)
    rownames(rates) <- paste("Chain",seq(1,x$nchains)," ")
    colnames(rates) <- paste("Acceptance rate:")
    rates <- as.table(rates)
    
    FF <- apply(x$Theta,2,cbind)
    # get posterior sd and mean
    overall <- rbind(apply(FF,2,mean),apply(FF,2,sd))
    rownames(overall) <- c("Post. mean","Post. sd")
    colnames(overall) <- paste("theta",seq(1,x$dim)," (",
                               x$specs[seq(1,x$dim)],")",sep="")
    all <- as.table(overall)
    
    test.summary[i] <- all
  }
  return(test.summary)
}
