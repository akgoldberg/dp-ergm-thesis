################################################################################
##                                 Testing Noise                              ##
################################################################################

# run tests comparing noise of smooth and restricted
run.noise.tests <- function(samples, dp.epsilons, dp.delta=1e-6, stop=1000, node=FALSE) {
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
          if (node) {
            tic()
            df <- rbindlist(list(df, run.noise.tests.node(sample, dp.epsilon, dp.delta)),
                            use.names=TRUE, fill=TRUE)
            toc()
          } else {
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

get.node.projection.bias <- function(sample) {
   # test (min, median, max, 1.5*max)
  dp.ks <- c(min(sample$deg), median(sample$deg), max(sample$deg), 1.5*max(sample$deg))
  dp.klevels <- c("min", "median", "max", "conservative")
 
  # collect data frames for different k
  df.rows <- data.frame()

  for (iter in 1:length(dp.ks)) {
    dp.k <- dp.ks[iter]
    
    # set up row to add to dataframe
    row <- list()
    row$n <- sample$n
    row$dp.k <- dp.k
    row$dp.klevel <- dp.klevels[iter]

    N <- length(sample)
    # iterate over nws in the sample
    for (i in 1:N) {
      nw <- sample$sample[[i]]
      true.stats <- summary(nw  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
      # if over-estimated degree, then don't need to get expected bias
      if (dp.k >= sample$deg[i]) {
        bias <- 0
        d.hat <- 0
      } 
      else {
        out <- projection.node.LP(nw, dp.k)
        nw.proj <- out$y
        d.hat <- out$d.hat
        proj.stats <- summary(nw.proj  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
        bias <- (proj.stats - true.stats)
      }
    }

    row$bias <- bias
    row$d.hat <- d.hat
    df.row <- data.frame(row)
    setDT(df.row, keep.rownames = TRUE)
    setnames(df.row, 1, "stat.name") 
    df.rows <- rbindlist(list(df.row, df.rows))
  }

  return(df.rows)
}

# visualize results of test
visualize.noise.tests <- function(tests.id,
                                  dp.klevels=c('min','median','conservative'),
                                  stat.names=c('edges', 'altkstar.0.5', 'gwesp.fixed.0.5', 'gwdsp.fixed.0.5')) {
  df.samples <- data.table(read.table(file=sprintf("obj/df/df.samples%d.txt",tests.id), sep = ",", header=TRUE))
  
  # refer to smooth as 'private local'
  df.samples$type = revalue(df.samples$type, c('smooth' = 'private local'))
  
  # get epsilon values in dataframe
  #eps.values <- unique(df.samples$eps)
  eps.values <- c(1.0)
  
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
           geom_line(size=1.25) +
           guides(shape=FALSE) +
           scale_x_continuous(name='n', breaks=seq(min(df.data$n),max(df.data$n), 200)) +
           scale_y_continuous(name='RMSE') +
           #scale_color_brewer('', type='seq', palette = 'YlGnBu', direction=0.5) + 
           scale_color_brewer('', palette = "Set1") +
           scale_linetype_manual('', values = ltps) +
           ggtitle(bquote(paste(epsilon, "=",.(this.eps/4.)))) +
           theme_gray() + 
           theme(plot.title = element_text(size=16, face='plain', hjust = 0.5), axis.text=element_text(size=14, angle=45, hjust=1, vjust=1), axis.title=element_text(size=14))
    }
    df.statvalues <- unique(df.samples[stat.name==this.stat.name, n, stat.value])
    # make plot of stat value
    plts[[length(eps.values)+1]] <- 
      ggplot(data=df.data, mapping=aes(x =n, y=stat.value)) +
        geom_point() +
        geom_line(size=1.25) +
        scale_x_continuous(name='n', breaks=seq(min(df.data$n),max(df.data$n), 200)) +
        scale_y_continuous(name='Value') +
        ggtitle("True Statistic") +
        theme_gray() +
        theme(plot.title = element_text(size=16, face='plain', hjust = 0.5), axis.text=element_text(size=14, angle=45, hjust=1, vjust=1), axis.title=element_text(size=14))
    
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
    plot.name <- sprintf("Model %d - %s", sample.map.id(tests.id), get.statname(this.stat.name))
    #title <- ggdraw() + draw_label(plot.name, fontface='bold', size = 12)
    prow <- plot_grid(plotlist=plts)
    #final.plt <- plot_grid(title, prow, legend, ncol=1, rel_heights = c(.1, 1., .15))
    #final.plt <- plot_grid(title, prow, ncol=1, rel_heights = c(.1, 1.))
    ggsave(filename = sprintf("plots/noise.to.use/brief%s.png", plot.name), width=10, height=5, plot = prow)
  }
}

################################################################################
################################################################################
##                                 Testing Inference                          ##
################################################################################
################################################################################

run.inference.tests <- function(samp.id, n, formula.rhs, dp.epsilon = 1.0, num.tests = 10,
                               burn.in=10000, main.iters=10000, sigma.epsilon=NULL, parallel=TRUE, non.private=TRUE,
                               method='restr', labels.priv=FALSE, attrs=NULL) {
  
  if (samp.id <= 5) {
    samples <- load.samples(samp.id, start=n, stop=n)[[1]]
    nw <- samples$sample[[25]]
    dp.k <- median(samples$deg)
    true.theta <- samples$theta
  } else {
    data('faux.mesa.high')
    nw <- faux.mesa.high
    dp.k <- 15
    true.theta <- NULL
  }
  
  formula <- nonsimp.update.formula(formula.rhs, nw ~ ., from.new=TRUE)
  print(sprintf("Running %d tests with sample=%d, n=%d, k=%d, eps=%g, method=%s", num.tests, samp.id, n, dp.k, sum(dp.epsilon), method))

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
                  sigma.epsilon, parallel, method, labels.priv, attrs)
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
                         labels.priv,
                         attrs,
                         print.out = 1000,
                         nchains = 3
                         ) {
  
  print(sprintf("Test: %d", t))
  if (!parallel) {
    tic()
  }
  
  if (method == 'nonpriv') {
    out <- bergm.orig(formula,
                      burn.in=burn.in, main.iters = main.iters, aux.iters = 0.2*choose(n,2),
                      sigma.epsilon = sigma.epsilon,
                      print.out=2500, nchains = 3)
    return(out)
  }
  if (method == "smooth") {
    dp.delta <- 1e-6
    nw.private <- make.private.smooth(formula, dp.epsilon, dp.delta)
  } 
  if (method == "rr") {
     nw.private <- make.private.rr(formula, dp.epsilon)
  }
  if (method == "restr") {
     nw.private <- make.private.restr(formula, dp.epsilon, dp.k, privacy.type = "edge", labels.priv=labels.priv, attrs=attrs)
  }
  
  private.out <- bergm.modified.private(nw.private$formula, nw.private$noise,
                                        burn.in=burn.in, main.iters = main.iters, aux.iters = 0.2*choose(n,2),
                                        sigma.epsilon = sigma.epsilon,
                                        print.out=print.out, nchains = 3)
  
  if (!parallel) toc()

  return(private.out)
}

load.inference.tests <- function(i, method, dp.epsilon) {
  if (method=="nonprivate") {
    fname = sprintf("obj/inference.tests/nonprivate%d", i)
    load(fname)
    return(nonprivate)
  } else {
    fname = sprintf("obj/inference.tests/inference.tests%d%s-eps%g", i, method, dp.epsilon)
    load(fname)
    return(inference.tests)
  } 
}


###############################################################################
##                      Node Level Noise Tests                               ##
###############################################################################

run.noise.node.tests <- function(samples, dp.epsilons=c(1.0, 2.0, 3.0), dp.delta=1e-6, stop=1000) {
  # dataframe to output
  df <- data.frame()
  
  # iterate over all samples (over synthetic networks of size n=100,200,...)
  for (sample in samples) {
    if (sample$n > stop) {
      break
    }
    print(sprintf("Testing samples of size %d", sample$n))
    df <- rbindlist(list(df, run.sample.noise.test.node(sample, dp.epsilons, dp.delta)),
                    use.names=TRUE, fill=TRUE)
    }
  return(data.table(df))
}

run.sample.noise.test.node <- function(sample, dp.epsilons, dp.delta) {
    df <- data.frame()

    # setup potential ks to test
    dp.ks <- c(min(sample$deg), median(sample$deg), max(sample$deg), 1.5*max(sample$deg))
    dp.klevels <- c("min", "median", "max", "conservative")
    
    N <- length(sample)
    
    for (projection.type in c('LP','trunc')) {
      tic(sprintf("%s Projection", projection.type))
      # iterate over values of cutoff k
      for (k in 1:length(dp.ks)) {
        dp.k <- dp.ks[k]
        dp.klevel <- dp.klevels[k]
        # iterate over nws in the sample
        for (i in 1:N) {
          nw <- sample$sample[[i]]
          # get projection
          if (projection.type == 'LP')  {
            if (sample$deg[[i]] <= dp.k) {
              proj.out <- list("y"=nw, "d.hat"=0)
            } 
            else {
              proj.out <- projection.node.LP(nw, dp.k)
            }
            proj.aux.info <- proj.out$d.hat
          }
          if (projection.type == 'trunc') {
            proj.out <- projection.node.trunc(nw, dp.k)
            proj.aux.info <- proj.out$Cs
          }
          formula <- formula(nw  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
          true.stats <- summary(formula)
          proj.stats <- summary(proj.out$y  ~ edges + altkstar(0.5, fixed=TRUE) + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE))
          terms <- ergm.getmodel(formula, nw)$terms
          
          for (noise.type in c('lap', 'cauchy')) {
            for (dp.epsilon in dp.epsilons) {
              row <- list("n"=sample$n, "i"=i, "dp.klevel"=dp.klevel, "proj.type"=projection.type, "noise.type"=noise.type,
                           "dp.epsilon"=dp.epsilon/4., "dp.k"=dp.k, "true.stats"=true.stats, "bias"=(proj.stats-true.stats))
              # draw laplace noise
              noise <- draw.lap.noise.restricted(terms, dp.epsilon, dp.k, privacy.type='node',
                          dp.deltaTot=dp.delta, proj.type=projection.type, noise.type=noise.type, proj.aux.info=proj.aux.info)
              row$smooth.sens <- noise$C
              row$noise.scale <- noise$level

              # make into a df and add on to current df
              df.row <- data.frame(row)
              setDT(df.row, keep.rownames = TRUE)
              setnames(df.row, 1, "stat.name") 
              df <- rbindlist(list(df, df.row), use.names=TRUE, fill=TRUE)
            }
          }
        }
      }
    toc()
    }
    return(df)
}

load.noise.test.node <- function(tests.id, trunc.only=FALSE) {
  if (trunc.only){
    add = "_node_trunc"
  }
  else {
    add = "_node"
  }
  df.samples <- data.table(read.table(file=sprintf("obj/df/df.samples%d%s.txt",tests.id,add), sep = ",", header=TRUE))
  return(df.samples)
}

plot.noise.tests <- function(tests.id, legend.pos='none') {
  df.samples <- load.noise.test.node(tests.id, trunc.only=FALSE)
  draw.noise <- function(noise.scale, noise.type) {
    if (noise.type == 'lap') {
      return(abs(rlaplace(n=1000, scale=noise.scale)))
    }
    if (noise.type == 'cauchy') {
      return(abs(rcauchy(n=1000, scale=noise.scale)))
    }
  }
  df.quantiles <- df.samples[,list("mean.stat.value"=mean(true.stats),
                                   "mean.scale"=mean(noise.scale),
                                   "median.err"=median(bias+draw.noise(noise.scale, noise.type)),
                                   "p25.err"=quantile(bias+draw.noise(noise.scale, noise.type),p=0.25),
                                   "p75.err"=quantile(bias+draw.noise(noise.scale, noise.type),p=0.75)),
                             by=list(n, stat.name, dp.klevel, dp.epsilon, proj.type, noise.type)]
  df.quantiles$stat.name <- sapply(df.quantiles$stat.name, get.statname)
  View(df.quantiles[dp.klevel == 'max' & dp.epsilon==0.5 & noise.type=='cauchy' & proj.type=='trunc',])
  ggplot(data=df.quantiles[dp.epsilon==0.5 & noise.type=='cauchy' & n>=400,], mapping = aes(x=n, y=median.err/mean.stat.value, color=dp.klevel)) +
    facet_wrap(stat.name~proj.type, nrow=4, ncol=2, scales = "free", labeller = label_wrap_gen(multi_line=FALSE)) +
    geom_line() +
    scale_color_discrete("Degree Cutoff:") +
    scale_x_continuous('n') +
    theme(axis.text=element_text(size=10), strip.text = element_text(size=8), title = element_text(size=12), legend.position=legend.pos) +
    scale_y_continuous('Relative Median Absolute Error') 
}

### ALREADY FIXED TRUNC TESTS ###
fix.noise.node.tests <- function(tests.id, trunc.only) {
  if (trunc.only){ add = "_node_trunc"}
  else { add = "_node" }
  fname <- sprintf("obj/df/df.samples%d%s.txt",tests.id,add)
  df.samples <- data.table(read.table(file=fname, sep = ",", header=TRUE))
  # fix factor of 3 in alt-k-star
  df.samples[stat.name == 'altkstar.0.5','noise.scale'] <- df.samples[stat.name == 'altkstar.0.5', noise.scale/3.]
  # fix lack of eps in alt k-two-path
  df.samples[stat.name == 'gwdsp.fixed.0.5','noise.scale'] <- df.samples[stat.name == 'gwdsp.fixed.0.5',noise.scale/dp.epsilon]
  # fix factors in front of smooth sens
  df.samples[noise.type == 'lap', 'smooth.sens'] <- df.samples[noise.type == 'lap', 2*smooth.sens]
  df.samples[noise.type == 'lap', 'noise.scale'] <- df.samples[noise.type == 'lap', 2*noise.scale]
  df.samples[noise.type == 'cauchy', 'smooth.sens'] <- df.samples[noise.type == 'cauchy', sqrt(2)*smooth.sens]
  df.samples[noise.type == 'cauchy', 'noise.scale'] <- df.samples[noise.type == 'cauchy', sqrt(2)*noise.scale]
  write.table(df.samples, file = fname, sep = ",", col.names = colnames(df.samples))
}