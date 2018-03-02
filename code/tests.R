
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
  #dp.ks <- c(min(sample$deg), median(sample$deg), max(sample$deg), 1.5*max(sample$deg))
 # dp.klevels <- c("min", "median", "max", "conservative")
 dp.ks <- c(0.75*min(sample$deg))
 dp.klevels <- c("aggressive")
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


# # 
# run.inference.tests <- function(formula, n, true.theta, dp.k, dp.epsilon = 1.0, num.tests = 10) {
#   print(sprintf("Running tests with n=%d, k=%d, eps=%g", n, dp.k, dp.epsilon))
  
#   tic("Runtime")
  
#  # print("Non-Private")
#  # nonprivate.out <- bergm.orig(formula,
#   #                   burn.in=5000, main.iters = 7500, aux.iters = 0.1*choose(n,2),
#    #                  sigma.epsilon = diag(c(0.001, 0.001, 0.0001)),
#     #                 print.out=2500, nchains = 1)
  
  
#   # run tests
#   outs <- list()
#   for (t in 1:num.tests) {
#     print(sprintf("Test: %d", t))
#     nw.private <- make.private.restr(formula, dp.epsilon, dp.k, privacy.type="edge")
#     private.out <- bergm.modified.private(nw.private$formula, nw.private$noise,
#                                           burn.in=7500, main.iters = 7500, aux.iters = 0.2*choose(n,2),
#                                           sigma.epsilon = diag(c(0.001, 0.001, 0.0001)),
#                                           print.out=1000, nchains = 1)
#     outs[[t]] <- private.out
#   }
#   toc()
  
#   return(list("private"=outs, "nonprivate"=nonprivate.out, "true"=true.theta))
# }

# tests0.higheps <- run.inference.tests(nw ~ edges + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE), n, true.theta, dp.k, dp.epsilon = 1.0)
# tests0.loweps <- run.inference.tests(nw ~ edges + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE), n, true.theta, dp.k, dp.epsilon = 0.5)