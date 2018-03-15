# make a data frame with multiple methods and epsilons from inference tests
make.inference.df <- function(model.id, methods=c("smooth", "restr", "rr"), dp.epsilons=c(1, 3)) {
  df.out <- data.frame()
  
  # get nonprivate df
  nonpriv.tests <- load.inference.tests(model.id, "nonprivate")
  nonpriv.df <- make.inference.testsdf.row(nonpriv.tests, nonpriv.tests$true, model.id, 0)
  df.out <- rbindlist(list(df.out, nonpriv.df), use.names=TRUE, fill=TRUE)
  
  test.id.start <- 0
  for (method in methods) {
    for (dp.epsilon in dp.epsilons) {
      fname <- sprintf("obj/inference.tests/inference.tests%d%s-eps%g", model.id, method, dp.epsilon)
      if (file.exists(fname)) {
        df.out <- rbindlist(list(df.out, get.summary.inference.tests(model.id, method, dp.epsilon, test.id.start)),
                            use.names=TRUE, fill=TRUE)
        test.id.start <- max(df.out$test.num)
      } else {
        print(sprintf("Could not load data for model %d, method %s, eps %g", model.id, method, dp.epsilon))
      }
    }
  }
  return(df.out)
}

# Make raw test data into dataframe
get.summary.inference.tests <- function(model.id, method, dp.epsilon, test.id.start) {
  # load the tests for this model, method and epsilon value
  tests <- load.inference.tests(model.id, method, dp.epsilon)
  df.out <- data.frame("model.id"=numeric(), "method"=character(), "stat.name"=character(),
                   "eps"=numeric(), "delta"=numeric(),
                    "stat.value"=numeric(),  "AR"=numeric(),
                   "post.mean"=numeric(), "post.se"=character(),
                   "noise.draw"=numeric(), "noise.level"=character(), "KL"=numeric())

  for (t in 1:length(tests$private)) {
    row.df <- make.inference.testsdf.row(tests$private[[t]], tests$true, model.id, t+test.id.start)
    df.out <- rbindlist(list(df.out, row.df), use.names=TRUE, fill=TRUE)
  }
  
  return(df.out)
}

# make one row of test df
make.inference.testsdf.row <- function(x, true.theta, model.id, test.num) {
    # get rid of 0's
    true.theta <- true.theta[true.theta != 0]
  
    row <- list()
    row$model.id <- model.id
    row$test.num <- test.num
    
    statnames <- names(x$stats)
    row$post.mean <- apply(x$Theta,c(2),mean)
    row$post.se <- getSE(x$Theta)
    names(row$post.mean) <- statnames
    names(row$post.se) <- statnames
    row$stat.value <- x$stats
    row$true.param.value <- true.theta
    row$AR <- mean(x$AR)
    
    row$KL <- computeKL(x$formula, true.theta, row$post.mean) 
    
    # non-private test inference run
    if (test.num == 0) {
      row$method <- "nonprivate"
    } 
    # private test inference run
    else {
      row$method <- x$noise$method
      row$eps <- x$noise$dp.epsilon
      row$delta <- x$noise$dp.delta
      row$noise.draw <- x$noise$draw
      row$noise.level <- x$noise$level
    }
    df.row <- data.frame(row)
    setDT(df.row, keep.rownames = TRUE)
    setnames(df.row, 1, "stat.name") 
    return(df.row)
}

# compute standard error
getSE <- function(Theta.out) {
  sd <- apply(Theta.out, 2, sd)
  eff.size <- apply(apply(Theta.out, 3, effectiveSize), 1, sum)
  return(sd/eff.size)
}

computeKL <- function(formula, true.theta, theta.other) {
  llr <- ergm.bridge.llr(formula,
                         from=true.theta, to=theta.other,
                         #control=control.ergm.bridge(nsteps=250, MCMC.burnin=1e5),
                         llronly=TRUE, verbose = FALSE)
  
  return(llr)
}

extract.nonprivate <- function(i) {
  tests <- load.inference.tests(i, "restr", 1.)
  #nonprivate <- load.inference.tests(i, "nonprivate", 1)
  nonprivate <- tests$nonprivate
  nonprivate$true <- tests$true
  fname = sprintf("obj/inference.tests/nonprivate%i", i)
  save(nonprivate, file=fname)
}
