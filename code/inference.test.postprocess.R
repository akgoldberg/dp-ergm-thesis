# make a data frame with multiple methods and epsilons from inference tests
make.inference.df <- function(model.id, methods=c("smooth", "restr", "rr"), dp.epsilons=c(1, 2, 3, 4)) {
  df.out <- data.frame()
  
  # get nonprivate df
  nonpriv.tests <- load.inference.tests(model.id, "nonprivate")
    if (model.id >= 6) {
    true.theta <- apply(nonpriv.tests$Theta,c(2),mean)
  } else {
    true.theta <- nonpriv.tests$true
  }
  nonpriv.df <- make.inference.testsdf.row(nonpriv.tests, true.theta, model.id, 0)
  df.out <- rbindlist(list(df.out, nonpriv.df), use.names=TRUE, fill=TRUE)
  
  test.id.start <- 0
  for (method in methods) {
    for (dp.epsilon in dp.epsilons) {
      fname <- sprintf("obj/inference.tests/inference.tests%d%s-eps%g", model.id, method, dp.epsilon)
      if (file.exists(fname)) {
        df.out <- rbindlist(list(df.out, get.summary.inference.tests(model.id, method, dp.epsilon, test.id.start, true.theta)),
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
get.summary.inference.tests <- function(model.id, method, dp.epsilon, test.id.start, true.theta) {
  # load the tests for this model, method and epsilon value
  tests <- load.inference.tests(model.id, method, dp.epsilon)
  df.out <- data.frame("model.id"=numeric(), "method"=character(), "stat.name"=character(),
                   "eps"=numeric(), "delta"=numeric(),
                    "stat.value"=numeric(),  "AR"=numeric(),
                   "post.mean"=numeric(), "post.mean.med"=numeric(), "post.se"=numeric(), "post.sd"=numeric(),
                   "noise.draw"=numeric(), "noise.level"=character(), "KL"=numeric())
  if (is.null(true.theta)) {
    true.theta <- tests$true
  }
  for (t in 1:length(tests$private)) {
    row.df <- make.inference.testsdf.row(tests$private[[t]], true.theta, model.id, t+test.id.start)
    df.out <- rbindlist(list(df.out, row.df), use.names=TRUE, fill=TRUE)
  }
  
  return(df.out)
}

# make one row of test df
make.inference.testsdf.row <- function(x, true.theta, model.id, test.num) {
     row <- list()
   
    row$model.id <- model.id
    row$test.num <- test.num
    
    statnames <- names(x$stats)
    row$post.mean <- apply(x$Theta,c(2),mean)
    row$post.mean.med <- apply(apply(x$Theta,c(2,3),mean), 1, median)
    row$post.se <- getSE(x$Theta)
    row$post.sd <- apply(Theta.out, 2, sd)
    names(row$post.mean) <- statnames
    names(row$post.se) <- statnames
    row$stat.value <- x$stats
    row$true.param.value <- true.theta
    row$AR <- mean(x$AR)

  
    # get rid of 0's
    true.theta <- true.theta[true.theta != 0]
    row$KL <- computeKL(x$formula, true.theta, row$post.mean) 
    
    # non-private test inference run
    if (test.num == 0) {
      row$method <- "nonprivate"
    } 
    # private test inference run
    else {
      row$method <- x$noise$method
      row$eps <- sum(x$noise$dp.epsilon)
      row$delta <- sum(x$noise$dp.delta)
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

# estimate KL divergence btw paramter distrs
computeKL <- function(formula, true.theta, theta.other) {
  llr <- ergm.bridge.llr(formula,
                         from=true.theta, to=theta.other,
                         #control=control.ergm.bridge(nsteps=250, MCMC.burnin=1e5),
                         llronly=TRUE, verbose = FALSE)
  
  return(llr)
}


# extract non-private tests from raw test data
extract.nonprivate <- function(i) {
  tests <- load.inference.tests(i, "restr", 1.)
  #nonprivate <- load.inference.tests(i, "nonprivate", 1)
  nonprivate <- tests$nonprivate
  nonprivate$true <- tests$true
  fname = sprintf("obj/inference.tests/nonprivate%i", i)
  save(nonprivate, file=fname)
}

# Load dfs for inference tests
load.df.inference.tests <- function(i) {
  df.tests <- read.table(file = sprintf("obj/df/df.inference.tests%d.txt",i), sep = ",", header=TRUE)
}

# cut out cases with worst KL
truncate.KL <- function(df.graph) {
  for (eps in unique(df.graph$eps)) {
    for (method in unique(df.graph$method)) {
      if (method != 'nonprivate') {
        KL.vec <- sort(df.graph[df.graph$eps == eps & df.graph$method == method, 'KL'])
        df.graph <- df.graph[!(df.graph$eps == eps & df.graph$method == method
                               & ((df.graph$KL %in% head(KL.vec, 1)) | df.graph$KL %in% tail(KL.vec, 2))),]
      }
    }
  }
  return(df.graph)
}

# get rows of tests with median KLs for their method and eps level
get.median.KLs <- function(df.tests) {
  df.tests$KL <- round(abs(df.tests$KL),3)
  df.tests <- data.table(df.tests)
  df.meds <- df.tests[stat.name == 'edges', list(quantile(KL, prob=c(0.25))), by=list(method, eps)]
  for (i in 1:dim(df.meds)[1]) {
      row <- df.meds[i,]
      df.tests <- df.tests[eps != row$eps | method != row$method | KL == row$V1,]
  }
  return(df.tests)
}

visualize.inference.tests <- function(tests.id) {
  df.tests <- load.df.inference.tests(tests.id)
  df.tests$KL <- abs(df.tests$KL)
  df.tests <- truncate.KL(df.tests)
  df.tests$method = revalue(df.tests$method, c('smooth' = 'private local', 'restr' = 'restricted'))
  
  df.graph <- df.tests[df.tests$method != 'nonprivate' & df.tests$stat.name == 'edges', ]
  df.graph$eps <- as.factor(df.graph$eps)
  print(sprintf("KL of non-private is %g", mean(df.tests[df.tests$method == 'nonprivate', 'KL'])))
  plt.KL <- ggplot(data = df.graph, mapping=aes(x=eps, y=KL, color=method)) +
    geom_boxplot(outlier.alpha = 0.5, width = 0.5) +
    scale_y_log10("KL Divergence from Ground Truth (Log Scale)", limits=c(1,NA), breaks=c(1,10,100,1000,10000)) +
    scale_x_discrete(expression(epsilon)) +
    theme(axis.text=element_text(size=8), axis.title = element_text(size=10))
  
  ggsave(filename = sprintf('plots/inference/KLplot%d.png', sample.map.id(tests.id)), plot = plt.KL)
  # MAE plots
  
  df.tests$squared.error <- (df.tests$post.mean - df.tests$true.param.value)**2
  df.tests <- data.table(df.tests)
  rel.mae.df <- unique(df.tests[stat.value != 0, (mean(sqrt(squared.error)/(abs(true.param.value)))), by=list(method, eps, stat.name)])
  nonpriv <- rel.mae.df[method=='nonprivate', list(stat.name, V1)]
  View(nonpriv)
  rel.mae.df <- rel.mae.df[method != 'nonprivate',]
  View(rel.mae.df)
  rel.mae.df$stat.name <- sapply(as.character(rel.mae.df$stat.name), get.statname, USE.NAMES = FALSE)
  if (tests.id == 6) {
    rel.mae.df = rel.mae.df[method %in% c('restricted', 'rr')]
  }
  # Labeler for facet grid
  plt.MAE <- ggplot(data = rel.mae.df, mapping=aes(x=stat.name, y=V1, fill=factor(method))) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(~ eps, ncol = 2, labeller = label_bquote(cols = epsilon == .(eps))) +
    scale_fill_discrete('Method') +
    scale_x_discrete('Statistic') +
    theme(axis.text=element_text(size=8, angle=45, hjust=1), axis.title = element_text(size=10)) +
    scale_y_continuous('Relative MAE', breaks=pretty_breaks(n=10)) 
  ggsave(filename = sprintf('plots/inference/MAEplot%d.png', sample.map.id(tests.id)), plot = plt.MAE)
}