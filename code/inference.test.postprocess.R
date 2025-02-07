# make a data frame with multiple methods and epsilons from inference tests
make.inference.df <- function(model.id, methods=c("nonpriv", "smooth", "restr", "rr"), dp.epsilons=c(1, 2, 3, 4)) {
  df.out <- data.frame()
  
  # get nonprivate df
  #nonpriv.tests <- load.inference.tests(model.id, "nonprivate")
  #  if (model.id >= 6) {
  #  true.theta <- apply(nonpriv.tests$Theta,c(2),mean)
  #} else {
  #  true.theta <- nonpriv.tests$true
  #}
  #nonpriv.df <- make.inference.testsdf.row(nonpriv.tests, true.theta, model.id, 0)
  #df.out <- rbindlist(list(df.out, nonpriv.df), use.names=TRUE, fill=TRUE)
  true.theta <- NULL
  test.id.start <- 0
  for (method in methods) {
    for (dp.epsilon in dp.epsilons) {
      print(sprintf("Method: %s, Epsilon: %d", method, dp.epsilon))
      fname <- sprintf("obj/inference.tests.new/inference.tests%d%s-eps%g", model.id, method, dp.epsilon)
      if (file.exists(fname) | (method == 'nonpriv' && dp.epsilon==1)) {
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
  tests <- load.inference.tests.new(model.id, method, dp.epsilon)
  df.out <- data.frame("model.id"=numeric(), "method"=character(), "stat.name"=character(),
                   "eps"=numeric(), "delta"=numeric(),
                    "stat.value"=numeric(),  "AR"=numeric(),
                   "post.mean"=numeric(), "post.mean.med"=numeric(), "post.se"=numeric(), "post.sd"=numeric(),
                   "noise.draw"=numeric(), "noise.level"=character())
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
    row$post.sd <- apply(x$Theta, 2, sd)
    names(row$post.mean) <- statnames
    names(row$post.se) <- statnames
    names(row$post.mean.med) <- statnames
    names(row$post.sd) <- statnames
    row$stat.value <- x$stats
    true.theta <- true.theta[true.theta != 0]
    row$true.param.value <- true.theta
    row$AR <- mean(x$AR)

  
    # get rid of 0's
    true.theta <- true.theta[true.theta != 0]
    #row$KL <- computeKL(x$formula, true.theta, row$post.mean.med) 
    
    # non-private test inference run
    if (test.num == 0) {
      row$method <- "nonprivate"
    } 
    # private test inference run
    else {
      row$method <- x$noise$method
      if (is.null(row$method)) {
        row$method <- 'nonpriv'
      }
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
load.df.inference.tests <- function(i, new=TRUE) {
  if (new) {
    read.table(file = sprintf("obj/df/df.inference.tests_new%d.txt",i), sep = ",", header=TRUE)
  } else {
    read.table(file = sprintf("obj/df/df.inference.tests%d.txt",i), sep = ",", header=TRUE)
  }
}

# cut out cases with worst KL
truncate.KL <- function(df.graph) {
  for (eps in unique(df.graph$eps)) {
    for (method in unique(df.graph$method)) {
      if (method != 'nonprivate') {
        KL.vec <- sort(df.graph[df.graph$eps == eps & df.graph$method == method, 'KL'])
        df.graph <- df.graph[!(df.graph$eps == eps & df.graph$method == method
                               & ((df.graph$KL %in% head(KL.vec, 0)) | df.graph$KL %in% tail(KL.vec, 5))),]
      }
    }
  }
  return(df.graph)
}

# get rows of tests with median KLs for their method and eps level
get.median.KLs <- function(df.tests, pr=0.5) {
  df.tests$KL <- round(abs(df.tests$KL),3)
  df.tests <- data.table(df.tests)
  df.meds <- df.tests[stat.name == 'edges', list(quantile(KL, prob=c(pr))), by=list(method, eps)]
  for (i in 1:dim(df.meds)[1]) {
      row <- df.meds[i,]
      df.tests <- df.tests[eps != row$eps | method != row$method | KL == row$V1,]
  }
  return(df.tests)
}

visualize.inference.tests <- function(tests.id) {
  options(scipen=1e5)
  if (tests.id >= 7) {
    legend.position = 'bottom'
  }  else {
    legend.position = 'none'
  }
  df.tests <- load.df.inference.tests(tests.id)
  df.tests$KL <- abs(df.tests$KL/10.)
  df.tests <- truncate.KL(df.tests)
  df.tests$method = revalue(df.tests$method, c('smooth' = 'private local', 'restr' = 'restricted'))
  
  df.graph <- df.tests[df.tests$method != 'nonprivate' & df.tests$stat.name == 'edges', ]
  df.graph$eps <- as.factor(df.graph$eps)
  print(sprintf("KL of non-private is %g", mean(df.tests[df.tests$method == 'nonprivate', 'KL'])))
  plt.KL <- ggplot(data = df.graph, mapping=aes(x=eps, y=KL, color=method)) +
    geom_boxplot(outlier.alpha = 0.5, width = 0.5) +
    scale_y_log10("KL Divergence (Log Scale)", breaks=c(0.1,1,10,100,1000,10000)) +
    scale_x_discrete(expression(epsilon)) +
    scale_color_discrete("Method:") +
    theme(axis.text=element_text(size=14), axis.title.y = element_text(size=16), axis.title.x = element_text(size=32), legend.position=legend.position)

  ggsave(filename = sprintf('plots/inference/KLplot%d.png', sample.map.id(tests.id)), plot = plt.KL)
  # MAE plots
  make_stat_plot(df.tests, 'RMSE', tests.id, legend.pos = legend.position)
  make_stat_plot(df.tests, 'MAE', tests.id, legend.pos = legend.position)
  make_stat_plot(df.tests, 'MSE', tests.id, legend.pos = legend.position)
}


make_stat_plot <- function(df.tests, stat.use, tests.id, use.rel=TRUE, legend.pos) {
  df.tests <- data.table(df.tests)
  
   if(tests.id == 7) {
    #df.tests <- df.tests[method != 'private local', ]
    use.rel <- FALSE
  } 

  if (stat.use == "MSE") {
      df.tests$squared.error <- (df.tests$post.mean.med - df.tests$true.param.value)**2
      use.rel <- FALSE
      rel.df <- unique(df.tests[stat.value != 0, mean(squared.error), by=list(method, eps, stat.name)])  
  }
  if (stat.use == "RMSE") {
      df.tests$squared.error <- (df.tests$post.mean.med - df.tests$true.param.value)**2
      if (use.rel) {
        rel.df <- unique(df.tests[stat.value != 0, (sqrt(mean(squared.error/true.param.value**2))), by=list(method, eps, stat.name)])
      } else {
        rel.df <- unique(df.tests[stat.value != 0, (sqrt(mean(squared.error))), by=list(method, eps, stat.name)])
      }
  }
  if (stat.use == "MAE") {
      df.tests$squared.error <- (df.tests$post.mean.med - df.tests$true.param.value)**2
      if (use.rel) {
        rel.df <- unique(df.tests[stat.value != 0, (mean(sqrt(squared.error/true.param.value**2))), by=list(method, eps, stat.name)])
      } else {
        rel.df <- unique(df.tests[stat.value != 0, (mean(sqrt(squared.error))), by=list(method, eps, stat.name)])
      }
  }
  
  rel.df <- rel.df[method != 'nonprivate',]
  rel.df$stat.name <- sapply(as.character(rel.df$stat.name), get.statname, USE.NAMES = FALSE)
  
  if(use.rel) {
    label.str <- sprintf('Relative %s', stat.use)
  } else {
    label.str <- sprintf('%s', stat.use)
  }
  
  plt <- ggplot(data = rel.df, mapping=aes(x=stat.name, y=V1, fill=factor(method))) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(~ eps, ncol = 2, labeller = label_bquote(cols = epsilon == .(eps)), scales="free") +
    scale_fill_discrete("Method:") +
    scale_x_discrete('Statistic') +
    theme(axis.text=element_text(size=14, angle=45, hjust=1, vjust=1), axis.title = element_text(size=16), legend.position=legend.pos) +
    scale_y_continuous(label.str, breaks=pretty_breaks(n=10)) 
  ggsave(filename = sprintf('plots/inference/%splot%d.png', stat.use, sample.map.id(tests.id)), plot = plt)
}

make.all.plots <- function() {
  for (i in c(1,3,5,7)) {
    visualize.inference.tests(i)
  }
}

view.ARs <- function(tests.id) {
  df.tests <- data.table(load.df.inference.tests(tests.id))
  df.tests$AR <- 100*df.tests$AR
  df.AR <- df.tests[stat.name == 'edges', list("min"=min(AR), "25p"=quantile(AR, .25),"median"=median(AR), "75p"=quantile(AR, .75), "max"=max(AR), "mean"=mean(AR)), by=list(method, eps)]
  View(df.AR)
}

trunc.points <- function(df.tests) {
  df.param.min <- df.tests[,list('param.est'=min(param.est)), by=list(method,eps,stat.name)]
  df.param.max <- df.tests[,list('param.est'=max(param.est)), by=list(method,eps,stat.name)]
  df.tests <- df.tests[!df.param.min, on=list(method, eps, stat.name, param.est)]
  df.tests <- df.tests[!df.param.max, on=list(method, eps, stat.name, param.est)]
  return(df.tests)
}

# PLOT BOXPLOTS OF PARAMS
plot.tests.params <- function(tests.id, legend.position="bottom", nonpriv.r = 0.025, dp.eps=c(1,3)) {
  df.tests <- data.table(load.df.inference.tests(tests.id))
  df.tests$stat.name <- sapply(df.tests$stat.name, get.statname)
  df.tests$param.est <- df.tests$post.mean.med
  df.tests$method = revalue(df.tests$method, c('restr' = 'restricted')) #'smooth' = 'private local', 
  df.tests <- trunc.points(df.tests[eps %in% dp.eps | method=='nonpriv',])
  df.tests.nonpriv <- trunc.points(trunc.points(df.tests[method != 'nonpriv',]))
  if (tests.id < 6) {
    df.tests.true <- df.tests[,list("val"=mean(true.param.value)), by=list(stat.name, method, eps)]
  }
  else {
    df.tests.true <- df.tests[method=="nonpriv", list("val"=mean(post.mean)), by=list(stat.name)]
    df.tests.true <- merge(df.tests, df.tests.true, by='stat.name', all=TRUE)
  } 
  plt.priv <- ggplot(df.tests.nonpriv, aes(x=method,y=param.est, fill=method)) +
    facet_wrap(stat.name~eps, ncol=length(dp.eps), scales='free', labeller = label_bquote(rows = .(stat.name)*",  "*epsilon == .(eps))) +
    geom_boxplot(outlier.alpha=0.5, alpha=0.8) +
    geom_hline(data=df.tests.true[method!='nonpriv',], aes(yintercept = val), alpha=0.8, linetype=2,color='black', size=1.1) +
    scale_y_continuous("Parameter Estimate", breaks=pretty_breaks(n=8)) +
    scale_fill_discrete("Method:") +
    scale_x_discrete("",breaks=NULL) +
    theme_grey() +
    theme(strip.text=element_text(size=8, margin=margin(0.01, 0, 0.01, 0, "cm")),
          axis.text.y=element_text(size=8),
          axis.text.x =element_text(size=0),
          axis.title.y = element_text(size=8),
          legend.position=legend.position,
          legend.margin=margin(0,0,0,0, "cm"))
  plt.legend <- get_legend(plt.priv)
  plt.priv <- plt.priv + theme(legend.position='none')
  
  data_dummylow <- df.tests.true[method == 'nonpriv', list('param.est'=(1-nonpriv.r)*val), by=list(method, stat.name)]
  data_dummyhigh <- df.tests.true[method == 'nonpriv', list('param.est'=(1+nonpriv.r)*val), by=list(method, stat.name)]
  plt.nonpriv <- ggplot(df.tests[method=='nonpriv',], aes(x=method,y=param.est)) +
    facet_wrap(~stat.name, ncol=1, scales='free', labeller = label_bquote(rows = .(stat.name)*",  Non-Private")) +
    geom_boxplot(outlier.alpha=0.5, alpha=0.8) +
    geom_hline(data=df.tests.true[method=='nonpriv',], aes(yintercept = val), alpha=0.8, linetype=2,color='black', size=1.1) +
    scale_y_continuous("Parameter Estimate", breaks=pretty_breaks(n=8)) +
    geom_blank(data=data_dummylow) +
    geom_blank(data=data_dummyhigh) +
    scale_x_discrete("",breaks=NULL) +
    theme_grey() +
    theme(strip.text=element_text(size=8, margin=margin(0.01, 0, 0.01, 0, "cm")),
          axis.text.y=element_text(size=8),
          axis.text.x =element_text(size=0),
          axis.title.y = element_text(size=0),
          plot.margin=margin(5.5,5.5,5.5,-5.5, "pt")) 
  
  plts <- plot_grid(plt.priv, plt.nonpriv, rel_widths = c(2, 1))
  plts
  plts
  #plts.legend <- plot_grid(plt.legend, plts, nrow=2, rel_heights = c(0.25, length(unique(df.tests.nonpriv$stat.name))))
  #ggsave(filename=sprintf('plots/inference/paramplot%d.png', sample.map.id(tests.id)), plot = plts)
  #ggsave(filename='plots/inference/legend.png', plot = plt.legend, height = 0.25, width = 4)
}

get.table.params <- function(tests.id, my.methods = c('nonpriv 0', 'restr 1', 'restr 2')) {
  df.tests <- data.table(load.df.inference.tests(tests.id))
  df.tests$stat.name <- sapply(df.tests$stat.name, get.statname)
  df.tests$param.est <- df.tests$post.mean.med
  df.tests$method <- paste(df.tests$method, df.tests$eps)
  df.tests <- trunc.points(df.tests[method %in% my.methods,])
  df.tests.nonpriv <- trunc.points(trunc.points(df.tests[method == 'nonpriv 0',]))
  df.tests[, eps:=NULL]
  df.tests.nonpriv[, eps:=NULL]
  if (tests.id < 6) {
    df.tests.true <- df.tests[,list("val"=mean(true.param.value)), by=list(stat.name, method)]
  }
  else {
    df.tests.true <- df.tests[method=="nonpriv 0", list("val"=mean(post.mean)), by=list(stat.name)]
    df.tests.true <- merge(df.tests, df.tests.true, by='stat.name', all=TRUE)
  } 
  df.summary <- df.tests.true[, list('avg.param'=mean(param.est), 'mse'=mean((param.est-val)**2), 'mae'=mean(abs(param.est-val))),
                                by=list(stat.name, method)]
  
  View(df.summary)
  get.df <- function (my.method) {
    if (my.method == 'nonpriv 0') {
      return(df.summary[method == my.method,-c('method','stat.name','mae', 'mse'),with=FALSE])
    } else {
      return(df.summary[method == my.method,-c('method','stat.name'),with=FALSE])
    }
  }
  
  df.out <- cbind.data.frame(unique(df.summary$stat.name),lapply(my.methods,get.df))
  View(df.out)
  print(xtable(df.out), include.rownames=FALSE)
}

get.table.med <- function(tests.id, my.methods = c('nonpriv 0', 'restr 1', 'restr 2')) {
  df.tests <- data.table(load.df.inference.tests(tests.id))
  df.tests$stat.name <- sapply(df.tests$stat.name, get.statname)
  df.tests$param.est <- df.tests$post.mean.med
  df.tests$method <- paste(df.tests$method, df.tests$eps)
  #df.tests <- trunc.points(df.tests[method %in% my.methods,])
  df.tests.nonpriv <- df.tests[method == 'nonpriv 0',]
  df.tests[, eps:=NULL]
  df.tests.nonpriv[, eps:=NULL]
  
  df.tests.true <- df.tests[method=="nonpriv 0", list("val"=mean(param.est), "val.sd"=mean(post.sd)), by=list(stat.name)]
  df.tests.true <- merge(df.tests, df.tests.true, by='stat.name', all=TRUE)
  df.tests.true$err <- df.tests.true[, abs(val - param.est)/val]
  df.totalerrs <- df.tests.true[, list('totalerr'=sum(err)), by=list(test.num, method)]
  df.mederrs <- df.totalerrs[, list('totalerr'=min(totalerr)), by=list(method)]
  df.totalerrs <- df.totalerrs[df.mederrs, on=c('totalerr')]
  
  df.tests.true <- df.tests.true[df.totalerrs, on=c('test.num')]
  View(df.tests.true)

  get.df <- function (my.method) {
    return(df.tests.true[method==my.method,c('param.est', 'post.sd'), with=FALSE])
  }
  
  df.out <- cbind.data.frame(unique(df.tests.true$stat.name),lapply(my.methods,get.df))
  View(df.out)
  print(xtable(df.out), include.rownames=FALSE)
}


plot.labeltests.params.nonpriv <- function(tests.id) {
  df.tests <- data.table(load.df.inference.tests(tests.id))[method != 'nonprivate',]
  df.tests$stat.name <- sapply(df.tests$stat.name, get.statname)
  df.tests$param.est <- df.tests$post.mean.med
  df.tests$method = revalue(df.tests$method, c('smooth' = 'private local', 'restr' = 'restricted'))

  plt.nonpriv <- ggplot(df.tests[method=='nonpriv',], aes(x=method,y=param.est)) +
    facet_wrap(~stat.name, ncol=1, scales='free', labeller = label_bquote(rows = .(stat.name)*",  Non-Private")) +
    geom_boxplot(outlier.alpha=0.5, alpha=0.8) +
    scale_y_continuous("Parameter Estimate", breaks=pretty_breaks(n=8)) +
    scale_x_discrete("",breaks=NULL) +
    theme_grey() +
    theme(strip.text=element_text(size=8, margin=margin(0.01, 0, 0.01, 0, "cm")),
          axis.text.y=element_text(size=8),
          axis.text.x =element_text(size=0),
          axis.title.y = element_text(size=0),
          plot.margin=margin(5.5,5.5,5.5,-5.5, "pt")) 
  plt.nonpriv
}


