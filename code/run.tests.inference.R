 args <- (commandArgs(TRUE))
 for(k in 1:length(args)){
      eval(parse(text=args[[k]]))
}

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
#setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')

if (i==1) {
    sigma.epsilon=diag(c(0.0005, 0.00025))
    inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE),
                                            dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                            non.private=FALSE, num.tests=25)
}
if (i==2) {
    sigma.epsilon=diag(c(0.0005, 0.00025))
    inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE),
                                            dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                            non.private=FALSE, num.tests=25)
}
if (i==3) {
    sigma.epsilon = diag(c(0.0005, 0.00025, 0.00001))
    inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE)+gwdsp(0.5,fixed=TRUE),
                                            dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                            non.private=FALSE, num.tests=25)
} 
if (i==5) {
    sigma.epsilon = diag(c(0.0001, 0.00005, 0.00005, 0.000005))
    inference.tests <- run.inference.tests(i, 300, ~edges+altkstar(0.5, fixed=TRUE)+gwesp(0.5, fixed=TRUE)+gwdsp(0.5,fixed=TRUE),
                                            dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                            non.private=FALSE, num.tests=25)
}

if (i==6) {
  form.rhs <- formula(~edges + nodefactor("Race") + nodefactor("Sex")
    + nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE) + altkstar(1.0,fixed=TRUE)
    + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE))
  
  sigma.epsilon = diag(c(0.00005, rep(0.00001, 11), 0.00001, 0.00001, 0.000001))
  #sigma.epsilon = NULL
  if (method != 'rr') {
    dp.epsilon <- c(rep(dp.epsilon/10., 5), rep(dp.epsilon/6., 3))
  }
  inference.tests <- run.inference.tests(i, 205, form.rhs,
                                         dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                         non.private=TRUE, num.tests=25, parallel = TRUE,
                                         burn.in=5000, main.iters=5000)
}
    
save(inference.tests, file=sprintf("inference.tests%d%s-eps%g", i, method, dp.epsilon))