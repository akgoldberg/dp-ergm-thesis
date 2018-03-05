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
    sigma.epsilon = diag(c(0.0005, 0.00025, 0.00025, 0.00001))
    inference.tests <- run.inference.tests(i, 300, ~edges+altkstar(0.5, fixed=TRUE)+gwesp(0.5, fixed=TRUE)+gwdsp(0.5,fixed=TRUE),
                                            dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
                                            non.private=TRUE, num.tests=25)
} 
    
save(inference.tests, file=sprintf("inference.tests%d%s-eps%g", i, method, dp.epsilon))