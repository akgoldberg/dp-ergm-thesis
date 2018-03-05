# get commands

 args <- (commandArgs(TRUE))
 for(k in 1:length(args)){
      eval(parse(text=args[[k]]))
}

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
#setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')

# print(i)
samples <- load.samples(i)
df.samples <- run.noise.tests(samples, c(0.1, 0.5, 1.), dp.delta=1e-6, stop=1000) 
write.table(df.samples, file = sprintf("df.samples%d.txt",i), sep = ",", col.names = colnames(df.samples))


# if (i==1) {
#     sigma.epsilon=diag(c(0.0005, 0.00025))
#     inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE),
#                                             dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
#                                             non.private=FALSE, num.tests=25)
# }
# if (i==2) {
#     sigma.epsilon=diag(c(0.0005, 0.00025))
#     inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE),
#                                             dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
#                                             non.private=FALSE, num.tests=25)
# }
# if (i==3) {
#     sigma.epsilon = diag(c(0.0005, 0.00025, 0.00001))
#     inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE)+gwdsp(0.5,fixed=TRUE),
#                                             dp.epsilon, method=method, sigma.epsilon=sigma.epsilon,
#                                             non.private=FALSE, num.tests=25)
# } 
    
# save(inference.tests, file=sprintf("inference.tests%d%s-eps%g", i, method, dp.epsilon))
