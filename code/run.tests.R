# get commands

 args <- (commandArgs(TRUE))
 for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
}

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
#setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')
#i <- myargument
# print(i)
# samples <- load.samples(i)
# df.samples <- run.noise.tests(samples, c(0.1, 0.5, 1.), dp.delta=1e-6, stop=1000) 
# write.table(df.samples, file = sprintf("df.samples%d.txt",i), sep = ",", col.names = colnames(df.samples))

#dp.epsilon=1.
inference.tests <- run.inference.tests(i, 300, ~edges+gwesp(0.5, fixed=TRUE), dp.epsilon, sigma.epsilon=diag(c(0.025, 0.01)))
save(inference.tests, file=sprintf("inference.tests%d-eps%g", i, dp.epsilon))