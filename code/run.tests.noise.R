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



