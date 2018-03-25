# get commands
args <- (commandArgs(TRUE))
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
#setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')

samples <- load.samples(i)
df.samples <- run.noise.node.tests(samples, stop=1000) 
write.table(df.samples, file = sprintf("df.samples%d_node_trunc.txt",i), sep = ",", col.names = colnames(df.samples))



