# get commands
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- as.numeric(sub("-","",myargument))

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
source('libraries.R')
i <- myargument
print(i)
samples <- load.samples(i)
df.samples <- run.noise.tests(samples, c(0.1, 0.5, 1.), dp.delta=1e-6, stop=1000) 
write.table(df.samples, file = sprintf("df.samples%d.txt",i), sep = ",", col.names = colnames(df.samples))
