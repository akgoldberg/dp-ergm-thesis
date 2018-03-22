# get commands
args <- (commandArgs(TRUE))
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
#setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')

df.tests <- make.inference.df(i)
write.table(df.tests, file = sprintf("df.inference.tests%d.txt",i), sep = ",", col.names = colnames(df.tests))
