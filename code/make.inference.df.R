# get commands
args <- (commandArgs(TRUE))
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

#setwd("/n/home06/akgoldberg/dp-ergm-thesis/code")
setwd("/Users/alexandergoldberg/Documents/Harvard/Senior Year/Thesis/project/code")
source('libraries.R')

# if (!file.exists(sprintf("obj/inference.tests/nonprivate%d", i))) {
#   extract.nonprivate(i)
# }
df.tests <- make.inference.df(i)
write.table(df.tests, file = sprintf("obj/df/df.inference.tests_new%d.txt",i), sep = ",", col.names = colnames(df.tests))
