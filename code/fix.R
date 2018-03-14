i <- 1
df.samples <- data.table(read.table(file=sprintf("obj/df/df.samples%d.txt",i), sep = ",", header=TRUE))
df.more.samples <- data.table(read.table(file=sprintf("obj/df/df.aggr.samples%d.txt",i), sep = ",", header=TRUE))
df.more.samples <- df.more.samples[type == 'restricted']
df.samples <- rbindlist(list(df.samples, df.more.samples),use.names=TRUE, fill=TRUE)
#View(df.samples)
#vec.eps <- c(rep("min",4), rep("median",4), rep("max",4), rep("conservative",4), rep(NA,4))
#df.samples <- df.samples[order(n, eps, dp.k)]
#df.samples[,"dp.klevel"] <- rep(vec.eps, 30)

write.table(df.samples, file = sprintf("obj/df/df.samples%d.txt",i), sep = ",", col.names = colnames(df.samples))
