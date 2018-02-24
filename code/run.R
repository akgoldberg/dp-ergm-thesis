# load code
setwd("~/Documents/Harvard/Senior Year/Thesis/project/code")
source("libraries.R", chdir=TRUE)
source("inference.R", chdir=TRUE)

# load data
data(package = 'ergm')
data(faux.mesa.high)

nw <- faux.mesa.high
dp.k <- 15
dp.epsilon <- 10

# initialize proposal distribution to be higher for theta for gwesp?
sigma.epsilon <- diag(c(0.1, 0.1))
sigma.prior <- diag(c(30, 30))
m.prior <- c(0,0)

dp.model <- dp.bergm(nw ~ edges + altkstar(0.5, fixed=TRUE) 
                      + gwesp(0.5, fixed=TRUE), 
                      dp.epsilon = dp.epsilon,
                      dp.k = dp.k,
                      privacy.type='edge',
                      dp.delta = NULL,
                      burn.in=30000,
                      main.iters=10000,
                      aux.iters=25000, 
                      m.prior = c(0,0,0), 
                      sigma.prior = sigma.prior, 
                      nchains = 1, 
                      sigma.epsilon = sigma.epsilon)


tic("Bergm training")
nonprivate.model <- bergm(nw ~ edges
                  + gwesp(0.5, fixed=TRUE), 
                  burn.in=5000,
                  main.iters=7500,
                  aux.iters = 25000,
                  m.prior = m.prior, 
                  sigma.prior = sigma.prior, 
                  nchains = 1, 
                  sigma.epsilon = sigma.epsilon)
toc()

# fit 300 edge network
dp.epsilon <- 1.
dp.k <- round(1.1*max(samples5.300$deg))
nw.private <- make.private(nw ~ edges + gwesp(0.5, fixed=TRUE) 
                               + gwdsp(0.5, fixed=TRUE),
                                dp.epsilon, dp.k)
fit.private <- bergm.orig.private(nw.private$formula,
                   nw.private$noise,
                   burn.in = 2500,
                   main.iters = 5000,
                   aux.iters = 0.25*choose(300,2), 
                   nchains = 1,
                   sigma.epsilon = diag(c(0.001, 0.001, 0.0001)),
                   print.out = 100)
bergm.output(dp.model, lag.max=100)