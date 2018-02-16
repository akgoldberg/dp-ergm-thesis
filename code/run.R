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
sigma.epsilon <- diag(c(0.05, 0.05, 0.05))
sigma.prior <- diag(c(30, 30, 30))
m.prior <- c(0,0,0)

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


nonprivate.model <- bergm(nw ~ edges + altkstar(0.5, fixed=TRUE) 
                  + gwesp(0.5, fixed=TRUE), 
                  burn.in=30000,
                  main.iters=10000,
                  aux.iters = 20000,
                  m.prior = 0, 
                  sigma.prior = sigma.prior, 
                  nchains = 1, 
                  sigma.epsilon = sigma.epsilon)

bergm.output(dp.model, lag.max=500)

dp.standard.model <- standard.bergm(nw ~ edges + altkstar(0.5, fixed=TRUE) 
                                    + gwesp(0.5, fixed=TRUE), 
                                    dp.epsilon = dp.epsilon,
                                    dp.k = dp.k,
                                    privacy.type='edge',
                                    dp.delta = NULL,
                                    burn.in=30000,
                                    main.iters=10000,
                                    aux.iters=20000, 
                                    m.prior = m.prior, 
                                    sigma.prior = sigma.prior, 
                                    nchains = 1, 
                                    sigma.epsilon = sigma.epsilon)


  