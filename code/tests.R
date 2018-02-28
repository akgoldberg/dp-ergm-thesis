run.inference.tests <- function(formula, n, true.theta, dp.k, dp.epsilon = 1.0, num.tests = 10) {
  print(sprintf("Running tests with n=%d, k=%d, eps=%g", n, dp.k, dp.epsilon))
  
  tic("Runtime")
  
 # print("Non-Private")
 # nonprivate.out <- bergm.orig(formula,
  #                   burn.in=5000, main.iters = 7500, aux.iters = 0.1*choose(n,2),
   #                  sigma.epsilon = diag(c(0.001, 0.001, 0.0001)),
    #                 print.out=2500, nchains = 1)
  
  
  # run tests
  outs <- list()
  for (t in 1:num.tests) {
    print(sprintf("Test: %d", t))
    nw.private <- make.private.restr(formula, dp.epsilon, dp.k, privacy.type="edge")
    private.out <- bergm.modified.private(nw.private$formula, nw.private$noise,
                                          burn.in=7500, main.iters = 7500, aux.iters = 0.2*choose(n,2),
                                          sigma.epsilon = diag(c(0.001, 0.001, 0.0001)),
                                          print.out=1000, nchains = 1)
    outs[[t]] <- private.out
  }
  toc()
  
  return(list("private"=outs, "nonprivate"=nonprivate.out, "true"=true.theta))
}

tests0.higheps <- run.inference.tests(nw ~ edges + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE), n, true.theta, dp.k, dp.epsilon = 1.0)
tests0.loweps <- run.inference.tests(nw ~ edges + gwesp(0.5, fixed=TRUE) + gwdsp(0.5, fixed=TRUE), n, true.theta, dp.k, dp.epsilon = 0.5)