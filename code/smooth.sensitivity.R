################################################################################
##                                Smooth Sensitivity                          ##
################################################################################

# Add noise to suff stats
make.private.smooth <- function (formula, dp.epsilon, dp.delta, max.deg = NULL, max.sp = NULL) {
  # get original network
  y <- ergm.getnetwork(formula)

  if (is.null(max.deg))   max.deg <- network.maxdegree(y)
  if (is.null(max.sp)) max.sp <- network.maxsharedpartners(y)
 
  model <- ergm.getmodel(formula, y)
  # draw laplace noise
  noise <- draw.lap.noise.smooth(model$terms, dp.epsilon, dp.delta, max.sp, max.deg)
  noise$dp.epsilon <- dp.epsilon
  noise$dp.delta <- dp.delta
  noise$method <- "smooth"
  return(list("formula" = formula, "noise" = noise))
}

# draw Laplace noise with smooth sensitivity for multiple terms
draw.lap.noise.smooth <- function (terms, dp.epsilonTot, dp.deltaTot, max.sp, max.deg) {
    # output vector of level of laplace noise to add
   noise.level <- rep(NA,  length(unlist(lapply(terms, function(x) {x$coef.names}))))

  # split epsilon evenly between terms if it is given as scalar
  if (length(dp.epsilonTot == 1)) {
    dp.epsilon <- rep(dp.epsilonTot/length(terms), length(terms))
  } else {
    if (length(dp.epsilonTot) != length(terms)) {
      print("Epsilon must be either scalar or specified for all terms.")
      return()
    }
    dp.epsilon <- dp.epsilonTot
  }
  
  # split delta evenly between terms if it is given as scalar
  if (length(dp.deltaTot == 1)) {
    dp.delta <- rep(dp.deltaTot/length(terms), length(terms))
  } else {
    if (length(dp.deltaTot) != length(terms)) {
      print("Delta must be either scalar or specified for all terms.")
      return()
    }
    dp.delta <- dp.deltaTot
  }

  # calculate amount of noise to add for each term
  i <- 1
  for (t in 1L:length(terms)) {
    name <- terms[[t]]$name
    param <- tail(terms[[t]]$inputs, 1)
    num.terms <- length(terms[[t]]$coef.names)
    print(num.terms)

    if (name == 'edges') {
    noise.level[[i]] <- 1/dp.epsilon[t]
    }
    if (name == 'altkstar') {
    noise.level[[i]] <- (2*(1./param))/dp.epsilon[t]
    }
    if (name == 'gwesp') {
    noise.level[[i]] <- draw.lap.noise.smooth.term(name, param, dp.epsilon[t], dp.delta[t], max.sp=max.sp)$level
    }
    if (name == 'gwdsp') {
    noise.level[[i]] <- draw.lap.noise.smooth.term(name, param, dp.epsilon[t], dp.delta[t], max.deg=max.deg)$level
    }

    # label-dependant terms
    if (name == 'nodematch') {
        noise.level[i:(i+num.terms)] <- 1/dp.epsilon[t]
    } 
    if (name == 'nodefactor') {
        noise.level[i:(i+num.terms)] <- 2/dp.epsilon[t]
    }

     i <- i+num.terms
  }
  print(noise.level)
   # draw Laplace noise to use
  noise.draw <- rlaplace(n = length(noise.level), scale = noise.level)
  return(list("level" = noise.level, "draw" = noise.draw))
}

# draw Laplace noise with smooth sensitivity for one term
draw.lap.noise.smooth.term <- function(term, param, dp.epsilon, dp.delta, max.sp=NULL, max.deg=NULL) {

    # compute epsilon and delta for algorithm
    dp.epsilon <- dp.epsilon/2.
    # split delta between terms requiring smooth sensitivity
    dp.delta <- 2*(dp.delta)/(exp(dp.epsilon))

    a <- log(1./dp.delta)/dp.epsilon

    # calculate amount of noise to add for each term
    if (term == 'gwesp') {
        ls.x <- 1./param + 2*max.sp
        gs.ls <- 2
        noise.level <- max(1e-6,(ls.x + rlaplace(n=1,s=gs.ls/dp.epsilon) + a*gs.ls)/dp.epsilon)
    }
    if (term == 'gwdsp') {
        ls.x <- 2*max.deg
        gs.ls <- 2
        noise.level <- max(1e-6,(ls.x + rlaplace(n=1,s=gs.ls/dp.epsilon) + a*gs.ls)/dp.epsilon)
    }
    # draw Laplace noise to use
    noise.draw <- rlaplace(n = 1, scale = noise.level)
    return(list("level" = noise.level, "draw" = noise.draw))
}

