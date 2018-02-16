################################################################################
##                                Smooth Sensitivity                          ##
################################################################################

# draw Laplace noise under restricted sensitivty
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
        noise.level <- (ls.x + rlaplace(n=1,s=gs.ls/dp.epsilon) + a*gs.ls)/dp.epsilon
    }
    if (term == 'gwdsp') {
        ls.x <- 2*max.deg
        gs.ls <- 2
        noise.level <- (ls.x + rlaplace(n=1,s=gs.ls/dp.epsilon) + a*gs.ls)/dp.epsilon
    }
    # draw Laplace noise to use
    noise.draw <- rlaplace(n = 1, s = noise.level)
    return(list("level" = noise.level, "draw" = noise.draw))
}

