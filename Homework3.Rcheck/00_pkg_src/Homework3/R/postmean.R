# importance sampling
## Don't think this is quite right but waited until last second to do bonus
## and can't find anything wrong
postmean <- function(beta, lower, upper) {
    ## importance sampling function
    importance <- function(beta, sigma) {
        n <- length(beta)
        ## sample from trial density
        hnorm <- function(N, sigma) {
            replicate(N, qnorm(0.5*(1 + runif(1)), 0, sqrt(pi/2)*sigma))
        }
        ## log of trial density
        ldhnorm <- function(x, sigma) {
            d <- ifelse(x<0, -Inf, log(2) + dnorm(x, 0, sqrt(pi/2)*sigma,
                                                  log=TRUE))
            return(d)
        }
        ## importance function
        lw <- function(x) beta - ldhnorm(beta, sigma)
        # sample from trial density with 1000 monte carlo samples
        U <- hnorm(1000, sigma)

        ## list of importance function values
        lps <- lw(U)
        lps <- lps - max(lps)

        ## Estimate of posterior mean
        I <- mean(exp(lps)*U) / mean(exp(lps))
        return(I)
    }
    sigma <- seq(lower, upper, length=50)

    postmeans <- sapply(sigma, function(x) importance(beta, x))
    plot(sigma, postmeans, xlab="Sigma", ylab="Posterior mean estimate", 
         type="l", main="Importance sampling")
}
