#Rejections sampling
postsample <- function(y, N, sigma) {
    ## function for drawing from half normal (using truncated normal)
    hnorm <- function(N, sigma) {
        replicate(N, qnorm(0.5*(1 + runif(1)), 0, sqrt(pi/2)*sigma))
    }

    ## use data/parameters from assignment if no arguments passed
    if(missing(y)) y <- c(20.100306, 2.272066, 3.796734, 2.265275, 3.480183)
    if(missing(sigma)) sigma=0.5
    if(missing(N)) N=1000
    mle <- 1/mean(y)
    ## since sample only 5, don't bother putting on log scale.
    lik.max <- prod(pexp(y, mle))
    post <- NULL
    while(length(post) < N) {
        u <- runif(1)
        beta <- hnorm(1, sigma)
        lik <- prod(pexp(y, beta))
        if( u <= lik/lik.max ) post <- append(post, beta)
    }
    return(post)
}
