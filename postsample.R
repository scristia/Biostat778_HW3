#Rejections sampling
postsample <- function(y, N, sigma) {
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

## function for drawing from half normal
hnorm <- function(N, sigma) {
    replicate(N, qnorm(0.5*(1 + runif(1)), 0, sigma))
}


y <- c(20.100306, 2.272066, 3.796734, 2.265275, 3.480183)
sigma <- 0.5
N <- 10000

post <- postsample(y, N, sigma)
hist(post, breaks=100, col="gray", border="gray", freq=FALSE)
