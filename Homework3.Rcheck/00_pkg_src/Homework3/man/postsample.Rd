\name{postsample}
\alias{postsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Rejection sampling
%%  ~~function to do ... ~~
}
\description{
    Rejection sampling
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\author{
    Stephen Cristiano
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, N, sigma) 
{
    hnorm <- function(N, sigma) {
        replicate(N, qnorm(0.5 * (1 + runif(1)), 0, sqrt(pi/2)/sigma))
    }
    if (missing(y)) 
        y <- c(20.100306, 2.272066, 3.796734, 2.265275, 3.480183)
    if (missing(sigma)) 
        sigma = 0.5
    if (missing(N)) 
        N = 1000
    mle <- 1/mean(y)
    lik.max <- prod(pexp(y, mle))
    post <- NULL
    while (length(post) < N) {
        u <- runif(1)
        beta <- hnorm(1, sigma)
        lik <- prod(pexp(y, beta))
        if (u <= lik/lik.max) 
            post <- append(post, beta)
    }
    return(post)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
