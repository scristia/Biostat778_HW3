\name{postmean}
\alias{postmean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    importance sampling
%%  ~~function to do ... ~~
}
\description{
    does importance sampling
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
\author{
    Stephen Cristiano
%%  ~~who you are~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (beta, lower, upper) 
{
    importance <- function(beta, sigma) {
        n <- length(beta)
        hnorm <- function(N, sigma) {
            replicate(N, qnorm(0.5 * (1 + runif(1)), 0, sqrt(pi/2) * 
                sigma))
        }
        ldhnorm <- function(x, sigma) {
            d <- ifelse(x < 0, -Inf, log(2) + dnorm(x, 0, sqrt(pi/2) * 
                sigma, log = TRUE))
            return(d)
        }
        lw <- function(x) beta - ldhnorm(beta, sigma)
        U <- hnorm(1000, sigma)
        lps <- lw(U)
        lps <- lps - max(lps)
        I <- mean(exp(lps) * U)/mean(exp(lps))
        return(I)
    }
    sigma <- seq(lower, upper, length = 50)
    postmeans <- sapply(sigma, function(x) importance(beta, x))
    plot(sigma, postmeans, xlab = "Sigma", ylab = "Posterior mean estimate", 
        type = "l", main = "Importance sampling")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
