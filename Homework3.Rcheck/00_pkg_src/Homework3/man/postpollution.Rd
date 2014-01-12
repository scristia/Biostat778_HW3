\name{postpollution}
\alias{postpollution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Hierarchical gibbs sampling with metropolis step
%%  ~~function to do ... ~~
}
\description{
    hybrid gibbs
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
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
function (y, x, g, N = 10000, burn = 1000) 
{
    tn <- max(g)
    lpost.beta <- function(y, x, alpha, beta, mu, tau2) {
        mu.i <- exp(alpha + beta * x)
        -sum(mu.i) + sum(y * log(mu.i)) - 1/(2 * tau2) * (beta - 
            mu)^2
    }
    lpost.alpha <- function(y, x, alpha, beta, sigma2) {
        mu.i <- exp(alpha + beta * x)
        -sum(mu.i) + sum(y * log(mu.i)) - 1/(2 * sigma2) * alpha^2
    }
    a <- b <- c <- d <- 1
    A <- 1
    beta <- rep(0, tn)
    alpha <- rep(0, tn)
    mu <- 0
    tau2 <- sigma2 <- 1
    MU <- TAU <- SIGMA <- rep(NA, N - burn)
    delta = 1.5
    for (i in 1:N) {
        for (j in 1:tn) {
            x.j <- x[g == j]
            y.j <- y[g == j]
            mu.j <- exp(alpha[j] + beta[j] * x.j)
            alpha.star <- rnorm(1, alpha[j], sqrt(delta))
            log.a <- lpost.alpha(y.j, x.j, alpha.star, beta[j], 
                sigma2)
            -lpost.alpha(y.j, x.j, alpha[j], beta[j], sigma2)
            if (log(runif(1)) < log.a) {
                alpha[j] <- alpha.star
            }
            beta.star <- rnorm(1, beta[j], sqrt(delta))
            log.b <- lpost.beta(y.j, x.j, alpha[j], beta.star, 
                mu, tau2)
            -lpost.beta(y.j, x.j, alpha[j], beta[j], mu, tau2)
            if (log(runif(1)) < log.b) {
                beta[j] <- beta.star
            }
        }
        beta.h <- mean(beta)
        alpha.h <- mean(alpha)
        mu <- rnorm(1, tn * beta.h/(1/A + tn), sqrt(1/(1/A + 
            tn)))
        tau2 <- 1/rgamma(1, a + tn/2, b + sum((beta - beta.h)^2)/2 + 
            tn/(A * (1/A + tn)) * beta.h^2/2)
        sigma2 <- 1/rgamma(1, c + tn/2, d + sum(alpha^2)/2)
        if (i > burn) {
            MU[i - burn] <- mu
            TAU[i - burn] <- tau2
            SIGMA[i - burn] <- sigma2
        }
    }
    return(list(MU = MU, TAU = TAU, SIGMA = SIGMA))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
