## Poisson hierarchical regression
## Takes a long time to run and I'm too lazy to optimize.
## Converges and creates very nice trace plots.
## Using 1 for all hyperparameters.
postpollution <- function(y, x, g, N=10000, burn=1000) {
    ## number of cities
    tn <- max(g)
    ## functions for the log posteriors in metropolis step
    lpost.beta <- function(y, x, alpha, beta, mu, tau2) {
        mu.i <- exp(alpha + beta*x)
        -sum(mu.i) + sum(y*log(mu.i)) - 1/(2*tau2) * (beta - mu)^2
    }
    lpost.alpha <- function(y, x, alpha, beta, sigma2) {
        mu.i <- exp(alpha + beta*x)
        -sum(mu.i) + sum(y*log(mu.i)) - 1/(2*sigma2) * alpha^2
    }

    ## just make all hyperpriors 1
    a <- b <- c <- d <- 1
    A <- 1
    ## initalize beta to all zeros?
    beta  <- rep(0, tn)
    alpha <- rep(0, tn)

    ## initialize mu, sigma, tau
    mu <- 0
    tau2 <- sigma2 <- 1

    ## Store
    MU <- TAU <- SIGMA <- rep(NA, N-burn)

    ## MCMC
    delta = 1.5 ## symmetry param
    for(i in 1:N) {
        ## city level parameters to update
        #########################
        #### METROPOLIS STEP ####
        #########################
        ## slow loop, could probably be made into vectorized components
        for (j in 1:tn) {
            x.j  <- x[g==j]
            y.j  <- y[g==j]
            ## city level means
            mu.j <- exp(alpha[j] + beta[j]*x.j)

            ## update alpha for each city
            alpha.star <- rnorm(1, alpha[j], sqrt(delta))
            log.a <- lpost.alpha(y.j, x.j, alpha.star, beta[j], sigma2)
            - lpost.alpha(y.j, x.j, alpha[j], beta[j], sigma2)
            if (log(runif(1)) < log.a) {
                alpha[j] <- alpha.star
            }

            ## update beta for each city
            beta.star <- rnorm(1, beta[j], sqrt(delta))
            log.b <- lpost.beta(y.j, x.j, alpha[j], beta.star, mu, tau2)
            - lpost.beta(y.j, x.j, alpha[j], beta[j], mu, tau2)
            if (log(runif(1)) < log.b) {
                beta[j] <- beta.star
            }
        }
        ####################
        #### GIBBS STEP ####
        ####################
        beta.h  <- mean(beta)
        alpha.h <- mean(alpha)

        ## update mu, tau, sigma
        mu <- rnorm(1, tn*beta.h/(1/A + tn), sqrt(1/(1/A + tn)))
        tau2 <- 1/rgamma(1, a + tn/2, b + sum((beta - beta.h)^2)/2 +
                         tn/(A*(1/A + tn)) * beta.h^2/2)
        sigma2 <- 1/rgamma(1, c + tn/2, d + sum(alpha^2)/2)

        ## store mu and tau
        if(i > burn) {
            MU[i-burn] <- mu
            TAU[i-burn] <- tau2
            SIGMA[i-burn] <- sigma2
        }

    }
    ## check acceptance probs
    return(list('MU'=MU, 'TAU'=TAU, 'SIGMA'=SIGMA))
}

dat <- read.csv(file="hw3_data.csv")
y <- dat$y
x <- dat$x
g <- dat$g

post <- postpollution(y, x, g, N=100, burn=10)
library('mcmcplots')
mcmcplot(matrix(c('mu'=post$MU, 'tau'=post$TAU, 'sigma'=post$SIGMA)))
