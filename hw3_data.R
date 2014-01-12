## Generate data for HW3 hierarchical model

set.seed(100)
N <- 100
nt <- 60 * 10
sigma <- 2
alpha <- rnorm(N, 0, sigma)
mu <- log(1.0004)
tau <- .0002
b <- rnorm(N, mu, tau)
x <- rnorm(N * nt, 0, 9)
g <- gl(N, nt)
y <- unlist(mapply(function(x, alpha, b) {
        log.mu <- alpha + b * x
        y <- rpois(nt, exp(log.mu))
}, split(x, g), alpha, b, SIMPLIFY = FALSE))

data <- data.frame(y = y, group = g, x = x)

write.csv(data, file = "hw3_data.csv", row.names = FALSE)

