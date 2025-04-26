#_____________________________________
# Statistical Inference II - Script 3
#_____________________________________

######################################
######## GIBBS SAMPLING ##############
######################################

# True value of parameters:
mu.true <- 10
sigma2.true <- 5
n <- 100

# Data generation:
set.seed(42)
y <- rnorm(n, mu.true, sqrt(sigma2.true))
ybar <- mean(y)

# Hyperparameters:
mu0 <- 0
sigma2.0 <- 1000

alpha <- 10
beta <- 10

# Prior expectation:
beta/(alpha-1)
# Prior variance:
(beta^2)/((alpha-1)*(alpha-2))

B <- 1000
mu.chain <- numeric(B)
sigma2.chain <- numeric(B)


# Initialization:

# We initialize each element by drawing a value
# from the corresponding prior:

set.seed(42)
mu.chain[1] <- rnorm(1, mu0, sqrt(sigma2.0))
sigma2.chain[1] <- 1/rgamma(1, alpha, rate = beta)

set.seed(42)
for(b in 2:B){
  # Draw mu from the F.C.
  mu.chain[b] <-
    rnorm(1,
          (n*ybar*sigma2.0 + sigma2.chain[b-1])/(n*sigma2.0+sigma2.chain[b-1]),
          sqrt((sigma2.chain[b-1]*sigma2.0)/(n*sigma2.0+sigma2.chain[b-1])))
  
  # Draw sigma2 from the F.C.
  sigma2.chain[b] <-
    1/rgamma(1, alpha + n/2,
             rate = beta + .5*sum((y-mu.chain[b])^2))
  
}


plot(mu.chain, sigma2.chain, pch=20)
points(mu.true, sigma2.true, col="red", pch=20)

# Traceplots:
par(mfrow=c(2,1))
plot(mu.chain[-1], pch=20, type="l")
plot(sigma2.chain[-1], pch=20, type="l")
par(mfrow=c(1,1))


mu.means <- numeric(B)
sigma2.means <- numeric(B)

mu.means[1] <- mu.chain[1]
sigma2.means[1] <- sigma2.chain[1]

for(b in 2:B){
  mu.means[b] <- mean(mu.chain[2:b])
  sigma2.means[b] <- mean(sigma2.chain[2:b])
}

par(mfrow=c(2,1))
plot(mu.means[-1], pch=20, type="l")
abline(h=mu.true, col="red")
plot(sigma2.means[-1], pch=20, type="l")
abline(h=sigma2.true, col="red")
par(mfrow=c(1,1))

par(mfrow=c(2,1))
acf(mu.chain)
acf(sigma2.chain)
par(mfrow=c(1,1))



# Removing the warm-up:
warm_perc <- .5

mu.new <- mu.chain[round(B*warm_perc+1):B]
sigma2.new <- sigma2.chain[round(B*warm_perc+1):B]

par(mfrow=c(2,1))
plot(mu.new, pch=20, type="l")
plot(sigma2.new, pch=20, type="l")
par(mfrow=c(1,1))


# Computing estimates:

mean(mu.new)
quantile(mu.new, probs = c(.025, .975))

mean(sigma2.new)
quantile(sigma2.new, probs = c(.025, .975))

mean(mu.new > 10)

# Histogram and kernel estimate:
par(mfrow=c(2,1))
hist(mu.new, prob=T, xlab=expression(mu),
     main="Posterior distribution of mu")
lines(density(mu.new), col="red")
hist(sigma2.new, prob=T, xlab=expression(sigma2),
     main="Posterior distribution of sigma2")
lines(density(sigma2.new), col="red")
par(mfrow=c(1,1))

plot(mu.new, sigma2.new, pch = 20)





###################################
####### Defining a function #######
###################################

normal_GS <- function(y, B = 5000, 
                      mu0 = 0, sigma2.0 = 1000, 
                      alpha = 10, beta = 10, 
                      warm_perc = .5, seed=42){
  
  mu.chain <- numeric(B)
  sigma2.chain <- numeric(B)
  ybar <- mean(y)
  n <- length(y)
  
  # Initialization:
  set.seed(seed)
  mu.chain[1] <- rnorm(1, mu0, sqrt(sigma2.0))
  sigma2.chain[1] <- 1/rgamma(1, alpha, rate = beta)
  
  for(b in 2:B){
    # Draw mu from the F.C.
    mu.chain[b] <-
      rnorm(1,
            (n*ybar*sigma2.0 + sigma2.chain[b-1])/(n*sigma2.0+sigma2.chain[b-1]),
            sqrt((sigma2.chain[b-1]*sigma2.0)/(n*sigma2.0+sigma2.chain[b-1])))
    
    # Draw sigma2 from the F.C.
    sigma2.chain[b] <-
      1/rgamma(1, alpha + n/2,
               rate = beta + .5*sum((y-mu.chain[b])^2))
  }
  
  mu.new <- mu.chain[round(B*warm_perc+1):B]
  sigma2.new <- sigma2.chain[round(B*warm_perc+1):B]
  
  return(cbind(mu.chain = mu.new, sigma2.chain = sigma2.new))
}


library(faraway)
data(gala)
str(gala)
help(gala)

y <- gala$Species
gala_GS <- normal_GS(y, B = 10000)

str(gala_GS)
head(gala_GS)


plot(gala_GS, pch = 20)

hist(gala_GS[,1], prob = T, 
     main="Posterior distribution of mu")
lines(density(gala_GS[,1]), col="red")

hist(gala_GS[,2], prob = T, 
     main="Posterior distribution of sigma2")
lines(density(gala_GS[,2]), col="red")

colMeans(gala_GS)
t(apply(gala_GS, 2, function(x) quantile(x, probs=c(.025, .975))))

# Computing probabilities:
mean(gala_GS[,1] > 70)
mean(gala_GS[,1] > 70 & gala_GS[,2] < 5500)



##################################
###### LINEAR REGRESSION #########

rm(list=ls())

# True values:
beta <- c(-2,5,3)
sigma2 <- 6

# Generating data:
n <- 100
set.seed(42)
X <- matrix(rnorm(2*n, 0, 50), ncol=2)
X <- cbind(rep(1,n), X)

y <- as.numeric(X%*%beta + rnorm(n, 0, sqrt(sigma2)))

summ <- summary(lm(y~X[,2]+X[,3])); summ



B <- 5000
beta.chain <- matrix(NA, ncol=3, nrow=B)
sigma2.chain <- numeric(B)

beta.chain[1,] <- rep(0,3)
sigma2.chain[1] <- 1

beta0 <- rep(0,3)
Sigma0 <- 100*diag(3)

a0 <- 10
b0 <- 10

# Prior Expectation for sigma2
b0/(a0-1)
# Prior Variance for sigma2
(b0^2)/((a0-1)*(a0-2))

library(MASS)
for(b in 2:B){
  
  Sigma.n <- solve(solve(Sigma0) + (t(X)%*%X)/sigma2.chain[b-1])
  beta.n <- Sigma.n %*% ((solve(Sigma0)%*%beta0) + (t(X)%*%y)/sigma2.chain[b-1])
  
  beta.chain[b,] <- mvrnorm(n=1, mu=beta.n, Sigma=Sigma.n)
  
  sigma2.chain[b] <- 
    1/rgamma(1, a0 + .5*n,
             rate = b0 + 
               0.5*(t(y-X%*%beta.chain[b,])%*%(y-X%*%beta.chain[b,])))
}




# Traceplots:
par(mfrow=c(2,2))
plot(beta.chain[,1], pch=20, type="l");abline(h=beta[1], col="red")
plot(beta.chain[,2], pch=20, type="l");abline(h=beta[2], col="red")
plot(beta.chain[,3], pch=20, type="l");abline(h=beta[3], col="red")
plot(sigma2.chain, pch=20, type="l");abline(h=sigma2, col="red")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(beta.chain[-1,1], pch=20, type="l");abline(h=beta[1], col="red")
plot(beta.chain[-1,2], pch=20, type="l");abline(h=beta[2], col="red")
plot(beta.chain[-1,3], pch=20, type="l");abline(h=beta[3], col="red")
plot(sigma2.chain[-1], pch=20, type="l");abline(h=sigma2, col="red")
par(mfrow=c(1,1))


beta_mean <- matrix(NA, ncol = 3, nrow = nrow(beta.chain))
sigma2.means <- numeric(length(sigma2.chain))

beta_mean[1,] <- beta.chain[1,]
sigma2.means[1] <- sigma2.chain[1]

for(b in 2:nrow(beta.chain)){
  beta_mean[b,1] <- mean(beta.chain[2:b,1])
  beta_mean[b,2] <- mean(beta.chain[2:b,2])
  beta_mean[b,3] <- mean(beta.chain[2:b,3])
  sigma2.means[b] <- mean(sigma2.chain[2:b])
}

par(mfrow=c(2,2))
plot(beta_mean[-1,1], pch=20, type="l")
plot(beta_mean[-1,2], pch=20, type="l")
plot(beta_mean[-1,3], pch=20, type="l")
plot(sigma2.means[-1], pch=20, type="l")
par(mfrow=c(1,1))

acf(beta_mean[,1])
acf(beta_mean[,2])
acf(beta_mean[,3])
acf(sigma2.means)

# Removing the warmu-up:
warm_perc <- .5

beta.new <- beta.chain[round(B*warm_perc+1):B,]
sigma2.new <- sigma2.chain[round(B*warm_perc+1):B]

# Computing estimates:
colMeans(beta.new)
mean(sigma2.new)


# Histogram and kernel estimate:
par(mfrow=c(2,2))
hist(beta.new[,1], prob=T, xlab=expression(beta0),
     main="Posterior distribution of beta0")
lines(density(beta.new[,1]), col="red")
##############################################
hist(beta.new[,2], prob=T, xlab=expression(beta1),
     main="Posterior distribution of beta1")
lines(density(beta.new[,2]), col="red")
##############################################
hist(beta.new[,3], prob=T, xlab=expression(beta2),
     main="Posterior distribution of beta2")
lines(density(beta.new[,3]), col="red")
##############################################
hist(sigma2.new, prob=T, xlab=expression(sigma2),
     main="Posterior distribution of sigma2")
lines(density(sigma2.new), col="red")
par(mfrow=c(1,1))





round(var(beta.new), 5)
round((2.483^2)*summ$cov.unscaled,5)





###################################
####### Defining a function #######
###################################
rm(list=ls())

LM_GS <- function(y, X, B = 5000, 
                  beta0 = rep(0, ncol(X)), 
                  Sigma0 = diag(ncol(X)), 
                  a0 = 10, b0 = 10, warm_perc = .5, seed=42){
  
  beta.chain <- matrix(NA, ncol=ncol(X), nrow=B)
  sigma2.chain <- numeric(B)
  n <- length(y)
  
  # Initialization:
  beta.chain[1,] <- rep(0, ncol(X))
  sigma2.chain[1] <- 1
  
  library(MASS)
  for(b in 2:B){
    Sigma.n <- solve(solve(Sigma0) + (t(X)%*%X)/sigma2.chain[b-1])
    beta.n <- Sigma.n %*% ((solve(Sigma0)%*%beta0) + (t(X)%*%y)/sigma2.chain[b-1])
    
    beta.chain[b,] <- mvrnorm(n=1, mu=beta.n, Sigma=Sigma.n)
    
    sigma2.chain[b] <- 
      1/rgamma(1, a0 + .5*n,
               rate = b0 + 
                 0.5*(t(y-X%*%beta.chain[b,])%*%(y-X%*%beta.chain[b,])))
  }
  
  beta.new <- beta.chain[round(B*warm_perc+1):B,]
  sigma2.new <- sigma2.chain[round(B*warm_perc+1):B]
  
  return(list(beta = beta.new, sigma2 = sigma2.new))
}

# Definition of y and X:
y <- gala$Species
X <- model.matrix(Species ~ ., data = gala[,-2])

# Fitting the model:
lm_gala <- LM_GS(y, X)
str(lm_gala)

# Extracting the elements of the chain:
betas <- lm_gala$beta
colnames(betas) <- colnames(X)
sigma2 <- lm_gala$sigma2

colMeans(betas)
t(apply(betas, 2, function(x) quantile(x, probs=c(.025, .975))))

mean(sigma2)


cov(betas)
