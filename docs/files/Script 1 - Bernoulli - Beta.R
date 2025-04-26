#_____________________________________
# Statistical Inference II - Script 1
#_____________________________________

setwd("G:/Il mio Drive/Didattica/PhD - Statistical Inference II/My_Scripts")

library(readr)
cardiac <- read.csv("cardiac.csv", header=T, sep=";")
str(cardiac)

y <- cardiac$Chd




par(mfrow=c(2,2))
curve(dbeta(x, 1, 1), main="Scenario A", xlab=expression(theta), ylab=expression(pi(theta)))
curve(dbeta(x, 10, 10), main="Scenario B", xlab=expression(theta), ylab=expression(pi(theta)))
curve(dbeta(x, 10, 5), main="Scenario C", xlab=expression(theta), ylab=expression(pi(theta)))
curve(dbeta(x, 5, 10), main="Scenario D", xlab=expression(theta), ylab=expression(pi(theta)))
par(mfrow=c(1,1))


L_Binom <- function(y, theta){
  s <- sum(y)
  n <- length(y)
  L <- theta^s * (1-theta)^(n-s)
  return(L)
}

L_Binom(y, theta=0.001)


# Scenario A:
a1 <- 1; b1 <- 1

# Prior distribution:
curve(dbeta(x, a1, b1), main="Prior A", ylim=c(0,1.4),
      xlab=expression(theta), lwd=2)
abline(v=a1/(a1+b1), col="red", lty="dashed")

# Likelihood function:
curve(L_Binom(y, theta=x), main="Likelihood",
      xlab=expression(theta), lwd=2)
abline(v=mean(y), col="blue", lty="dashed")


# Posterior distribution:
n <- length(y)
a1.post <- a1 + sum(y)
b1.post <- b1 + n - sum(y)

curve(dbeta(x, a1.post, b1.post), main="Scenario A",
      xlab=expression(theta), lwd=2)
curve(dbeta(x, a1, b1), lty="dashed", add=T, lwd=2)
abline(v=a1/(a1+b1), col="red", lty="dashed", lwd=2)
abline(v=mean(y), col="green", lty="dotted", lwd=2)
abline(v=a1.post/(a1.post+b1.post), col="blue", lty="dashed", lwd=2)
legend(.8,8, c("Prior", "Posterior", "Prior Mean", "MLE", "Post. Mean"), lwd=2,
       col=c("black", "black", "red", "green", "blue"), lty=c(2,1,2,3,2), bty="n")


# The chosen prior distribution does not favor 
# any value of theta.
# Thus, we selected a non-informative prior.
# The posterior distribution is centered 
# around the MLE.


# Scenario B:
a2 <- 10; b2 <- 10

# Prior distribution:
curve(dbeta(x, a2, b2), main="Prior B",
      xlab=expression(theta), lwd=2)
abline(v=a2/(a2+b2), col="red", lty="dashed")

# Likelihood function:
curve(L_Binom(y, theta=x), main="Likelihood",
      xlab=expression(theta), lwd=2)
abline(v=mean(y), col="blue", lty="dashed")


# Posterior distribution:
n <- length(y)
a2.post <- a2 + sum(y)
b2.post <- b2 + n - sum(y)

curve(dbeta(x, a2.post, b2.post), main="Scenario B",
      xlab=expression(theta), lwd=2)
curve(dbeta(x, a2, b2), lty="dashed", add=T, lwd=2)
abline(v=a2/(a2+b2), col="red", lty="dashed", lwd=2)
abline(v=mean(y), col="green", lty="dotted", lwd=2)
abline(v=a2.post/(a2.post+b2.post), col="blue", lty="dashed", lwd=2)
legend(.8,8, c("Prior", "Posterior", "Prior Mean", "MLE", "Post. Mean"), lwd=2,
       col=c("black", "black", "red", "green", "blue"), lty=c(2,1,2,3,2), bty="n")

# Now we are still considering a symmetric prior
# with prior mean = 0.5, but it does not
# give same density/probability to each 
# value/interval of theta.
# Now the posterior mean is a weighted average 
# between prior mean and MLE.





# Scenario C:
a3 <- 10; b3 <- 5

# Posterior distribution:
n <- length(y)
a3.post <- a3 + sum(y)
b3.post <- b3 + n - sum(y)

curve(dbeta(x, a3.post, b3.post), main="Scenario C",
      xlab=expression(theta), lwd=2)
curve(dbeta(x, a3, b3), lty="dashed", add=T, lwd=2)
abline(v=a3/(a3+b3), col="red", lty="dashed", lwd=2)
abline(v=mean(y), col="green", lty="dotted", lwd=2)
abline(v=a3.post/(a3.post+b3.post), col="blue", lty="dashed", lwd=2)
legend(.8,8, c("Prior", "Posterior", "Prior Mean", "MLE", "Post. Mean"), lwd=2,
       col=c("black", "black", "red", "green", "blue"), lty=c(2,1,2,3,2), bty="n")





# Scenario D:
a4 <- 5; b4 <- 10

# Posterior distribution:
n <- length(y)
a4.post <- a4 + sum(y)
b4.post <- b4 + n - sum(y)

curve(dbeta(x, a4.post, b4.post), main="Scenario D",
      xlab=expression(theta), lwd=2)
curve(dbeta(x, a4, b4), lty="dashed", add=T, lwd=2)
abline(v=a4/(a4+b4), col="red", lty="dashed", lwd=2)
abline(v=mean(y), col="green", lty="dotted", lwd=2)
abline(v=a4.post/(a4.post+b4.post), col="blue", lty="dashed", lwd=2)
legend(.8,8, c("Prior", "Posterior", "Prior Mean", "MLE", "Post. Mean"), lwd=2,
       col=c("black", "black", "red", "green", "blue"), lty=c(2,1,2,3,2), bty="n")






# We may summarize the (prior and) posterior distribution
# by some measures:
Beta.Exp <- function(a, b) a/(a+b)
Beta.Var <- function(a, b)(a*b)/((a+b)^2*(a+b+1))
Beta.Mode <- function(a, b){
  ifelse(a>1 & b>1, (a-1)/(a+b-2), NA)
}
Beta.Median <- function(a, b) qbeta(.5, a, b)


measures <- matrix(NA, ncol=6, nrow=4)
rownames(measures) <- c("A", "B", "C", "D")
colnames(measures) <- c("a.post", "b.post", "Exp", "Mode", "Median", "Var")

measures[,1] <- c(a1.post,a2.post,a3.post,a4.post)
measures[,2] <- c(b1.post,b2.post,b3.post,b4.post)
for(i in 1:4){
  a <- measures[i,1]
  b <- measures[i,2]
  measures[i,3] <- Beta.Exp(a, b)
  measures[i,4] <- Beta.Mode(a, b)
  measures[i,5] <- Beta.Median(a, b)
  measures[i,6] <- Beta.Var(a, b)
}

round(measures,4)

# In this simple scenario, even considering 
# different priors, we do not have posteriors 
# leading to very different point estimates.

# This is because the empirical evidence 
# (i.e., our data) is much
# stronger than our prior belief. 

# Let us consider an additional scenario where 
# we strongly believe that most people 
# do not have cardiovascular disease.

a5 <- 100
b5 <- 1000

# Prior distribution:
curve(dbeta(x, a5, b5), main="Prior E",
      xlab=expression(theta))
abline(v=a5/(a5+b5), col="red", lty="dashed")

# Posterior distribution:
n <- length(y)
a5.post <- a5 + sum(y)
b5.post <- b5 + n - sum(y)

curve(dbeta(x, a5.post, b5.post), main="Scenario E",
      xlab=expression(theta), lwd=2)
curve(dbeta(x, a5, b5), lty="dashed", add=T, lwd=2)
abline(v=a5/(a5+b5), col="red", lty="dashed", lwd=2)
abline(v=mean(y), col="green", lty="dotted", lwd=2)
abline(v=a5.post/(a5.post+b5.post), col="blue", lty="dashed", lwd=2)
legend(.8,40, c("Prior", "Posterior", "Prior Mean", "MLE", "Post. Mean"), lwd=2,
       col=c("black", "black", "red", "green", "blue"), lty=c(2,1,2,3,2), bty="n")


measures <- matrix(NA, ncol=6, nrow=5)
rownames(measures) <- c("A", "B", "C", "D", "E")
colnames(measures) <- c("a.post", "b.post", "Exp", "Mode", "Median", "Var")

measures[,1] <- c(a1.post,a2.post,a3.post,a4.post,a5.post)
measures[,2] <- c(b1.post,b2.post,b3.post,b4.post,b5.post)
for(i in 1:5){
  a <- measures[i,1]
  b <- measures[i,2]
  measures[i,3] <- Beta.Exp(a, b)
  measures[i,4] <- Beta.Mode(a, b)
  measures[i,5] <- Beta.Median(a, b)
  measures[i,6] <- Beta.Var(a, b)
}

round(measures,4)




# Interval Estimates:

# Credible Sets (CS).
# The CS are very easy to compute. Indeed, we only have to compute
# quantiles of the posterior distribution.


# Let us consider scenario D:
a <- a4.post
b <- b4.post

# Scenario D:
curve(dbeta(x, a, b), main="Scenario D",
      xlab=expression(theta), lwd=2)

# CS of level 0.95
CS <- qbeta(c(0.025,0.975), a, b);CS
abline(v=CS,lty="dashed", col="blue")
legend(0.8,8,c("Posterior", "Credible Set"),
       lty=c(1,2), col=c("black", "blue"), bty="n")


# Highest posterior density (HPD):

# Let's start by considering a random value for h
h <- 2
curve(dbeta(x, a, b), 
      ylab=expression(paste(pi,"(",theta,"|x)")),
      xlab=expression(theta))
abline(h=h,lty=2)

# We need to find the values of theta such that p(theta|y) = h

translated <- function(x, a, b) dbeta(x,a, b)-h 
# We decrease the density curve by h.
# Now, the values we are looking for are such that translated = 0

curve(translated(x, a, b), add=T,lty=3,lwd=2)
abline(h=0,lty=3)

hpd1 <- uniroot(translated,c(.2, .4),a,b)$root;hpd1
hpd2 <- uniroot(translated,c(.45, .5),a,b)$root;hpd2
integrate(dbeta, lower=hpd1, upper=hpd2, shape1=a, shape2=b)

# We have a probability equal to 0.9152 --> Decrease h

# Iterative way:
h.grid <- seq(1, 2, by = 0.01)   

res <- matrix(NA, ncol=4, nrow=length(h.grid))
colnames(res) <- c("HPD1", "HPD2", "level", "h")

for(i in 1:length(h.grid)){        
  translated <- function(x, a, b) dbeta(x,a, b)-h.grid[i] 
  
  hpd1<-uniroot(translated,c(.2,.4),a,b)$root
  hpd2<-uniroot(translated,c(.43,.55),a,b)$root
  
  I <- integrate(dbeta, lower=hpd1, upper=hpd2, shape1=a, shape2=b)$value
  
  iter <- i
  res[i,] <- c(hpd1, hpd2, I, h.grid[i])
  if(I <= 0.95) break
}

res[1:iter,]
h.grid[iter]

HPD <- res[iter-1,-c(3,4)]; HPD
CS


# Graphically:
curve(dbeta(x, a, b), main="Scenario D",
      xlab=expression(theta), lwd=2)

abline(v=CS,lty="dashed", col="blue")
abline(v=HPD,lty="dashed", col="red")

legend(0.8,8,c("Posterior", "Credible Set", "HPD"),
       lty=c(1,2,2), col=c("black", "blue", "red"), bty="n")


# HPD and CS are very similar, because the posterior is 
# unimodal and almost symmetrical.



# Hypothesis testing:

# Let consider H0: theta <= 0.3 vs H1: theta > 0.3
pbeta(.3, a4, b4)
pbeta(.3, a4.post, b4.post)

ODDS.prior <- pbeta(.3, a4, b4)/(1-pbeta(.3, a4, b4)); ODDS.prior
ODDS.post <- pbeta(.3, a4.post, b4.post)/(1-pbeta(.3, a4.post, b4.post)); ODDS.post

Bayes_Factor <- ODDS.post/ODDS.prior; Bayes_Factor

# Data do not support H0: by adding data, our confidence on
# H0 strongly decreases



# HOMEWORK: consider a sample of size n = 10 composed of
# Bernoulli trials, for which we observed s = 3 successes.
# Consider that a priori theta ~ Beta(2, 8)

# 1. Which is the posterior distribution of theta?
# 2. Compute the MLE and the prior and posterior means. 
#    Are these results reasonable?
# 3. Compute and compare the CS and the HPD.
# 4. Let H0 : theta in (0.3, 0.5)
#        H1 : theta in (0, 0.3) U (0.5, 1)
#    Which hypothesis can you support?


