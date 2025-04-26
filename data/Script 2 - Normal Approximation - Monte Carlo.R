#_____________________________________
# Statistical Inference II - Script 2
#_____________________________________

setwd("G:/Il mio Drive/Didattica/PhD - Statistical Inference II")

#######################################

# Normal Approximation for the Bernoulli-Beta posterior:

# Let's consider results of a general election, where
# party A collected 150'000 votes out of the 240'000 total votes
sum_y <- 150000
n <- 240000

# Hyperparameters:
# Results of a previous general election showed that party A
# collected 49'000 votes out of 270'000 total votes:
a <- 49000
b <- 270000

Beta.Mode <- function(a, b){
  ifelse(a>1 & b>1, (a-1)/(a+b-2), NA)
}

# Posterior distribution:
a.post <- a + sum_y
b.post <- b + n - sum_y

post.mode <- Beta.Mode(a.post, b.post)
sigma_approx <- ((a.post-1)/(post.mode^2) + (b.post-1)/(1-post.mode)^2)^(-1)



curve(dbeta(x, a.post, b.post), xlim=c(.3,.4))
curve(dnorm(x, post.mode, sqrt(sigma_approx)), add=T, col=2)


# Iterative way:
h.grid <- seq(100, 105, by = 0.001)

res <- matrix(NA, ncol=4, nrow=length(h.grid))
colnames(res) <- c("HPD1", "HPD2", "level", "h")

for(i in 1:length(h.grid)){
  translated <- function(x, a, b) dbeta(x,a, b)-h.grid[i]

  hpd1<-uniroot(translated,c(.34,post.mode-.0005),a.post,b.post)$root
  hpd2<-uniroot(translated,c(post.mode+.0005,.36),a.post,b.post)$root

  I <- integrate(dbeta, lower=hpd1, upper=hpd2, shape1=a.post, shape2=b.post)$value

  iter <- i
  res[i,] <- c(hpd1, hpd2, I, h.grid[i])
  if(I <= 0.95) break
}


res[iter,]

HPD <- res[iter-1,-c(3,4)]; HPD
qnorm(c(.025,.975), post.mode, sqrt(sigma_approx))





# MONTE CARLO

# A biostatistician wants to compare two treatments:
# 1 - a new vaccine for COVID-19
# 2 - a placebo

# He/She collected a sample of 200 patients and divided
# them into 2 groups of equal size by randomization.

n1 <- 100
n2 <- 100

# The first group will receive the new vaccine, whereas
# the second group will receive the placebo.

# Then, the biostatistician registers how many patients
# died in each group (Y1 and Y2).

# Y1 ~ Bern(theta1)
# Y2 ~ Bern(theta2)

# theta represents the probability of death in a group.

# The biostatistician chooses to use non-informative priors:

a1 <- 1
b1 <- 1

a2 <- 1
b2 <- 1

# Once the priors have been selected, He/She looks
# at the data:

s1 <- 20
s2 <- 80

# Update of the hyperparameters:
a1.post <- a1 + s1
b1.post <- b1 + n1 - s1

a2.post <- a2 + s2
b2.post <- b2 + n2 - s2


B <- 100000


#Posterior densities:
curve(dbeta(x ,a1.post,b1.post), 0,1, lty=1, ylim=c(0,10),
      xlab=expression(theta),ylab=expression(pi),lwd=2)
curve(dbeta(x ,a2.post,b2.post),add=T,lty=3,lwd=2)
legend(0.4,10,c("Posterior 1","Posterior 2"),lty=c(1,3),cex=0.8,lwd=2,bty="n")

# These posteriors seem to indicate that theta1
# and theta2 are concentrated around different
# values.

set.seed(42)
theta1 <- rbeta(B, a1.post, b1.post)
theta2 <- rbeta(B, a2.post, b2.post)


eta <- log(theta1/(1-theta1)*(1-theta2)/theta2)

mean(eta)
var(eta)

CS <- quantile(eta, p=c(.025,.975)); CS


# Histogram and kernel estimate:
hist(eta, prob=T, xlab=expression(eta),
     main="Posterior distribution of log OR")
lines(density(eta), col="red")

# The biostatistician concludes that, with a probability
# of 95%, the log OR is included in (-3.453943, -2.069884),
# meaning that the odds of death in the "Vaccine" group
# is lower than the odds of death in the "Placebo" group.




