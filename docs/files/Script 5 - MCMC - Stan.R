#_____________________________________
# Statistical Inference II - Script 5
#_____________________________________

#################################
######## HMC with Stan ##########
#################################

# https://mc-stan.org/users/interfaces/
# https://mc-stan.org/users/interfaces/rstan

# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

rm(list=ls())
setwd("G:/Il mio Drive/Didattica/PhD - Statistical Inference II/My_Scripts")

library(tidyverse)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
library(loo)
library(car)

LM_Model <- rstan::stan_model(file="LM_Model.stan")

data("Prestige")
help("Prestige")

Prestige2 <- Prestige
Prestige2 <- Prestige2[complete.cases(Prestige2),-5]
Prestige2$income <- log(Prestige2$income)


y <- Prestige2$prestige
mod <- lm(prestige ~ ., data=Prestige2)

X <- model.matrix(mod)
n <- nrow(X)
K <- ncol(X)

data.stan <- list(
  y = y,
  X = X,
  n = n,
  K = K,
  a0 = 10,
  b0 = 10,
  beta0 = rep(0, K),
  s2_0 = rep(10, K)
) 

n.iter <- 10000
nchain <- 2

LM_stan <- rstan::sampling(
  object = LM_Model,
  data = data.stan,
  warmup = 0.5*n.iter, iter = n.iter,
  thin=1, chains = nchain
  #, pars=c("beta", "sigma2")
  )

summ <- summary(LM_stan, pars = c("beta", "sigma2"))$summary
rownames(summ)[1:K] <- colnames(X)
round(summ, 3)

rstan::traceplot(LM_stan, pars = c('beta', 'sigma2'), inc_warmup = TRUE)
rstan::traceplot(LM_stan, pars = c('beta', 'sigma2'), inc_warmup = FALSE)


post_chain <- rstan::extract(LM_stan, c('beta', 'sigma2')) 
betas <- post_chain$beta
colnames(betas) <- colnames(X)

sigma2 <- post_chain$sigma2

df_plot <- data.frame(value = c(betas[,1], betas[,2], 
                                betas[,3], betas[,4], 
                                betas[,5], betas[,6],
                                sigma2), 
                      param = rep(c(colnames(X), 'sigma2'), 
                                  each = nrow(post_chain$beta)))


df_plot %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  theme_minimal() + 
  theme(legend.position="false")


acf(betas[,1])
acf(betas[,2])
acf(betas[,3])
acf(betas[,4])
acf(betas[,5])
acf(betas[,5])
acf(betas[,6])
acf(sigma2)


mean(betas[,2] > 4 & betas[,4] > -.5)

###################################################à


######################################
######## Logistic Regression: ########
######################################
rm(list = ls())
cardiac <- read.csv("cardiac.csv", header=T, sep=";")

y <- cardiac$Chd
X <- model.matrix(lm(Chd ~ ., data=cardiac))
n <- nrow(X)
K <- ncol(X)

Logistic_Model <- rstan::stan_model(file="Logistic_Model.stan")

data.stan_logis <- list(
  y = y,
  X = X,
  n = n,
  K = K,
  beta0 = rep(0, K),
  s2_0 = rep(10, K)
) 


n.iter <- 10000
nchain <- 2

Logistic_stan <- rstan::sampling(
  object = Logistic_Model,
  data = data.stan_logis,
  warmup = 0.5*n.iter, iter = n.iter,
  thin = 1, chains = nchain
)

summ <- summary(Logistic_stan, pars = "beta")$summary
rownames(summ)[1:K] <- colnames(X)
round(summ, 3)

rstan::traceplot(Logistic_stan, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Logistic_stan, pars = "beta", inc_warmup = FALSE)


post_chain <- rstan::extract(Logistic_stan, "beta") 
betas <- post_chain$beta
colnames(betas) <- colnames(X)


df_plot <- data.frame(value = c(betas[,1], betas[,2]), 
                      param = rep(colnames(X), 
                                  each = nrow(post_chain$beta)))

df_plot %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  theme_minimal() + 
  theme(legend.position="false")



# Evaluating probabilities:

probs_stan <- rstan::extract(Logistic_stan, "p")$p
plot(cardiac$Age, colMeans(probs_stan), pch=20, type="l")
plot(cardiac$Age, colMeans(probs_stan), pch=20)
