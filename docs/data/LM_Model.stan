data{
	int<lower = 0> n;  // number of obs
	int<lower = 0> K;  // number of covariates (including the intercept)
	
	vector[n] y;     // response
	matrix[n, K] X;	 // covariates

  real a0; 
  real b0;
	vector[K] beta0;	
	vector[K] s2_0;
}

parameters{
	real<lower = 0> sigma2;
	vector[K] beta;
}

transformed parameters {
	vector[n] mu;
  mu = X * beta;
}

model{
	// Prior:
	sigma2 ~ inv_gamma(a0, b0);
	
	for(k in 1:K){
		beta[k] ~ normal(beta0[k], sqrt(s2_0[k]));	
	}

	// Likelihood:
	y ~ normal(mu, sqrt(sigma2));

}

generated quantities{
  
  	vector[n] log_lik;
  	
  	for (j in 1:n){
    		log_lik[j] = normal_lpdf(y[j] | mu[j], sqrt(sigma2));
  	}
}
