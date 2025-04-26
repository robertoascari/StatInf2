///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> p;   // number of covariates
	
	int<lower = 0> Y[N];  	// response vector
	matrix[N, p] X;   	// design matrix 
  
  vector[p] beta0;
  vector[p] s20;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p] beta;        	// regression coefficients 
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	for(i in 1:N)
	{
	  mu[i] = exp(row(X, i) * beta);
	}
}

////////////////// MODEL ////////////////////////
model {
	// Likelihood     
	for(n in 1:N)
	{
		Y[n] ~ poisson(mu[n]);  
	} 

	// Prior
	for(j in 1:p) 
	{
	 	beta[j] ~ normal(beta0[j], sqrt(s20[j]));
	}
}

////////////////// GENERATED QUANTITIES ////////////////////////
generated quantities 
{
  vector[N] log_lik;
  for(j in 1:N) 
	{
    log_lik[j] = poisson_lpmf(Y[j] | mu[j]);
  }
}
