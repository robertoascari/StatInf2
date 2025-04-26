///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 1> n;       // number of data
	int<lower = 1> K;       // number of covariates (including the intercept)
	int<lower = 0, upper = 1> y[n];   // response vector
	matrix[n, K] X;   		  // design matrix
  
  vector[K] beta0;
  vector[K] s2_0;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[K] beta;        // regression coefficients 
}

transformed parameters {
	vector[n] eta;
	vector<lower=0>[n] Odds;
	vector<lower=0,upper=1>[n] p;
	
  eta = X * beta;
  
  for(i in 1:n){
    Odds[i] = exp(eta[i]);
    p[i] = Odds[i]/(1+Odds[i]);
  }
  
}

////////////////// MODEL ////////////////////////
model {
  // Prior
	for (j in 1:K){
	 	beta[j] ~ normal(beta0, sqrt(s2_0));
	}
	
  // Likelihood     
	for (s in 1:n){
		y[s] ~ bernoulli(p[s]);  
	} 
}

