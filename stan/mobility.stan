data{
  int<lower=0> N;
  int<lower=0> deaths[N];
  real mob[N];
  
  //variant specific attributes
  // h: dbution for interval from case to death.
  int<lower=0>  L_h_wildtype;
  real<lower=0> h_wildtype[L_h_wildtype];
  int<lower=0>  L_h_alpha;
  real<lower=0> h_alpha[L_h_alpha];
  // w: serial interval distribution
  int<lower=0>  L_w_wildtype;
  real<lower=0> w_wildtype[L_w_wildtype];
  int<lower=0>  L_w_alpha;
  real<lower=0> w_alpha[L_w_alpha];

  // proportion phi of alpha-attributable deaths
  real<lower=0, upper=1> phi[N];
}

transformed data{

  // attribute deaths to variants
  real D_wildtype[N];
  real D_alpha[N];
 
  // dot product function does not really exist for real * integer....
  for (n in 1:N){
    D_alpha[n] =  deaths[n] * phi[n];
    D_alpha[n] =  deaths[n] * (1-phi[n]);
  }

  //pre-calculate values of the gamma pdf
  // int deaths0 = 0; // assume 0 dvectoreaths on day 0
  // real omega[N];
  // real h[N];
  // real omega0; 
  // real h0;
  // for (n in 1:N){
    // convert mean and sd to alpha and beta
    // h[n] = exp(gamma_lpdf(n | (18.8)^2 / (8.46)^2, 18.8 / (8.46)^2)); 
    // omega[n] = exp(gamma_lpdf(n | (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2)); 
  // }
  // h0 = exp(gamma_lpdf(0 | (18.8)^2 / (8.46)^2, 18.8 / (8.46)^2)); 
  // omega0 = exp(gamma_lpdf(0 | (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2));
}

parameters{
  real R0_wildtype;
  real R0_alpha;
  real beta_wildtype;
  real beta_alpha;
  real<lower=0> delta;
}

transformed parameters{
  real R_wildtype[N];
  real R_alpha[N];
  real RD_wildtype[N];
  real RD_alpha[N];

  for (n in 1:N){
    R_wildtype[n] = exp(log(R0_wildtype) - beta_wildtype * (1-mob[n]));
    R_alpha[n] = exp(log(R0_alpha) - beta_alpha * (1-mob[n]));
  }


  // what happens when n=1 below...
  for (n in 1:N)
  {
    if(n < L_h_wildtype)
    {
      RD_wildtype[n] = R0_wildtype * h_wildtype[n]; 
    } else {
      RD_wildtype[n] = 0;
    }

    if(n>1)
    {
      for (s in 1:min(L_h_wildtype, n-1)){ 

        RD_wildtype[n] += R_wildtype[n-s] * h_wildtype[s];
      }
    }
    RD_wildtype[n] += 0; //R_wildtype[n] * h0; 
  }

  // what happens when n=1 below...
  for (n in 1:N)
  {
    if(n < L_h_alpha)
    {
      RD_alpha[n] = R0_alpha * h_alpha[n]; 
    } else {
      RD_alpha[n] = 0;
    }

    if(n>1)
    {
      for (s in 1:min(L_h_alpha, n-1)){ 
        RD_alpha[n] += R_alpha[n-s] * h_alpha[s];
      }
    }
    RD_alpha[n] += 0; //R_alpha[n] * h0; 
  }
  // for (n in 1:N){
  //   // s = 0
  //   RD_wildtype[n] = R0_wildtype * h_wildtype[n]; 
  //   RD_alpha[n] = R0_alpha * h_alpha[n]; 
  //   // s = 1, ..., n-1
  //     RD_alpha[n] += R_alpha[s] * h_alpha[n-s];
  //   }
  //   // I think I am assuming the below is 0...
  //   // s = n
  //   RD_wildtype[n] += 0; //R_wildtype[n] * h0; 
  //   RD_alpha[n] += 0; // R_alpha[n] * h0; 
  // }
}

model{
  // really don t like the name, instead conv_D_w.
  real conv_phiDw_wildtype;
  real conv_phiDw_alpha;
  int idx_wildtype;
  int idx_alpha;
  
  // Priors: 
  R0_wildtype ~ uniform(0,5);
  R0_alpha ~ uniform(0,5);
  beta_wildtype ~ uniform(-100,100);
  beta_alpha ~ uniform(-100,100);
  delta ~ exponential(1);
  
  // Model:
  // Do not consider the first deaths, there are no past deaths which 
  //   the model can use to justify their existence.
  // deaths[1] ~ neg_binomial_2(RD[1] * (deaths0 * omega[1] + deaths[1] * omega0), delta);
  //
  for (n in 2:N){
    // sumD_omega = deaths0 * omega[n]; // s = 0
    // initialise convolutions values
    conv_phiDw_wildtype = 0;
    conv_phiDw_alpha = 0;
    idx_wildtype = min(L_w_wildtype,  n-1 );
    idx_alpha = min(L_w_alpha,  n-1 );
    
    for (s in 1:idx_wildtype){ // s = 1, ..., n-1
      conv_phiDw_wildtype += D_wildtype[n-s] * w_wildtype[s];
    }

    for (s in 1:idx_alpha){ // s = 1, ..., n-1
      conv_phiDw_alpha += D_alpha[n-s] * w_alpha[s];
    }

    // no no no: w[0] should be 0
    // sumD_omega += deaths[n] * omega0; // s = n
    // deaths[n] ~ neg_binomial_2(RD_wildtype[n] * conv_phiDw_wildtype + RD_alpha[n] * conv_phiDw_alpha , delta);

  }
}

generated quantities {


}
