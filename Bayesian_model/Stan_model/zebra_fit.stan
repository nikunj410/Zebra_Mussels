functions { 
  real lq_port(int n, int Time, int intro_port, real N_knot, real h, matrix Q_tran, int nonzero, array[] int v,array[] int u, real lambda, real deff, vector p_sampling, array[,] int survey){
    real lq = 0;
    vector[n] pop, growth, p_detection;                                                                                          
    vector[nonzero] w = csr_extract_w(add_diag(deff*Q_tran,1));                                                             // non-zero elements of P. Used to create a sparce matrix of P.
    for (j in 1:Time){
      pop =  j==1 ? N_knot*one_hot_vector(n,intro_port) : csr_matrix_times_vector(n,n,w,v,u,rows_dot_product(growth,pop));  // simulating population dynamics using sparse matrix. N_t+1 = P'*G(N_t)*N_t
      growth = pow(lambda,1-pop);                                                                                           // Rietkerk Growth model [λ^(1-N)]
      p_detection = -expm1(-pop*h);                                                                                         // detection probability. θ_detection = 1-exp(-N*h).
      lq += bernoulli_lpmf(survey[j,] | rows_dot_product(p_sampling,p_detection));}                                         // observation model y~bernoulli(θ_sampling,θ_detection)
    return(lq);}                                                                                                            // returns log likelihood for observing zebra mussel for a given port of introduction 
  
  matrix Leslie_Matrix(int n_ages, vector sigma, vector f){                                                                 // Leslie Matrix elements. See Casagrandi et al Freshwater Biology (2007)
    matrix[n_ages, n_ages] Leslie = rep_matrix(0, n_ages, n_ages);
    Leslie[1, 2:n_ages] = sigma[1] * f' * 1e6;                                                                              // age-dependent fecundity
    for (i in 1:n_ages) if(i != n_ages) Leslie[i+1,i] = sigma[i+1]; else Leslie[i,i] = sigma[i+1];                          // transition probability from one class to nexr
    return(Leslie);}
}

data {
  int <lower = 1> n;                                    // number of ports = 160
  int <lower = 0> Time;                                 // time frame = 40 years (Jan 1980 to Dec 2019)
  int n_ages;                                           // number of age classes = 4
  int <lower = 1, upper = n> n_port;                    // number of possible ports of introduction
  int <lower = 0, upper = 1> prior_sampling;            // prior_sampling == 1 ? sample from prior : sample from posterior
  int <lower = 0, upper = 1> posterior_checks;          // posterior_checks == 1 ? simulate observations from posterior : skip
  array [n_port] int <lower = 1, upper = n> index_port; // ID's of possible ports of introduction
  array [Time,n] int <lower = 0, upper = 1> survey;     // presence data 
  int <lower = 0, upper = n*n> nonzero;                 // number of non-zero elements of P matrix
  matrix <lower = 0, upper = 1> [n,n] A;                // relative weighted adjacency matrix. Sum(A) = 1
  matrix[n,n] Dist;                                     // (standardized) distance between ports
  vector[n_ages + 1] sigma_mu, sigma_sd;                // mean and sd of priors for elements of the Leslie Matrix (see Casagrandi et al Freshwater Biology (2007)
  vector[n_ages - 1] f_mu, f_sd;                        // mean and sd of priors for elements of the Leslie Matrix (see Casagrandi et al Freshwater Biology (2007)
} 

transformed data{
  matrix <lower = -1, upper = 1> [n,n] Q_tran = (A - diag_matrix(A*ones_vector(n)))'; // Transpose of the rate transition matrix
  matrix[n,n] Dist_sq = square(Dist);                                                 // Square of standardized distance b/w ports
  real dmax = inv(max(-diagonal(Q_tran)));                                            // Throretically maximum value of De
  array[nonzero] int v = csr_extract_v(Q_tran);                                       // column indices for non-zero values P'
  array[n+1] int u = csr_extract_u(Q_tran);                                           // row start indices for non-zero values P'
  int time_sim = 2 * Time;                                                            // Time forward simulation = 80 years
  real dt = inv(3.0);                                                                 // One unit of time in simulations correspond to 4 months => dt = 1/3 years
  real start_year = 1980;                                                             // Time zero in simulation corresponds to year 1980
}

parameters {
  // parameters for the Life-history dynamics model
  vector[n_ages + 1] sigma_std;         // transition probability from one age class to next (centered)
  vector[n_ages - 1] f_std;             // age dependent fecundity (centered)
  
  // parameters of the spatial population dynamics model (aka process model)
  real <upper = 0> log_No;              // order of magnitude of density in year 1980 at port of introduction
  real <lower = 0, upper = dmax> deff;  // Effective dispersal; the lower limit ensures P is a stochastic matrix
  
  // parameters of the observation model
  real <lower = 0> h;                   // population as a fraction of carrying capacity at which the detection probability of the zebra mussel is 0.63; h = K*pi
  real beta_mu_obs;                     // inv_logit(beta_mu_obs) is the mean sampling probability
  vector [n] beta_obs;                  // inv_logit(beta_obs[i]) is the sampling probability at port i
  vector <lower = 0> [3] K_hyper;       // hyperparameters of the covariance matrix——eta, rho, and sigma
} 

transformed parameters {
  // transformed parameters to model Life-history dynamics
  vector <lower = 0, upper = 1> [n_ages + 1] sigma = sigma_mu + rows_dot_product(sigma_sd, sigma_std);   // transition probability from one age class to next (non-centered).
  vector <lower = 0> [n_ages - 1] f = f_mu + rows_dot_product(f_sd, f_std);                              // age dependent fecundity (non-centered)
  matrix <lower = 0> [n_ages, n_ages] Leslie = Leslie_Matrix(n_ages, sigma, f);                          // Leslie Matrix
  
  // transformed parameters to model spatial population dynamics
  vector[n_port] lq = log(uniform_simplex(n_port));                                                      // uniform discrete prior for the port of introduction
  real <lower = 1> lambda = pow(max(abs(get_real(eigenvalues(Leslie)))),dt);                             // λ is the dominant eigenvalue of the Leslie Matrix. Rescaled to simulation time
  real <lower = 0, upper = 1> N_knot = pow(10,log_No);                                                   // Density of zebra mussel at the port of introduction in 1980
  
  // transformed parameters to model the observation process
  vector[n] p_sampling = inv_logit(beta_mu_obs + cholesky_decompose(add_diag(square(K_hyper[1])*         // probability of sampling ports
                         exp(-0.5*Dist_sq*inv_square(K_hyper[2])),square(K_hyper[3])))*beta_obs); 
  
  //marginalizing over port of introduction. Execute only if  prior_sampling = FALSE
  if (prior_sampling == 0) for(i in 1:n_port) lq[i] += lq_port(n,Time,index_port[i],N_knot,h,Q_tran,nonzero,v,u,lambda,deff,p_sampling,survey);
}

model {
  target += n_port > 1 ? log_sum_exp(lq) : lq[1];                  // summing up the likelihood from all possible ports of introduction
  deff ~ exponential(1);                                           // wealky informative prior. De has to be << dmax ~ 10
  log_No ~ normal(-5,2.5);                                         // wealky informative prior
  h ~ exponential(1);                                              // wealky informative prior
  beta_obs ~ std_normal();                                         // wealky informative prior
  K_hyper ~ exponential(2);                                        // wealky informative prior
  beta_mu_obs ~ normal(0,1.5);                                     // wealky informative prior
  sigma_std ~ std_normal();                                        // informative prior based on data from Casagrandi et al Freshwater Biology (2007)
  f_std ~ std_normal();                                            // informative prior based on data from Casagrandi et al Freshwater Biology (2007)
} 

generated quantities{ 
  vector[n_port] p_intro_port = n_port > 1 ? exp(lq-log_sum_exp(lq)) : rep_vector(1,1); // probability that zebra mussels were introduced at a particular port
  int intro_port = index_port[categorical_rng(p_intro_port)];                           // sampling p_intro_port using a categorical random variable
  real<lower = 0, upper = 1> p_mus = inv_logit(beta_mu_obs);                            // mean probability of sampling ports
  real<lower = 0, upper = 1> p_Kd = -expm1(-h);                                         // probability of detection when the port is at carrying capacity
  real Deff = inv(dt)*deff;                                                             // scaled De to real time
  real r = inv(dt)*log(lambda);                                                         // r = ln(λ); where r is the per-capita growth rate scaled to real time
  real ln_N_star = uniform_rng(-log(1.02e9*22.321),-log(1.5e8*22.321));                 // -ln(zebra mussel carrying capacity)
  real intro_year = (ln_N_star-log(N_knot))*inv(r)+start_year;                          // year when the zebra mussel density was 1e-9K
  vector[posterior_checks == 1 ? n : 0] Pop;                                            // dummy variable to forward simulate zebra mussel population 
  array[posterior_checks == 1 ? time_sim : 0, posterior_checks == 1 ? n : 0] int Survey;// dummy variable to forward simulate observation process 
  
  if (posterior_checks == 1){
    for (j in 1 : time_sim){                                                            // simulating observations from 1980 to 2060
      Pop = j==1 ?  N_knot*one_hot_vector(n,intro_port) : add_diag(deff*Q_tran,1)*rows_dot_product(pow(lambda,1-Pop),Pop);
      Survey[j, : ] = bernoulli_rng(rows_dot_product(p_sampling, -expm1(-Pop * h)));}}
}
