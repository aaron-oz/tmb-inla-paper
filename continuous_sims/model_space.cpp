// /////////////////////////////////////////////////////////////////////////////
// RB AOZ
// 2019
// Template file for fitting space only models
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 3. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// /////////////////////////////////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// helper function for pc prior on precision of multi-var normal
// type 2 gumbel
template<class Type>
Type dPCPriPrec(Type tau, Type u, Type a, int give_log=0)
{
  Type logres;
  Type lambda = -log(a) / u; //) P( 1/sqrt(tau) > u ) = a 
  logres = log(lambda/Type(2.0)) - Type(3.0/2.0)*log(tau) - lambda*pow(tau, -.5);
  if(give_log)return logres; else return exp(logres);
}

// helper function for pc prior on log(precision) of multi-var normal
// internally, INLA places the prior on log(precision)
template<class Type>
Type dPCPriLogPrec(Type log_prec, Type u, Type a, int give_log=0)
{
  Type logres;
  Type lambda = -log(a) / u; //) P( 1/sqrt(tau) > u ) = a 
  logres = log(lambda/Type(2.0)) - Type(1.0/2.0)*log_prec - lambda*exp(-log_prec/2.0);
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type dPCPriSPDE(Type logtau, Type logkappa,
                Type matern_par_a, Type matern_par_b, 
                Type matern_par_c, Type matern_par_d,
                //vector<Type> matern_pri(4),
                int give_log=0)
{
  
  // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  
  Type penalty; // prior contribution to jnll
  
  Type d = 2.;  // dimension
  Type lambda1 = -log(matern_par_b) * pow(matern_par_a, d/2.);
  Type lambda2 = -log(matern_par_d) / matern_par_c;
  Type range   = sqrt(8.0) / exp(logkappa);
  Type sigma   = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  
  penalty = (-d/2. - 1.) * log(range) - lambda1 * pow(range, -d/2.) - lambda2 * sigma;
  // Note: (rho, sigma) --> (x=log kappa, y=log tau) -->
  //  transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
  //  --> Jacobian: |J| propto e^(-y -2x)
  Type jacobian = - logtau - 2.0*logkappa;
  penalty += jacobian;
  
  if(give_log)return penalty; else return exp(penalty);
}

///////////////////////////
// our main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ~~~~~~~~~------------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------
  // ~~~~~~~~~------------------------------------------------------~~

  // normalization flag
  DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( num_i );   // Number of data points in space
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh

  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_i );   // Num occurrences (deaths) per binomial experiment at point i (clust)
  DATA_VECTOR( n_i );   // Trials per cluster
  DATA_MATRIX( X_alpha );  // Covariate 'design matrix' for just intercept (i.e. 1col matrix of all 1s)
  DATA_MATRIX( X_betas );   // Covariate design matrix excluding intercept columm

  // SPDE objects
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( Aproj );   // Used to project spatial mesh to data locations

  // Options
  DATA_VECTOR( options ); // TODO make this a named list!
  // options[0] == 1 : turn on adreport
  // options[1] == 1 : use priors
  // options[2] == 1 : fit with intercept
  // options[3] == 1 : fit with cov effects
  // options[4] == 1 : fit with cluster RE
  // options[5] == 0 : use normal data lik
  // options[5] == 1 : use binom  data lik
  // options[6] == 1 : use normalization trick

  // Prior specifications
  DATA_VECTOR( norm_prec_pri );
  DATA_VECTOR( clust_prec_pri );
  DATA_VECTOR( alphaj_pri );
  // DATA_VECTOR( logtau_pri );
  // DATA_VECTOR( logkappa_pri );
  DATA_VECTOR( matern_pri);   // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit:    rho0
  Type matern_par_b = matern_pri[1]; // range prob:     alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit: sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob:  alpha_sigma

  // Fixed effects
  PARAMETER( alpha );             // Intercept
  PARAMETER_VECTOR( betas );      // Covariate coefficients
  // PARAMETER( log_gauss_sigma );  // if using normal likelihood, log sd of single normal obs
  PARAMETER( log_gauss_prec );      // if using normal likelihood, log prec of single normal obs
  PARAMETER( log_tau );           // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( log_kappa );         // Log of INLA kappa (related to spatial correlation and range)
  // PARAMETER( log_clust_sigma );// Log of SD for cluster RE
  PARAMETER( log_clust_prec );    // Log of prec for cluster RE
  PARAMETER_VECTOR( clust_i );    // random effect estimate for each cluster observation
  
  // Random effects
  PARAMETER_ARRAY( Epsilon_s );  // Random effect for each spatial mesh location. Currently a 1d array of num_s

  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  
  // objective function -- joint negative log-likelihood
  Type jnll = 0;
  //parallel_accumulator<Type> jnll(this); // parallelize jnll NOTE: seems to break with omp>1 and mkl >=1

  // print parallel info
  // max_parallel_regions = omp_get_max_threads();
  max_parallel_regions = 1;
  // printf("This is thread %d\n", max_parallel_regions);

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);

  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * log_tau) * exp(2.0 * log_kappa));
  // Type log_clust_prec = -2.0 * log_clust_sigma; // needed if using non-pc prior
  // Type log_gauss_prec = -2.0 * log_gauss_sigma;   // needed if using non-pc prior
  // Type clust_sigma = exp(log_clust_sigma);       // used when passing in log(sigma) instead of log(prec)
  // Type gauss_sigma = exp(log_gauss_sigma);         // used when passing in log(sigma) instead of log(prec)
  // Type clust_prec = exp(-2.0 * log_clust_sigma); // used when passing in log(sigma) instead of log(prec)
  // Type gauss_prec = exp(-2.0 * log_gauss_sigma);   // used when passing in log(sigma) instead of log(prec)
  Type clust_sigma = exp(-.5 * log_clust_prec);
  Type gauss_sigma = exp(-.5 * log_gauss_prec);
  Type clust_prec  = exp(log_clust_prec);
  Type gauss_prec  = exp(log_gauss_prec);
  
    
  // Define objects for derived values
  vector<Type> fe_i(num_i);              // main effect alpha + X_betas %*% t(betas)
  vector<Type> latent_field_i(num_i);    // Logit estimated prob for each point i
  vector<Type> epsilon_s(num_s);         // Epsilon_s (array) unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);     // value of gmrf at data points
  
  // evaluate fixed effects for intercept and covs if applicable
  fe_i = X_alpha * Type(0.0); // initialize
  if(options[3] == 1){
    fe_i = X_betas * betas.matrix(); // init w/ covariate effects if using
  }
  
  if(options[2] == 1){
    for (int i = 0; i < num_i; i++){
      fe_i[i] = fe_i[i] + alpha; // add on intercept if using
    }
  }

  // Transform GMRFs and make vector form
  for(int s = 0; s < num_s; s++){
    epsilon_s[s] = Epsilon_s(s); // This is probably unnecssary in space-only...
  }

  // Project GP approx from mesh points to data points 
  projepsilon_i = Aproj * epsilon_s.matrix();

  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) priors
  // 2) GP field
  // 3) data
  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-

  ///////// 
  // (1) // 
  /////////
  // the random effects. we do this so to do the normalization outside of every optimization step
  // 'GP' field contribution (i.e. log-lik of Gaussian-Markov random fields, GMRFs)
  // NOTE: likelihoods from namespace 'density' already return NEGATIVE log-liks so we add
  //       other likelihoods return positive log-liks
  if(options[6] == 1){
    // then we are not calculating the normalizing constant in the inner opt.
    // that norm constant means taking an expensive determinant of Q_ss
    jnll += GMRF(Q_ss, false)(epsilon_s);
  }else{
    jnll += GMRF(Q_ss)(epsilon_s);
  }
  
  // cluster contribution to the likelihood
  if(options[4] == 1 ){
    for (int i = 0; i < num_i; i++){
      jnll -= dnorm(clust_i(i), Type(0.0), clust_sigma, true);
    }
  }
  
  if(options[6] == 1){
    if (flag == 0) return jnll; // return without data ll contrib to avoid unneccesary log(det(Q)) calcs
  }
  

  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood (if options[1]==1)
  if(options[1] == 1) {
    
    // add in priors for spde gp
    // jnll -= dnorm(log_tau,   logtau_pri[0], sqrt(1/logtau_pri[1]), true);     // N(mean, sd) prior for logtau
    // jnll -= dnorm(log_kappa, logkappa_pri[0], sqrt(1/logkappa_pri[1]), true); // N(mean, sd) prior for logkappa
    jnll -= dPCPriSPDE(log_tau, log_kappa,
                       matern_par_a, matern_par_b, matern_par_c, matern_par_d,
                       true);
    
    // prior for intercept
    if(options[2] == 1){
      jnll -= dnorm(alpha, alphaj_pri[0], alphaj_pri[1], true); // N(mean, sd)
    }
    
    // prior for covariate coefs
    if(options[3] == 1){
      for( int j = 0; j < betas.size(); j++){
        jnll -= dnorm(betas(j), alphaj_pri[0], alphaj_pri[1], true); // N(mean, sd)
      }
    }
    
    // prior for log(cluster RE prec)
    if(options[4] == 1){

      // pc.prior on precision
      // jnll -= dPCPriPrec(lclust_prec, clust_prec_pri[0], clust_prec_pri[1], true); //type2gumbel. P(1/sqrt(prec)>u)=a
      jnll -= dPCPriLogPrec(log_clust_prec, clust_prec_pri[0], clust_prec_pri[1], true); // P(1/sqrt(prec)>u)=a
      
      // in tmb, log(X)~logGamma(shape, scale) where X~Gamma(shape, scale)
      // almost like INLA, except INLA uses params (shape, inv-scale) which is *_prec_pri is defined
      //jnll -= dlgamma(log_clust_prec, clust_prec_pri[0], 1.0 / clust_prec_pri[1], true); // tmb takes (shape, scale)
      // don't forget the jacobian!
      // optimizer is passing in log(sigma), but the prior is defined on log(tau)
      // pi(log(sigma)) = pi(log(tau))*|dflog(sigma)/dlog(sigma)|, where log(tau) = f(log(sigma))
      // so, f(x) = -2x --> jacoabian = 2
      //jnll -= log(2.0);
    }
    
    // prior for log(obs prec) if using normal data lik
    if(options[5] == 0){
      
      // pc.prior on precision
      // jnll -= dPCPriPrec(gauss_prec, norm_prec_pri[0], norm_prec_pri[1], true); //type2gumbel. P(1/sqrt(prec)>u)=a
      jnll -= dPCPriLogPrec(log_gauss_prec, norm_prec_pri[0], norm_prec_pri[1], true); // P(1/sqrt(prec)>u)=a
      
      // in tmb, log(X)~logGamma(shape, scale) where X~Gamma(shape, scale)
      // almost like INLA, except INLA uses params (shape, inv-scale) which is *_prec_pri is defined
      //jnll -= dlgamma(log_gauss_prec, norm_prec_pri[0], 1.0 / norm_prec_pri[1], true); // tmb takes (shape, scale)
      // don't forget the jacobian!
      // optimizer is passing in log(sigma), but the prior is defined on log(tau)
      // pi(log(sigma)) = pi(log(tau))*|dflog(sigma)/dlog(sigma)|, where log(tau) = f(log(sigma))
      // so, f(x) = -2x --> jacoabian = 2
      //jnll -= log(2.0);
    }
    
  } 

  /////////
  // (3) //
  /////////
  // Likelihood contribution from each datapoint i

  for (int i = 0; i < num_i; i++){

    // latent field estimate at each obs
    latent_field_i(i) = fe_i(i) + projepsilon_i(i);
    // add on cluster if using
    if(options[4] == 1){ 
      latent_field_i(i) = latent_field_i(i) + clust_i(i);
    }

    // and add data contribution to jnll
    if(!isNA(y_i(i))){

      // if normal
      if(options[5] ==  0){
      	jnll -= dnorm( y_i(i), latent_field_i(i), gauss_sigma/sqrt(n_i(i)), true );
      }

      // if binom
      if(options[5] ==  1){
	      // Uses the dbinom_robust function, which takes the logit probability
      	jnll -= dbinom_robust( y_i(i), n_i(i), latent_field_i(i), true );
      }
      
    } // !isNA

  } // for( i )
  

  // ~~~~~~~~~~~
  // ADREPORT
  // ~~~~~~~~~~~
  if(options[0] == 1){
    // ADREPORT(betas);
    // ADREPORT(Epsilon_s);
    ADREPORT(clust_prec);
    ADREPORT(gauss_prec);
    ADREPORT(sp_range);
    ADREPORT(sp_sigma);
  }

  return jnll;
}




