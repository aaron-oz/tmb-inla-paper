// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from:
//     R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// helper function to use the same penalized complexity prior on
//  matern params that is used in INLA

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
// the main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------

  // normalization flag - used for speed-up
  DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll

  // Indices
  DATA_INTEGER( num_i );   // Number of data points in space
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh

  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_i );   // Num occurrences (deaths) per binomial experiment at point i (clust)
  DATA_VECTOR( n_i );   // Trials per cluster
  DATA_MATRIX( X_alpha );  // Covariate 'design matrix' for just intercept (i.e. 1col matrix of all 1s)

  // SPDE objects
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( Aproj );   // Used to project spatial mesh to data locations

  // Options
  DATA_VECTOR( options ); // TODO make this a named list!
  // options[0] == 1 : use normalization trick

  // Prior specifications
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( matern_pri);   // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit:    rho0
  Type matern_par_b = matern_pri[1]; // range prob:     alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit: sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob:  alpha_sigma

  // Fixed effects
  PARAMETER( alpha );             // Intercept
  PARAMETER( log_tau );           // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( log_kappa );         // Log of INLA kappa (related to spatial correlation and range)

  // Random effects
  PARAMETER_ARRAY( Epsilon_s );  // Random effect for each spatial mesh vertex

  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~

  // objective function -- joint negative log-likelihood
  Type jnll = 0;

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);

  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * log_tau) * exp(2.0 * log_kappa));

  // Define objects for derived values
  vector<Type> fe_i(num_i);              // main effect alpha + X_betas %*% t(betas)
  vector<Type> latent_field_i(num_i);    // Logit estimated prob for each cluster i
  vector<Type> epsilon_s(num_s);         // Epsilon_s (array) unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);     // value of gmrf at data points

  // fixed effects is just alpha in this example
  fe_i = X_alpha * Type(1.0); // initialize
  //for (int i = 0; i < num_i; i++){
  //  fe_i[i] = fe_i[i] + alpha; // add on intercept if using
  //}

  // Transform GMRFs and make vector form
  //for(int s = 0; s < num_s; s++){
  //  epsilon_s[s] = Epsilon_s(s); // This is probably unnecssary in space-only...
  //}

  // Project GP approx from mesh points to data points
  projepsilon_i = Aproj * epsilon_s.matrix();

  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) priors
  // 2) GP field
  // 3) data
  // ~~~~~~~~~------------------------------------------------~~-

  /////////
  // (1) //
  /////////
  // the random effects. we do this first so to do the normalization outside of every optimization step
  // 'GP' field contribution (i.e. log-lik of Gaussian-Markov random fields, GMRFs)
  // NOTE: likelihoods from namespace 'density' already return NEGATIVE log-liks so we add
  //       other likelihoods return positive log-liks
  if(options[0] == 1){
    // then we are not calculating the normalizing constant in the inner opt.
    // that norm constant means taking an expensive determinant of Q_ss
    jnll += GMRF(Q_ss, false)(epsilon_s);
    if (flag == 0) return jnll; // return without data ll contrib to avoid unneccesary log(det(Q)) calcs
  }else{
    jnll += GMRF(Q_ss)(epsilon_s);
  }

  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood (if options[1]==1)

  // add in priors for spde gp
  jnll -= dPCPriSPDE(log_tau, log_kappa,
                     matern_par_a, matern_par_b, matern_par_c, matern_par_d,
                     true);

  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri[0], alpha_pri[1], true); // N(mean, sd)

  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i

  for (int i = 0; i < num_i; i++){

    // latent field estimate at each obs
    latent_field_i(i) = fe_i(i) + projepsilon_i(i);

    // and add data contribution to jnll
    if(!isNA(y_i(i))){

     // Uses the dbinom_robust function, which takes the logit probability
      	jnll -= dbinom_robust( y_i(i), n_i(i), latent_field_i(i), true );

    } // !isNA

  } // for( i )


  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options[1]==1){
    ADREPORT(sp_range);
    ADREPORT(sp_sigma);
  }

  return jnll;

}
