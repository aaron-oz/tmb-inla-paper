// /////////////////////////////////////////////////////////////////////////////
// AOZ
// 2020
// Template file for fitting discrete space-only BYM2 models
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// /////////////////////////////////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for logit beta prior on theta = logit(probability)
//   s.t. probability = p ~ Beta(a, b)
template<class Type>
Type dlogitbeta(Type logit_p, Type a, Type b, int give_log=0)
{
  Type part1 = lgamma(a + b) - lgamma(a)  - lgamma(b);
  Type part2 = (a - 1) * (logit_p - log(1 + exp(logit_p)));
  Type part3 = (b - 1) * log( 1 - exp(logit_p)/(1 + exp(logit_p)));
  Type part4 =  logit_p - 2 * log( 1 + exp(logit_p));

  Type logres = part1 + part2 + part3 + part4;

  if(give_log)return logres; else return exp(logres);
}

// helper function for pc prior on log(precision) of multi-var normal
// where stdev sigma ~ N(u, s) & sigma > 0
// internally, INLA places the prior on log(precision)
template<class Type>
Type dlogtgaussian(Type log_prec, Type u, Type s, int give_log=0)
{
  Type part1 = -0.5 * log(8 * M_PI) - log(s);
  Type part2 = -0.5 * 1 / pow(s, 2) * pow( (exp(-log_prec / 2) - u ), 2) - log_prec / 2;

  Type logres = part1 + part2;

  if(give_log)return logres; else return exp(logres);
}


// make a structure to house our named boolean option flags
template<class Type>
struct option_list {

  int adreport_on;
  int normal_trick;
  int verbose;

  option_list(SEXP x){
    adreport_on  = asVector<int>(getListElement(x,"adreport_on"))[0];
    normal_trick = asVector<int>(getListElement(x,"normal_trick"))[0];
    verbose = asVector<int>(getListElement(x,"verbose"))[0];
  }
};


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
  // ~~~~~~~~~------------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------------~~

  // normalization flag
  DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll

  // Indices
  DATA_INTEGER( N );   // Number of spatial groups

  // Data
  DATA_VECTOR( y_i );      // Num occurrences (eg deaths) of poisson outcome per group
  DATA_VECTOR( e_i );      // populations per group
  DATA_MATRIX( X );  // Covariate 'design matrix', including column of 1s for intercept

  // Options
  DATA_STRUCT(options, option_list);  // boolean vector of options to be used to select different models/modelling options: see above

  // Prior specifications
  DATA_VECTOR( coef_pri );  //TODO specifying u and prec for alpha                      with alpha~Norm(u,prec)
  // DATA_VECTOR( g_var_pri ); //TODO specifying a and b for logitbeta       on logit(phi) with phi~Beta(a,b)
  // DATA_VECTOR( o_var_pri ); //TODO specifying u and prec for logtgaussian on log(tau)   with sigma~N(u,prec)

  if(options.verbose == 1){
    std::cout << "Data Loaded.\n";
  }

  // Fixed effects
  PARAMETER_VECTOR( betas ); // Intercept + betas
  // PARAMETER( log_grp_sigma ); // if using normal likelihood, log sd of single normal obs
  // PARAMETER( log_obs_sigma ); // logit of BYM2 mixing parameter

  // Random effects
  // PARAMETER_ARRAY( Epsilon_s );    // effects for each spatial
  // region. length 2N. first half is
  // total, seconf half is structured

  if(options.verbose == 1){
    std::cout << "Params Loaded.\n";
  }

  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~

  // objective function -- joint negative log-likelihood
  vector<Type> jnll_comp(3);
  Type jnll = 0; // sum of jnll_comp

  // set/print parallel info
  // max_parallel_regions = omp_get_max_threads();
  // printf("This is thread %d\n", max_parallel_regions);
  max_parallel_regions = 1;

  // Transform some of our parameters
  // Type bym2_phi = exp(logit_phi) / (1 + exp(logit_phi));
  // Type bym2_prec = exp(log_tau);
  // Type bym2_sigma = exp(-.5 * log_tau);

  // Define objects for derived values
  // vector<Type> fe_i(num_i); // main effect alpha + X_betas %*% t(betas)
  // vector<Type> total_group_i( N );  // latent field for each region
  vector<Type> risk_i( N );          // exp(total_group) risk for each region
  vector<Type> fe_i( N );          // exp(total_group) risk for each region

  // if(options.verbose == 1){
  //   std::cout << "Transforms initialized.\n";
  //   std::cout << printf("logit_phi is: %f\n", logit_phi);
  //   std::cout << printf("phi is: %f\n", bym2_phi);
  //   std::cout << printf("log_tau is: %f\n", log_tau);
  //   std::cout << printf("tau is: %f\n", bym2_prec);
  // }

  // evaluate fixed effects for intercept and covs if applicable
  fe_i = X * betas.matrix(); // init w/ covariate effects if using

  // get the poisson risk
  for(int i=0; i < N; i++){

    risk_i[i] = exp(fe_i[i]);

    if(options.verbose == 1){
      std::cout << printf("group mean %i is :%f\n", i, fe_i[i]);
      std::cout << printf("risk %i is: %f\n", i, risk_i[i]);
    }
  }

  // // Transform GMRFs and make vector form
  // for(int s = 0; s < N; s++){
  //   epsilon_s[s] = Epsilon_s(s);
  //   epsilon_s[N + s] = Epsilon_s(N + s);
  //   // unnecessary to start with the array and unlist in space only,
  //   // but helpful in sXt
  // }

  // if(options.verbose == 1){
  //   std::cout << "constructed epsilon_s.\n";
  // }

  // constrain the structured component of the GMRF to sum-to-zero
  // the constrained vector is calculated as:
  // x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  // for constraint Ax=e, and for GMRF x with precision Q
  // matrix<Type> Q_inv_s = invertSparseMatrix(Q_s);

  // if(options.verbose == 1){
  //   std::cout << "constructed Q_inv_s.\n";
  // }
  // matrix<Type> Q_inv(full_N, full_N);
  // for(int i=0; i < full_N; i++){
  //   for(int j=0; j < full_N; j++){
  //     Q_inv(i, j) = Q_inv_s.coeffRef(i, j);
  //   }
  // }

  // // column vector
  // vector<Type> A(full_N); // sum all components
  // for(int i = 0; i < N; i++){
  //   A( i )     = 0; // no constraints on total effects
  //   A( N + i ) = 1; // sum-to-0 on structured effects
  // }

  // if(options.verbose == 1){
  //   std::cout << "Constructed A.\n";
  // }

  // // row vector
  // matrix<Type> A_t = A.matrix().transpose();
  // if(options.verbose == 1){
  //   std::cout << "Constructed A_t.\n";
  // }

  // matrix<Type> QinvA_mat = Q_inv_s * A.matrix();

  // if(options.verbose == 1){
  //   std::cout << "Constructed QinvA_mat \n";
  //   //std::cout << printf("QinvA_mat dims are: %i x %i\n", QinvA_mat.rows(), QinvA_mat.cols());
  // }


  // vector<Type> QinvA(full_N);
  // for(int i=0;i<full_N;i++){
  //   QinvA(i) = QinvA_mat(i,0);
  // }

  // if(options.verbose == 1){
  //   std::cout << "Constructed QinvA \n";
  // }

  // Type e = 0;


  // Type AQinvA; // b/c of this particular A, we can find this
  //                  // value, (AQ^{-1}A')^{-1}, with summation and
  //                  // division
  // for(int i=N; i < full_N; i++){
  //   for(int j=N; j < full_N; j++){
  //     AQinvA += Q_inv(N+i, N+j);
  //   }
  // }

  // matrix<Type> AQinvA;
  // AQinvA = A_t * QinvA_mat;

  // Type AQinvA_inv = 1 / AQinvA(0,0);

  // if(options.verbose == 1){
  //   std::cout << "Constructed AQinvA_inv.\n";
  // }

  // // Type Ax = 0;
  // // for(int i=N; i < full_N; i++){
  // //   Ax += epsilon_s(i);
  // // }
  // Type Ax = (A * epsilon_s).sum();

  // if(options.verbose == 1){
  //   std::cout << "Constructed Ax.\n";
  // }

  // if(options.verbose == 1){
  //   std::cout << printf("QinvA: %f\n", QinvA);
  //   std::cout << printf("AQinvA_inv: %f\n", AQinvA_inv);
  //   std::cout << printf("Ax: %f\n", Ax);
  // }

  // with all the pieces, we can now constrain the vector
  //vector<Type> epsilon_s_c = epsilon_s - QinvA * (AQinvA_inv * Ax);
  // vector<Type> epsilon_s_c = epsilon_s - QinvA * AQinvA_inv * Ax;

  // if(options.verbose == 1){
  //   std::cout << "Constructed epsilon_s_c.\n";
  // }

  // //
  // for(int s = 0; s < N; s++){
  //   total_latent_i[s] = epsilon_s_c(s);     // first have is total
  //   struc_latent_i[s] = epsilon_s_c(N + s); // second half is structured effect
  // }

  // if(options.verbose == 1){
  //   std::cout << "Transforms filled in.\n";
  // }


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
  //       other likelihoods return positive log-liks so we subtract

  // With the constrained epsilon, their density is:
  // pi(epsilon | A*epsilon) = pi(A*epsilon|epsilon) * pi(epsilon) / pi (A*epsilon)
  // where
  // pi(epsilon)  is the usual GMRF with precision Q_s
  // pi(A*epsilon|epsilon) is a 1d gaussian with mean A*0 = 0, and covariance AQinvA^T
  // pi(Ax) is 0 if the constraint is not met, or the constant |AA^T|^{-.5} if it is metabolic
  //        with A = (0, 0, ..., 0, 1, 1, ..., 1), |AA^T| = 1, and log(1) = 0, so this contributes nothing

  // if(options.normal_trick == 1){
  //   // then we are not calculating the normalizing constant in the inner opt.
  //   // that norm constant means taking an expensive determinant of Q_s
  //   jnll_comp[0] -= dnorm(Ax, Type(0.0), Type(sqrt(AQinvA(0,0))), true); // N(mean, sd)
  //   jnll_comp[0] += GMRF(Q_s, false)(epsilon_s);
  // }else{
  //   jnll_comp[0] -= dnorm(Ax, Type(0.0), Type(sqrt(AQinvA(0,0))), true); // N(mean, sd)
  //   jnll_comp[0] += GMRF(Q_s)(epsilon_s);
  // }

  // // add other random effects here to get the 'flag' option working correctly

  // if(options.normal_trick == 1){
  //   if (flag == 0) return jnll_comp[0]; // return without data ll contrib to avoid unneccesary log(det(Q)) calcs
  // }

  // if(options.verbose == 1){
  //   std::cout << printf("\nGMRF contrib added: %f\n", jnll_comp[0]);
  // }

  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood

  // prior for intercept
  for(int i=0; i < N; i++){
    jnll_comp[1] -= dnorm(betas[i],
                          coef_pri[0], 1/sqrt(coef_pri[1]),
                          true);
  }

  // // add in hyperpriors for bym2
  // jnll_comp[1] -= dlogitbeta(logit_phi,
  //                            bym2_phi_pri[0], bym2_phi_pri[1],
  //                            true);

  // jnll_comp[1] -= dlogtgaussian(log_tau,
  //                               bym2_tau_pri[0], 1/sqrt(bym2_tau_pri[1]),
  //                               true);

  if(options.verbose == 1){
    std::cout << printf("\nPrior contrib added: %f\n", jnll_comp[1]);
  }

  /////////
  // (3) //
  /////////
  // Likelihood contribution from each datapoint i

  for (int i = 0; i < N; i++){

    // link to the poisson rate per person
    //    risk_i(i) = exp(alpha + total_latent_i(i));

    // and add data contribution to jnll
    jnll_comp[2] -= dpois( y_i(i), risk_i(i) * e_i(i), true );

  } // for( i )

  if(options.verbose == 1){
    std::cout << printf("\nData contrib added: %f\n", jnll_comp[2]);
  }

  // combine all parts of contribs to the jnll
  jnll += jnll_comp.sum();


  // ~~~~~~~~~~~
  // ADREPORT
  // get means and prec/cov of transformed params
  // include (untransformed) params if you want joint prec between params and transforms
  // ~~~~~~~~~~~

  //REPORT(QinvA_mat);
  //REPORT(AQinvA_inv);
  //REPORT(Ax);

  if(options.adreport_on == 1){

  }

  if(options.verbose == 1){
    std::cout << printf("Post adreport pre return jnll: %f\n", jnll);
  }

  return jnll;
}
