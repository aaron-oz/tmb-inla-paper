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
  DATA_VECTOR( y_i );  // Num occurrences (eg deaths) of poisson outcome per group
  DATA_VECTOR( e_i );  // populations per group

  // Options
  DATA_STRUCT(options, option_list);  // boolean vector of options to be used to select different models/modelling options: see above

  // Prior specifications
  DATA_VECTOR( alpha_pri );  //TODO specifying u and prec for alpha                      with alpha~Norm(u,prec)
  DATA_VECTOR( g_tau_pri ); //TODO specifying a and b for logitbeta       on logit(phi) with phi~Beta(a,b)

  if(options.verbose == 1){
    std::cout << "Data Loaded.\n";
  }

  // Fixed effects
  PARAMETER( alpha );     // Intercept
  PARAMETER( log_g_tau ); // if using normal likelihood, log sd of single normal obs

  // Random effects
  PARAMETER_VECTOR( g_re_i );    // effects for each group

  if(options.verbose == 1){
    std::cout << "Params Loaded.\n";
  }

  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  // ~~~~~~~~~------------------------------------------------~~

  // objective function -- joint negative log-likelihood
  vector<Type> jnll_comp(2); // (0) is prior contrib, (1) is data contrib
  Type jnll = 0; // sum of jnll_comp

  // set/print parallel info
  // max_parallel_regions = omp_get_max_threads();
  max_parallel_regions = 1;

  // Transform some of our parameters
  Type g_prec  = exp(log_g_tau);
  Type g_sigma = exp(-.5 * log_g_tau);

  // Define objects for derived values
  vector<Type> g_mean_i( N ); // (int + g_re_i) linear predictor for each group
  vector<Type> risk_i( N );   // exp(g_mean_i)  risk for each grp

  if(options.verbose == 1){
    std::cout << "Transforms initialized.\n";
    std::cout << printf("g_sigma is: %f\n", g_sigma);
  }

  // get the total group mean and the poisson risk
  for(int i=0; i < N; i++){

    g_mean_i[i] = alpha + g_re_i[i];

    risk_i[i] = exp(g_mean_i[i]);

    if(options.verbose == 1){
      std::cout << printf("group mean %i is :%f\n", i, g_mean_i[i]);
      std::cout << printf("risk %i is: %f\n", i, risk_i[i]);
    }
  }

  if(options.verbose == 1){
    std::cout << "constructed risks .\n";
  }

  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) priors
  // 2) data
  // ~~~~~~~~~------------------------------------------------~~-
  // ~~~~~~~~~------------------------------------------------~~-


  /////////
  // (1) // Prior contributions to joint likelihood
  /////////

  // prior for intercept
  jnll_comp[0] -= dnorm(alpha,
                       alpha_pri[0], 1/sqrt(alpha_pri[1]),
                       true);

  // RE distribution/prior
  jnll_comp[0] -= sum(dnorm(g_re_i,
                        0.0, g_sigma,
                            true));

  // prior for group precision
  jnll_comp[0] -= dlogtgaussian(log_g_tau,
                                g_tau_pri[0], 1/sqrt(g_tau_pri[1]),
                                true);

  if(options.verbose == 1){
    std::cout << printf("\nPrior contrib added: %f\n", jnll_comp[0]);
  }

  /////////
  // (2) // Likelihood contribution from each datapoint i
  /////////


  for (int i = 0; i < N; i++){

    // and add data contribution to jnll
    jnll_comp[1] -= dpois( y_i(i), risk_i(i) * e_i(i), true );

  } // for( i )

  if(options.verbose == 1){
    std::cout << printf("\nData contrib added: %f\n", jnll_comp[1]);
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
    ADREPORT(g_prec);
  }

  if(options.verbose == 1){
    std::cout << printf("Post adreport pre return jnll: %f\n", jnll);
  }

  return jnll;
}
