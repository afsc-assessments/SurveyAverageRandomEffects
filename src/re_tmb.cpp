// Template Model Builder version of the RE model. Started 8/2020
// Cole

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_IVECTOR(yrs);
  DATA_IVECTOR(yrs_srv);
  DATA_IVECTOR(yrs_srv_ind);
  DATA_VECTOR(srv_est);
  DATA_VECTOR(srv_cv);

  PARAMETER(logSdLam);
  PARAMETER_VECTOR(biom);

  int nyrs=yrs.size();
  int nobs=yrs_srv.size();
  vector<Type> srv_sd(srv_cv.size());
  srv_sd = srv_cv.array()*srv_cv.array();
  srv_sd = sqrt(log(1+srv_sd));

  Type jnll=0;
  // The random effect likelihood contribution
  for(int i=1; i<nyrs; i++){
    jnll-=dnorm(biom(i-1),biom(i), exp(logSdLam), true);
  }
  // The observational likelihood
  for(int i=0; i<nobs; i++){
    jnll-=dnorm(biom(yrs_srv_ind(i)), log(srv_est(i)), srv_sd(i), true);
  }
  jnll += pow(logSdLam+1.5,2);// Modest penalty to keep process error from getting too large...

  REPORT(biom);
  REPORT(srv_sd);
  ADREPORT(biom)
  return(jnll);
}
