// Template Model Builder version of the RE model. Started 8/2020
// Cole

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(styr);
  DATA_INTEGER(endyr);
  DATA_IVECTOR(yrs_srv);
  DATA_IVECTOR(yrs_srv_ind);
  DATA_VECTOR(srv_est);
  DATA_VECTOR(srv_cv);
  int nobs=yrs_srv.size();

  PARAMETER(logSdLam);
  PARAMETER_VECTOR(biom);

  // Gotta be a better way to do this
  int nyrs=endyr-styr+1; //yrs.size();
  vector<int> yrs(nyrs);
  for(int i=0; i<nyrs; i++) yrs(i)=styr+i;
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

  REPORT(biom);
  ADREPORT(biom)
  return(jnll);
}
