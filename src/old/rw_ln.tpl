 // random walk model to test simulations from Paul Spencer
DATA_SECTION
  init_int styr
  init_int endyr
  // !! endyr = endyr+5;
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);
  init_int nobs
  init_ivector yrs_srv(1,nobs)
  init_vector srv_est(1,nobs)
  init_vector srv_cv(1,nobs)
  vector srv_sd(1,nobs)
  number meany
  vector yvar(1,nobs)
  vector yconst(1,nobs)
  // !! srv_sd = elem_prod(srv_cv,srv_est);
  !! srv_sd = log(1.+elem_prod(srv_cv,srv_cv));
  !! yvar   = srv_sd;
  !! srv_sd = sqrt(srv_sd);
  !! yconst = log(2.0*M_PI*yvar);
  // !! cout <<srv_sd<<endl;

PARAMETER_SECTION
  init_bounded_number var_p(0.0001,2)
  sdreport_vector biomsd(styr,endyr);
  random_effects_vector biom(styr,endyr);
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0.0;
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1),biom(i),var_p);
  }

  for(int i=1; i<=nobs; ++i)
  {
    obs(biom(yrs_srv(i)),i);
  }
  if (sd_phase()) biomsd = exp(biom-mean(yvar)/2.);

FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& var_p)
  jnll+=0.5*(log(2.0*M_PI*var_p)+square(biom2-biom1)/var_p);
  // cout<<biom2<<" "<<biom1<<" "<<jnll<<endl;

FUNCTION void obs(const dvariable& biom, int i)
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));
  // jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)-yvar(i)*0.5))/yvar(i));

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(3000);

REPORT_SECTION
  biomsd = exp(biom-mean(yvar)/2.);
  report <<sqrt(var_p)<<" " << biomsd <<endl;
  // biomsd = meany*biom;
  // report << "yr"  <<endl;
  // report <<  yr   <<endl;
  // report << "est" <<endl;
GLOBALS_SECTION
  #include <admodel.h>
	#undef REPORT
	#define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  adstring sppname;

FINAL_SECTION
  dvar_vector UCI = elem_prod(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector LCI = elem_div(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  //  =Y3*EXP(2*SQRT(LN(1+Z3^2/Y3^2)))
  // dvar_vector c2 = biomsd+1.96*biomsd.sd;
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs);
  write_R(LCI);
  write_R(biom);
  write_R(UCI);
  mysum.close();
 /*
 */
