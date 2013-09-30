<<<<<<< HEAD
 // Random walk model for survey averaging using log-scale process and observation errors
DATA_SECTION
  init_int styr
  init_int endyr
=======
<<<<<<< HEAD
 // Random walk model for survey averaging using log-scale process and observation errors
DATA_SECTION
  init_int styr
  init_int endyr
=======
 // random walk model to test simulations from Paul Spencer
DATA_SECTION
  init_int styr
  init_int endyr
  // !! endyr = endyr+5;
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);
  init_int nobs
  init_ivector yrs_srv(1,nobs)
  init_vector srv_est(1,nobs)
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  // init_vector cv(1,nobs)
  init_vector srv_cv(1,nobs)
  vector srv_sd(1,nobs)
  !! if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est);
  !! srv_sd = elem_prod(srv_cv,srv_cv);
  !! srv_sd += 1;
  !! srv_sd = sqrt(log(srv_sd));
  number meany
  vector yvar(1,nobs)
  vector yconst(1,nobs)
 
  !! yvar=elem_prod(srv_sd,srv_sd);
  !! yconst= log(2.0*M_PI*yvar);
 
PARAMETER_SECTION
  init_number logSdLam
  sdreport_vector biomsd(styr,endyr);
  sdreport_vector biomA(styr,endyr);
<<<<<<< HEAD
=======
=======
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
  init_bounded_number logSdLam(-10,2)
  sdreport_vector biomsd(styr,endyr);
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  random_effects_vector biom(styr,endyr);
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0.0;
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1),biom(i),logSdLam);
  }

  for(int i=1; i<=nobs; ++i)
  {
    obs(biom(yrs_srv(i)),i);
  }
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }

SEPARABLE_FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);

SEPARABLE_FUNCTION void obs(const dvariable& biom, int i)
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));
<<<<<<< HEAD
=======
=======
  if (sd_phase()) biomsd = exp(biom);

FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  // cout<<biom2<<" "<<biom1<<" "<<jnll<<endl;

FUNCTION void obs(const dvariable& biom, int i)
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));
  // jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)-yvar(i)*0.5))/yvar(i));
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(3000);

REPORT_SECTION
<<<<<<< HEAD
  biomsd = biom;
  report << biomsd <<endl;
=======
<<<<<<< HEAD
  biomsd = biom;
  report << biomsd <<endl;
=======
  // biomsd = biom;
  report <<exp(logSdLam)<<" "<< exp(biom) <<endl;
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  // biomsd = meany*biom;
  // report << "yr"  <<endl;
  // report <<  yr   <<endl;
  // report << "est" <<endl;
GLOBALS_SECTION
  #include <admodel.h>
<<<<<<< HEAD
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
=======
<<<<<<< HEAD
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
=======
	#undef REPORT
	#define write_R(object) mysum << #object "\n" << object << endl;
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  ofstream mysum("rwout.rep");
  adstring sppname;

FINAL_SECTION
<<<<<<< HEAD
  /*
=======
<<<<<<< HEAD
  /*
=======
  // Aritmetic scale Confidence intervals (biomsd is on aritmetic scale)
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  dvar_vector UCI = elem_prod(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector LCI = elem_div(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector upp90th =  elem_prod(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  dvar_vector low90th =  elem_div(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  */
  dvar_vector UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_vector LCI = exp(biomsd-1.96*biomsd.sd);
  dvar_vector upp90th = exp(biomsd+1.645*biomsd.sd);
  dvar_vector low90th = exp(biomsd-1.645*biomsd.sd);

<<<<<<< HEAD
=======
=======
  /* 
  // Lognormal
  dvar_vector UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_vector LCI = exp(biomsd-1.96*biomsd.sd);
  */
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs);
  write_R(LCI);
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
  write_R(biomA);
  write_R(UCI);
  write_R(low90th);
  write_R(upp90th);
  write_R(biomsd);
  write_R(biomsd.sd);

  mysum.close();
<<<<<<< HEAD
=======

  mysum.close();
=======
  write_R(biom);
  write_R(UCI);
  write_R(upp90th);
  mysum.close();
 /*
 */
>>>>>>> 3d6e0e3c71d5eaf492b181a815dfc59d977a38af
>>>>>>> 7ee022b1356c0e7b9bfc99a1366a220192512353
