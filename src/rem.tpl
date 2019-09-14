 // Random walk model for survey averaging using log-scale process and observation errors
DATA_SECTION

  !!CLASS ofstream evalout("evalout.prj");

  // Start/end years
  init_int styr;
  init_int endyr;
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);

  // Define number of process error parameters
  init_int n_PE;

  // Define number of catchability parameters
  init_int n_q;  

  // Define bottom trawl survey data
  init_int num_indx;
  init_ivector PE_vec(1,num_indx);
  init_int nobs;
  init_number BTS_wt; 
  init_ivector yrs_srv(1,nobs);
  init_matrix srv_est(1,nobs,1,num_indx);
  init_matrix srv_cv(1,nobs,1,num_indx);
  matrix srv_sd(1,nobs,1,num_indx);
  matrix yvar(1,nobs,1,num_indx);
  matrix yconst(1,nobs,1,num_indx);
 LOCAL_CALCS
  logdat( styr);
  logdat( endyr);
  // Define number of process error parameters
  logdat( n_PE);

  // Define number of catchability parameters
  logdat( n_q);  

  // Define bottom trawl survey data
  logdat( num_indx);
  logdat( PE_vec);
  logdat( nobs);
  logdat( BTS_wt); 
  logdat( yrs_srv);
  logdat( srv_est);
  logdat( srv_cv);
  if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est+0.0001);
  srv_sd = elem_prod(srv_cv,srv_cv)+1;
  srv_sd = sqrt(log(srv_sd));
  yvar = elem_prod(srv_sd,srv_sd);
  yconst = log(2.0*M_PI*yvar);
 END_CALCS

  // Define srv2_est survey data
  init_int num_indx_srv2;
  init_ivector PE_vec_srv2(1,num_indx_srv2);
  init_int nobs_srv2;
  init_number srv2_wt;
  init_ivector yrs_srv2(1,nobs_srv2);
  init_matrix srv_est_srv2(1,nobs_srv2,1,num_indx_srv2);
  init_matrix srv_cv_srv2(1,nobs_srv2,1,num_indx_srv2);
  matrix srv_sd_srv2(1,nobs_srv2,1,num_indx_srv2);
 LOCAL_CALCS
  logdat( num_indx_srv2);
  logdat( PE_vec_srv2);
  logdat( nobs_srv2);
  logdat( srv2_wt); 
  logdat( yrs_srv2);
  logdat( srv_est_srv2);
  logdat( srv_cv_srv2);
  if (mean(srv_cv_srv2)>5) srv_cv_srv2 = elem_div(srv_cv_srv2,srv_est_srv2+0.0001);
  srv_sd_srv2 = elem_prod(srv_cv_srv2,srv_cv_srv2) + 1;
  srv_sd_srv2 = sqrt(log(srv_sd_srv2));
 END_CALCS
  matrix yvar_srv2(1,nobs_srv2,1,num_indx_srv2);
  matrix yconst_srv2(1,nobs_srv2,1,num_indx_srv2);
  !! yvar_srv2 = elem_prod(srv_sd_srv2,srv_sd_srv2);
  !! yconst_srv2 = log(2.0*M_PI*yvar_srv2);

PARAMETER_SECTION
  //init_bounded_number logSdq_srv2(0,2);
  init_bounded_vector logSdLam(1,n_PE,-5,2,1);
  init_vector log_q_srv2(1,n_q,2); 
  //init_number log_q_srv2(2);
  sdreport_matrix biomsd(styr,endyr,1,num_indx);
  sdreport_matrix biomA(styr,endyr,1,num_indx);
  sdreport_matrix srv2_est(styr,endyr,1,num_indx_srv2);
  //random_effects_vector log_q_srv2(1,n_q,2);  
  random_effects_matrix biom(styr,endyr,1,num_indx);
  sdreport_number Like;
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0.0;

  for(int j=1; j<=num_indx; ++j)
  {
    for(int i=styr+1; i<=endyr; ++i)
    {
      step(biom(i-1,j),biom(i,j),logSdLam(PE_vec(j)));
    }
    for(int i=1; i<=nobs; ++i)
    {
      if(srv_est(i,j)>-1) obs(biom(yrs_srv(i),j),i,j);
    }
  }

  //For now - hard wire srv2 survey estimates
  for(int i=styr; i<=endyr; ++i)
  {
    srv2_est(i) = mfexp(log_q_srv2 + biom(i));
    //srv2_est(i,2) = exp(log_q_srv2(2))*exp(biom(i,2));
    //srv2_est(i,3) = exp(log_q_srv2(3))*exp(biom(i,3));
  }

  for(int j=1; j<=num_indx_srv2; ++j)
  {
    for(int i=1; i<=nobs_srv2; ++i)
    {
      obs_srv2(srv2_est(yrs_srv2(i),j),i,j);
    }
  }

  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }

  if (mceval_phase()){
   evalout<<logSdLam<<" "<<jnll<<endl;}


SEPARABLE_FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);

SEPARABLE_FUNCTION void obs(const dvariable& biom, int i, int j)
  jnll+=BTS_wt*0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));

SEPARABLE_FUNCTION void obs_srv2(const dvariable& biom, int i, int j)
  jnll+=srv2_wt*0.5*(yconst_srv2(i,j) + square(log(biom)-log(srv_est_srv2(i,j)))/yvar_srv2(i,j));

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(777000);
  arrmblsize = 777000;

REPORT_SECTION
  biomsd = biom;
  Like = jnll;
  //report << biom << endl;
  report << srv2_est << endl;

GLOBALS_SECTION
  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  ofstream checkin("log_read.rep");
  #undef logdat
  #define logdat(object) checkin << #object "\n" << object << endl;
  adstring sppname;

FINAL_SECTION
  dvar_vector srv_est_TOT     = rowsum(srv_est);
  dvar_vector srv_est_TOT_srv2  = rowsum(srv_est_srv2);
  dvar_vector biom_TOT        = rowsum(biomA);
  dvar_vector SD_numer        = rowsum(elem_prod(exp(2*biomsd+square(biomsd.sd)),(exp(square(biomsd.sd))-1)));
  dvar_vector SD_denom        = square(rowsum(exp(biomsd+0.5*square(biomsd.sd))));
  dvar_vector SD_biom_TOT     = sqrt(log(elem_div(SD_numer,SD_denom)+1));
  dvar_vector biom_TOT_UCI    = exp(log(biom_TOT)+1.96*SD_biom_TOT);
  dvar_vector biom_TOT_LCI    = exp(log(biom_TOT)-1.96*SD_biom_TOT);
  dvar_vector biom_TOT_srv2     = rowsum(srv2_est);
  dvar_vector SD_numer_srv2     = rowsum(elem_prod(exp(2*log(srv2_est)+square(biomsd.sd)),(exp(square(biomsd.sd))-1)));
  dvar_vector SD_denom_srv2     = square(rowsum(exp(log(srv2_est)+0.5*square(biomsd.sd))));
  dvar_vector SD_biom_TOT_srv2  = sqrt(log(elem_div(SD_numer_srv2,SD_denom_srv2)+1));
  dvar_vector biom_TOT_UCI_srv2 = exp(log(biom_TOT_srv2)+1.96*SD_biom_TOT_srv2);
  dvar_vector biom_TOT_LCI_srv2 = exp(log(biom_TOT_srv2)-1.96*SD_biom_TOT_srv2);
  dvar_matrix UCI             = exp(biomsd+1.96*biomsd.sd);
  dvar_matrix LCI             = exp(biomsd-1.96*biomsd.sd);

  write_R(yrs_srv);
  write_R(srv_est_TOT);
  write_R(yrs_srv2);
  write_R(srv_est_TOT_srv2);
  write_R(log_q_srv2);
  write_R(yrs);
  write_R(biom_TOT);
  write_R(SD_biom_TOT);
  write_R(biom_TOT_UCI);
  write_R(biom_TOT_LCI);
  write_R(biom_TOT_srv2);
  write_R(SD_biom_TOT_srv2);
  write_R(biom_TOT_UCI_srv2);
  write_R(biom_TOT_LCI_srv2);
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs_srv2);
  write_R(srv_est_srv2);
  write_R(srv_sd_srv2);
  write_R(yrs);
  write_R(LCI);
  write_R(biomA);
  write_R(UCI);
  write_R(biomsd);
  write_R(biomsd.sd);

  mysum.close();
