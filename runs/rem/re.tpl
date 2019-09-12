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

  // Define LL survey data
  init_int num_indx_LL;
  init_ivector PE_vec_LL(1,num_indx_LL);
  init_int nobs_LL;
  init_number LL_wt;
  init_ivector yrs_srv_LL(1,nobs_LL);
  init_matrix srv_est_LL(1,nobs_LL,1,num_indx_LL);
  init_matrix srv_cv_LL(1,nobs_LL,1,num_indx_LL);
  matrix srv_sd_LL(1,nobs_LL,1,num_indx_LL);
  !! if (mean(srv_cv_LL)>5) srv_cv_LL = elem_div(srv_cv_LL,srv_est_LL+0.0001);
  !! srv_sd_LL = elem_prod(srv_cv_LL,srv_cv_LL) + 1;
  !! srv_sd_LL = sqrt(log(srv_sd_LL));
  matrix yvar_LL(1,nobs_LL,1,num_indx_LL);
  matrix yconst_LL(1,nobs_LL,1,num_indx_LL);
  !! yvar_LL = elem_prod(srv_sd_LL,srv_sd_LL);
  !! yconst_LL = log(2.0*M_PI*yvar_LL);

PARAMETER_SECTION
  //init_bounded_number logSdq_LL(0,2);
  init_bounded_vector logSdLam(1,n_PE,-5,2,1);
  init_vector log_q_LL(1,n_q,2); 
  //init_number log_q_LL(2);
  sdreport_matrix biomsd(styr,endyr,1,num_indx);
  sdreport_matrix biomA(styr,endyr,1,num_indx);
  sdreport_matrix LL_est(styr,endyr,1,num_indx_LL);
  //random_effects_vector log_q_LL(1,n_q,2);  
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

  //For now - hard wire LL survey estimates
  for(int i=styr; i<=endyr; ++i)
  {
    //LL_est(i,1) = exp(log_q_LL)*exp(biom(i,1));
    //LL_est(i,2) = exp(log_q_LL)*exp(biom(i,2));
    //LL_est(i,3) = exp(log_q_LL)*exp(biom(i,3));
    LL_est(i,1) = exp(log_q_LL(1))*exp(biom(i,1));
    LL_est(i,2) = exp(log_q_LL(2))*exp(biom(i,2));
    LL_est(i,3) = exp(log_q_LL(3))*exp(biom(i,3));
  }

  for(int j=1; j<=num_indx_LL; ++j)
  {
  for(int i=1; i<=nobs_LL; ++i)
  {
    obs_LL(LL_est(yrs_srv_LL(i),j),i,j);
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

SEPARABLE_FUNCTION void obs_LL(const dvariable& biom, int i, int j)
  jnll+=LL_wt*0.5*(yconst_LL(i,j) + square(log(biom)-log(srv_est_LL(i,j)))/yvar_LL(i,j));

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(777000);
  arrmblsize = 777000;

REPORT_SECTION
  biomsd = biom;
  Like = jnll;
  //report << biom << endl;
  report << LL_est << endl;

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
  dvar_vector srv_est_TOT = rowsum(srv_est);
  dvar_vector srv_est_TOT_LL = rowsum(srv_est_LL);
  dvar_vector biom_TOT = rowsum(biomA);
  dvar_vector SD_numer = rowsum(elem_prod(exp(2*biomsd+square(biomsd.sd)),(exp(square(biomsd.sd))-1)));
  dvar_vector SD_denom = square(rowsum(exp(biomsd+0.5*square(biomsd.sd))));
  dvar_vector SD_biom_TOT = sqrt(log(elem_div(SD_numer,SD_denom)+1));
  dvar_vector biom_TOT_UCI = exp(log(biom_TOT)+1.96*SD_biom_TOT);
  dvar_vector biom_TOT_LCI = exp(log(biom_TOT)-1.96*SD_biom_TOT);
  dvar_vector biom_TOT_LL = rowsum(LL_est);
  dvar_vector SD_numer_LL = rowsum(elem_prod(exp(2*log(LL_est)+square(biomsd.sd)),(exp(square(biomsd.sd))-1)));
  dvar_vector SD_denom_LL = square(rowsum(exp(log(LL_est)+0.5*square(biomsd.sd))));
  dvar_vector SD_biom_TOT_LL = sqrt(log(elem_div(SD_numer_LL,SD_denom_LL)+1));
  dvar_vector biom_TOT_UCI_LL = exp(log(biom_TOT_LL)+1.96*SD_biom_TOT_LL);
  dvar_vector biom_TOT_LCI_LL = exp(log(biom_TOT_LL)-1.96*SD_biom_TOT_LL);
  dvar_matrix UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_matrix LCI = exp(biomsd-1.96*biomsd.sd);

  write_R(yrs_srv);
  write_R(srv_est_TOT);
  write_R(yrs_srv_LL);
  write_R(srv_est_TOT_LL);
  write_R(log_q_LL);
  write_R(yrs);
  write_R(biom_TOT);
  write_R(SD_biom_TOT);
  write_R(biom_TOT_UCI);
  write_R(biom_TOT_LCI);
  write_R(biom_TOT_LL);
  write_R(SD_biom_TOT_LL);
  write_R(biom_TOT_UCI_LL);
  write_R(biom_TOT_LCI_LL);
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs_srv_LL);
  write_R(srv_est_LL);
  write_R(srv_sd_LL);
  write_R(yrs);
  write_R(LCI);
  write_R(biomA);
  write_R(UCI);
  write_R(biomsd);
  write_R(biomsd.sd);

  mysum.close();
