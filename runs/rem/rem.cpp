#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  ofstream checkin("log_read.rep");
  #undef logdat
  #define logdat(object) checkin << #object "\n" << object << endl;
  adstring sppname;
#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <rem.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
  styr.allocate("styr");
  endyr.allocate("endyr");
  yrs.allocate(styr,endyr);
 yrs.fill_seqadd(styr,1);
  n_PE.allocate("n_PE");
  n_q.allocate("n_q");
  num_indx.allocate("num_indx");
  PE_vec.allocate(1,num_indx,"PE_vec");
  nobs.allocate("nobs");
  BTS_wt.allocate("BTS_wt");
  yrs_srv.allocate(1,nobs,"yrs_srv");
  srv_est.allocate(1,nobs,1,num_indx,"srv_est");
  srv_cv.allocate(1,nobs,1,num_indx,"srv_cv");
  srv_sd.allocate(1,nobs,1,num_indx);
  yvar.allocate(1,nobs,1,num_indx);
  yconst.allocate(1,nobs,1,num_indx);
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
  num_indx_srv2.allocate("num_indx_srv2");
  PE_vec_srv2.allocate(1,num_indx_srv2,"PE_vec_srv2");
  nobs_srv2.allocate("nobs_srv2");
  srv2_wt.allocate("srv2_wt");
  yrs_srv2.allocate(1,nobs_srv2,"yrs_srv2");
  srv_est_srv2.allocate(1,nobs_srv2,1,num_indx_srv2,"srv_est_srv2");
  srv_cv_srv2.allocate(1,nobs_srv2,1,num_indx_srv2,"srv_cv_srv2");
  srv_sd_srv2.allocate(1,nobs_srv2,1,num_indx_srv2);
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
  yvar_srv2.allocate(1,nobs_srv2,1,num_indx_srv2);
  yconst_srv2.allocate(1,nobs_srv2,1,num_indx_srv2);
 yvar_srv2 = elem_prod(srv_sd_srv2,srv_sd_srv2);
 yconst_srv2 = log(2.0*M_PI*yvar_srv2);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  logSdLam.allocate(1,n_PE,-5,2,1,"logSdLam");
  log_q_srv2.allocate(1,n_q,2,"log_q_srv2");
  biomsd.allocate(styr,endyr,1,num_indx,"biomsd");
  biomA.allocate(styr,endyr,1,num_indx,"biomA");
  srv2_est.allocate(styr,endyr,1,num_indx_srv2,"srv2_est");
  biom.allocate(styr,endyr,1,num_indx,"biom");
  Like.allocate("Like");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  jnll =0.0;
  ofstream& evalout= *pad_evalout;
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
}

void SEPFUN1  model_parameters::step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs(const dvariable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=BTS_wt*0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs_srv2(const dvariable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=srv2_wt*0.5*(yconst_srv2(i,j) + square(log(biom)-log(srv_est_srv2(i,j)))/yvar_srv2(i,j));
  end_df1b2_funnel();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  biomsd = biom;
  Like = jnll;
  //report << biom << endl;
  report << srv2_est << endl;
}

void model_parameters::final_calcs()
{
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
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(777000);
  arrmblsize = 777000;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
df1b2variable::pool = new adpool();
initial_df1b2params::varsptr = new P_INITIAL_DF1B2PARAMS[1000];
{
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
    mp.computations(argc,argv);
}
delete [] init_df1b2variable::list;
init_df1b2variable::list = NULL;
delete [] initial_df1b2params::varsptr;
initial_df1b2params::varsptr = NULL;
delete df1b2variable::pool;
df1b2variable::pool = NULL;
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

void model_parameters::preliminary_calculations(void){
  #if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

  #endif

}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  jnll =0.0;
  ofstream& evalout= *pad_evalout;
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
}

void   df1b2_pre_parameters::step(const funnel_init_df1b2variable& biom1, const funnel_init_df1b2variable& biom2, const funnel_init_df1b2variable& logSdLam)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  df1b2variable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs(const funnel_init_df1b2variable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=BTS_wt*0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs_srv2(const funnel_init_df1b2variable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=srv2_wt*0.5*(yconst_srv2(i,j) + square(log(biom)-log(srv_est_srv2(i,j)))/yvar_srv2(i,j));
  end_df1b2_funnel();
}

void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
df1b2_gradlist::set_no_derivatives(); 
quadratic_prior::in_qp_calculations=1; 
}  

void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
(*re_objective_function_value::pobjfun)=0; 
other_separable_stuff_begin(); 
f1b2gradlist->reset();  
if (!quadratic_prior::in_qp_calculations) 
{ 
df1b2_gradlist::set_yes_derivatives();  
} 
funnel_init_var::allocate_all();  
}  

void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
lapprox->do_separable_stuff(); 
other_separable_stuff_end(); 
funnel_init_var::deallocate_all(); 
} 

void model_parameters::begin_df1b2_funnel(void) 
{ 
if (lapprox)  
{  
{  
begin_funnel_stuff();  
}  
}  
}  

void model_parameters::end_df1b2_funnel(void) 
{  
if (lapprox)  
{  
end_df1b2_funnel_stuff();  
}  
} 
void df1b2_parameters::deallocate() 
{
  logSdLam.deallocate();
  log_q_srv2.deallocate();
  biomsd.deallocate();
  biomA.deallocate();
  srv2_est.deallocate();
  biom.deallocate();
  Like.deallocate();
  prior_function_value.deallocate();
  likelihood_function_value.deallocate();
  jnll.deallocate();
} 
void df1b2_parameters::allocate(void) 
{
  logSdLam.allocate(1,n_PE,-5,2,1,"logSdLam");
  log_q_srv2.allocate(1,n_q,2,"log_q_srv2");
  biomsd.allocate(styr,endyr,1,num_indx,"biomsd");
  biomA.allocate(styr,endyr,1,num_indx,"biomA");
  srv2_est.allocate(styr,endyr,1,num_indx_srv2,"srv2_est");
  biom.allocate(styr,endyr,1,num_indx,"biom");
  Like.allocate("Like");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
