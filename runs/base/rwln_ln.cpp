  #include <admodel.h>
	#undef REPORT
	#define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  adstring sppname;
#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <rwln_ln.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  styr.allocate("styr");
  endyr.allocate("endyr");
  yrs.allocate(styr,endyr);
 yrs.fill_seqadd(styr,1);
  nobs.allocate("nobs");
  yrs_srv.allocate(1,nobs,"yrs_srv");
  srv_est.allocate(1,nobs,"srv_est");
  srv_cv.allocate(1,nobs,"srv_cv");
  srv_sd.allocate(1,nobs);
  yvar.allocate(1,nobs);
  yconst.allocate(1,nobs);
 srv_sd = sqrt(log(1+elem_prod(srv_cv,srv_cv)));
 yvar=elem_prod(srv_sd,srv_sd);
 yconst= log(2.0*M_PI*yvar);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  logSdLam.allocate("logSdLam");
  mu.allocate("mu");
  biomsd.allocate(styr,endyr,"biomsd");
  u.allocate(styr,endyr,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  jnll =0.0;
  dvariable sig = exp(logSdLam);
  dvariable var = sig*sig ;
  jnll = .1* norm2(u);
  for(int i=styr+1; i<=endyr; ++i)
  {
    // step(biom(i-1),biom(i),logSdLam);
    // step(u(i-1),u(i),mu,logSdLam);
    // jnll          += 0.5*(log(2.0*M_PI*var)+square(sig*(u(i-1) - u(i)))/var);
    jnll += 0.5*square(u(i-1) - u(i));
    // cout<<jnll<<endl;
  }
  for(int i=1; i<=nobs; ++i)
  {
    // obs(u(yrs_srv(i)),mu, logSdLam, i);
    jnll+=0.5*(yconst(i) + square( log(mu * mfexp(sig*u(yrs_srv(i)))) - log(srv_est(i))) /yvar(i));
  }
  if (sd_phase()) biomsd = mu * mfexp(sig*u);
  /*
  SEPARABLE_FUNCTION void step(const dvariable& u1, const dvariable& u2, const dvariable& mu, const dvariable& logSdLam)
 dvariable sig = exp(logSdLam);
 dvariable var = sig*sig;
 jnll          += 0.5*(log(2.0*M_PI*var)+square(sig*(u1-u2))/var);
 
  SEPARABLE_FUNCTION void obs(const dvariable& u,const dvariable& mu,const dvariable& logSdLam, int i)
  dvariable sig = exp(logSdLam);
  jnll+=0.5*(yconst(i) + square(mu+sig*u-srv_est(i))/yvar(i));
  */
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << mu*mfexp(exp(logSdLam)*u) <<endl;
  // biomsd = meany*biom;
  // report << "yr"  <<endl;
  // report <<  yr   <<endl;
  // report << "est" <<endl;
}

void model_parameters::final_calcs()
{
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
  write_R(mu*u);
  write_R(UCI);
  mysum.close();
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(3000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
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
{}

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
  df1b2variable sig = exp(logSdLam);
  df1b2variable var = sig*sig ;
  jnll = .1* norm2(u);
  for(int i=styr+1; i<=endyr; ++i)
  {
    // step(biom(i-1),biom(i),logSdLam);
    // step(u(i-1),u(i),mu,logSdLam);
    // jnll          += 0.5*(log(2.0*M_PI*var)+square(sig*(u(i-1) - u(i)))/var);
    jnll += 0.5*square(u(i-1) - u(i));
    // cout<<jnll<<endl;
  }
  for(int i=1; i<=nobs; ++i)
  {
    // obs(u(yrs_srv(i)),mu, logSdLam, i);
    jnll+=0.5*(yconst(i) + square( log(mu * mfexp(sig*u(yrs_srv(i)))) - log(srv_est(i))) /yvar(i));
  }
  if (sd_phase()) biomsd = mu * mfexp(sig*u);
  /*
  SEPARABLE_FUNCTION void step(const df1b2variable& u1, const df1b2variable& u2, const df1b2variable& mu, const df1b2variable& logSdLam)
 df1b2variable sig = exp(logSdLam);
 df1b2variable var = sig*sig;
 jnll          += 0.5*(log(2.0*M_PI*var)+square(sig*(u1-u2))/var);
 
  SEPARABLE_FUNCTION void obs(const df1b2variable& u,const df1b2variable& mu,const df1b2variable& logSdLam, int i)
  df1b2variable sig = exp(logSdLam);
  jnll+=0.5*(yconst(i) + square(mu+sig*u-srv_est(i))/yvar(i));
  */
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

void df1b2_parameters::allocate(void) 
{
  logSdLam.allocate("logSdLam");
  mu.allocate("mu");
  biomsd.allocate(styr,endyr,"biomsd");
  u.allocate(styr,endyr,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
