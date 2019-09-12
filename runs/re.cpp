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
  adstring sppname;
#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <re.htp>

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
 if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est);
 srv_sd = elem_prod(srv_cv,srv_cv);
 srv_sd += 1;
 srv_sd = sqrt(log(srv_sd));
  yvar.allocate(1,nobs);
  yconst.allocate(1,nobs);
 yvar=elem_prod(srv_sd,srv_sd);
 yconst= log(2.0*M_PI*yvar);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  logSdLam.allocate(1,"logSdLam");
  biomsd.allocate(styr,endyr,"biomsd");
  biomA.allocate(styr,endyr,"biomA");
  biom.allocate(styr,endyr,"biom");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  jnll =0.0;
  jnll=0.0;
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1),biom(i),logSdLam);
  }
  for(int i=1; i<=nobs; ++i)
  {
    obs(biom(yrs_srv(i)),i);
  }
  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }
	if (mceval_phase())
  {
    biomA = exp(biom);
    biomsd = biom;
	  cout <<biomA<<endl;
	}
	jnll += square(logSdLam+1.5);
}

void SEPFUN1  model_parameters::step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
{
  begin_df1b2_funnel();
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs(const dvariable& biom, int i)
{
  begin_df1b2_funnel();
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));
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
  report << biomsd <<endl;
  // biomsd = meany*biom;
  // report << "yr"  <<endl;
  // report <<  yr   <<endl;
  // report << "est" <<endl;
}

void model_parameters::final_calcs()
{
  /*
  dvar_vector UCI = elem_prod(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector LCI = elem_div(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector upp90th =  elem_prod(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector low90th =  elem_div(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  */
  dvar_vector UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_vector LCI = exp(biomsd-1.96*biomsd.sd);
  dvar_vector upp90th = exp(biomsd+1.645*biomsd.sd);
  dvar_vector low90th = exp(biomsd-1.645*biomsd.sd);
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs);
  write_R(LCI);
  write_R(biomA);
  write_R(UCI);
  write_R(low90th);
  write_R(upp90th);
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
  gradient_structure::set_MAX_NVAR_OFFSET(3000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
    mp.computations(argc,argv);
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
  jnll=0.0;
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1),biom(i),logSdLam);
  }
  for(int i=1; i<=nobs; ++i)
  {
    obs(biom(yrs_srv(i)),i);
  }
  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }
	if (mceval_phase())
  {
    biomA = exp(biom);
    biomsd = biom;
	  cout <<biomA<<endl;
	}
	jnll += square(logSdLam+1.5);
}

void   df1b2_pre_parameters::step(const funnel_init_df1b2variable& biom1, const funnel_init_df1b2variable& biom2, const funnel_init_df1b2variable& logSdLam)
{
  begin_df1b2_funnel();
  df1b2variable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs(const funnel_init_df1b2variable& biom, int i)
{
  begin_df1b2_funnel();
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));
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
  biomsd.deallocate();
  biomA.deallocate();
  biom.deallocate();
  prior_function_value.deallocate();
  likelihood_function_value.deallocate();
  jnll.deallocate();
} 
void df1b2_parameters::allocate(void) 
{
  logSdLam.allocate(1,"logSdLam");
  biomsd.allocate(styr,endyr,"biomsd");
  biomA.allocate(styr,endyr,"biomA");
  biom.allocate(styr,endyr,"biom");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
