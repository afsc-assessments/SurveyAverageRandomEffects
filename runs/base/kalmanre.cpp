#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <kalmanre.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  yrs.allocate("yrs");
  yrs_index.allocate(1,yrs,"yrs_index");
  index.allocate(1,yrs,"index");
  index_CV.allocate(1,yrs,"index_CV");
  index_sd.allocate(1,yrs,"index_sd");
 cout<<index_sd<<endl;
}

void model_parameters::initializationfunction(void)
{
  log_sigw.set_initial_value(0.064375);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  log_sigw.allocate(1,"log_sigw");
  var.allocate("var");
  #ifndef NO_AD_INITIALIZE
  var.initialize();
  #endif
  kalman_index.allocate(1,yrs,"kalman_index");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  obfun.allocate("obfun");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  obfun =0.0;
 /* var = square(mfexp(log_sigw));
 // Year 1
   prior_index(1)=B_0;
   prior_var(1)=var;
   MSE(1)=prior_var(1)+square(index_sd(1));
   e_t(1)=index(1)-prior_index(1);
   kalman_index(1)=prior_index(1)+(e_t(1)*prior_var(1))/MSE(1);
   post_var(1)=prior_var(1)-(square(prior_var(1)))/MSE(1);
 // Rest of years
   for(i=2; i<=yrs; i++){
   prior_index(i)=kalman_index(i-1);
   prior_var(i)=post_var(i-1)+var;
   MSE(i)=prior_var(i)+square(index_sd(i));
   e_t(i)=index(i)-prior_index(i);
   kalman_index(i)=prior_index(i)+(e_t(i)*prior_var(i))/MSE(i);
   post_var(i)=prior_var(i)-(square(prior_var(i)))/MSE(i);
   obfun +=0.5*(log(MSE(i))+(square(e_t(i))/MSE(i)));
       
  } */
  dvariable var;
  var = mfexp(log_sigw);
  cout<<var<<" var"<<endl;
  for(int i=2; i<=yrs; ++i){
    obfun += log(var) +0.5*square(kalman_index(i)-kalman_index(i-1))/(var*var);
  }
  for(int i=1; i<=yrs; ++i){
    obfun += log(index_sd(i)) + 0.5*square(kalman_index(i)-index(i))/(index_sd(i)*index_sd(i));
  }
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
  report<<"Index values"<<endl<<index<<endl;
  report<<"Kalman values"<<endl<<kalman_index<<endl;
  
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

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  obfun =0.0;
 /* var = square(mfexp(log_sigw));
 // Year 1
   prior_index(1)=B_0;
   prior_var(1)=var;
   MSE(1)=prior_var(1)+square(index_sd(1));
   e_t(1)=index(1)-prior_index(1);
   kalman_index(1)=prior_index(1)+(e_t(1)*prior_var(1))/MSE(1);
   post_var(1)=prior_var(1)-(square(prior_var(1)))/MSE(1);
 // Rest of years
   for(i=2; i<=yrs; i++){
   prior_index(i)=kalman_index(i-1);
   prior_var(i)=post_var(i-1)+var;
   MSE(i)=prior_var(i)+square(index_sd(i));
   e_t(i)=index(i)-prior_index(i);
   kalman_index(i)=prior_index(i)+(e_t(i)*prior_var(i))/MSE(i);
   post_var(i)=prior_var(i)-(square(prior_var(i)))/MSE(i);
   obfun +=0.5*(log(MSE(i))+(square(e_t(i))/MSE(i)));
       
  } */
  df1b2variable var;
  var = mfexp(log_sigw);
  cout<<var<<" var"<<endl;
  for(int i=2; i<=yrs; ++i){
    obfun += log(var) +0.5*square(kalman_index(i)-kalman_index(i-1))/(var*var);
  }
  for(int i=1; i<=yrs; ++i){
    obfun += log(index_sd(i)) + 0.5*square(kalman_index(i)-index(i))/(index_sd(i)*index_sd(i));
  }
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
  log_sigw.allocate(1,"log_sigw");
  var.allocate("var");
  #ifndef NO_AD_INITIALIZE
  var.initialize();
  #endif
  kalman_index.allocate(1,yrs,"kalman_index");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  obfun.allocate("obfun");  /* ADOBJECTIVEFUNCTION */
}
