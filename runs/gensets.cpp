  #include <admodel.h>
  adstring simname;
	#undef write_sim
	#define write_sim(object) simout << "# " #object "\n" << object << endl;
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <gensets.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  nyrs.allocate(1,n,"nyrs");
  yrs.allocate(1,n,1,nyrs,"yrs");
  biom.allocate(1,n,1,nyrs,"biom");
  cv.allocate(1,n,1,nyrs,"cv");
  char buffer [33];
  for (int i=1;i<=n;i++)
  {
    simname = "sims/sim_"+ adstring(itoa(i,buffer,10)) + ".dat";
    ofstream simout(simname);
    simout<< 1960<<endl<<2013<<endl;
    write_sim(nyrs(i));
    write_sim(yrs(i));
    write_sim(biom(i));
    write_sim(  cv(i));
    simout.close();
  }
  exit(1);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  b1.allocate("b1");
  yfit.allocate(1,n,"yfit");
  #ifndef NO_AD_INITIALIZE
    yfit.initialize();
  #endif
  RSS.allocate("RSS");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  RSS =0.0;
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
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
