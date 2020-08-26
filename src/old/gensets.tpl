DATA_SECTION
  init_int n
  init_ivector nyrs(1,n)
  init_matrix yrs(1,n,1,nyrs)
  init_matrix biom(1,n,1,nyrs)
  init_matrix cv(1,n,1,nyrs)
 LOCAL_CALCS
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
 END_CALCS

PARAMETER_SECTION
  init_number b1
  vector yfit(1,n)
  objective_function_value RSS

PROCEDURE_SECTION

GLOBALS_SECTION
  #include <admodel.h>
  adstring simname;
	#undef write_sim
	#define write_sim(object) simout << "# " #object "\n" << object << endl;
