DATA_SECTION
  int i
  init_int yrs
  init_ivector yrs_index(1,yrs)
  init_vector index(1,yrs)
  init_vector index_CV(1,yrs)
  init_vector index_sd(1,yrs)
  !! cout<<index_sd<<endl;
INITIALIZATION_SECTION
  log_sigw 0.064375

PARAMETER_SECTION
  init_number log_sigw(1);
 // init_number B_0(1);
  number var;
//  init_number logSdy
  random_effects_vector kalman_index(1,yrs);
//  vector kalman_index(1,yrs);
 // vector prior_index(1,yrs);
//  vector prior_var(1,yrs);
 // vector e_t(1,yrs);
 // vector MSE(1,yrs);
//  vector post_index(1,yrs);
 // vector post_var(1,yrs);
//  vector loglike(1,yrs);
  
  // init_vector lam(1,N);
  objective_function_value obfun;

PROCEDURE_SECTION
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
REPORT_SECTION
  report<<"Index values"<<endl<<index<<endl;
  report<<"Kalman values"<<endl<<kalman_index<<endl;
//  report<<"prior variance"<<endl<<prior_var<<endl;
//  report<<"posterior variance"<<endl<<post_var<<endl;
//  report<<"error values"<<endl<<e_t<<endl;
//  report<<"MSE values"<<endl<<MSE<<endl;
//  report<<"prior index"<<prior_index<<endl;
  
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(3000);
