double prob_dens_fpt__absorb_unit_sphere_r0_eq_0(const double t)
{
 

  if (t==0.0) return 0.0;


  if (t>=0.25)   // series for large times
  {
   const double expofac=-((double) 9.86960440108935861883449)*t;    
                                    //pi_sq

   if (t<1.5)
   {
       if (t<0.6)
       {
           if (t<0.3)
           {
             return ((double) 19.73920880217871723766898)*(exp(expofac)-4.0*exp(4.0*expofac)+ 9.0*exp(9.0*expofac)-16.0*exp(16.0*expofac));
           }                //two_pi_sq
           return ((double) 19.73920880217871723766898)*(exp(expofac)-4.0*exp(4.0*expofac)+ 9.0*exp(9.0*expofac));
       }                  //two_pi_sq
       return ((double) 19.73920880217871723766898)*(exp(expofac)-4.0*exp(4.0*expofac));
   }                     //two_pi_sq
   return ((double) 19.73920880217871723766898)*exp(expofac);
                        //two_pi_sq
  }
  else  // series for short times
  {
   const double expofac=-0.25/t;
   const double two_t=2.0*t;
   
   if (t>0.04)
   {
       if (t>0.125)
       {
         return ((double) 0.2820947917738781434740)/(sqrt(t)*t*t)* (exp(expofac)*(1.0-two_t )+exp(9.0*expofac)*(9.0-two_t)+ exp(25.0*expofac)*(25.0-two_t));
       }              //one_d_two_sqrt_pi
       return ((double) 0.2820947917738781434740)/(sqrt(t)*t*t)* (exp(expofac)*(1.0-two_t )+exp(9.0*expofac)*(9.0-two_t));
   }                       //one_d_two_sqrt_pi
   return ((double) 0.2820947917738781434740)/(sqrt(t)*t*t)* exp(expofac)*(1.0-two_t );
  }                 //one_d_two_sqrt_pi
}