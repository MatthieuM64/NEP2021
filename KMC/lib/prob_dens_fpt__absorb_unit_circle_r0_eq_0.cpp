double prob_dens_fpt__absorb_unit_circle_r0_eq_0(const double t)
{

 static const double dens_prob_coeffs_sphere[]={9.264517580270208621681950383807e+00, -3.244577754672615939504108912007e+01, 6.375873709815058757416930950523e+01, -1.014500816746580385484944343021e+02, 1.445768632842918714938201859900e+02, -1.925230830350328337516009072838e+02, 2.448449156625870628374044408466e+02, -3.012026962876648103640785315704e+02, 3.613255335640216472347336785456e+02, -4.249907913096698091621445022049e+02, 4.920112476660101617495594704727e+02 }; //, -5.622265875052649993464184392143e+02, 6.354975101421602830517045661089e+02, -7.117015006454910030535809683286e+02, 7.907297051860868471946656265683e+02, -8.724845651458942553213638439223e+02, 9.568779882389944466431767682286e+02, -1.043829909339556121385804604407e+03, 1.133267140297022455140089811903e+03, -1.225122438127859947226674968774e+03, 1.319333740986874200141152503656e+03, -1.415843534955639517044350899650e+03, 1.514598324178592170136144268227e+03, -1.615548183616207643985087979894e+03, 1.718646378552378851435925761730e+03, -1.823849038564806943486852838330e+03, 1.931114876325525079420299471286e+03, -2.040404943603000351002107676925e+03};



 static const double dens_prob_expo_sphere[]={ -5.783185962946784521175995758456e+00, -3.047126234366208639907816317502e+01, -7.488700679069518344488904131014e+01, -1.390402844264598490015914235159e+02, -2.229323036176341569545620905221e+02, -3.265633529323284561822452014657e+02, -4.499335285180355267273004172701e+02, -5.930428696559552880436708062030e+02, -7.558913947839329622048319078704e+02, -9.384791134756942104843853843088e+02, -1.140806031099644870886625197458e+03}; //, -1.362872150854104466433176789779e+03, -1.604677474740231538001886216051e+03, -1.866222004061852903637839162150e+03, -2.147505739697844129854411310411e+03, -2.448528682258052414632583408480e+03, -2.769290832176358603290670812634e+03, -3.109792189768248983066127026002e+03, -3.470032755267558487311631648519e+03, -3.850012528850570343054190150127e+03, -4.249731510652209489687193849814e+03, -4.669189700777160279792635783751e+03, -5.108387099307648160175013108735e+03, -5.567323706308981979916711528707e+03, -6.045999521833564131820777022782e+03, -6.544414545923834084238845651162e+03, -7.062568778614457577863273598391e+03, -7.600462219933974562935697753222e+03};


 static const  double a=(( double) -5136.28518882094164887105639867650725371433537596128749041572721646282806552949);
 static const  double b=(( double) 1.0372029103575564047778792543082567798937510383250835584e5);
 static const  double c=(( double) -1.97865296397587452777963659995478965537583984827369745991571544e6);
 static const  double d=(( double) 2.230991541757848153613294707003107916673465607387497169211683e7);
   
   
 static const  double a2=(( double) -5143.32780379979613086595618594509226184890443955305840563653);
 static const  double b2=(( double) 1.065849721702697425503098134660093497169734879010511101633246160e5);
 static const  double c2=(( double) -2.3960413474540472298307632642251478263987980116927876246900624751e6);
 static const  double d2=(( double) 4.976974564956018268089028209045047146684625360553612229911e7);  
 static const  double e2=(( double) -7.99923736113403906007321258151484609528066907261643243416861e8);
 static const  double f2=(( double) 8.2309302022708892399388555231636277076946298373592058274883e9);
 static const  double g2=(( double) -3.92767642844452193733951085658753518280527227524902459966624158583e10);  

 if (t<=0.001) return 0.0;  // dens(t)<10^-100  for 0<t<0.001


 if(t<=0.015)
 {
   double one_over_t=1.0/t;
   
   double t_sq=t*t;
   double t3=t_sq*t;
   double t4=t_sq*t_sq;
  
   //std::cout<<exp(-0.25*one_over_t)*(   (0.5*one_over_t-0.5)*one_over_t  -0.5*t +2.0*t_sq -22.5*t3+296.0*t4+ a*t3*t_sq + b*t4*t_sq +c*t3*t4 +  d*t4*t4 )<<std::endl;
   
   return exp(-0.25*one_over_t)*(   (0.5*one_over_t-0.5)*one_over_t  -0.5*t +2.0*t_sq -22.5*t3+296.0*t4+ a*t3*t_sq + b*t4*t_sq +c*t3*t4 +  d*t4*t4 );
 }
 
 if(t<=0.032)
 {
   double one_over_t=1.0/t;
   
   double t_sq=t*t;
   double t3=t_sq*t;
   double t4=t_sq*t_sq;
   double t5=t3*t_sq;
   double t6=t3*t3;
   
   return exp(-0.25*one_over_t)*((0.5*one_over_t-0.5)*one_over_t  -0.5*t +2.0*t_sq -22.5*t3+296.0*t4+ a2*t5 + b2*t6 +c2*t3*t4 +  d2*t4*t4 + e2* t4*t5 +f2*t5*t5 +g2*t5*t6);
 }
 
   
   


   if (t<0.087)   // 0.032<t<0.087
   {
     if (t<0.044)   // 0.032<t<0.044
     {
        if (t<0.037)    // 0.032<t<0.037
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t) + dens_prob_coeffs_sphere[6]*exp(dens_prob_expo_sphere[6]*t) + dens_prob_coeffs_sphere[7]*exp(dens_prob_expo_sphere[7]*t) + dens_prob_coeffs_sphere[8]*exp(dens_prob_expo_sphere[8]*t) + dens_prob_coeffs_sphere[9]*exp(dens_prob_expo_sphere[9]*t) + dens_prob_coeffs_sphere[10]*exp(dens_prob_expo_sphere[10]*t);
        }
                        //0.037<t<0.044
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t) + dens_prob_coeffs_sphere[6]*exp(dens_prob_expo_sphere[6]*t) + dens_prob_coeffs_sphere[7]*exp(dens_prob_expo_sphere[7]*t) + dens_prob_coeffs_sphere[8]*exp(dens_prob_expo_sphere[8]*t) + dens_prob_coeffs_sphere[9]*exp(dens_prob_expo_sphere[9]*t);
     }
      // 0.044<t<0.087
        if (t<0.053)   // 0.044<t<0.053
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t) + dens_prob_coeffs_sphere[6]*exp(dens_prob_expo_sphere[6]*t) + dens_prob_coeffs_sphere[7]*exp(dens_prob_expo_sphere[7]*t) + dens_prob_coeffs_sphere[8]*exp(dens_prob_expo_sphere[8]*t);
        }
        if (t<0.067)   // 0.053<t<0.067
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t) + dens_prob_coeffs_sphere[6]*exp(dens_prob_expo_sphere[6]*t) + dens_prob_coeffs_sphere[7]*exp(dens_prob_expo_sphere[7]*t);
        }
                      // 0.067<t<0.087
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t) + dens_prob_coeffs_sphere[6]*exp(dens_prob_expo_sphere[6]*t);
     
   }
   //t>0.087

     if (t<0.278)   // 0.087<t<0.278
     {
        if (t<0.118)   // 0.087<t<0.118
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t) + dens_prob_coeffs_sphere[5]*exp(dens_prob_expo_sphere[5]*t);
        }
        if (t<0.172)   // 0.118<t<0.172
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t) + dens_prob_coeffs_sphere[4]*exp(dens_prob_expo_sphere[4]*t);
        }
                       // 0.172<t<0.278
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t) + dens_prob_coeffs_sphere[3]*exp(dens_prob_expo_sphere[3]*t);
     }
                  
                  // 0.278<t<infinity    

        if (t<0.53)   // 0.278<t<0.53
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t) +                 dens_prob_coeffs_sphere[2]*exp(dens_prob_expo_sphere[2]*t);
        }
        if (t<1.45)   // 0.53<t<1.45
        {
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t) + dens_prob_coeffs_sphere[1]*exp(dens_prob_expo_sphere[1]*t);
        }
                       // 1.45<t<\infinity
           return  dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t);

   abort();

   //double sum= dens_prob_coeffs_sphere[0]*exp(dens_prob_expo_sphere[0]*t);  
   //for(int loop=1;loop<=10;loop++) {sum+= dens_prob_coeffs_sphere[loop]*exp(dens_prob_expo_sphere[loop]*t);};
   //return sum;
}
