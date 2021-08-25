/*C++ CODE - "NARROW ESCAPE PROBLEM IN TWO-SHELL SPHERICAL DOMAINS" by M. MANGEAT AND H. RIEGER
  3D NUMERICAL SIMULATIONS WITH KINETIC MONTE CARLO METHOD - WRITTEN BY M. MANGEAT (2021)*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//Personal libraries.
#include "lib/random.cpp"
#include "lib/special_functions.cpp"

//Kinetic Monte Carlo (KMC) library created by Karsten Schwarz.
#include "lib/prob_dens_fpt__absorb_unit_sphere_r0_eq_0.cpp"
#include "lib/rand_fpt__absorb_sphere_r0_eq_0.cpp"

////////////////////////////////////////
///// CLASS FOR BROWNIAN PARTICLES /////
////////////////////////////////////////

class particle
{
	long double x,y,z;

	public:

	particle(const double &chi, const double &V1, const double &V2);
	particle();
	void moveBC(const double &dt, const double &theta, const double &phi, const double &chi, const double &D1, const double &D2, const double &V1, const double &V2);
	double move(bool &continue_motion, const double &epsilon, const double &chi, const double &D1, const double &D2, const double &V1, const double &V2);
};

//Initialization of the particle to calculate the GMFPT, according to Ps(X).
particle::particle(const double &chi, const double &V1, const double &V2)
{
	double u=ran();
	double theta=acos(2*ran()-1);
	double phi=2*M_PI*ran();
	
	//partition function.
	static const double Z=exp(-V1)+(exp(-V2)-exp(-V1))*cube(chi);	
	
	double r=0;
	//in the inner-shell.
	if (u<exp(-V2)*cube(chi)/Z)
	{
		r=cbrt(Z*u/exp(-V2));
	}
	//in the outer-shell.
	else
	{
		r=cbrt(Z*u/exp(-V1)+(1-exp(-V2)/exp(-V1))*cube(chi));
	}
		
	x=r*cos(theta);
	y=r*sin(theta)*cos(phi);
	z=r*sin(theta)*sin(phi);
}

//Initialization of the particle to calculate t(0,0,0).
particle::particle()
{
	x=0.;
	y=0.;
	z=0.;
}

//Motion close to the interface between the two layers.
void particle::moveBC(const double &dt, const double &theta, const double &phi, const double &chi, const double &D1, const double &D2, const double &V1, const double &V2)
{
	static const double alpha=sqrt(D2)*exp(-V2)/(sqrt(D1)*exp(-V1)+sqrt(D2)*exp(-V2));
	double dx=0.,dy=0.,dz=0.;
	//motion throught the inner-shell(2).
	if (ran()<alpha)
	{
		do
		{
			double dr=-fabs(gsl_ran_ugaussian(GSL_r))*sqrt(2*D2*dt);
			double dtheta=gsl_ran_ugaussian(GSL_r)*sqrt(2*D2*dt);
			double dphi=gsl_ran_ugaussian(GSL_r)*sqrt(2*D2*dt);			
			dx=cos(theta)*dr          - sin(theta)*dtheta;
			dy=sin(theta)*cos(phi)*dr + cos(theta)*cos(phi)*dtheta - sin(theta)*sin(phi)*dphi;
			dz=sin(theta)*sin(phi)*dr + cos(theta)*sin(phi)*dtheta + sin(theta)*cos(phi)*dphi;
		}
		while(square(x+dx)+square(y+dy)+square(z+dz)>chi*chi);
	}
	//motion throught the outer-shell(1).
	else
	{
		do
		{
			double dr=fabs(gsl_ran_ugaussian(GSL_r))*sqrt(2*D1*dt);
			double dtheta=gsl_ran_ugaussian(GSL_r)*sqrt(2*D1*dt);
			double dphi=gsl_ran_ugaussian(GSL_r)*sqrt(2*D1*dt);			
			dx=cos(theta)*dr          - sin(theta)*dtheta;
			dy=sin(theta)*cos(phi)*dr + cos(theta)*cos(phi)*dtheta - sin(theta)*sin(phi)*dphi;
			dz=sin(theta)*sin(phi)*dr + cos(theta)*sin(phi)*dtheta + sin(theta)*cos(phi)*dphi;
		}
		while(square(x+dx)+square(y+dy)+square(z+dz)<chi*chi);
	}
	x+=dx;
	y+=dy;
	z+=dz;
}

//Motion of the particle
double particle::move(bool &continue_motion, const double &epsilon, const double &chi, const double &D1, const double &D2, const double &V1, const double &V2)
{
	static const double eta=1e-4; //minimal distance to boundaries.
	static const double min_rprotect=1e-3; //minimal protection radius for KMC method (close to external boundary).
	static const double dt=eta*eta/2./max(D1,D2); //time-step chosen wrt eta.
	
	//Position of the particle in spherical coordinates.
	const double r=sqrt(x*x+y*y+z*z);
	const double theta=acos(x/r);
	const double phi=atan2(z,y);
	
	//Distances to outer and inner boundaries.
	const double dist_outer=1.-r;
	const double dist_inner=fabs(chi-r);	
	
	//Absorption close to the escape region -> stop the motion.
	if (dist_outer<eta and fabs(theta)-epsilon/2.<eta)
	{
		continue_motion=false;
		return 0.;
	}
	//Protection + motion close to the outer boundary -> return the FPT according to KMC method (sphere sector).
	else if (dist_outer<eta)
	{
		//new position on a circular spherical arc.
		double r_protect=min(min_rprotect,2*sin(theta/2.-epsilon/4.));
		const double A=r_protect/2.;
		const double theta_rand=acos(A+(1-A)*ran());
		const double phi_rand=2*M_PI*ran();
		
		const double dr=-r_protect*cos(theta_rand);
		const double dtheta=r_protect*sin(theta_rand)*cos(phi_rand);
		const double dphi=r_protect*sin(theta_rand)*sin(phi_rand);

		x=cos(theta)          + cos(theta)*dr          - sin(theta)*dtheta;
		y=sin(theta)*cos(phi) + sin(theta)*cos(phi)*dr + cos(theta)*cos(phi)*dtheta - sin(theta)*sin(phi)*dphi;
		z=sin(theta)*sin(phi) + sin(theta)*sin(phi)*dr + cos(theta)*sin(phi)*dtheta + sin(theta)*cos(phi)*dphi;		
		
		//FPT according to the KMC library.
		return rand_fpt__absorb_sphere_r0_eq_0(D1,r_protect);	
	}
	//Motion close to the interface between the two shells (inner boundary) -> Langevin equation for the dynamics.
	else if (dist_inner<eta)
	{
		double dx=0.,dy=0.,dz=0.;
		double t=0.;
		//motion in inner-shell(2).
		if (r<chi)
		{
			dx=gsl_ran_ugaussian(GSL_r)*sqrt(2*D2*dt);
			dy=gsl_ran_ugaussian(GSL_r)*sqrt(2*D2*dt);
			dz=gsl_ran_ugaussian(GSL_r)*sqrt(2*D2*dt);
			//if crosses the interface.
			if (square(x+dx)+square(y+dy)+square(z+dz)>chi*chi)
			{
				//mean time to reach the interface when the distance is very small (and reached in one hop).
				t=square(r-chi)/(6*D2);
				//move the particle on the interface.
				x=chi*cos(theta);
				y=chi*sin(theta)*cos(phi);
				z=chi*sin(theta)*sin(phi);
				//motion starting on the interface (over time dt).
				moveBC(dt,theta,phi,chi,D1,D2,V1,V2);				
			}
			//else apply the dynamic.
			else
			{
				x+=dx;
				y+=dy;
				z+=dz;
			}
		
		}
		//motion in outer-shell(1).
		else
		{
			dx=gsl_ran_ugaussian(GSL_r)*sqrt(2*D1*dt);
			dy=gsl_ran_ugaussian(GSL_r)*sqrt(2*D1*dt);
			dz=gsl_ran_ugaussian(GSL_r)*sqrt(2*D1*dt);
			//if crosses the interface.
			if (square(x+dx)+square(y+dy)+square(z+dz)<chi*chi)
			{
				//mean time to reach the interface when the distance is very small (and reached in one hop).
				t=square(r-chi)/(6*D2);
				//move the particle on the interface.
				x=chi*cos(theta);
				y=chi*sin(theta)*cos(phi);
				z=chi*sin(theta)*sin(phi);
				//motion starting on the interface (over time dt).
				moveBC(dt,theta,phi,chi,D1,D2,V1,V2);
			}
			//else apply the dynamic.s
			else
			{
				x+=dx;
				y+=dy;
				z+=dz;
			}
		}
		//time of this hop: dt + MFPT before crossing the interface (=0 if no interface cross).
		return t+dt;
		
	}
	//Protection + movement far from boundaries -> return the FPT according to KMC method (disk).
	else
	{
		//new position on a sphere.
		const double r_protect=min(dist_outer,dist_inner);		
		const double theta_rand=acos(2*ran()-1);
		const double phi_rand=2*M_PI*ran();
		
		x+=r_protect*cos(theta_rand);
		y+=r_protect*sin(theta_rand)*cos(phi_rand);
		z+=r_protect*sin(theta_rand)*sin(phi_rand);
		
		//FPT according to the KMC library.
		if (r>chi)
		{
			return rand_fpt__absorb_sphere_r0_eq_0(D1,r_protect);
		}
		else
		{
			return rand_fpt__absorb_sphere_r0_eq_0(D2,r_protect);
		}
	}
}

///////////////////////////////////////
///// READ COMMAND LINE ARGUMENTS /////
///////////////////////////////////////

void ReadCommandLine(int argc, char** argv, int &Npart, double &epsilon, double &Delta, double &D1, double &D2, double &V1, double &V2, int &ran)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-Npart=" ))
		{
			Npart=atoi(argv[i]+7);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-Delta=" ))
		{
			Delta=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-D1=" ))
		{
			D1=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-D2=" ))
		{
			D2=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-V1=" ))
		{
			V1=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-V2=" ))
		{
			V2=atof(argv[i]+4);
		}		
		else if (strstr(argv[i], "-ran=" ))
		{
			ran=atoi(argv[i]+5);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: epsilon=escape angle, Delta=cortex width, D1=outer diffusion, D2=inner diffusion, V1=outer potential, V2=inner potential.
	double epsilon=0.5, Delta=0.25, D1=5., D2=1., V1=0., V2=2.;
	
	//Numerical parameter: Npart=number of particles, ran=random number generator index.
	int Npart=100000, ran=0;

	//Import parameters in command line.
	ReadCommandLine(argc,argv,Npart,epsilon,Delta,D1,D2,V1,V2,ran);

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index (initial time) = " << ran << endl;
	gsl_rng_set(GSL_r,ran);

	//Create the datafile.
	stringstream ss;
	int sys=system("mkdir -p ./data_NEP3d_KMC/");
	ss << "./data_NEP3d_KMC/NEP3d_KMC_epsilon=" << epsilon << "_Delta=" << Delta << "_D=" << D1/D2 << "_V=" << V2-V1 << "_ran=" << ran << ".txt";
	string nameData = ss.str();
	cout << nameData << endl;
	
	//Create the average value of GMFPT and CMFPT.
	double CMFPT=0., GMFPT=0.;

	for(int k=1;k<=Npart;k++)
	{
		//create the particle to calculate the GMFPT.
		particle partGMFPT(1-Delta,V1,V2);
		bool continue_motion=true;
		double gfpt=0;
		while (continue_motion)
		{
			gfpt+=partGMFPT.move(continue_motion,epsilon,1-Delta,D1,D2,V1,V2);
		}		
		GMFPT+=gfpt;				
		
		//create the particle to calculate the CMFPT.
		particle searcherCMFPT;
		continue_motion=true;
		double cfpt=0;
		while (continue_motion)
		{
			cfpt+=searcherCMFPT.move(continue_motion,epsilon,1-Delta,D1,D2,V1,V2);
		}
		CMFPT+=cfpt;

		//when 0.1% of the particles have reached the escape region.
		if (k%(Npart/1000)==0 or k==Npart)
		{
			cout << 100*double(k)/Npart << " % of MFPT obtained (" << k << " samples) -GMFPT=" << D1*GMFPT/k << " -CMFPT=" << D1*CMFPT/k << running_time.TimeRun(" ") << endl;

			ofstream fileData(nameData.c_str(),ios::trunc);
			fileData.precision(8);
			fileData << "#Npart\tGMFPT\tCMFPT" << endl;
			fileData << k << "\t" << D1*GMFPT/k << "\t" << D1*CMFPT/k << endl;
			fileData.close();
		}
	}
	return 0;
}
