#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Randomly arange gas molecule around the center particle without overlapping.
  	- The velocity is picked from  the Maxwell-Boltzumann distribution
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_vapor(void) {
	double dis=15;	/*	minimum vapor-vapor, vapor-ion, vapor-gas distance */
	int Nsofar=vars->Molecules.size();

  	// Maxwell-Boltzumann distribution generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-vars->domainL*0.5,vars->domainL*0.5);

  	// Sampling vapor molecules
	int i=0;
	do {
 		Molecule a;
		double minDis=10000.0;
   		// sample random position (x, y, z)
		a.qx=r(mt), a.qy=r(mt), a.qz=r(mt);

		// calculate minimum distance from existing atoms (min_dis)
		// min distance from existing vapor molecules
		for(auto &b : vars->Molecules) {
			if(loop==0){
				for(auto &c : b.inAtoms) {
					double dx = a.qx - b.qx;
					double dy = a.qy - b.qy;
					double dz = a.qz - b.qz;
					adjust_periodic(dx, dy, dz, vars->domainL);
					double d=sqrt(dx*dx+dy*dy+dz*dz);
					if(d<minDis) minDis=d; // minimum vapor-vapor distance
				}
			}
			// min distance from existing gas & vapor molecules
			else{
				double dx=a.qx-b.qx;
				double dy=a.qy-b.qy;
				double dz=a.qz-b.qz;
				adjust_periodic(dx, dy, dz, vars->domainL);
				double d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<minDis) minDis=d; // minimum gas-gas distance
			}
			loop++;
		}

 		// if min_dis is less than criteria, add the sampled molecule (a)
		if(minDis>dis){
			a.px=distvapor(engine)*1e-5;
			a.py=distvapor(engine)*1e-5;
			a.pz=distvapor(engine)*1e-5;
			a.fx=a.fy=a.fz=0;
			a.mass=pp->Mvapor;
			a.type=2;
			a.id=i;
			a.bonds=vars->bonds_v;
			a.angles=vars->angles_v;
			a.dihedrals=vars->dihedrals_v;
			vars->Molecules.push_back(a);
			vars->MolID[2].push_back(i+Nsofar);
			i++;
		}
		//collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}
