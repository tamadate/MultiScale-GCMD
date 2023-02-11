#include "../md.hpp"

/*########################################################################################

-----Initialization-----
intialization:
This intialize the calculation. Reading initial positions and setting initial velocities.
reintialization:
This makes connection between thermal relaxation and main diffusion coeficient calculation. It reset position, time, pair list and margine length.

#######################################################################################*/

/////////////////////////////////////////////////////////////////////
/*
	- Randomly aranging gas molecule around an ion with avoiding the
	overlapping.
	- The velocity is picked from  the Maxwell-Boltzumann
	distribution
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_gas(void) {
	double dis=15;	/*	minimum gas-gas, gas-ion distance */
	int Nsofar=vars->position.size();
	int Ngas=vars->positionGas.size();
	random_device seed;
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-vars->domainL*0.5,vars->domainL*0.5);

  	// Sampling gas molecules
	int i=0;
	do {
		int Ntotal=Nsofar+i*Ngas;
		double min_dis=10000.0;
		// sample random position (x, y, z)
		double x[3];
		x[0]=r(mt), x[1]=r(mt), x[2]=r(mt);

		// calculate minimum distance from existing atoms (min_dis)
		for(auto &c : vars->position) {
			double dx=x[0]-c[0];
			double dy=x[1]-c[1];
			double dz=x[2]-c[2];
			adjust_periodic(dx, dy, dz, vars->domainL);
			double d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_dis) min_dis=d; // minimum gas-ion distance
		}

		// if min_dis is less than criteria, add the sampled molecule (a)
		if(min_dis>dis){
			// sample velocities from Maxwell-Boltzumann distribution
			// set mass and id
			for (int j=0;j<Ngas;j++){
				// Maxwell-Boltzumann distribution generator
				normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
				std::array<double,3> v;
				std::array<double,3> xx;
				std::array<double,3> f;
				xx[0]=vars->positionGas[j][0]+x[0];
				xx[1]=vars->positionGas[j][1]+x[1];
				xx[2]=vars->positionGas[j][2]+x[2];
				v[0]=distgas(mt)*1e-5;
				v[1]=distgas(mt)*1e-5;
				v[2]=distgas(mt)*1e-5;
				f[0]=f[1]=f[2]=0;
				int type=0;
				
				vars->position.push_back(xx);
				vars->velocity.push_back(v);
				vars->force.push_back(f);
				vars->charge.push_back(0);
				vars->mass.push_back(vars->atypes[type].mass);
				vars->type.push_back(type);
				vars->MolID[1].push_back(Ntotal+j);
				vars->region.push_back(AA);
				vars->weight.push_back(1);
				std::vector<int> pdum;
				vars->ignorePairs.push_back(pdum);
				vars->pairs.push_back(pdum);
				
			}

			for(auto b:vars->bonds_g) {
				b.atom1+=Ntotal;
				b.atom2+=Ntotal;
				vars->bonds.push_back(b);
				vars->ignorePairs[b.atom1].push_back(b.atom2);
				vars->ignorePairs[b.atom2].push_back(b.atom1);
			}
			i++;
		}
	} while(i<Nof_around_gas);
}
