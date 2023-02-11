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
	int Nsofar=vars->position.size();
	int Nvapor=vars->positionVapor.size();

  	// Maxwell-Boltzumann distribution generator
	random_device seed;
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-vars->domainL*0.5,vars->domainL*0.5);

  	// Sampling vapor molecules
	int i=0;
	do {
		double minDis=10000.0;
   		// sample random position (x, y, z)
		double x[3];
		x[0]=r(mt), x[1]=r(mt), x[2]=r(mt);
		int Ntotal=Nsofar+i*Nvapor;

		// calculate minimum distance from existing atoms (min_dis)
		// min distance from existing vapor molecule
		for(auto &c : vars->position) {
			double dx=x[0]-c[0];
			double dy=x[1]-c[1];
			double dz=x[2]-c[2];
			adjust_periodic(dx, dy, dz, vars->domainL);
			double d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<minDis) minDis=d; // minimum gas-ion distance
		}

 		// if min_dis is less than criteria, add the sampled molecule (a)
		if(minDis>dis){
			for (int j=0;j<Nvapor;j++){
				// Maxwell-Boltzumann distribution generator
				normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
				std::array<double,3> v;
				std::array<double,3> xx;
				std::array<double,3> f;
				xx[0]=vars->positionVapor[j][0]+x[0];
				xx[1]=vars->positionVapor[j][1]+x[1];
				xx[2]=vars->positionVapor[j][2]+x[2];
				v[0]=distvapor(mt)*1e-5;
				v[1]=distvapor(mt)*1e-5;
				v[2]=distvapor(mt)*1e-5;
				f[0]=f[1]=f[2]=0;
				int type=vars->typeVapor[j];
				vars->position.push_back(xx);
				vars->velocity.push_back(v);
				vars->force.push_back(f);
				vars->charge.push_back(vars->chargeVapor[j]);
				vars->mass.push_back(vars->atypes[type].mass);
				vars->type.push_back(type);
				vars->MolID[2].push_back(Ntotal+j);
				vars->weight.push_back(1);
				vars->region.push_back(AA);
				std::vector<int> pdum;
				vars->ignorePairs.push_back(pdum);
				vars->pairs.push_back(pdum);
			}

			for(auto b:vars->bonds_v) {
				b.atom1+=Ntotal;
				b.atom2+=Ntotal;
				vars->bonds.push_back(b);
				vars->ignorePairs[b.atom1].push_back(b.atom2);
				vars->ignorePairs[b.atom2].push_back(b.atom1);
			}
			for(auto b:vars->angles_v) {
				b.atom1+=Ntotal;
				b.atom2+=Ntotal;
				b.atom3+=Ntotal;
				vars->angles.push_back(b);
				vars->ignorePairs[b.atom1].push_back(b.atom3);
				vars->ignorePairs[b.atom3].push_back(b.atom1);
			}
			for(auto b:vars->dihedrals_v) {
				b.atom1+=Ntotal;
				b.atom2+=Ntotal;
				b.atom3+=Ntotal;
				b.atom4+=Ntotal;
				vars->dihedrals.push_back(b);
				vars->ignorePairs[b.atom1].push_back(b.atom4);
				vars->ignorePairs[b.atom4].push_back(b.atom1);
			}
			i++;
		}
		//collisionFlagVapor.push_back(0);
	} while(i<Nof_around_vapor);
}
