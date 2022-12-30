#include "thermostat.hpp"

void
Thermostat::computeTnow(void){
	double K=0;
	double Nout=0;
	int Nmol=vars->Molecules.size();
	Molecule *mols=vars->Molecules.data();

	for(int i=1; i<Nmol; i++){
		if(vars->Region[i]==CG){
			K += mols[i].px * mols[i].px * mols[i].mass;
			K += mols[i].py * mols[i].py * mols[i].mass;
			K += mols[i].pz * mols[i].pz * mols[i].mass;
			Nout++;
		}
	}
	K*= (0.5 * real_to_kcalmol);
	Tnow=K/double(Nout)*kb_real_inv/1.5;
};
