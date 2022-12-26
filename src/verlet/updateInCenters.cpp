#include "../md.hpp"

void
MD:: updateInCenters(void){
	Molecule *mols = vars->Molecules.data();
	int Nmol=vars->Molecules.size();
	for (int i=0;i<Nmol;++i){
		if(vars->Region[i]==CG) continue;
		mols[i].qx=mols[i].qy=mols[i].qz=0;
		mols[i].px=mols[i].py=mols[i].pz=0;
		for (auto &c : mols[i].inAtoms){
			double massRatio=c.mass/mols[i].mass;
			mols[i].qx += c.qx*massRatio;
			mols[i].qy += c.qy*massRatio;
			mols[i].qz += c.qz*massRatio;
			mols[i].px += c.px*massRatio;
			mols[i].py += c.py*massRatio;
			mols[i].pz += c.pz*massRatio;
		}
	}
}
