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

void
MD:: updateIonCenter(void){
	Molecule *mols = vars->Molecules.data();
	mols[0].qx=mols[0].qy=mols[0].qz=0;
	mols[0].px=mols[0].py=mols[0].pz=0;
	for (auto &c : mols[0].inAtoms){
		double massRatio=c.mass/mols[0].mass;
		mols[0].qx += c.qx*massRatio;
		mols[0].qy += c.qy*massRatio;
		mols[0].qz += c.qz*massRatio;
		mols[0].px += c.px*massRatio;
		mols[0].py += c.py*massRatio;
		mols[0].pz += c.pz*massRatio;
	}
}
