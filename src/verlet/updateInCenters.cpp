#include "../md.hpp"

void
MD:: updateInCenters(Molecule &mol){
	//Molecule *mols = vars->Molecules.data();
	mol.qx=mol.qy=mol.qz=0;
	mol.px=mol.py=mol.pz=0;
	for (auto &c : mol.inAtoms){
		double massRatio=c.mass/mol.mass;
		mol.qx += c.qx*massRatio;
		mol.qy += c.qy*massRatio;
		mol.qz += c.qz*massRatio;
		mol.px += c.px*massRatio;
		mol.py += c.py*massRatio;
		mol.pz += c.pz*massRatio;
	}
}

void
MD:: updateInCentersAll(void){
	Molecule *mols = vars->Molecules.data();
	int Nmol=vars->Molecules.size();
	for (int i=0;i<Nmol;++i){
		if(vars->Region[i]==CG) continue;
		updateInCenters(mols[i]);
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
