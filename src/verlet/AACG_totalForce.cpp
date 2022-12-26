#include "../md.hpp"

//	Calculate total at AA-CG overlapping region

void
MD::AACG_totalForce(void) {
	Molecule *mols = vars->Molecules.data();
	int Nmol=vars->Molecules.size();
	for(int i=0;i<Nmol;i++){
		if(vars->Region[i]==AACG){
			for(auto &a : mols[i].inAtoms){
				double massRatio=a.mass/mols[i].mass;
				a.fx+=mols[i].fx*massRatio;
				a.fy+=mols[i].fy*massRatio;
				a.fz+=mols[i].fz*massRatio;
			}
		}
	}
}
