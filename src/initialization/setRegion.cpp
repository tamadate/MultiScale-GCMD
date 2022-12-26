#include "../md.hpp"

void
MD::setRegion(void) {
  Molecule *mols=vars->Molecules.data();
  int Nmol=vars->Molecules.size();
  vars->Region.resize(Nmol);
  // center particle is always AA
  vars->Region[0]=AA;
  // set gas molecule region
  for (auto i : vars->MolID[1]){
    double dx=mols[i].qx-mols[0].qx;
    double dy=mols[i].qy-mols[0].qy;
    double dz=mols[i].qz-mols[0].qz;
    double dr2=dx*dx+dy*dy+dz*dz;
		if(dr2<RO2) {
      vars->Region[i]=AA;
      makeDiatomicProp_in(mols[i]);
    }
		else {
			if(dr2<RI2) {
				vars->Region[i]=AACG;
				makeDiatomicProp_in(mols[i]);
			}
			else vars->Region[i]=CG;
		}
  }
  // set vapor molecule region
  for (auto i : vars->MolID[2]){
    double dx=mols[i].qx-mols[0].qx;
    double dy=mols[i].qy-mols[0].qy;
    double dz=mols[i].qz-mols[0].qz;
    double dr2=dx*dx+dy*dy+dz*dz;
		if(dr2<RO2) {
      vars->Region[i]=AA;
      makePolyatomicProp_in(mols[i]);
    }
		else {
			if(dr2<RI2) {
				vars->Region[i]=AACG;
				makePolyatomicProp_in(mols[i]);
			}
			else vars->Region[i]=CG;
		}
  }
}