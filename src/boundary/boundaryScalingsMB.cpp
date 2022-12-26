#include "../md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>

void
MD::boundary_scaling_ion_move(void){
	double HL=vars->domainL*0.5;
	int flag, flagx, flagy, flagz;
	random_device seed;
	mt19937 mt(seed());
	Molecule *mols = vars -> Molecules.data();

	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	for (auto i : vars->MolID[1]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		if (dx < -HL) mols[i].qx += vars->domainL, flagx--, flag++;
		if (dy < -HL) mols[i].qy += vars->domainL, flagy--, flag++;
		if (dz < -HL) mols[i].qz += vars->domainL, flagz--, flag++;
		if (dx > HL) mols[i].qx -= vars->domainL, flagx++, flag++;
		if (dy > HL) mols[i].qy -= vars->domainL, flagy++, flag++;
		if (dz > HL) mols[i].qz -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			mols[i].px = distgas(mt) *1e-5;
			mols[i].py = distgas(mt) *1e-5;
			mols[i].pz = distgas(mt) *1e-5;
		}
	}

	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	for (auto i : vars->MolID[2]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		if (mols[i].qx < mols[0].qx-HL) mols[i].qx += vars->domainL, flag++;
		if (mols[i].qy < mols[0].qy-HL) mols[i].qy += vars->domainL, flag++;
		if (mols[i].qz < mols[0].qz-HL) mols[i].qz += vars->domainL, flag++;
		if (mols[i].qx > mols[0].qx+HL) mols[i].qx -= vars->domainL, flag++;
		if (mols[i].qy > mols[0].qy+HL) mols[i].qy -= vars->domainL, flag++;
		if (mols[i].qz > mols[0].qz+HL) mols[i].qz -= vars->domainL, flag++;
		if (flag>0) {
			mols[i].px = distvapor(mt) *1e-5;
			mols[i].py = distvapor(mt) *1e-5;
			mols[i].pz = distvapor(mt) *1e-5;
		}
	}

}
