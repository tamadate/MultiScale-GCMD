#include "../md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>

void
MD::boundary_scaling_gas_move(void){
	double HL=vars->domainL*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	Molecule *mols = vars -> Molecules.data();

	for (auto i : vars->MolID[1]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mols[i].w=weightFunc(dr2);

		if (dx < -HL) mols[i].qx += vars->domainL, flagx--, flag++;
		if (dy < -HL) mols[i].qy += vars->domainL, flagy--, flag++;
		if (dz < -HL) mols[i].qz += vars->domainL, flagz--, flag++;
		if (dx > HL) mols[i].qx -= vars->domainL, flagx++, flag++;
		if (dy > HL) mols[i].qy -= vars->domainL, flagy++, flag++;
		if (dz > HL) mols[i].qz -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
			vx = mols[i].px;
			vy = mols[i].py;
			vz = mols[i].pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdist->vflux[mbdist->number]*1e-5;
			mod_factor= vMB/v;
			mols[i].px = vx * mod_factor;
			mols[i].py = vy * mod_factor;
			mols[i].pz = vz * mod_factor;
			mbdist->number++;
		}
	}
}

void
MD::boundary_scaling_vapor_move(void){
	double HL=vars->domainL*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	Molecule *mols = vars -> Molecules.data();

	for (auto i : vars->MolID[2]){
		flag=flagx=flagy=flagz=0;
		double dx=mols[i].qx-mols[0].qx;
		double dy=mols[i].qy-mols[0].qy;
		double dz=mols[i].qz-mols[0].qz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mols[i].w=weightFunc(dr2);

		if (dx < -HL) mols[i].qx += vars->domainL, flagx--, flag++;
		if (dy < -HL) mols[i].qy += vars->domainL, flagy--, flag++;
		if (dz < -HL) mols[i].qz += vars->domainL, flagz--, flag++;
		if (dx > HL) mols[i].qx -= vars->domainL, flagx++, flag++;
		if (dy > HL) mols[i].qy -= vars->domainL, flagy++, flag++;
		if (dz > HL) mols[i].qz -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
			vx=mols[i].px;
			vy=mols[i].py;
			vz=mols[i].pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			mod_factor= vMB/v;
			mols[i].px = vx * mod_factor;
			mols[i].py = vy * mod_factor;
			mols[i].pz = vz * mod_factor;
			mbdistV->number++;
		}
	}
}
