#include "../md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>

void
MD::boundary_scaling_gas_move(void){
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	#pragma omp parallel for
	for (auto i : vars->MolID[1]){
		double ionx=pre_ion[0];
		double iony=pre_ion[1];
		double ionz=pre_ion[2];
		int flag, flagx, flagy, flagz;
		double HL=vars->domainL*0.5;
		Molecule &mol=vars->Molecules[i];
		flag=flagx=flagy=flagz=0;
		double dx=mol.qx-ionx;
		double dy=mol.qy-iony;
		double dz=mol.qz-ionz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mol.w=weightFunc(dr2);

		if (dx < -HL) mol.qx += vars->domainL, flagx--, flag++;
		if (dy < -HL) mol.qy += vars->domainL, flagy--, flag++;
		if (dz < -HL) mol.qz += vars->domainL, flagz--, flag++;
		if (dx > HL) mol.qx -= vars->domainL, flagx++, flag++;
		if (dy > HL) mol.qy -= vars->domainL, flagy++, flag++;
		if (dz > HL) mol.qz -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
			vx = mol.px;
			vy = mol.py;
			vz = mol.pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdist->vflux[mbdist->number]*1e-5;
			mod_factor= vMB/v;
			mol.px = vx * mod_factor;
			mol.py = vy * mod_factor;
			mol.pz = vz * mod_factor;
			mbdist->number++;
		}
	}
}

void
MD::boundary_scaling_vapor_move(void){
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	#pragma omp parallel for
	for (auto i : vars->MolID[2]){
		double ionx=pre_ion[0];
		double iony=pre_ion[1];
		double ionz=pre_ion[2];
		int flag, flagx, flagy, flagz;
		double HL=vars->domainL*0.5;
		Molecule &mol=vars->Molecules[i];
		flag=flagx=flagy=flagz=0;
		double dx=mol.qx-ionx;
		double dy=mol.qy-iony;
		double dz=mol.qz-ionz;
		double dr2=dx*dx+dy*dy+dz*dz;
		mol.w=weightFunc(dr2);

		if (dx < -HL) mol.qx += vars->domainL, flagx--, flag++;
		if (dy < -HL) mol.qy += vars->domainL, flagy--, flag++;
		if (dz < -HL) mol.qz += vars->domainL, flagz--, flag++;
		if (dx > HL) mol.qx -= vars->domainL, flagx++, flag++;
		if (dy > HL) mol.qy -= vars->domainL, flagy++, flag++;
		if (dz > HL) mol.qz -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
			vx=mol.px;
			vy=mol.py;
			vz=mol.pz;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			mod_factor= vMB/v;
			mol.px = vx * mod_factor;
			mol.py = vy * mod_factor;
			mol.pz = vz * mod_factor;
			mbdistV->number++;
		}
	}
}
