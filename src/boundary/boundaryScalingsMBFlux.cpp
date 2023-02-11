#include "../md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>

void
MD::boundary_scaling_gas_move(void){
	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *v=vars->velocity.data();
	
	for (auto i : vars->MolID[1]){
		int flag, flagx, flagy, flagz;
		double HL=vars->domainL*0.5;
		flag=flagx=flagy=flagz=0;
		double dx=x[i][0]-pre_ion[0];
		double dy=x[i][1]-pre_ion[1];
		double dz=x[i][2]-pre_ion[2];
		double dr2=dx*dx+dy*dy+dz*dz;

		if (dx < -HL) x[i][0] += vars->domainL, flagx--, flag++;
		if (dy < -HL) x[i][1] += vars->domainL, flagy--, flag++;
		if (dz < -HL) x[i][2] += vars->domainL, flagz--, flag++;
		if (dx > HL) x[i][0] -= vars->domainL, flagx++, flag++;
		if (dy > HL) x[i][1] -= vars->domainL, flagy++, flag++;
		if (dz > HL) x[i][2] -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
			double vx = v[i][0];
			double vy = v[i][1];
			double vz = v[i][2];
			double v2 = vx*vx+vy*vy+vz*vz;
			double vnew = sqrt(v2);
			double vMB = mbdist->vflux[mbdist->number]*1e-5;
			double mod_factor= vMB/vnew;
			v[i][0] = vx * mod_factor;
			v[i][1] = vy * mod_factor;
			v[i][2] = vz * mod_factor;
			mbdist->number++;
		}
	}
}

void
MD::boundary_scaling_vapor_move(void){
	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *v=vars->velocity.data();
	for (auto i : vars->MolID[2]){
		int flag, flagx, flagy, flagz;
		double HL=vars->domainL*0.5;
		flag=flagx=flagy=flagz=0;
		double dx=x[i][0]-pre_ion[0];
		double dy=x[i][1]-pre_ion[1];
		double dz=x[i][2]-pre_ion[2];
		double dr2=dx*dx+dy*dy+dz*dz;

		if (dx < -HL) x[i][0] += vars->domainL, flagx--, flag++;
		if (dy < -HL) x[i][1] += vars->domainL, flagy--, flag++;
		if (dz < -HL) x[i][2] += vars->domainL, flagz--, flag++;
		if (dx > HL) x[i][0] -= vars->domainL, flagx++, flag++;
		if (dy > HL) x[i][1] -= vars->domainL, flagy++, flag++;
		if (dz > HL) x[i][2] -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cgas,pp->mgas,T);}
			double vx = v[i][0];
			double vy = v[i][1];
			double vz = v[i][2];
			double v2 = vx*vx+vy*vy+vz*vz;
			double vnew = sqrt(v2);
			double vMB = mbdistV->vflux[mbdistV->number]*1e-5;
			double mod_factor= vMB/vnew;
			v[i][0] = vx * mod_factor;
			v[i][1] = vy * mod_factor;
			v[i][2] = vz * mod_factor;
			mbdistV->number++;
		}
	}
}
