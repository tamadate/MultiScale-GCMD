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

	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *v=vars->velocity.data();
	std::array<double,3> *xmol=vars->molPosition.data();

	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	for (auto i : vars->group[1]){
		flag=flagx=flagy=flagz=0;
		double dx=xmol[i][0]-origin[0];
		double dy=xmol[i][1]-origin[1];
		double dz=xmol[i][2]-origin[2];
		if (dx < -HL) xmol[i][0] += vars->domainL, flag++;
		if (dy < -HL) xmol[i][1] += vars->domainL, flag++;
		if (dz < -HL) xmol[i][2] += vars->domainL, flag++;
		if (dx > HL) xmol[i][0] -= vars->domainL, flag++;
		if (dy > HL) xmol[i][1] -= vars->domainL, flag++;
		if (dz > HL) xmol[i][2] -= vars->domainL, flag++;
		if (flag>0) {
			v[i][0] = distgas(mt) *1e-5;
			v[i][1] = distgas(mt) *1e-5;
			v[i][2] = distgas(mt) *1e-5;
		}
	}

	for (auto i : vars->MolID[2]){
		flag=flagx=flagy=flagz=0;
		double dx=xmol[i][0]-origin[0];
		double dy=xmol[i][1]-origin[1];
		double dz=xmol[i][2]-origin[2];
		if (dx < -HL) xmol[i][0] += vars->domainL, flagx--, flag++;
		if (dy < -HL) xmol[i][1] += vars->domainL, flagy--, flag++;
		if (dz < -HL) xmol[i][2] += vars->domainL, flagz--, flag++;
		if (dx > HL) xmol[i][0] -= vars->domainL, flagx++, flag++;
		if (dy > HL) xmol[i][1] -= vars->domainL, flagy++, flag++;
		if (dz > HL) xmol[i][2] -= vars->domainL, flagz++, flag++;
		if (flag>0) {
			v[i][0] = distgas(mt) *1e-5;
			v[i][1] = distgas(mt) *1e-5;
			v[i][2] = distgas(mt) *1e-5;
		}
	}
}
