//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------
void
MD::velocity_scaling(void) {
	obs->computeProps(vars,0);
	double Tp;
	Tp=300;
//	Tp=2500-2.2e-4*vars->time;
	for (auto &a : vars->Molecules[0].inAtoms){
		double ratio=sqrt(Tp/obs->Tin[0]);
		a.px*=ratio;
		a.py*=ratio;
		a.pz*=ratio;
	}
}

void
MD::setNVTion(double temp){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
	pp->Tnh_ion=temp;
}
