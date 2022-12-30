//------------------------------------------------------------------------
#include "../md.hpp"
#include "thermostat.hpp"
//------------------------------------------------------------------------

void
MD::nosehoover_zeta(void){
	obs->computeProps(vars,0);
	int g=vars->Molecules[0].inAtoms.size()*3;
	double Q_inv = 0.0001;
	vars->zeta_ion += (obs->Tin[0] - pp->Tnh_ion)*g*kb_real*Q_inv*dt;
}


void
MD::nosehoover_ion(void){
	double Coeff=exp(-vars->zeta_ion*0.5*dt);
	for (auto &a : vars->Molecules[0].inAtoms) {
		a.px *= Coeff;
    a.py *= Coeff;
    a.pz *= Coeff;
	}
}

void
NoseHoover::NH_zeta(void){
	if(interval==loop){
		computeTnow();
		int g=1;
		for(auto i : vars->MolID[1]) if(vars->Region[i]==CG) g++;
		for(auto i : vars->MolID[2]) if(vars->Region[i]==CG) g++;
		g*=3;
		zeta += (Tnow - T)*g*kb_real*Q_inv*dt;
	}
}

void 
NoseHoover::Tcontrol(int LOOP){
	if(interval==loop){
		double Coeff=exp(-zeta*0.5*dt);
		for (int j=1;j<3;j++){
			for (auto i : vars->MolID[j]) {
				if(vars->Region[i]==CG){
					vars->Molecules[i].px *= Coeff;
					vars->Molecules[i].py *= Coeff;
					vars->Molecules[i].pz *= Coeff;
				}
			}
		}
		if(LOOP==1) loop=0;
	}
}

void
MD::setNVTion(double temp){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
	pp->Tnh_ion=temp;
}


