//------------------------------------------------------------------------
#include "../md.hpp"
#include "thermostat.hpp"
//------------------------------------------------------------------------

void
MD::nosehoover_zeta(void){
	obs->computeProps(vars,0);
	int g=vars->MolID[0].size()*3;
	double Q_inv = 0.0001;
	vars->zeta_ion += (obs->Tin[0] - pp->Tnh_ion)*g*kb_real*Q_inv*dt;
}


void
MD::nosehoover_ion(void){
	double Coeff=exp(-vars->zeta_ion*0.5*dt);
	std::array<double,3> *v=vars->velocity.data();
	for (auto i : vars->MolID[0]) {
		v[i][0] *= Coeff;
		v[i][1] *= Coeff;
		v[i][2] *= Coeff;
	}
}

void
NoseHoover::NH_zeta(void){
	if(interval==loop){
		computeTnow();
		int g=vars->position.size()+1;
		g*=3;
		zeta += (Tnow - T)*g*kb_real*Q_inv*dt;
	}
}

void 
NoseHoover::Tcontrol(int LOOP){
	if(interval==loop){
		double Coeff=exp(-zeta*0.5*dt);
		std::array<double,3> *v=vars->velocity.data();
		for(auto v :vars->velocity){
			v[0] *= Coeff;
			v[1] *= Coeff;
			v[2] *= Coeff;
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


