
#include "thermostat.hpp"

void
velocityScaling::Tcontrol(int LOOP) {
	if(interval==loop){
		obs->computeProps(vars,1);
		double Coeff=sqrt(T/obs->Tout[1]);
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


