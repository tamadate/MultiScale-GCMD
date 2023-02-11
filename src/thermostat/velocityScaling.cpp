#include "thermostat.hpp"

void
velocityScaling::Tcontrol(int LOOP) {
	if(interval==loop){
		computeTnow();
		double Coeff=sqrt(T/Tnow);
		std::array<double,3> *v=vars->velocity.data();
		for(auto v :vars->velocity){
			v[0] *= Coeff;
			v[1] *= Coeff;
			v[2] *= Coeff;
		}
		if(LOOP==1) loop=0;
	}
}


