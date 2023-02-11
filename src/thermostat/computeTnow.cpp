#include "thermostat.hpp"

void
Thermostat::computeTnow(void){
	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *v=vars->velocity.data();
	double *m=vars->mass.data();
	int Natom=vars->position.size();
	double K=0;
	for(int i=1; i<Natom; i++){
		K += v[i][0] * v[i][0] * m[i];
		K += v[i][1] * v[i][1] * m[i];
		K += v[i][2] * v[i][2] * m[i];
	}
	K*= (0.5 * real_to_kcalmol);
	Tnow=K/double(Natom)*kb_real_inv/1.5;
};
