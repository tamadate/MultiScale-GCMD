#include "potential.hpp"

//	Force from electric field

void
PotentialEfield::compute(Variables *vars, FLAG *flags) {
	int Natom=vars->position.size();
    for (int i=0;i<Natom;i++) {
		double coeff=5.55e-5*vars->charge[i];
		vars->force[i][0]+=coeff*Ecoeff[0];
		vars->force[i][1]+=coeff*Ecoeff[1];
		vars->force[i][2]+=coeff*Ecoeff[2];
	}

}
