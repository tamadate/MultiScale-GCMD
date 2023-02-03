#include "potential.hpp"

//	Force from electric field

void
PotentialEfield::compute(Variables *vars, FLAG *flags) {
    for (auto &a : vars->Molecules[0].inAtoms) {
		a.fx+=5.55e-5*a.charge*Ecoeff[0];
		a.fy+=5.55e-5*a.charge*Ecoeff[1];
		a.fz+=5.55e-5*a.charge*Ecoeff[2];
	}

}
