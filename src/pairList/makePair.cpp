#include "../md.hpp"

/////////////////////////////////////////////////////////////////////
/*
	make pair list
	- always gas-ion pair list was calculated
	- if gas-gas interaction is ON, make pair_gasgas
	- get max velocity of gas molecule vmax
	- set update loop = margine_length/(vmax*dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair(void){
	periodicAtoms();
	/*boundary_scaling_gas_move();
	boundary_scaling_vapor_move();
	pre_ion[0]=origin[0];
	pre_ion[1]=origin[1];
	pre_ion[2]=origin[2];
	boundary_scaling_ion_move();*/

	setCellIndex();
	make_pairLJ();
	updateRegion();
	// set number of steps for next pair list update
	double vmax2 = 0.0;
	for (auto &v : vars->velocity) {
		double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
		if (vmax2 < v2) vmax2 = v2 ;
	}
	double vmax = sqrt(vmax2);

	loop_update=MARGIN/(vmax*dt)*0.5;
	loop=0;
}
