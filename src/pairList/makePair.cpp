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
	boundary_scaling_gas_move();
	boundary_scaling_vapor_move();
	pre_ion[0]=vars->Molecules[0].qx;
	pre_ion[1]=vars->Molecules[0].qy;
	pre_ion[2]=vars->Molecules[0].qz;
	boundary_scaling_ion_move();

	updateInOut();

	vars->setCellIndex(vars->domainL*0.5);
	vars->pairsLJ.clear();
	vars->pairsLJCoul.clear();
	vars->pairsLJHybrid.clear();
	vars->pairsLJCoulHybrid.clear();
	make_pairLJ();
 	make_pairLJCoul();
	make_pairLJHybrid();
	make_pairsLJCoulHybrid();
	// set number of steps for next pair list update
	double vmax2 = 0.0;
	for (auto i : vars->MolID[1]) {
		double px=vars->Molecules[i].px;
		double py=vars->Molecules[i].py;
		double pz=vars->Molecules[i].pz;
		double v2 = px*px + py*py + pz*pz;
		if (vmax2 < v2) vmax2 = v2 ;
	}
	//double vion2=vars->Molecules[0][0].px*vars->Molecules[0][0].px+vars->Molecules[0][0].py*vars->Molecules[0][0].py+vars->Molecules[0][0].pz*vars->Molecules[0][0].pz;
	//if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2);

	loop_update=MARGIN/(vmax*dt)*0.5;
	loop=0;
}
