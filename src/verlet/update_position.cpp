#include "../md.hpp"

//  Update position 
void
MD::update_position(void) {
	vars->times.tpos-=omp_get_wtime();
    int Nmol=vars->Molecules.size();
    Molecule *mols = vars -> Molecules.data();

	for (int i=0;i<Nmol;i++){
		if(vars->Region[i]==CG){
			mols[i].qx += mols[i].px * dt;
			mols[i].qy += mols[i].py * dt;
			mols[i].qz += mols[i].pz * dt;
		}
		else{
			for (auto &a : mols[i].inAtoms){
				a.qx += a.px * dt;
				a.qy += a.py * dt;
				a.qz += a.pz * dt;
				a.fx=a.fy=a.fz=0.0;
			}
		}
		mols[i].fx=mols[i].fy=mols[i].fz=0.0;
	}
	updateInCenters();
	boundary_scaling_gas_move();
	boundary_scaling_vapor_move();
    vars->times.tpos+=omp_get_wtime();
}
