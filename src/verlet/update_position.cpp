#include "../md.hpp"

//  Update position 
void
MD::update_position(void) {
	vars->times.tpos-=omp_get_wtime();
    int Nmol=vars->Molecules.size();
	//omp_set_num_threads(omp_get_num_threads());
	#pragma omp parallel for
	for (int i=0;i<Nmol;i++){
		Molecule &mol=vars->Molecules[i];
		if(vars->Region[i]==CG){
			mol.qx += mol.px * dt;
			mol.qy += mol.py * dt;
			mol.qz += mol.pz * dt;
		}
		else{
			for (auto &a : mol.inAtoms){
				a.qx += a.px * dt;
				a.qy += a.py * dt;
				a.qz += a.pz * dt;
				a.fx=a.fy=a.fz=0.0;
			}
			updateInCenters(mol);
		}
		mol.fx=mol.fy=mol.fz=0.0;
	}
	//boundary_scaling_gas_move();
	//boundary_scaling_vapor_move();
	vars->times.tpos+=omp_get_wtime();
}
