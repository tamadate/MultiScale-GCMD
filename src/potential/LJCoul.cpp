//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-vapor (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialLJCoul::compute(Variables *vars, FLAG *flags) {
	Molecule *mols = vars->Molecules.data();
	Atom *ions = vars->Molecules[0].inAtoms.data();
	vars->times.tvi-=omp_get_wtime();
	for(auto &p : vars->pairsLJCoul){
		int i=p.i;
		int j=p.j;
		for (auto &at1 : mols[i].inAtoms){
			for (auto &at2 : mols[j].inAtoms){
				double dx = at1.qx - at2.qx;
				double dy = at1.qy - at2.qy;
				double dz = at1.qz - at2.qz;
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=at1.type;
				int type2=at2.type;
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * at1.charge * at2.charge * sqrt(r2inv);
				double force_pair = (force_lj + force_coul)*r2inv;
				if(force_pair>100) {
					double xx=0;
				}
				at1.fx += force_pair * dx;
				at1.fy += force_pair * dy;
				at1.fz += force_pair * dz;
				at2.fx -= force_pair * dx;
				at2.fy -= force_pair * dy;
				at2.fz -= force_pair * dz;
				if(flags->eflag) {
					vars->Utotal.Uvi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					vars->Utotal.Uvi+=force_coul;
				}
			}
		}
	}
	vars->times.tvi+=omp_get_wtime();
}
