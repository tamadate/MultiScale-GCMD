//------------------------------------------------------------------------
#include "potential.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-ion (Coulombic+LJ)
*/
/////////////////////////////////////////////////////////////////////
void
Potential::compute(Variables *vars, FLAG *flags){
	/*
	Atom *ions = vars->Molecules[0].inAtoms.data();
	vars->times.tion-=omp_get_wtime();
	int ipsize=vars->ion_pairs.size();
	for (int ip=0;ip<ipsize;ip++){
		int i=vars->ion_pairs[ip].i;
		int j=vars->ion_pairs[ip].j;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=ions[i].type;
		int type2=ions[j].type;
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		double force_coul = qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		ions[i].fx += force_pair * dx;
		ions[i].fy += force_pair * dy;
		ions[i].fz += force_pair * dz;
		ions[j].fx -= force_pair * dx;
		ions[j].fy -= force_pair * dy;
		ions[j].fz -= force_pair * dz;
		if(flags->eflag) {
			vars->Utotal.Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->Utotal.Uion+=force_coul;
		}
		if(flags->vflag){
			vars->totalVirial+=force_lj+force_coul;
		}
	}
	vars->times.tion+=omp_get_wtime();*/
}

