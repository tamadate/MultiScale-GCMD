#include "potentialAMBER.hpp"


void
PotentialAMBER::computeLong(Variables *vars, FLAG *flags) {
	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *f=vars->force.data();
	int lpsize=longPair.size();
	for (int ip=0;ip<lpsize;ip++) {
		int i=longPair[ip].i;
		int j=longPair[ip].j;
		double dx = x[i][0] - x[j][0];
		double dy = x[i][1] - x[j][1];
		double dz = x[i][2] - x[j][2];
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=vars->type[i];
		int type2=vars->type[j];
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		double force_coul = qqrd2e * vars->charge[i] * vars->charge[j] * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		f[i][0] += force_pair * dx;
		f[i][1] += force_pair * dy;
		f[i][2] += force_pair * dz;
		f[j][0] -= force_pair * dx;
		f[j][1] -= force_pair * dy;
		f[j][2] -= force_pair * dz;
		if(flags->eflag) {
			vars->Utotal.Uion+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			vars->Utotal.Uion+=force_coul;
		}
	}
}
