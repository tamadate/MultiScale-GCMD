#include "potentialAMBER.hpp"

void
PotentialAMBER::computeBond(Variables *vars, FLAG *flags) {
	std::array<double,3> *x = vars->position.data();
	std::array<double,3> *f = vars->force.data();
	Bond_type *btypes = vars->btypes.data();
	int bsize=vars->bonds.size();
	for (int ib=0;ib<bsize;ib++) {
		int i=vars->bonds[ib].atom1;
		int j=vars->bonds[ib].atom2;
		int type=vars->bonds[ib].type;
		double dx = x[i][0] - x[j][0];
		double dy = x[i][1] - x[j][1];
		double dz = x[i][2] - x[j][2];
		adjust_periodic(dx, dy, dz, vars->domainL);
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r = sqrt(rsq);
		double dr = (r-btypes[type].coeff[1]);
		double rk = btypes[type].coeff[0] * dr;
		double force_bond_harmonic;
		force_bond_harmonic = -2.0*rk/r;
		f[i][0] += force_bond_harmonic * dx;
		f[i][1] += force_bond_harmonic * dy;
		f[i][2] += force_bond_harmonic * dz;
		f[j][0] -= force_bond_harmonic * dx;
		f[j][1] -= force_bond_harmonic * dy;
		f[j][2] -= force_bond_harmonic * dz;
		if(flags->eflag) vars->Utotal.Uion+=rk*dr;
	}
}
