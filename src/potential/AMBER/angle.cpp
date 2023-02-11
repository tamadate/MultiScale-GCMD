#include "potentialAMBER.hpp"

void
PotentialAMBER::computeAngle(Variables *vars, FLAG *flags) {
	std::array<double,3> *x = vars->position.data();
	std::array<double,3> *f = vars->force.data();
	Angle_type *ctypes = vars->ctypes.data();
	Angle *angles=vars->angles.data();
	int asize=vars->angles.size();
	for (int ian=0;ian<asize;ian++) {
		double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];
		int i=angles[ian].atom1;
		int j=angles[ian].atom2;
		int k=angles[ian].atom3;
		int type=angles[ian].type;
		dx1 = x[i][0] - x[j][0];
		dy1 = x[i][1] - x[j][1];
		dz1 = x[i][2] - x[j][2];
		adjust_periodic(dx1, dy1, dz1, vars->domainL);
		rsq1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
		r1 = sqrt(rsq1);
		dx2 = x[k][0] - x[j][0];
		dy2 = x[k][1] - x[j][1];
		dz2 = x[k][2] - x[j][2];
		adjust_periodic(dx2, dy2, dz2, vars->domainL);
		rsq2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
		r2 = sqrt(rsq2);
		C = dx1*dx2 + dy1*dy2 + dz1*dz2;
		C /= r1*r2;	// cos(theta)
		Cs = 1/(sqrt(1.0-C*C));	// 1/sin(theta)
		dtheta = acos(C) - ctypes[type].coeff[1]; // q-q0
		tk = ctypes[type].coeff[0] * dtheta; //k(q-q0)
		a = -2.0 * tk * Cs;
		a11 = a*C / rsq1;
		a12 = -a / (r1*r2);
		a22 = a*C / rsq2;
		f1[0] = a11*dx1 + a12*dx2;
		f1[1] = a11*dy1 + a12*dy2;
		f1[2] = a11*dz1 + a12*dz2;
		f3[0] = a22*dx2 + a12*dx1;
		f3[1] = a22*dy2 + a12*dy1;
		f3[2] = a22*dz2 + a12*dz1;
		f[i][0] += f1[0];
		f[i][1] += f1[1];
		f[i][2] += f1[2];
  		f[j][0] -= f1[0] + f3[0];
		f[j][1] -= f1[1] + f3[1];
		f[j][2] -= f1[2] + f3[2];
		f[k][0] += f3[0];
		f[k][1] += f3[1];
		f[k][2] += f3[2];
    if (flags->eflag) vars->Utotal.Uion+= tk*dtheta;
	}
}
