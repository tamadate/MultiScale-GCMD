//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::computeProps(Variables *vars,int molID){
	std::array<double,3> *v=vars->velocity.data();
	Kin[molID]=0;
	int Natom=vars->MolID[molID].size();
	for(auto i : vars->MolID[molID]){
		Kin[molID] += v[i][0] * v[i][0] * vars->mass[i];
		Kin[molID] += v[i][1] * v[i][1] * vars->mass[i];
		Kin[molID] += v[i][2] * v[i][2] * vars->mass[i];
	}
	Kin[molID]*= (0.5 * real_to_kcalmol);
	Tin[molID]=Kin[molID]/double(Natom)*coeff;
};



double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T) {
	double phi = 0.0;
	/*const int ps = pairs.size();
	Gas *gases = vars->gases.data();
	for (int k = 0; k < ps; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = gases[j].qx - gases[i].qx;
		double dy = gases[j].qy - gases[i].qy;
		double dz = gases[j].qz - gases[i].qz;
		adjust_periodic(dx, dy, dz);
		double r2 = (dx * dx + dy * dy + dz * dz);
		double r2inv= 1/r2;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (r2 < CL2){
			double r6inv = r2inv * r2inv * r2inv;
			phi += r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		}
	}*/
	phi = phi * Cpress;
	return  p/T*Treal + phi;
}
