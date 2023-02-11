#include "potential.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on vapor-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialLJHybrid::compute(Variables *vars, FLAG *flags) {
	vars->times.tvv-=omp_get_wtime();
	int Natom=vars->position.size();
	std::array<double,3> *x=vars->position.data();
	std::array<double,3> *f=vars->force.data();
	for(int i=0;i<Natom;i++){
		for(auto j : vars->pairs[i]){
			int UNION=vars->region[i]&vars->region[j];	// Region[i] union Region[j]
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			//adjust_periodic(dx, dy, dz, vars->domainL);
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			int type1=vars->type[i];
			int type2=vars->type[j];
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_coul = qqrd2e * vars->charge[i] * vars->charge[j] * sqrt(r2inv);
			double w=1;
			if(vars->region[i]==AACG || vars->region[i]==AACG){
				double w1=vars->weight[i];
				double w2=vars->weight[j];
				w=w1*w2;
			}
			double force_pair = (force_lj + force_coul)*r2inv*w;
			f[i][0] += force_pair * dx;
			f[i][1] += force_pair * dy;
			f[i][2] += force_pair * dz;
			f[j][0] -= force_pair * dx;
			f[j][1] -= force_pair * dy;
			f[j][0] -= force_pair * dz;
			if(force_pair>1000) {
				double xx=0;
			}
			if(flags->eflag) {
				vars->Utotal.Uvv+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0)*w;
			}
		}
	}
	vars->times.tvv+=omp_get_wtime();
}