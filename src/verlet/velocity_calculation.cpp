#include "../md.hpp"

//	- Update velocity (half of a time step, dt/2)

void
MD::velocity_calculation(void) {
	vars->times.tvel-=omp_get_wtime();
	double const Coeff=0.5*dt*4.184e-4;
	int i=0;
	for (auto &mol : vars->Molecules){
		if(vars->Region[i]==CG){
			double Coeff2=Coeff/mol.mass;
			mol.px += mol.fx *Coeff2;
			mol.py += mol.fy *Coeff2;
			mol.pz += mol.fz *Coeff2;
			const int th=1;
			if(mol.px*mol.px>th || mol.py*mol.py>th || mol.pz*mol.pz>th) {
				double xx=0;
			}
		}
		else{
			for (auto &a : mol.inAtoms){
				double Coeff2=Coeff/a.mass;
				a.px += a.fx *Coeff2;
				a.py += a.fy *Coeff2;
				a.pz += a.fz *Coeff2;
				const int th=10;
				if(a.px*a.px>th || a.py*a.py>th || a.pz*a.pz>th){
					double xx=0;
				}
			}
		}
		i++;
	}
	vars->times.tvel+=omp_get_wtime();
}
