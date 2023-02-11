#include "../md.hpp"

//	- Update velocity (half of a time step, dt/2)

void
MD::velocity_calculation(void) {
	vars->times.tvel-=omp_get_wtime();
	double const Coeff=0.5*dt*4.184e-4;
	const int Nmol=vars->position.size();
	//omp_set_num_threads(4);
	#pragma omp parallel for
	for (int i=0; i<Nmol; ++i){
		double Coeff2=Coeff/vars->mass[i];
		vars->velocity[i][0] += vars->force[i][0] *Coeff2;
		vars->velocity[i][1] += vars->force[i][1] *Coeff2;
		vars->velocity[i][2] += vars->force[i][2] *Coeff2;
	}
	vars->times.tvel+=omp_get_wtime();
}
