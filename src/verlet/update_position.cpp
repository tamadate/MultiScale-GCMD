#include "../md.hpp"

//  Update position 
void
MD::update_position(void) {
	vars->times.tpos-=omp_get_wtime();
    int Nmol=vars->position.size();
	//omp_set_num_threads(omp_get_num_threads());
	#pragma omp parallel for
	for (int i=0;i<Nmol;i++){
		vars->position[i][0] += vars->velocity[i][0] * dt;
		vars->position[i][1] += vars->velocity[i][1] * dt;
		vars->position[i][2] += vars->velocity[i][2] * dt;
		vars->force[i][0]=vars->force[i][1]=vars->force[i][2]=0.0;
	}
	setOrigin();
	vars->times.tpos+=omp_get_wtime();
}
