#include "../md.hpp"

//------------------------------------------------------------//
/*	check necessity of the pair list updating	*/
//------------------------------------------------------------//
void
MD::check_pairlist(void){
	vars->times.tpair-=omp_get_wtime();
	loop++;
	if(loop>loop_update){
		make_pair();
	}
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
	vars->times.tpair+=omp_get_wtime();
}
