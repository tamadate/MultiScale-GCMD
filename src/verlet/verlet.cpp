//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Compute a domain (NVE or NVT)
	- v(t) -> v(t+dt/2)
	- Calculate velocity and potition of ion's center of mass.
	- Temperature control (velocity scaling, this case is NVT)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatmic interaction
	- Calculate ion-gas interatominc interaction
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
void
MD::verlet(void) {
	velocity_calculation(); //	v(t) -> v(t+dt/2) using F(x(t))
	//if(flags->velocity_scaling==1)	velocity_scaling();

	update_position();

	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	check_pairlist();
	vars->totalVirial=0;
	for (auto &inter : InterInter) inter->compute(vars,flags);
	for (auto &intra : IntraInter) intra->compute(vars,flags);
	AACG_totalForce();

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	//(flags->nose_hoover_gas==1)	nosehoover_gas();
	flags->eflag=0;
}