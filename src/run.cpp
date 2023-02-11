//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	- Diffusion coefficient calculation
	- MD simulation performed for (step_relax) steps as a thermal
	relaxation, and main calculation run for (Noftimestep) steps.
	- Among of two calculations, a calculation is re-initialized.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
*/
/////////////////////////////////////////////////////////////////////
void
MD::run(char** argv) {
	/****Thermal relaxation****/
  	const int logger=10000;
	setPotential(flags,0);
	for (auto &a : IntraInter) {a->printName();}
	for (auto &a : InterInter) {a->printName();}
	for (itime=0; itime < step_relax; itime++) {
		if (itime%OBSERVE==0) {
			display(1);
			exportDump();
      		//exportDumpOut(); // secret command
			vars->time=(itime*dt);
		}
    	verlet();
	}

/*
****Reinitialization****
- Reset time (time -> 0).
- Re-initializating the ion and gas positions and velocities.
Position -> 0 elta, Traslational velocity -> 0
- Set ion's center of mass (maybe -> 0), make pair list for initial
step of simulation, reset the margine size.
*/
	vars->time=0;
	itime=0;
	setPotential(flags,1);
	make_pair();
	setNVE();

	/****Main simulaiton****/
	for (itime=0; itime<Noftimestep; itime++) {
		if(itime%logger==0){
		output();
		vars->time+=(logger*dt);
		}
		if (itime%OBSERVE==0) {
			display(0);
			exportDump();
			//exportDumpOut();
		}
		if ((itime+1)%OBSERVE==0) {
			flags->eflag=1;
			vars->Uzero();
		}
		verlet();
	}
}
