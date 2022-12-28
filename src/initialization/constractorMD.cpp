//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, int calcNumber) {
	startTime=omp_get_wtime();
	calculation_number =  calcNumber;

	vars = new Variables();
	obs = new Observer();
	pp = new Physical();
	flags = new FLAG();
	mbdist = new MBdist();
	mbdistV = new MBdist();

	// Default calculation parameters
	dt = 0.5;	/*	fs	*/
	CUTOFF = 20.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	OBSERVE=10000000;
	T=300;
	p=1e5;
	positionLogStep=0;
	totalVaporIn=1;

	setCondition(condfile);
	//thermo = new NoseHoover(vars,obs,dt,T,10000);
	thermo = new velocityScaling(vars,obs,T);
	output_initial();	
	read_initial();	// Read initial ion structure and vapor structure
	vars->ionInitialVelocity(T);
	setPotential(flags,1);

	// This function is currently not working
	if(flags->takeOver==1){
		takeOver();
		//cout<<takeOverFile<<endl;
	}

	initialization_gas();	//Set initial positions & velocities for gas
 	initialization_vapor();	//Set initial positions & velocities for vapor
	setRegion();
	updateInCenters();
	boundary_scaling_ion_move();
	make_pair();
	margin_length = MARGIN;
	vars->tzero();

	mbdist -> makeWeightedMB(pp->cgas,pp->mgas,T);
	mbdistV -> makeWeightedMB(pp->cvapor,pp->mvapor,T);

}
