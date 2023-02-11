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
	CUTOFF = 15.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	OBSERVE=10000000;
	T=300;
	p=1e5;

	readCondFile(condfile);	
	pp->readIonProp(atomFile);	// Get ion total mass & net charge
	pp->readVaporProp(vaporFile);	// Get vapor total mass & net charge
	pp->setPhysicalProp(gastype,T,p);	// Set physical properties (mass of gas molecule, gas viscosity, etc...)

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

	setOrigin();
	initialization_gas();	//Set initial positions & velocities for gas
 	initialization_vapor();	//Set initial positions & velocities for vapor
	mbdist -> makeWeightedMB(pp->cgas,pp->mgas,T);
	mbdistV -> makeWeightedMB(pp->cvapor,pp->mvapor,T);
	make_pair();
	vars->tzero();
}
