#include "../md.hpp"

void
MD::setCondition(char* condfile){
	readCondFile(condfile);	
	pp->readIonProp(atomFile);	// Get ion total mass & net charge
	pp->readVaporProp(vaporFile);	// Get vapor total mass & net charge
	pp->setPhysicalProp(gastype,T,p);	// Set physical properties (mass of gas molecule, gas viscosity, etc...)
}
