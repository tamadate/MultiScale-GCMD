#include "../md.hpp"

void
MD::setCondition(char* condfile){
	readCondFile(condfile);	
	pp->readIonProp(atomFile);	// Get ion total mass & net charge
	// Add ion molecule to Molecules vector
	Molecule mol;
	mol.qx=mol.qy=mol.qz=0;
	mol.px=mol.py=mol.pz=0;
	mol.fx=mol.fy=mol.fz=0;
	mol.type=-1;
	mol.mass=pp->Mion;
	vars->Molecules.push_back(mol);
	pp->readVaporProp(vaporFile);	// Get vapor total mass & net charge
	pp->setPhysicalProp(gastype,T,p);	// Set physical properties (mass of gas molecule, gas viscosity, etc...)
}
