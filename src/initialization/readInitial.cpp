#include "../md.hpp"


void
MD::read_initial(void) {
  // Add ion molecule to Molecules vector
	Molecule mol;
	mol.qx=mol.qy=mol.qz=0;
	mol.px=mol.py=mol.pz=0;
	mol.fx=mol.fy=mol.fz=0;
	mol.type=-1;
	mol.mass=pp->Mion;
	mol.w=1;
	vars->MolID[0].push_back(0);
	vars->Molecules.push_back(mol);
  
  vars->readIonFile(atomFile);
  vars->readVaporFile(vaporFile);
  vars->setBMHPotential();
  vars->setCrossPotentials(); // LJ potential parameters (BL rule)
  vars->ionRotation();
}
