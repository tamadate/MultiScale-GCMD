#include "../md.hpp"


void
MD::read_initial(void) {
	// Add ion molecule to Molecules vector
	vars->readIonFile(atomFile);
	vars->readVaporFile(vaporFile);
	vars->setCrossPotentials(); // LJ potential parameters (BL rule)
	vars->ionRotation();
}
