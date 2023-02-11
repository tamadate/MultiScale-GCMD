#include "../md.hpp"

/*########################################################################################

-----Periodic conditions-----

#######################################################################################*/

/**************************periodic**********************************/
void
adjust_periodic(double &dx, double &dy, double &dz, double d_size) {
	const double LH = d_size * 0.5;
	if (dx < -LH)dx += d_size;
	if (dx > LH) dx -= d_size;
	if (dy < -LH)dy += d_size;
	if (dy > LH) dy -= d_size;
	if (dz < -LH)dz += d_size;
	if (dz > LH) dz -= d_size;
}

/************************periodec condition for in gas***************************/


#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>

void
MD::periodicAtoms(void){
	double HL=vars->domainL*0.5;
	int flag;
	std::array<double,3> *x=vars->position.data();
	int Natom=vars->position.size();
	for (int i=0;i<Natom;i++){
		flag=0;
		double dx=x[i][0]-origin[0];
		double dy=x[i][1]-origin[1];
		double dz=x[i][2]-origin[2];
		if (dx < -HL) x[i][0] += vars->domainL;
		if (dy < -HL) x[i][1] += vars->domainL;
		if (dz < -HL) x[i][2] += vars->domainL;
		if (dx > HL) x[i][0] -= vars->domainL;
		if (dy > HL) x[i][1] -= vars->domainL;
		if (dz > HL) x[i][2] -= vars->domainL;
	}
}
