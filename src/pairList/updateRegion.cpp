#include "../md.hpp"

//------------------------------------------------------------//
/*	Make update regions  */
//------------------------------------------------------------//
void
MD::updateRegion(void){
	std::array<double,3> *x=vars->position.data();
	int Natom=vars->position.size();
	for (int i=0;i<Natom;i++){
		double dx = x[i][0] - origin[0];
		double dy = x[i][1] - origin[1];
		double dz = x[i][2] - origin[2];
		double dr2=dx*dx+dy*dy+dz*dz;
		if(dr2<RO2) vars->region[i]=AA;
		else {
			if(dr2<RI2) vars->region[i]=AACG;
			else vars->region[i]=CG;
		}
	}
}
