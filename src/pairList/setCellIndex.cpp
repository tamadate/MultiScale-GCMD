#include "../md.hpp"

//------------------------------------------------------------//
/*	Make ion-gas interaction pair list  */
//------------------------------------------------------------//

void
MD::setCellIndex(void){
	for(auto &a : vars->cellIndex){
		for(auto &b : a){
			for(auto &c : b){
				c.clear();
			}
		}
	}
	double dcellInv=1/vars->dcell;
	double HL=vars->domainL*0.5;
	for (int molID=1;molID<2;molID++){
		for (auto i : vars->MolID[1]){
			double dx=vars->position[i][0]-origin[0]+HL;
			double dy=vars->position[i][1]-origin[1]+HL;
			double dz=vars->position[i][2]-origin[2]+HL;
			int ix=int(dx*dcellInv);
			int iy=int(dy*dcellInv);
			int iz=int(dz*dcellInv);
			vars->cellIndex[ix][iy][iz].push_back(i);
		}
	}
}
