#include "../variables.hpp"

//------------------------------------------------------------//
/*	Make ion-gas interaction pair list  */
//------------------------------------------------------------//

void
Variables::setCellIndex(double HL){
	for(auto &a : cellIndex){
		for(auto &b : a){
			for(auto &c : b){
				c.clear();
			}
		}
	}
	double dcellInv=1/dcell;
	for (int molID=1;molID<2;molID++){
		for (auto i : MolID[1]){
			double dx=Molecules[i].qx-Molecules[0].qx+HL;
			double dy=Molecules[i].qy-Molecules[0].qy+HL;
			double dz=Molecules[i].qz-Molecules[0].qz+HL;
			int ix=int(dx*dcellInv);
			int iy=int(dy*dcellInv);
			int iz=int(dz*dcellInv);
			cellIndex[ix][iy][iz].push_back(i);
		}
	}
}
