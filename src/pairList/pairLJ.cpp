#include "../md.hpp"

//------------------------------------------------------------//
/*	Make ion-gas interaction pair list  */
//------------------------------------------------------------//

void
MD::make_pairLJ(void){
	for(auto &pair:vars->pairs) pair.clear();
	std::array<double,3> *x=vars->position.data();
	// make gas-vapor pair list
	for (auto i : vars->MolID[1]){
		for (auto j : vars->MolID[2]){
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2) vars->pairs[i].push_back(j);
		}
	}
	// make gas-gas pair list
	for (int i=0;i<vars->cellIndexSize;i++){
		for (int j=0;j<vars->cellIndexSize;j++){
			for (int k=0;k<vars->cellIndexSize;k++){
				make_pairLJHybridSelf(i,j,k);

				make_pairLJHybridSub(i,j,k,i+1,j,k);
				make_pairLJHybridSub(i,j,k,i+1,j+1,k);
				make_pairLJHybridSub(i,j,k,i+1,j+1,k+1);
				make_pairLJHybridSub(i,j,k,i,j+1,k);
				make_pairLJHybridSub(i,j,k,i,j+1,k+1);
				make_pairLJHybridSub(i,j,k,i,j,k+1);
				make_pairLJHybridSub(i,j,k,i+1,j,k+1);

				make_pairLJHybridSub(i+1,j,k,i,j+1,k+1);
				make_pairLJHybridSub(i,j+1,k,i+1,j,k+1);
				make_pairLJHybridSub(i,j,k+1,i+1,j+1,k);

				make_pairLJHybridSub(i,j,k+1,i+1,j,k);
				make_pairLJHybridSub(i,j+1,k,i,j,k+1);
				make_pairLJHybridSub(i,j+1,k,i+1,j,k);
			}
		}
	}
	// make ion-gas pair list
	for (auto i : vars->MolID[1]){
		if(vars->region[i]==CG) continue;
		for(auto j : vars->MolID[0]){
			double dx=x[i][0]-x[j][0];
			double dy=x[i][1]-x[j][1];
			double dz=x[i][2]-x[j][2];
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2) {
				vars->pairs[i].push_back(j);
			}
		}
	}
	// make ion-vapor pair list
	for (auto i : vars->MolID[2]){
		if(vars->region[i]==CG) continue;
		for(auto j : vars->MolID[0]){
			vars->pairs[i].push_back(j);
		}
	}
	// make vapor-vapor pair list
	int vs=vars->MolID[2].size();
	for (int I=0; I<vs-1; I++){
		int i=vars->MolID[2][I];
		if(vars->region[i]==CG) continue;
		for (int J=I+1; J<vs; J++){
			int j=vars->MolID[2][J];
			if(vars->region[j]==CG) continue;
			if(checkOverlap(i,j)) vars->pairs[i].push_back(j);
		}
	}
	// make ion-ion pair list
	int is=vars->MolID[0].size();
	for (int I=0; I<is-1; I++){
		int i=vars->MolID[0][I];
		for (int J=I+1; J<is; J++){
			int j=vars->MolID[0][J];
			if(checkOverlap(i,j)) vars->pairs[i].push_back(j);
		}
	}
}


void
MD::make_pairLJHybridSelf(int i, int j, int k){
	std::array<double,3> *pos=vars->position.data();
	int Nmol=vars->cellIndex[i][j][k].size();
	if(Nmol>1){
		for (int i1=0;i1<Nmol-1;i1++){
			int I=vars->cellIndex[i][j][k][i1];
			for (int i2=i1+1;i2<Nmol;i2++){
				int J=vars->cellIndex[i][j][k][i2];
				double dx = pos[I][0] - pos[J][0];
				double dy = pos[I][1] - pos[J][1];
				double dz = pos[I][2] - pos[J][2];
				double r2 = (dx * dx + dy * dy + dz * dz);
				if (r2 < ML2){
					if(checkOverlap(I,J)) vars->pairs[I].push_back(J);
				}
			}
		}
	}
}

void
MD::make_pairLJHybridSub(int i1, int j1,int k1, int i2,int j2, int k2){
	std::array<double,3> *x=vars->position.data();
	if(i2==vars->cellIndexSize) i2=0;
	if(j2==vars->cellIndexSize) j2=0;
	if(k2==vars->cellIndexSize) k2=0;
	if(i1==vars->cellIndexSize) i1=0;
	if(j1==vars->cellIndexSize) j1=0;
	if(k1==vars->cellIndexSize) k1=0;
	for(auto I : vars->cellIndex[i1][j1][k1]){
		for(auto J : vars->cellIndex[i2][j2][k2]){
			double dx = x[I][0] - x[J][0];
			double dy = x[I][1] - x[J][1];
			double dz = x[I][2] - x[J][2];
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				if(checkOverlap(I,J)) vars->pairs[I].push_back(J);
			}
		}
	}
}

bool
MD::checkOverlap(int i, int j){
	bool flag=true;
	for (auto J : vars->ignorePairs[i]){
		if (J==j) {
			flag=false; 
			break;
		}
	}
	return flag;
}
