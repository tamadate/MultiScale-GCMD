
#include "../md.hpp"

//------------------------------------------------------------//
/*	Make gas-gas & gas-vapor interactions pair list */
//------------------------------------------------------------//
void
MD::make_pairLJHybrid(void){
	vars->pairsLJHybrid.clear();
	Molecule *mols=vars->Molecules.data();
	// make gas-vapor pair list
	for (auto i : vars->MolID[1]){
		for (auto j : vars->MolID[2]){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJHybrid.push_back(p);
			}
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

	//make_pairLJHybridTest();
}

void
MD::make_pairLJHybridSelf(int i, int j, int k){
	Molecule *mols=vars->Molecules.data();
	int Nmol=vars->cellIndex[i][j][k].size();
	if(Nmol>1){
		for (int i1=0;i1<Nmol-1;i1++){
			int I=vars->cellIndex[i][j][k][i1];
			for (int i2=i1+1;i2<Nmol;i2++){
				int J=vars->cellIndex[i][j][k][i2];
				double dx = mols[I].qx - mols[J].qx;
				double dy = mols[I].qy - mols[J].qy;
				double dz = mols[I].qz - mols[J].qz;
				double r2 = (dx * dx + dy * dy + dz * dz);
				if (r2 < ML2){
					Pair p;
					p.i=I;
					p.j=J;
					vars->pairsLJHybrid.push_back(p);
				}
			}
		}
	}
}

void
MD::make_pairLJHybridSub(int i1, int j1,int k1, int i2,int j2, int k2){
	Molecule *mols=vars->Molecules.data();
	if(i2==vars->cellIndexSize) i2=0;
	if(j2==vars->cellIndexSize) j2=0;
	if(k2==vars->cellIndexSize) k2=0;
	if(i1==vars->cellIndexSize) i1=0;
	if(j1==vars->cellIndexSize) j1=0;
	if(k1==vars->cellIndexSize) k1=0;
	for(auto I : vars->cellIndex[i1][j1][k1]){
		for(auto J : vars->cellIndex[i2][j2][k2]){
			double dx = mols[I].qx - mols[J].qx;
			double dy = mols[I].qy - mols[J].qy;
			double dz = mols[I].qz - mols[J].qz;
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=I;
				p.j=J;
				vars->pairsLJHybrid.push_back(p);
			}
		}
	}
}

void
MD::make_pairLJHybridTest(void){
	vars->pairsLJHybridTest.clear();
	Molecule *mols=vars->Molecules.data();

	// make gas-vapor pair list
	for (auto i : vars->MolID[1]){
		for (auto j : vars->MolID[2]){
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJHybridTest.push_back(p);
			}
		}
	}

	// make gas-gas pair list
	int gs=vars->MolID[1].size();
	for (int I=0; I<gs-1; I++){
		int i=vars->MolID[1][I];
		for (int J=I+1; J<gs; J++){
			int j=vars->MolID[1][J];
			double dx = mols[i].qx - mols[j].qx;
			double dy = mols[i].qy - mols[j].qy;
			double dz = mols[i].qz - mols[j].qz;
			adjust_periodic(dx, dy, dz, vars->domainL);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				vars->pairsLJHybridTest.push_back(p);
			}
		}
	}

	if(vars->pairsLJHybrid.size()!=vars->pairsLJHybridTest.size()){
		cout<<"pair list error"<<endl;
	}
	for(auto &p1 : vars->pairsLJHybrid){
		int i=p1.i;
		int j=p1.j;
		if(p1.i>p1.j) {
			i=p1.j;
			j=p1.i;
		}
		int flag=0;
		for(auto &p2 : vars->pairsLJHybridTest){
			if(p2.i==i&&p2.j==j){
				flag=1;
				break;
			}
		}
		if(flag==0){
			Molecule MOL1=vars->Molecules[p1.i];
			Molecule MOL2=vars->Molecules[p1.j];
			double dcellInv=1/vars->dcell;
			double HL=vars->domainL*0.5;
			double dx=vars->Molecules[p1.i].qx-vars->Molecules[0].qx+HL;
			double dy=vars->Molecules[p1.i].qy-vars->Molecules[0].qy+HL;
			double dz=vars->Molecules[p1.i].qz-vars->Molecules[0].qz+HL;
			int ix=int(dx*dcellInv);
			int iy=int(dy*dcellInv);
			int iz=int(dz*dcellInv);
			double dx2=vars->Molecules[p1.j].qx-vars->Molecules[0].qx+HL;
			double dy2=vars->Molecules[p1.j].qy-vars->Molecules[0].qy+HL;
			double dz2=vars->Molecules[p1.j].qz-vars->Molecules[0].qz+HL;
			int ix2=int(dx*dcellInv);
			int iy2=int(dy*dcellInv);
			int iz2=int(dz*dcellInv);
		}
		
	}
}
