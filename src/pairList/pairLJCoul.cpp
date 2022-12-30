#include "../md.hpp"

//------------------------------------------------------------//
/*	Make ion-gas interaction pair list  */
//------------------------------------------------------------//

void
MD::make_pairLJCoul(void){
	Molecule *mols=vars->Molecules.data();
	for (auto i : vars->MolID[2]){
		if(vars->Region[i]==CG) continue;
		Pair p;
		p.i=i;
		p.j=0;
		vars->pairsLJCoul.push_back(p);
	}
}
