#include "../md.hpp"

void
MD::setIgnorePairs(void){
	for (auto &b : vars->bonds) {
		int i=b.atom1;
		int j=b.atom2;
		vars->ignorePairs[i].push_back(j);
        vars->ignorePairs[j].push_back(i);
	}
	for (auto &b : vars->angles) {
		int i=b.atom1;
		int j=b.atom3;
		vars->ignorePairs[i].push_back(j);
        vars->ignorePairs[i].push_back(i);
	}
	for (auto &b : vars->dihedrals) {
		int i=b.atom1;
		int j=b.atom4;
		vars->ignorePairs[i].push_back(j);
        vars->ignorePairs[j].push_back(i);
	}
}
