#include "../variables.hpp"

void
Variables::ionRotation(void){
  random_device seed;
	double A,B,C;
	A=seed(),B=seed(),C=seed();
	for(auto i:MolID[0]) {
        double x=position[i][0];
        double y=position[i][1];
        double z=position[i][2];
        ROTATION(position[i][0],position[i][1],position[i][2],A,B,C,x,y,z);
    }
    int Nion=MolID[0].size();
    for(int i=0;i<Nion-1;i++) {
        for(int j=i+1;j<Nion;j++){
            int flag=0;
            for (auto &d : dihedrals) {
                int I=d.atom1;
                int J=d.atom2;
                int K=d.atom3;
                int L=d.atom4;
                if(i==I){if(j==J||j==K||j==L) flag=1;}
                if(i==J){if(j==I||j==K||j==L) flag=1;}
                if(i==K){if(j==J||j==I||j==L) flag=1;}
                if(i==L){if(j==J||j==K||j==I) flag=1;}
            }
            if (flag==0){
                Pair p;
                p.i=i;
                p.j=j;
                ion_pairs.push_back(p);
            }
        }
    }
}