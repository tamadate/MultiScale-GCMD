#include "../variables.hpp"

void
Variables::ionRotation(void){
  random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	for(auto &a : Molecules[0].inAtoms) {x=a.qx,y=a.qy,z=a.qz; ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z);}
    int is=Molecules[0].inAtoms.size();
    for(int i=0;i<is-1;i++) {
        for(int j=i+1;j<is;j++){
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