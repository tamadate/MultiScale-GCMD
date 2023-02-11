#include "../md.hpp"

void
MD::updateWeightFunc(void){
  int Natom=vars->position.size();
  std::array<double,3> *x=vars->position.data();
  for(int i=0;i<Natom;i++){
    if(vars->region[i]==CG) vars->weight[i]=0;
    if(vars->region[i]==AA) vars->weight[i]=1;
    if(vars->region[i]==AACG) {
      double dx=x[i][0]-origin[0];
      double dy=x[i][1]-origin[1];
      double dz=x[i][2]-origin[2];
      double dr2=dx*dx+dy*dy+dz*dz;
      double weight=0;
      if(dr2<RCG2) {
        weight=1;
      } else if(dr2>RAA2) {
        weight=0;
      } else {
        double dr=(sqrt(dr2)-RCG);
        double C=cos(M_PI/BoundL*0.25*dr);
        weight=C*C;
      }
      vars->weight[i]=weight;
    }
  }
}
