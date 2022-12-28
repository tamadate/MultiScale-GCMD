#include "../variables.hpp"

void
Variables::setCrossPotentials(void){
  int Natypes=atypes.size();
  pair_coeff.resize(Natypes);
  for (int i=0;i<Natypes;i++){
    pair_coeff[i].resize(Natypes);
    for (int j=0;j<Natypes;j++){
      pair_coeff[i][j].resize(2);
    }
  }
  for (int i=0;i<Natypes;i++){
    for (int j=0;j<Natypes;j++){
      double epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
      double sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
      //sigma=sqrt(atypes_v[i].coeff2*atypes_v[j].coeff2);
      pair_coeff[i][j][0]=48 * epu*pow(sigma,12.0);
      pair_coeff[i][j][1]=24 * epu*pow(sigma,6.0);
    }
  }

  double epu11=0.14397; // gas CG potential
  double sigma11=3.798; // gas CG potential
  double epu22=0.1521; // vapor CG potential
  double sigma22=3.15061; // vapor CG potential
  double epu12=sqrt(epu11*epu22);
  double sigma12=(sigma11+sigma11)*0.5;
  pair_coeff_CG[1][1][0]= 48 * epu11*pow(sigma11,12.0);
  pair_coeff_CG[1][1][1]= 24 * epu11*pow(sigma11,6.0);
  pair_coeff_CG[2][2][0]= 48 * epu22*pow(sigma22,12.0);
  pair_coeff_CG[2][2][1]= 24 * epu22*pow(sigma22,6.0);
  pair_coeff_CG[1][2][0]=pair_coeff_CG[2][1][0]= 48 * epu12*pow(sigma12,12.0);
  pair_coeff_CG[1][2][1]=pair_coeff_CG[2][1][1]= 24 * epu12*pow(sigma12,6.0);
}
