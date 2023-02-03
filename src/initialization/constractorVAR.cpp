#include "../variables.hpp"

Variables::Variables(void) {
  time = 0.0;
  Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
  MolID.resize(3);
  initialCellIndex();
}


void
Variables::initialCellIndex(void){
  cellIndexSize=5;
  cellIndex.resize(cellIndexSize);
  for(auto &a : cellIndex) {
    a.resize(cellIndexSize);
    for(auto &b : a) {
      b.resize(cellIndexSize);
    }
  }
}


void
Variables::setBMHPotential(void){
  //difine Born-Mayer-Huggins conefficients for "only" NaCl
  //https://doi.org/10.1063/1.1522375
  //https://doi.org/10.1016/0022-3697(64)90160-X
  //0:A/rho, 1:6C, 2:8D, 3:sigma, 4:1/rho
  //NaNa, A=25.4435kJ/mol, C=101.1719kJ/mol, D=48.1771kJ/mol, sigma=2.340A, 1/rho=3.1546A-1
  //NaCl, A=20.3548kJ/mol, C=674.4793kJ/mol, D=837.077kJ/mol, sigma=2.755A, 1/rho=3.1546A-1
  //ClCl, A=15.2661kJ/mol, C=6985.6786kJ/mol, D=14031.5785kJ/mol, sigma=3.170A, 1/rho=3.1546A-1
  	bornCoeff[0][0][0]=25.4435*3.1546/4.184;
  	bornCoeff[0][1][0]=bornCoeff[1][0][0]=20.3548*3.1546/4.184;
  	bornCoeff[1][1][0]=15.2661*3.1546/4.184;

  	bornCoeff[0][0][1]=6*101.1719/4.184;
  	bornCoeff[0][1][1]=bornCoeff[1][0][1]=6*674.4793/4.184;
  	bornCoeff[1][1][1]=6*6985.6786/4.184;

  	bornCoeff[0][0][2]=8*48.1771/4.184;
  	bornCoeff[0][1][2]=bornCoeff[1][0][2]=8*837.077/4.184;
  	bornCoeff[1][1][2]=8*14031.5785/4.184;

  	bornCoeff[0][0][3]=2.340;
  	bornCoeff[0][1][3]=bornCoeff[1][0][3]=2.755;
  	bornCoeff[1][1][3]=3.170;

  	bornCoeff[0][0][4]=3.1546;
  	bornCoeff[0][1][4]=bornCoeff[1][0][4]=3.1546;
  	bornCoeff[1][1][4]=3.1546;
}

void
Variables::ionInitialVelocity(double T) {
	random_device seed;
	default_random_engine engine(seed());
  for(auto &a : Molecules[0].inAtoms) {
    double matom=a.mass*1e-3/Nw;
    normal_distribution<> dist(0.0, sqrt(kb*T/matom));
    a.px=dist(engine)*1e-5;
    a.py=dist(engine)*1e-5;
    a.pz=dist(engine)*1e-5;
	}
}

