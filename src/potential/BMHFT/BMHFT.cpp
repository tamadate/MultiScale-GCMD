#include "potentialBMHFT.hpp"
/*########################################################################################

Compute ion intramolecular interaction _ Born-Mayer-Huggins-Fumi-Tosi
- Currently, NaCl is only avairable particle

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialBorn::compute(Variables *vars, FLAG *flags) {
	Atom *ions = vars->Molecules[0].inAtoms.data();
	const int is = vars->Molecules[0].inAtoms.size();
	vars->times.tion-=omp_get_wtime();
	for(int i=0; i<is-1; i++){
		for(int j=i+1; j<is; j++){
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			double r=sqrt(rsq);
			int type1=ions[i].type;
			int type2=ions[j].type;

			//vars->bornCoeff[type1][type2][0]=A
			//vars->bornCoeff[type1][type2][1]=6C
			//vars->bornCoeff[type1][type2][2]=8D
			//vars->bornCoeff[type1][type2][3]=sigma
			//vars->bornCoeff[type1][type2][4]=1/rho

			double rexp=exp(pairCoeff[type1][type2][4]*(pairCoeff[type1][type2][3]-r));
			double force1 = pairCoeff[type1][type2][4]*pairCoeff[type1][type2][0]*r*rexp;
			double force2 = -pairCoeff[type1][type2][1]*r6inv;
			double force3 = -pairCoeff[type1][type2][2]*r6inv*r2inv;
			double force_coul = qqrd2e*ions[i].charge*ions[j].charge*sqrt(r2inv);
			double force_pair = (force1+force2+force3+force_coul)*r2inv;

			ions[i].fx += force_pair * dx;
			ions[i].fy += force_pair * dy;
			ions[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Uion+=force_coul;
				vars->Utotal.Uion+=rexp*pairCoeff[type1][type2][0];
				vars->Utotal.Uion-=pairCoeff[type1][type2][1]/6.0*r6inv;
				vars->Utotal.Uion-=pairCoeff[type1][type2][2]/8.0*r6inv*r2inv;
			}
			//vars->totalVirial+=force_lj;
		}
	}
	vars->times.tion+=omp_get_wtime();
}

PotentialBorn::PotentialBorn(Variables *vars){
	//difine Born-Mayer-Huggins conefficients for "only" NaCl
	//https://doi.org/10.1063/1.1522375
	//https://doi.org/10.1016/0022-3697(64)90160-X
	//0:A/rho, 1:6C, 2:8D, 3:sigma, 4:1/rho
	//NaNa, A=25.4435kJ/mol, C=101.1719kJ/mol, D=48.1771kJ/mol, sigma=2.340A, 1/rho=3.1546A-1
	//NaCl, A=20.3548kJ/mol, C=674.4793kJ/mol, D=837.077kJ/mol, sigma=2.755A, 1/rho=3.1546A-1
	//ClCl, A=15.2661kJ/mol, C=6985.6786kJ/mol, D=14031.5785kJ/mol, sigma=3.170A, 1/rho=3.1546A-1
  	double NaNa0=25.4435*3.1546/4.184;
  	double NaCl0=20.3548*3.1546/4.184;
  	double ClCl0=15.2661*3.1546/4.184;

  	double NaNa1=6*101.1719/4.184;
  	double NaCl1=6*674.4793/4.184;
  	double ClCl1=6*6985.6786/4.184;

  	double NaNa2=8*48.1771/4.184;
  	double NaCl2=8*837.077/4.184;
  	double ClCl2=8*14031.5785/4.184;

  	double NaNa3=2.340;
  	double NaCl3=2.755;
  	double ClCl3=3.170;

  	double NaNa4=3.1546;
  	double NaCl4=3.1546;
  	double ClCl4=3.1546;

	int Natypes=vars->atypes.size();
	pairCoeff.resize(Natypes);
	for (int i=0;i<Natypes;i++){
		pairCoeff[i].resize(Natypes);
		for (int j=0;j<Natypes;j++){
			pairCoeff[i][j].resize(5);
		}
	}
	for (int i=0;i<Natypes;i++){
		for (int j=0;j<Natypes;j++){
			if(vars->atypes[i].name=="Na+" && vars->atypes[j].name=="Na+"){
				pairCoeff[i][j][0]=NaNa0;
				pairCoeff[i][j][1]=NaNa1;
				pairCoeff[i][j][2]=NaNa2;
				pairCoeff[i][j][3]=NaNa3;
				pairCoeff[i][j][4]=NaNa4;
			}
			if(vars->atypes[i].name=="Na+" && vars->atypes[j].name=="Cl-"){
				pairCoeff[i][j][0]=pairCoeff[j][i][0]=NaCl0;
				pairCoeff[i][j][1]=pairCoeff[j][i][1]=NaCl1;
				pairCoeff[i][j][2]=pairCoeff[j][i][2]=NaCl2;
				pairCoeff[i][j][3]=pairCoeff[j][i][3]=NaCl3;
				pairCoeff[i][j][4]=pairCoeff[j][i][4]=NaCl4;
			}
			if(vars->atypes[i].name=="Cl-" && vars->atypes[j].name=="Cl-"){
				pairCoeff[i][j][0]=ClCl0;
				pairCoeff[i][j][1]=ClCl1;
				pairCoeff[i][j][2]=ClCl2;
				pairCoeff[i][j][3]=ClCl3;
				pairCoeff[i][j][4]=ClCl4;
			}
		}
	}
	potName="Ion BMH";
}
