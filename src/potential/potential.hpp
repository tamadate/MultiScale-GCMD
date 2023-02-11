#pragma once
#include "basePotential.hpp"
#include "AMBER/potentialAMBER.hpp"
#include "StillingerWeber/potentialSW.hpp"
#include "Tersoff/potentialTersoff.hpp"
#include "BMHFT/potentialBMHFT.hpp"


class PotentialLJHybrid : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJHybrid(){potName="LJ-AA/CG(gas-gas, gas-vapor)";};
		~PotentialLJHybrid(){};
};




class PotentialEfield : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		double Ecoeff[3];
		void compute(Variables *vars, FLAG *flags);
		PotentialEfield(double Ex, double Ey, double Ez){
			Ecoeff[0]=Ex;
			Ecoeff[1]=Ey;
			Ecoeff[2]=Ez;
			potName="Ion E-filed";
		};
		~PotentialEfield(){};
};

class PotentialIonDipole : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		double alphagas;
		double zion;
		void compute(Variables *vars, FLAG *flags);
		PotentialIonDipole(){potName="Induced dipole ion-gas";};
		~PotentialIonDipole(){};
};

