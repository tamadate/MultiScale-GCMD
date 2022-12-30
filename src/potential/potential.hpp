#pragma once
#include "basePotential.hpp"
#include "AMBER/potentialAMBER.hpp"
#include "StillingerWeber/potentialSW.hpp"
#include "Tersoff/potentialTersoff.hpp"
#include "BMHFT/potentialBMHFT.hpp"

class PotentialLJ : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJ(){potName="LJ(gas-ion)";};
		~PotentialLJ(){};
};

class PotentialLJCoulHybrid : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJCoulHybrid(){potName="LJ+Coulombic-AA/CG(vapor-vapor)";};
		~PotentialLJCoulHybrid(){};
};

class PotentialLJHybrid : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJHybrid(){potName="LJ-AA/CG(gas-gas, gas-vapor)";};
		~PotentialLJHybrid(){};
};

class PotentialLJCoul : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialLJCoul(){potName="LJ+Coulombic(vapor-ion)";};
		~PotentialLJCoul(){};
};

class PotentialGasGas : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasGas(){potName="LJ gas-gas";};
		~PotentialGasGas(){};
};



class PotentialGasIntra : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasIntra(){potName="Gas intra";};
		~PotentialGasIntra(){};
};


class PotentialEfield : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		double Ecoeff[3];
		void compute(Variables *vars, FLAG *flags);
		PotentialEfield(){potName="Ion E-filed";};
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

class PotentialVaporIntra : public Potential {
	private:
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporIntra(){potName="AMBER vapor";};
		~PotentialVaporIntra(){};
};
