#pragma once
#include "constants.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	Variables(void);
	~Variables(void){};

	/*variables*/
  	int Nth;
	// length is number of atoms in domain
	std::vector<std::array<double,3>> position;
	std::vector<std::array<double,3>> velocity;
	std::vector<std::array<double,3>> force;
	std::vector<int> type;
	std::vector<double> charge;
	std::vector<double> mass;
	std::vector<double> weight;
	std::vector<std::vector<int>> pairs;
	std::vector<std::vector<int>> ignorePairs;
	std::vector<int> region;
	
	// length is number of molecules in domain
	std::vector<std::vector<int>> molecule;
	std::vector<std::array<double,3>> molPosition;

	// length is number of groups (typically 3: ion, gas, vapor)
	std::vector<std::vector<int>> group;
	std::vector<std::vector<int>> MolID;

	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;

	double time;
	double zeta_ion;
	double domainL;

	Potentials Utotal;
	Times times;

	void Uzero(void)	{
		Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
 	}

	void tzero(void)	{times.tion=times.tgas=times.tvap=times.tgi=times.tvv=times.tvg=times.tvi=times.tpair=0;}
	double Usum(void)	{return Utotal.Uion+Utotal.Ugas+Utotal.Uvap+Utotal.Ugi+Utotal.Ugg+Utotal.Uvg+Utotal.Uvi+Utotal.Uvv;}

	/*vectors for potential calculation*/
	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;

	std::vector<std::array<double,3>> positionVapor;
	std::vector<int> typeVapor;
	std::vector<double> massVapor;
	std::vector<double> chargeVapor;
	std::vector<Bond> bonds_v;
	std::vector<Angle> angles_v;
	std::vector<Dihedral> dihedrals_v;

	std::vector<std::array<double,3>> positionGas;
	std::vector<int> typeGas;
	std::vector<double> massGas;
	std::vector<double> chargeGas;
	std::vector<Bond> bonds_g;
	/*std::vector<Angle> angles_g;
	std::vector<Dihedral> dihedrals_g;*/
	void gasSetting(int gastype);

	std::vector<Pair> ion_pairs;
	std::vector<vector<vector<double>>> pair_coeff;
	std::vector<vector<vector<vector<int>>>> cellIndex;
	int cellIndexSize;
	double dcell;
	void initialCellIndex(double CLML);
	double bornCoeff[2][2][5];

	void setGasPotentials(void);
	void setBMHPotential(void);
	void setCrossPotentials(void);
	
	/*initialization and export to dump file*/
	void read_initial(char* ionFile, char* vaporFile);
	void readIonFile(char* infile);
	void readVaporFile(char* infile);
	void ionInitialVelocity(double T);
	void ionRotation(void);
	double totalPotential;
	double totalVirial;

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
	  X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
	  Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
	  Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
	}


private:
};
//------------------------------------------------------------------------
