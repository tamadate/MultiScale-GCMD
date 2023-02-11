#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "flags.hpp"
#include "boundary/MBdist.hpp"
#include "thermostat/thermostat.hpp"

//------------------------------------------------------------------------

class MD {
	private:


	public:

		double startTime;
		int Nth;
		int calculation_number;

		int gastype;	/*1:He, 2:Ar, 3:N2*/
		long int step_relax;
		long int step_repre;
		long int Noftimestep;
		double p;
		double T;
		int Nof_around_gas;
		int Nof_around_vapor;
		int OBSERVE;
		double Ecoeff[3];

		double dt;
		double CUTOFF;
		double MARGIN;
		double ML2;
		double CL2;
		double V;
		void setCondition(char* condfile);
		void readCondFile(char* condfile);
		void read_initial(void);

		long int itime;
		std::vector<long int> collisionFlagGas;
		std::vector<long int> collisionFlagVapor;
		std::vector<Potential*> InterInter;
		std::vector<Potential*> IntraInter;
		void setPotential(FLAG *flags,int mode);

		Variables *vars;
		Observer *obs;
		Physical *pp;
		FLAG *flags;
		MBdist *mbdist;
		MBdist *mbdistV;
		Thermostat *thermo;

	//	General functions
		bool errorCheck(void);

	//	velocity verlet
		void run(char** argv);
		void verlet(void);
		void update_position(void);
		void velocity_calculation(void);
		void forceCombine(void);
		void updateWeightFunc(void);

	//	pair list
		void setIgnorePairs(void);
		bool checkOverlap(int i, int j);
		void make_pair(void);
		void make_pairLJ(void);
		void check_pairlist(void);
		void make_pairLJHybridSub(int i1, int j1,int k1, int i2,int j2, int k2);
		void make_pairLJHybridSelf(int i, int j, int k);
		void setCellIndex(void);
		void updateRegion(void);
		void setOrigin(void);

	//	initialization
		void initialization_gas(void);
		void initialization_vapor(void);
		void takeOver(void);

	//	periodic
		void boundary_scaling_gas_move(void);
		void boundary_scaling_ion_move(void);
		void boundary_scaling_vapor_move(void);
		void periodicAtoms(void);
		int loop, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
		double pre_ion[3];

	//	analysis (calculating position and velocity of center of mass)
		double gyration;

	//	export
		void exportDump(void);

	/*other*/
		void output(void);
		void output_temp(double gastemp, double iontemp);
		void output_initial(void);
		void output_gas_collision(long int initime);
		void output_vapor_collision(long int initime);
		void display(int output_ONOFF);
		char filepath[100];
		char atomFile[100];
		char vaporFile[100];
		char takeOverFile[100];
		char vaporStickFile[100];
		char filepath_gyration[100];
		void gyration_initial(void);
		void gyration_out(MD *md2);
		string gyration_path;
		double crsq;
		double origin[3];

		void velocity_scaling(void);
		void nosehoover_ion(void);
		void nosehoover_zeta(void);
		double zeta;
		void setNVE(void);
		void setNVTion(double temp);


		double del2,CD2,rmin2;
		double totalPotential;


		MD(char* condfile,int calcNumber);
		~MD(void) {
			delete vars;
			delete obs;
			delete pp;
			delete flags;
			delete mbdist;
			delete mbdistV;
		}
};




//------------------------------------------------------------------------
