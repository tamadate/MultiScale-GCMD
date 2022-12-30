#include "../md.hpp"

void
MD::readCondFile(char* condfile){
	ifstream stream(condfile);
	string str;
	int iflag=0;
	cout<<"**************************************************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"************  Calculation conditions  ************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"**************************************************"<<endl;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		std::vector<string> readings;
		istringstream stream(str);
		string reading;
		while(getline(stream,reading,'\t')) {
			readings.push_back(reading);
		}
		// Input file name for ion molecule
		if(readings[0]=="Input"){
			strcpy(atomFile,readings[1].c_str());
			cout<<"Atom file -->\t\t"<<readings[1]<<endl;
		}
		// Input file name for vapor molecule
		if(readings[0]=="Vapor"){
			strcpy(vaporFile,readings[1].c_str());
			Nof_around_vapor=stoi(readings[2]);
			cout<<"Vapor file -->\t\t"<<readings[1]<<endl;
			cout<<"Number of vapors\t"<<Nof_around_vapor<<endl;
		}
		// Restart condition file name (option)
		if(readings[0]=="TakeOver"){
			strcpy(takeOverFile,readings[1].c_str());
			flags->takeOver=1;
			cout<<"Take over from -->\t"<<readings[1]<<endl;
		}
		// Vapor sticking position file name
		if(readings[0]=="VaporStickPositions"){
			strcpy(vaporStickFile,readings[1].c_str());
			positionLogStep=stoi(readings[2]);
			ifstream stream2(vaporStickFile);
			string str2;
			while(getline(stream2,str2)) {
				if(str2.length()==0) continue;
				stickPositionList.push_back(stoi(str2)-1);
			}
			cout<<"Vapor stick position file -->\t\t"<<readings[1]<<endl;
		}
		// N of gas molecule and specie 
		if(readings[0]=="Gas"){
			Atom_type at;
			if(readings[1]=="He"){
				gastype=1;
				at.mass=4.027;
				at.name="He(g)";
				at.coeff1=0.0203;
				at.coeff2=2.556;
			}
			if(readings[1]=="N2"){
				gastype=2;
				at.mass=14.01;
				at.name="N(g)";
				at.coeff1=0.1098;
				at.coeff2=3.27351824993;
			}
			if(readings[1]=="N2monoatomic"){
				gastype=3;
				at.mass=28.02;
				at.name="N2(g)";
				at.coeff1=0.14397;
				at.coeff2=3.798;
			}
			if(readings[1]=="Ar"){
				gastype=4;
				at.mass=28.02;
				at.name="Ar(g)";
				at.coeff1=0.14397;
				at.coeff2=3.798;
			}
			vars->atypes.push_back(at);
			Nof_around_gas=stoi(readings[2]);
			// Push gas molecule potential parameters
			vars->pair_coeff.resize(1);	// Resize only for gas molecule
			double epu=at.coeff1;
			double sigma=at.coeff2;
			vars->pair_coeff.resize(1);
			vars->pair_coeff[0].resize(1);
			vars->pair_coeff[0][0].resize(2);
			vars->pair_coeff[0][0][0]=48 * epu*pow(sigma,12.0);
			vars->pair_coeff[0][0][1]=24 * epu*pow(sigma,6.0);
			cout<<"Gastype\t\t\t"<<readings[1]<<endl;
			cout<<"Number of gases\t\t"<<Nof_around_gas<<endl;
		}
		// Set temperature
		if(readings[0]=="Temperature"){
			T=stod(readings[1]);
			cout<<"Temperature\t\t"<<T<<" K"<<endl;
		}
		// Set pressure
		if(readings[0]=="Pressure"){
			p=stod(readings[1]);
			cout<<"Pressure\t\t"<<p<<" Pa"<<endl;
		}
		// Set time step
		if (readings[0]=="dt") {
			dt=stod(readings[1]);
			cout<<"Time step\t\t"<<dt<<" fs"<<endl;
		}
		// Set total time steps
		if (readings[0]=="TotalSteps") {
			Noftimestep=stod(readings[1]);
			cout<<"Total steps\t\t"<<float(Noftimestep)<<endl;
		}
		// Set time steps for thermal relaxation
		if (readings[0]=="RelaxSteps") {
			step_relax=stod(readings[1]);
			cout<<"Relax steps\t\t"<<float(step_relax)<<endl;
		}
		// Set cutoff diameter
		if (readings[0]=="CutOff") {
			CUTOFF=stod(readings[1]);
			cout<<"Cutoff\t\t\t"<<CUTOFF<<" ang."<<endl;
		}
		// Set margin size
		if (readings[0]=="Margin") {
			MARGIN=stod(readings[1]);
			cout<<"Margin size\t\t"<<MARGIN<<" ang."<<endl;
		}
		// Set output file name
		if (readings[0]=="Output") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dump";
			string tmp2=ss.str();
			pp->dump_path=new char[tmp2.length()+1];
			strcpy(pp->dump_path,tmp2.c_str());
			OBSERVE=stoi(readings[2]);
			FILE*f=fopen(pp->dump_path, "w");
			fclose(f);
			cout<<"Dump file -->\t\t"<<tmp2<<endl;
		}
		// Set emsemble of ion
		if (readings[0]=="NVTion") {
			if (readings[1]=="OFF") {
				flags->nose_hoover_ion=0;
				cout<<"Nose-Hoover for ion --> OFF"<<endl;
			}
			else if(readings[1]=="scale") {
				flags->nose_hoover_ion=0;
				flags->velocity_scaling=1;
				cout<<"Nose-Hoover for ion --> OFF\nVelocity scaling for ion --> ON"<<endl;
			}
			else {
				flags->nose_hoover_ion=1;
				pp->Tnh_ion=stod(readings[1]);
				cout<<"Nose-Hoover for ion --> ON --> "<<pp->Tnh_ion<<" K"<<endl;
			}
		}
		// Set emsemble of gas (can not use this option now)
		if (readings[0]=="NVTgas") {
			if (readings[1]=="OFF") {
				flags->nose_hoover_gas=0;
				cout<<"Nose-Hoover for gas --> OFF"<<endl;
			}
			else {
				flags->nose_hoover_gas=1;
				pp->Tnh_gas=stod(readings[1]);
				cout<<"Nose-Hoover for gas --> ON --> "<<pp->Tnh_gas<<" K"<<endl;
			}
		}
		// Set interactions
		if (readings[0]=="Interactions") {continue;}
		// Gas-Gas interaction
		if (readings[1]=="gg") {
			if (readings[2]=="LJ") flags->inter_gg=1;
			else if (readings[2]=="OFF") flags->inter_gg=0;
			else printf("**************Uknown gas gas parameter was found**************\n");
		}
		// Gas-Ion interaction
		if (readings[1]=="gi"||readings[1]=="ig") {
			if (readings[2]=="LJ") flags->force_lj=1;
			else if (readings[2]=="ion dipole") flags->force_ion_dipole=1;
			else if (readings[2]=="OFF") {
				flags->force_ion_dipole=0;
				flags->force_lj=0;
			}
			else printf("**************Uknown gas ion parameter was found**************\n");
		}
		// Ion intramolecular interaction
		if (readings[1]=="ion") {
			if (readings[2]=="AMBER") flags->intra_AMBER=1;
			else if (readings[2]=="Stilinger-Weber") flags->force_sw=1;
			else if (readings[2]=="Tersoff") flags->force_ters=1;
			else if (readings[2]=="BMH") flags->force_born=1;
			else printf("**************Uknown ion parameter was found**************\n");
		}
		// Vapor-Ion interaction
		if (readings[1]=="vi"||readings[1]=="iv") {
			if (readings[2]=="LJcoul") flags->inter_vi=1;
			else if (readings[2]=="OFF") flags->inter_vi=0;
			else printf("**************Uknown vapor ion parameter was found**************\n");
		}
		// Vapor-Vapor interaction
		if (readings[1]=="vv"||readings[1]=="vv") {
			if (readings[2]=="LJcoul") flags->inter_vv=1;
			else if (readings[2]=="OFF") flags->inter_vv=0;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		// Gas-Vapor interaction
		if (readings[1]=="gv"||readings[1]=="vg") {
			if (readings[2]=="LJ") flags->inter_vg=1;
			else if (readings[2]=="OFF") flags->inter_vg=0;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		// Electric field
		if (readings[1]=="Efield") {
			flags->efield=1;
			Ecoeff[0]=stod(readings[2]);
			Ecoeff[1]=stod(readings[3]);
			Ecoeff[2]=stod(readings[4]);
		}
		// Gyration radius of ion
		if (readings[0]=="Gyration") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dat";
			string tmp=ss.str();
			pp->gyration_path=new char[tmp.length()+10];
			strcpy(pp->gyration_path,tmp.c_str());
			flags->gyration=1;
			FILE*f=fopen(pp->gyration_path, "w");
			fclose(f);
			cout<<"Gyration --> ON -->\t"<<tmp<<endl;
		}
	}
	vars->domainL=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;  // Calculation box size
	vars->dcell=vars->domainL/double(vars->cellIndexSize);
	V=vars->domainL*vars->domainL*vars->domainL;		// Calculation box volume
	CL2 = (CUTOFF)*(CUTOFF);
	ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
	cout<<"Cut off length\t\t"<<CUTOFF<<" ang."<<endl;
	cout<<"Margin length\t\t"<<MARGIN<<" ang."<<endl;
	cout<<"Domain size\t\t"<<vars->domainL<<" ang."<<endl;
}
