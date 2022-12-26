#include "../md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::output_initial(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "gas_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "K_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "U_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
}


void
MD::output(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *ions = vars->Molecules.data();
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, ions[0].qx, ions[0].qy, ions[0].qz, ions[0].px, ions[0].py, ions[0].pz);
	fclose(f);
}


void
MD::output_gas_collision(long int initime){
	sprintf(filepath, "gas_collision_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%e %e\n", initime*dt, itime*dt);
	fclose(f);
}

void
MD::output_vapor_collision(long int initime){
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%e %e\n", initime*dt, itime*dt);
	fclose(f);
}

void
MD::Ovin(int i){
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *mols=vars->Molecules.data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,mols[i].qx-mols[0].qx,mols[i].qy-mols[0].qy,mols[i].qz-mols[0].qz,\
	mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
	fclose(f);
}

void
MD::Ovout(int i){
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	Molecule *mols=vars->Molecules.data();
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,mols[i].qx-mols[0].qx,mols[i].qy-mols[0].qy,mols[i].qz-mols[0].qz,\
	mols[i].px-mols[0].px,mols[i].py-mols[0].py,mols[i].pz-mols[0].pz);
	fclose(f);
}


void
MD::exportDumpIn(void) {
	Molecule *mols = vars->Molecules.data();
	int Nmols=vars->Molecules.size();
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Natom=0;
	for (int i=0;i<Nmols;i++){
	//	if(vars->Region[i]!=AA) Natom++;
		if(vars->Region[i]!=CG) Natom+=mols[i].inAtoms.size();
	}
	double HL=vars->domainL*0.5;
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz transparency\n",
	count, Natom, -HL, HL, -HL, HL, -HL, HL);

	// Set center of calculation domain
	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->Molecules[0].qx;
		Y=vars->Molecules[0].qy;
		Z=vars->Molecules[0].qz;
	}

	int ID=0;
	for (int i=0;i<Nmols;i++){
	/*	if(vars->Region[i]!=AA) {
			fprintf(f, "%d %d %f %f %f %f %f %f %f\n", ID, mols[i].type, mols[i].qx-X, mols[i].qy-Y, mols[i].qz-Z, mols[i].px, mols[i].py, mols[i].pz, mols[i].w);
			ID++;
		}*/
		if(vars->Region[i]!=CG) {
			for(auto &b : mols[i].inAtoms){
				fprintf(f, "%d %s %f %f %f %f %f %f %f\n", ID, (vars->atypes[b.type].name).c_str(), b.qx-X, b.qy-Y, b.qz-Z, b.px, b.py, b.pz, (1-mols[i].w));
				ID++;
			}
		}
	}	
	fclose(f);
}

void
MD::exportDumpOut(void) {
	Molecule *mols = vars->Molecules.data();
	int Nmols=vars->Molecules.size();
	int count = vars->time;
	FILE*f=fopen("out.dump", "a");
	int Natom=0;
	for (int i=0;i<Nmols;i++){
		if(vars->Region[i]!=AA) Natom++;
	}

	double HL=vars->domainL*0.5;
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz transparency\n",
	count, Natom, -HL, HL, -HL, HL, -HL, HL);

	// Set center of calculation domain
	double X,Y,Z;
	//X=Y=Z=0;
	//if(flags->dump_fix==1) {
		X=vars->Molecules[0].qx;
		Y=vars->Molecules[0].qy;
		Z=vars->Molecules[0].qz;
	//}

	int ID=0;
	for (int i=0;i<Nmols;i++){
		if(vars->Region[i]!=AA) {
			fprintf(f, "%d %d %f %f %f %f %f %f %f\n", ID, mols[i].type, mols[i].qx-X, mols[i].qy-Y, mols[i].qz-Z, mols[i].px, mols[i].py, mols[i].pz, mols[i].w);
			ID++;
		}
	}	
	fclose(f);
}
