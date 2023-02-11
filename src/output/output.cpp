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
	double *mass=vars->mass.data();
	std::array<double,3> *v=vars->velocity.data();
	double Mass=0;
	double vx=0;
	double vy=0;
	double vz=0;
	for(auto i : vars->MolID[0]) Mass+=mass[i];
	for(auto i : vars->MolID[0]) {
		double massRatio=mass[i]/Mass;
		vx+=v[i][0]*massRatio;
		vy+=v[i][1]*massRatio;
		vz+=v[i][2]*massRatio;
	}
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, origin[0], origin[1], origin[2], vx, vy, vz);
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
MD::exportDump(void) {
	int Natom=vars->position.size();
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	double HL=vars->domainL*0.5;
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
	count, Natom, -HL, HL, -HL, HL, -HL, HL);

	// Set center of calculation domain
	double X,Y,Z;
	X=Y=Z=0;
	//if(flags->dump_fix==1) {
		X=origin[0];
		Y=origin[1];
		Z=origin[2];
	//}

	std::array<double,3> *v=vars->velocity.data();
	std::array<double,3> *x=vars->position.data();
	for (int i=0;i<Natom;i++){
		int aid=vars->type[i];
		fprintf(f, "%d %s %f %f %f %f %f %f\n", i, (vars->atypes[aid].name).c_str(), x[i][0]-X, x[i][1]-Y, x[i][2]-Z, v[i][0], v[i][1], v[i][2]);
	}	
	fclose(f);
}