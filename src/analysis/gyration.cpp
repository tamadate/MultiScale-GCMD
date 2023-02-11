#include "../md.hpp"

void
MD::gyration_out(MD *md2){
	FILE*f=fopen(pp->gyration_path, "a");
	double numerator=0;
	std::array<double,3> *x=vars->position.data();
	for(auto i : vars->MolID[0]){
		double dx=x[i][0]-origin[0];
		double dy=x[i][1]-origin[1];
		double dz=x[i][2]-origin[2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=vars->mass[i]*rsq;
	}
	fprintf(f, "%f\t%f\t",vars->time,sqrt(numerator/pp->Mion));
	fclose(f);
}
