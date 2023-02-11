#include "../md.hpp"

//  Update position 
void
MD::setOrigin(void) {
	//omp_set_num_threads(omp_get_num_threads());
	std::array<double,3> *x=vars->position.data();
    double *m=vars->mass.data();
    origin[0]=0;
    origin[1]=0;
    origin[2]=0;
	for (auto i : vars->MolID[0]){
        double massRatio=m[i]/pp->Mion;
		origin[0] += x[i][0]*massRatio;
		origin[1] += x[i][1]*massRatio;
		origin[2] += x[i][2]*massRatio;
	}
}
