#include "../md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/


void
MD::display(int output_ONOFF){
	obs->computeProps(vars,0);
	obs->computeProps(vars,1);
	obs->computeProps(vars,2);
    double virial=0;//ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
    double gaspress=(kb*Nof_around_gas*obs->Tout[1] + vars->totalVirial/3.0*6.95e-21)/(V*1e-30);
    double U = vars->Usum();
	double Kin=0;
	for(auto &kin : obs->Kin) Kin+=kin;
	double Kout=0;
	for(auto &kout : obs->Kout) Kout+=kout;

    std::cout << "-------------------TIME = " << vars->time/1000.0 << " ps----------------------" << endl;
	cout<<"*Inside propeties"<<endl;
	printf("  %-15s%-12s%-12s%-12s\n", "Prop", "Ion", "Gas", "Vapor");
	printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "K-Energy", obs->Kin[0], obs->Kin[1], obs->Kin[2]);
	printf("  %-15s%-12.2f%-12.2f%-12.2f\n", "Temperature", obs->Tin[0], obs->Tin[1], obs->Tin[2]);
	printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "Uintra", vars->Utotal.Uion, vars->Utotal.Ugas, vars->Utotal.Uvap);
	printf("  %-15s%-12.2e%-12.2e%-12.2e\n", "Uinter(Ion)", 0.0, vars->Utotal.Ugi, vars->Utotal.Uvi);
	printf("  %-15s%-12s%-12.2e%-12.2e\n", "Uinter(Gas)", "-", vars->Utotal.Ugg, vars->Utotal.Uvg);
	printf("  %-15s%-12s%-12s%-12.2e\n", "Uinter(Vapor)", "-", "-", vars->Utotal.Uvv);

	cout<<"*Out side propeties"<<endl;
	printf("  %-15s%-12s%-12s%-12s\n", "Prop", "Ion", "Gas", "Vapor");
	printf("  %-15s%-12s%-12.2e%-12.2e\n", "K-Energy", "-", obs->Kout[1], obs->Kout[2]);
	printf("  %-15s%-12s%-12.2f%-12.2f\n", "Temperature", "-", obs->Tout[1], obs->Tout[2]);

	cout<<"*System total propeties"<<endl;
	printf("  K = %1.2e	U = %1.2e	Press = %f\n",Kin+Kout, U, gaspress/101300.0);

	cout<<"*Times"<<endl;
	printf("  tion  = %1.1f s	tgas = %1.1f s 	tvap = %1.1f s\n",vars->times.tion,vars->times.tgas,vars->times.tvap);
	printf("  tvi   = %1.1f s	tgi  = %1.1f s	tvg  = %1.1f s	tvv  = %1.1f s\n",vars->times.tvi,vars->times.tgi,vars->times.tvg,vars->times.tvv);
	printf("  tpair = %1.1f s	tpos = %1.1f s	tvel = %1.1f s	tetc = %1.1f s\n",vars->times.tpair,vars->times.tpos,vars->times.tvel,vars->times.tetc);
	printf("  tpot  = %1.1f s	ttot = %1.1f s\n",(vars->times.tvi+vars->times.tgi+vars->times.tvv+vars->times.tvg+vars->times.tion+vars->times.tgas+vars->times.tvap), omp_get_wtime()-startTime);
	printf("  NCPU = %d\n",Nth);

	cout <<endl;

	sprintf(filepath, "K_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	fprintf(f,"%e %e %e %e %e %e\n",vars->time,obs->Kin[0],obs->Kin[1],obs->Kin[2],obs->Kout[1],obs->Kout[2]);
	fclose(f);

	sprintf(filepath, "U_%d.dat", int(calculation_number));
	f=fopen(filepath, "a");
	fprintf(f,"%e %e %e %e %e %e %e %e %e\n",vars->time,vars->Utotal.Uion,vars->Utotal.Ugas,vars->Utotal.Uvap,vars->Utotal.Ugi,vars->Utotal.Ugg,vars->Utotal.Uvg,vars->Utotal.Uvi,vars->Utotal.Uvv);
	fclose(f);
}
