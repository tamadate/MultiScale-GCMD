//------------------------------------------------------------------------
#include "../PhysicalProp.hpp"
//------------------------------------------------------------------------

void
Physical::setPhysicalProp(int gastype, double T, double p){
	if(gastype==1) {
		Mgas=MHe;
		myu=myuHe;
		alphagas=alphaHe;
	}
	if(gastype==2) {
		Mgas=MN2;
		myu=myuN2*pow(T/TrefN2,1.5)*(TrefN2+SN2)/(T+SN2);
		alphagas=alphaN2/2.0;
	}
	if(gastype==3) {
		Mgas=MN2;
		myu=myuN2*pow(T/TrefN2,1.5)*(TrefN2+SN2)/(T+SN2);
		alphagas=alphaN2;
	}
	if(gastype==4) {
		Mgas=MAr;
		myu=myuAr*pow(T/TrefAr,1.5)*(TrefAr+SAr)/(T+SAr);
		alphagas=alphaAr;
	}
	mvapor=Mvapor/Nw/1000.0;
	mgas=Mgas/Nw/1000.0;
	m=Mion/Nw/1000.0;
	m_gas=mgas*m/(mgas+m);
	m_N2=mN2*m/(mN2+m);
	c=sqrt(8*kb*T/M_PI/m);
	cgas=sqrt(8*kb*T/M_PI/mgas);
	cvapor=sqrt(8*kb*T/M_PI/mvapor);

	printf("Ion mass\t\t%f g/mol\nIon charges\t\t%f\n", Mion,z);
	printf("Vapor mass\t\t%f g/mol\n", Mvapor);
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
}
