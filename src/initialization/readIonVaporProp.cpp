//------------------------------------------------------------------------
#include "../PhysicalProp.hpp"
//------------------------------------------------------------------------

void
Physical::readIonProp(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	Mion=0;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="atoms") {iflag=2; continue;}
		if (str=="bonds") {break;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) mass=stod(tmp);
				loop++;
			}
			atype.push_back(mass);
		}
		if (iflag==2) {
			int loop=0;
			int type, id;
			double charge, mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) type=stoi(tmp);
				if (loop==2) charge=stod(tmp);
				loop++;
			}
			mass=atype[type-1];
			Mion+=mass;
			z+=charge;
		}
	}
}

void
Physical::readVaporProp(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	Mvapor=0;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="atoms") {iflag=2; continue;}
		if (str=="bonds") {break;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) mass=stod(tmp);
				loop++;
			}
			atype_v.push_back(mass);
		}
		if (iflag==2) {
			int loop=0;
			int type, id;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) type=stoi(tmp)-1;
				loop++;
			}
			double mass=atype_v[type];
			Mvapor+=mass;
		}
	}
}
