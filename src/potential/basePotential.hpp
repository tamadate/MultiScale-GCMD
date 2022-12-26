#pragma once
#include "../variables.hpp"

//------------------------------------------------------------------------
class Potential {
	private:
	public:
		string potName="potential";
		virtual void printName(void) {cout<<potName<<endl;}
		virtual void compute(Variables *vars, FLAG *flags);
		Potential(){};
		~Potential(){};
};
