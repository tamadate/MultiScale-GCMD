#include "../../variables.hpp"
#include "../basePotential.hpp"

class PotentialBorn : public Potential {
	private:
	public:
		std::vector<std::vector<std::vector<double>>> pairCoeff;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialBorn(Variables *vars);
		~PotentialBorn(){};
};
