
#include "../../variables.hpp"
#include "../basePotential.hpp"

class PotentialAMBER : public Potential {
	private:
		std::vector<Pair> longPair;
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void computeLong(Variables *vars, FLAG *flags);
		void computeBond(Variables *vars, FLAG *flags);
		void computeAngle(Variables *vars, FLAG *flags);
		void computeDihedral(Variables *vars, FLAG *flags);
		void initialAMBER(Variables *vars, FLAG *flags);
		PotentialAMBER(Variables *vars, FLAG *flags){
			initialAMBER(vars,flags);
			potName="AMBER ion";
		};
		~PotentialAMBER(){};
};
