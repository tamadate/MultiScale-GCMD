//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


bool
MD::errorCheck(void){
    bool returnValue=true;
    double AAsize=(Rinter+BoundL+BoundMergin)*2;
    if(AAsize>vars->domainL) {
        returnValue=false;
        cout<<"***Error: AA region is larger than domain length!\n"<<endl;
    }
    return returnValue;
}