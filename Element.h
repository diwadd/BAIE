#include <vector>
#include "ElementEnergyLine.h"

#ifndef Element_H

#define Element_H

class Element {

    private:
        int Z;
        double W;
        std::vector<ElementEnergyLine> energyLines;
        
    public:
        Element();
        Element(int Z, double W);
        
        int  getZ();
        double getW();
        std::vector<ElementEnergyLine> getEnergyLines();
        
        void setZ(int Z);
};

#endif