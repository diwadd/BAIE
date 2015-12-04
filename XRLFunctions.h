#ifndef XRLFunctions_H

#define XRLFunctions_H

class XRLFunctions {
    
    public:
        static double xrlCSTotal(int Z,double energy);
        static double xrlCSFluorLineKissel(int Z, int xrlint, double energy);
        static double xrlLineEnergy(int Z, int xrlint);

};

#endif