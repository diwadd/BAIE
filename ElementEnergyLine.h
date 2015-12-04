#ifndef ElementEnergyLine_H

#define ElementEnergyLine_H

class ElementEnergyLine{
    
    private:
        int Z;
        int xrlint;
        double energy;
        
    public:
        ElementEnergyLine();
        ElementEnergyLine(int Z, int xrlint);
        ElementEnergyLine(int Z, int xrlint, double energy);
        
        int getZ();
        int getXrlint();
        double getEnergy();
        
        void setZ(int Z);
        void setXrlint(int xrlint);
        void setEnergy(double energy);
};

#endif