#include "ElementEnergyLine.h"
#include "XRLFunctions.h"

ElementEnergyLine::ElementEnergyLine(){   
}

ElementEnergyLine::ElementEnergyLine(int Z, int xrlint){
    this->Z = Z;
    this->xrlint = xrlint;
    this->energy = XRLFunctions::xrlLineEnergy(Z,xrlint);
}

ElementEnergyLine::ElementEnergyLine(int Z, int xrlint, double energy){
    this->Z = Z;
    this->xrlint = xrlint;
    this->energy = energy;
}

int ElementEnergyLine::getZ(){
    return this->Z;
}

int ElementEnergyLine::getXrlint(){
    return this->xrlint;
}

double ElementEnergyLine::getEnergy(){
    return this->energy;
}
        
void ElementEnergyLine::setZ(int Z){
    this->Z = Z;   
}

void ElementEnergyLine::setXrlint(int xrlint){
    this->xrlint = xrlint;
}

void ElementEnergyLine::setEnergy(double energy){
    this->energy = energy;
}