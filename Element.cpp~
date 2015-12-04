#include "Element.h"
#include "ElementEnergyLine.h"
#include <vector>

Element::Element() {
}

Element::Element(int Z, double W){
    this->Z = Z;
    this->W = W;
    
    std::vector<ElementEnergyLine> energyLinesVector(4);
    
    for(int i = 0; i < 4; i++){
        energyLinesVector[i] = ElementEnergyLine(Z,i);
    }
    this->energyLines = energyLinesVector;
}

void Element::setZ(int Z){
    this->Z = Z;
}

int Element::getZ() {
    return this->Z;
}

double Element::getW(){
    return this->W;
}

std::vector<ElementEnergyLine> Element::getEnergyLines(){
    return  this->energyLines;   
}
