#include <vector>
#include "Compound.h"

Compound::Compound() {
}

Compound::Compound(std::vector<Element> compoundElements){
    this->compoundElements = compoundElements;
}
        
void Compound::addElement(Element ele){
    (this->compoundElements).push_back(ele);
}
        
std::vector<Element> Compound::getElements(){
    return this->compoundElements;
}
