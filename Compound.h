#include <vector>
#include "Element.h"

#ifndef Compound_H

#define Compound_H

class Compound {
  
    private:
       std::vector<Element> compoundElements; 
       
    public:
        Compound();
        Compound(std::vector<Element> compoundElements);
        
        void addElement(Element ele);
        
        std::vector<Element> getElements();
        
        
};

#endif
