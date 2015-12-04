#include "IncidentEnergyObj.h"

IncidentEnergyObj::IncidentEnergyObj(){
}

IncidentEnergyObj::IncidentEnergyObj(double incidentEnergy){
    this->incidentEnergy = incidentEnergy;
}
        
double IncidentEnergyObj::getIncidentEnergy(){
    return this->incidentEnergy;
}

void IncidentEnergyObj::setIncidentEnergy(double incidentEnergy){
    this->incidentEnergy = incidentEnergy;
}