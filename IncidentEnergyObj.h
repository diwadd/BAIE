#ifndef IncidentEnergyObj_H

#define IncidentEnergyObj_H

class IncidentEnergyObj {
  
    private:
        double incidentEnergy;
        
    public:
        IncidentEnergyObj();
        IncidentEnergyObj(double incidentEnergy);
        
        double getIncidentEnergy();
        void setIncidentEnergy(double incidentEnergy);
        
};

#endif