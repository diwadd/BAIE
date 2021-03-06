#ifndef CinFunctions_H

#define CinFunctions_H

#include "mex.h"
#include "windows.h"
#include "Element.h"
#include "XRLFunctions.h"
#include "Compound.h"
#include "IncidentEnergyObj.h"

double muTotal(Compound *compound,double energy);
double muStar(Compound *compound, double energy, double lambdaEnergy, double theta1, double theta2); // in the current paper version this is  known as mu bar

bool atomKroneckerDelta(Element *element1, Element *element2);

double xL(Compound *compound, double energy, double lambdaEnergy, double nuEnergy, double theta1, double theta2 );
double xLCurly(Compound *compound, double energy, double lambdaEnergy, double nuEnergy, double theta1, double theta2, double mT);

double xdL(Element *nElement, double energy, double nuEnergy, Compound *compound, double theta1);
double xSjnu( Element *iElement, Element *jElement,Compound* compound, double energy, ElementEnergyLine *lambda, ElementEnergyLine *nu, double theta1, double theta2);

double dividingFactor(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2);
double CinBA(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2);

double CinBAThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT);
double CinIE(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2);

double Cin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2);
double yPart(Element *iElement,double energy, ElementEnergyLine *lambda,Compound *compound, double theta1, double theta2);

double xSjnuThin( Element *iElement, Element *jElement,Compound* compound, double energy, ElementEnergyLine *lambda, ElementEnergyLine *nu, double theta1, double theta2, double mT);
double dividingFactorThin(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2,double mT);            

double xdLThin(Element *nElement, double energy,double lambdaEnergy ,double nuEnergy, Compound *compound, double theta1, double theta2, double mT); 
double CinIEThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT);

double CinThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT);
double yPartThin(Element *iElement,double energy, ElementEnergyLine *lambda,Compound *compound, double theta1, double theta2, double mT);

double dividingFactorOnePlus(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2);
double dividingFactorThinOnePlus(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2, double mT);

double CinIEThinHOnePlus(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT);
double CinBAThinHOnePlus(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT);

double E1Fun(double x);

#endif
