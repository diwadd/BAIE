#include "mex.h"
#include "matrix.h"
#include "windows.h"
#include "Element.h"
#include "XRLFunctions.h"
#include "Compound.h"
#include "IncidentEnergyObj.h"
#include "CinFunctioncs.h"
#include <vector>
#include <cmath>

#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/complex/acos.hpp>


#define ZEROForComp 0.1
#define zeroTolerance 10000000

bool atomKroneckerDelta(Element *element1, Element *element2){
    
    if( element1->getZ() == element2->getZ() )
        return true;
    else
        return false;
}

double E1Fun(double x){
    
    if( (ceilf(x * zeroTolerance) / zeroTolerance) == 0  ) 
        return 0.0;    
    
    // here we use the relation E1(x) = -Ei(-x)
    return -1.0*boost::math::expint(-1.0*x);
    
}

/*
 * the summed quantity in the Sharmen equation
 */

double xSjnu( Element *iElement, Element *jElement,Compound* compound, double energy, ElementEnergyLine *lambda, ElementEnergyLine *nu, double theta1, double theta2){
   
    double xS = (jElement->getW()) * XRLFunctions::xrlCSFluorLineKissel( iElement->getZ() , lambda->getXrlint(), nu->getEnergy()) * XRLFunctions::xrlCSFluorLineKissel( jElement->getZ() , nu->getXrlint(), energy) *xL(compound,energy, lambda->getEnergy(), nu->getEnergy(), theta1, theta2 );
    
    // because xL can be undetermined
    if (mxIsNaN(xS)==1){
        return 0.0;
    } else {
        return xS;
    }
}


/*
 * the dividing factor is defined here as tau_{i,\lambda}(E) + \sum_{j,\nu}S_{i,j,\lambda,\nu}
 * this is a somewhat diffrent definition as in the original formulas yet it is equaivalent
 */

double dividingFactor(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
   
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
            temp = temp + xSjnu(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2);
        }
    }
    
    // convenction with old dividing factor
    //double df = 1.0 + (temp/XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) );
    double df = XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) + temp;
    
    if (mxIsNaN(df)==1){
        return XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy);
    } else {
        return df;
    }
}


double dividingFactorOnePlus(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
   
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
            temp = temp + xSjnu(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2);
        }
    }
    
    // convenction with old dividing factor
    double df = 1.0 + (temp/XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) );
    //double df = XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) + temp;
    
    if (mxIsNaN(df)==1){
        return XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy);
    } else {
        return df;
    }
}

/*
 * Dividing factor for thin samples
 */

double dividingFactorThin(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2, double mT){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
   
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
            // xSjnuThin stands for the van Dyck, van Espen integral result
            temp = temp + xSjnuThin(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2, mT);
        }
    }

    // convenction with old dividing factor
    //double df = 1.0 + (temp/XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) );
    double muS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
    double df = XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy)*(1 - exp(-1.0*muS*mT)) + temp;

    if (mxIsNaN(df)==1){
        return XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy);
    } else {
        return df;
    }
}

/*
 * Dividing factor for thin samples - thin variant
 */

double dividingFactorThinOnePlus(Element *iElement, ElementEnergyLine *lambda, double energy, Compound* compound, double theta1, double theta2, double mT){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
   
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
            // xSjnuThin stands for the van Dyck, van Espen integral result
            temp = temp + xSjnuThin(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2, mT);
        }
    }

    // convenction with old dividing factor
    //double df = 1.0 + (temp/XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) );
    double muS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
    double df = 1.0 + temp/(XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy)*(1 - exp(-1.0*muS*mT)));

    if (mxIsNaN(df)==1){
        return XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy);
    } else {
        return df;
    }
}

/*
 * the mu star term
 */

double muStar(Compound *compound, double energy, double lambdaEnergy, double theta1, double theta2) {
    
    return ((muTotal(compound, energy)/cos(theta1)) + (muTotal(compound,lambdaEnergy)/cos(theta2)));
    
}

/*
 * the derivative of the L term
 */

double xdL(Element *nElement, double energy, double nuEnergy, Compound *compound, double theta1){
    
   // depending on the value of muTotal can be undetermined (for example when Ejnu=0) 
   double A =  XRLFunctions::xrlCSTotal(nElement->getZ(), energy)*(nElement->getW())*cos(theta1)/(2.0*muTotal(compound,energy)); 
   double B1 = 1.0 / (muTotal(compound,nuEnergy)*cos(theta1));
   double B2 = muTotal(compound,energy)*B1;
   double C1 = 1.0/muTotal(compound,energy);
   return A*(B1*(1.0/(1.0 + B2)    )  - C1*log(1.0 + B2)              );
    
}

/*
 * the L term of the Sherman equation
 */

double xL(Compound *compound, double energy, double lambdaEnergy, double nuEnergy, double theta1, double theta2 ){
    
    // the value of xL can be undetermined for example for Eilambda=0
    // this usually occurs for ligth elements
    return 0.5*( (cos(theta1)/ muTotal(compound, energy) )*log( 1.0 +  (muTotal(compound, energy)/( muTotal(compound, nuEnergy)*cos(theta1) )   ) )     +  (cos(theta2)/muTotal(compound,lambdaEnergy) ) * log( 1.0 +  (muTotal(compound, lambdaEnergy)/( muTotal(compound, nuEnergy)*cos(theta2) )   ) )  );
    
}

/*
 * total attenuation coefficient mu
 *
 */

double muTotal(Compound *compound, double energy){
    
    int numberOfElements = ((*compound).getElements()).size();

    double muTotal = 0.0;
    for(int i = 0; i < numberOfElements; i++) {
        muTotal = muTotal + ((*compound).getElements())[i].getW()*XRLFunctions::xrlCSTotal(((*compound).getElements())[i].getZ(), energy  );
    }
    
    return muTotal;
}



/*
 * makeElementVector
 */

std::vector<Element> makeElementVector(double *zArray, double* wArray, mwSize m) {
    
    std::vector<Element> elementVector;
    
    for(mwSize i = 0; i < m; i++){

        Element *tempElement = new Element((int)(zArray[i]), wArray[i]);
        elementVector.push_back( *tempElement);
    }
    return elementVector;
}


/*
 * beam attenuation part of Cin - thick variant
 */

double CinBA(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2){
    
    // convenction with old dividing factor
    //double oneOverDividingFactor = 1.0/dividingFactor(iElement, lambda, energy, compound, theta1, theta2);
    double tauOverDividingFactor = XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) / dividingFactor(iElement, lambda, energy, compound, theta1, theta2);
    double A = atomKroneckerDelta(iElement, nElement) - XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  muStar(compound, energy, lambda->getEnergy(), theta1, theta2)  ); 
    return tauOverDividingFactor*A;
}


/*
 * beam attenuation part of Cin - thin variant
 */

double CinBAThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT){
    
    double tauRound = ceilf(XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) * zeroTolerance) / zeroTolerance;
    if(  tauRound == 0.0 ){
        return 0.0;
    }
    
    // convenction with old dividing factor
    //double oneOverDividingFactor = 1.0/dividingFactor(iElement, lambda, energy, compound, theta1, theta2);
    double mS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
    double tauOverDividingFactor = (XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy)*(1 - exp(-1.0*mS*mT) ) ) / dividingFactorThin(iElement, lambda, energy, compound, theta1, theta2, mT);
    double oneMinusF = 1.0 - mS*mT*( exp(-mS*mT) / ( 1 - exp(-mS*mT) ) );
    double A = atomKroneckerDelta(iElement, nElement) - ( XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  mS  ) ) * oneMinusF; 
    return tauOverDividingFactor*A;
    
}


/*
 * CinBa--Thick-no h division
 */

double CinBAThinHOnePlus(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT){
    
    double tauRound = ceilf(XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) * zeroTolerance) / zeroTolerance;
    if(  tauRound == 0.0 ){
        return 0.0;
    }
    
    // convenction with old dividing factor
    //double oneOverDividingFactor = 1.0/dividingFactor(iElement, lambda, energy, compound, theta1, theta2);
    double mS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
    double tauOverDividingFactor = 1.0 / dividingFactorThinOnePlus(iElement, lambda, energy, compound, theta1, theta2, mT);
    double oneMinusF = 1.0 - mS*mT*( exp(-mS*mT) / ( 1 - exp(-mS*mT) ) );
    double A = atomKroneckerDelta(iElement, nElement) - ( XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  mS  ) ) * oneMinusF; 
    return tauOverDividingFactor*A;
    
}

/*
 * indirect excitation part of Cin
 */

double CinIE(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
    
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
            
            
            double xLv = xL(compound, energy, lambda->getEnergy(), (elementSet[i].getEnergyLines())[j].getEnergy(), theta1, theta2 );
            double xdLv = xdL(nElement,  energy,  (elementSet[i].getEnergyLines())[j].getEnergy(), compound, theta1);
            
            double part1 = atomKroneckerDelta( &elementSet[i], nElement) -  XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  muStar(compound, energy, lambda->getEnergy(), theta1, theta2)  );
            double part2 = xdLv/xLv;
            double S = xSjnu(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2) ;
            
            if(mxIsNaN(part1)==1) {continue;}
            if(mxIsNaN(part2)==1) {continue;}
            
            temp = temp +  S*(part1 + part2);
        }
    }
    
    // convenction with old dividing factor
    //double cieReturn = temp/(  XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) *dividingFactor(iElement, lambda, energy, compound, theta1, theta2)  );     
    double cieReturn = temp/( dividingFactor(iElement, lambda, energy, compound, theta1, theta2) );     
 
 
    if (mxIsNaN(cieReturn)==1){
        return 0.0;
    } else {
        return cieReturn;
    }
    
}



/*
 * Cin
 */

double Cin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2){

    /*
     * when \tau_{i,\lambda}(E) is zero then the formulas do not make any sense because there is no fluorescence
     * therefore there is no point in counting the c coefficients
     */
    
    double tauRound = ceilf(XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) * zeroTolerance) / zeroTolerance;
    if(  tauRound == 0.0 ){
        return 0.0;
    }
    
    double baPart = CinBA(iElement, nElement, energy, lambda, compound, theta1, theta2);
    double iePart = CinIE(iElement, nElement, energy, lambda, compound, theta1, theta2);
    return baPart + iePart;

}


/*
 * Cin Thin
 */

double CinThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT){

    /*
     * when \tau_{i,\lambda}(E) is zero then the formulas do not make any sense because there is no fluorescence
     * therefore there is no point in counting the c coefficients
     */
    
    double tauRound = ceilf(XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) * zeroTolerance) / zeroTolerance;
    if(  tauRound == 0.0 ){
        return 0.0;
    }
    
    double baPart = CinBAThin(iElement, nElement, energy, lambda, compound, theta1, theta2, mT);
    double iePart = CinIEThin(iElement, nElement, energy, lambda, compound, theta1, theta2, mT);
    return baPart + iePart;

}


/*
 * yPart
 */

double yPart(Element *iElement,double energy, ElementEnergyLine *lambda,Compound *compound, double theta1, double theta2){
    
        double dF = dividingFactor(iElement, lambda, energy, compound, theta1, theta2);
        // convenction with old dividing factor
        //double yP = XRLFunctions::xrlCSFluorLineKissel( iElement->getZ() , lambda->getXrlint(), energy) / muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
        double yP = 1.0 / muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
        
        
        return yP*dF;

}


/*
 * yPart van Dyck
 */

double yPartThin(Element *iElement,double energy, ElementEnergyLine *lambda,Compound *compound, double theta1, double theta2, double mT){
    
        double dF = dividingFactorThin(iElement, lambda, energy, compound, theta1, theta2, mT);
        // convenction with old dividing factor
        //double yP = XRLFunctions::xrlCSFluorLineKissel( iElement->getZ() , lambda->getXrlint(), energy) / muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
        double mS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
        double yP = ( 1.0 - exp(-mS*mT) )/ mS ;
        
        
        return yP*dF;

}



/*
 * L as given by van Dyck
 */

double xLCurly(Compound *compound, double energy, double lambdaEnergy, double nuEnergy, double theta1, double theta2, double mT){

    double mEil = muTotal(compound, lambdaEnergy);
    double mEjnu = muTotal(compound, nuEnergy);
    double mE = muTotal(compound, energy);    
    
    // when the sums in the IE term run over energy lines of the i atom 
    // we get infiniteies beacuse these sum elements should not exist
    // they are cancelled by the tau's 
    
    // F1, F2 ... F5 - in the future these have to be replaced with separate functions
    // so we don't repeat code
    double F1 = log( abs(1 - mE /(cos(theta1) * mEjnu) ) ) + E1Fun(  ( mEjnu - mE/cos(theta1) )*mT  );
    double F2 = log( abs(1 - mEil /(cos(theta2) * mEjnu) ) ) + E1Fun(  ( mEjnu - mEil/cos(theta2) )*mT  );
    double F3 = log(1 + mE /(cos(theta1) * mEjnu ) ) + E1Fun( ( mEjnu + mE/cos(theta1) )*mT );
    double F4 = log(1 + mEil /(cos(theta2) * mEjnu ) ) + E1Fun( ( mEjnu + mEil/cos(theta2) )*mT );
    double F5 = ( exp(-mEil*mT/cos(theta2)) + exp(-mE*mT/cos(theta1)) )*E1Fun(mEjnu*mT);
    
    double muS = muStar(compound, energy, lambdaEnergy, theta1, theta2);
    
    double sum1 = ( cos(theta1)/(2.0*mE) )*(exp(-muS*mT)*F1+F3-F5);
    double sum2 = ( cos(theta2)/(2.0*mEil) )*(exp(-muS*mT)*F2+F4-F5);
    
    return (sum1 + sum2);
    
}


/*
 *van Dyck integration result
 */

double xSjnuThin( Element *iElement, Element *jElement,Compound* compound, double energy, ElementEnergyLine *lambda, ElementEnergyLine *nu, double theta1, double theta2, double mT){
    
    
    double tauila = XRLFunctions::xrlCSFluorLineKissel( iElement->getZ() , lambda->getXrlint(), nu->getEnergy());
    double taujnu = XRLFunctions::xrlCSFluorLineKissel( jElement->getZ() , nu->getXrlint(), energy);
    
    if ( ( (ceilf(tauila * zeroTolerance) / zeroTolerance) == 0 )  || ( (ceilf(taujnu * zeroTolerance) / zeroTolerance) == 0 ) ) 
        return 0.0;
    
    double LCurly = xLCurly(compound, energy, lambda->getEnergy(), nu->getEnergy(), theta1, theta2, mT);
    
    double res = (jElement->getW())*tauila*taujnu*LCurly; 
    return res;
    
}


/*
derivative of van Dyck equation
this should reduce to xdL for T -> infinity
the absolute value from van Dyck is not present
 */

double xdLThin(Element *nElement, double energy,double lambdaEnergy ,double nuEnergy, Compound *compound, double theta1, double theta2, double mT){
    
    
    double mEil = muTotal(compound, lambdaEnergy);
    double mEjnu = muTotal(compound, nuEnergy );
    double mE = muTotal(compound, energy );
    double mEn = XRLFunctions::xrlCSTotal(nElement->getZ(), energy);
    double wn = nElement->getW();
    
    // F1, F2 ... F5 - in the future these have to be replaced with separate functions
    // so we don't repeat code
    double F1 = log( abs(1 - mE /(cos(theta1) * mEjnu) ) ) + E1Fun(  ( mEjnu - mE/cos(theta1) )*mT  );
    double F2 = log( abs(1 - mEil /(cos(theta2) * mEjnu) ) ) + E1Fun(  ( mEjnu - mEil/cos(theta2) )*mT  );
    double F3 = log(1 + mE /(cos(theta1) * mEjnu ) ) + E1Fun( ( mEjnu + mE/cos(theta1) )*mT );
    double F4 = log(1 + mEil /(cos(theta2) * mEjnu ) ) + E1Fun( ( mEjnu + mEil/cos(theta2) )*mT );
    double F5 = ( exp(-mEil*mT/cos(theta2)) + exp(-mE*mT/cos(theta1)) )*E1Fun(mEjnu*mT);
    
    double muS = muStar(compound, energy, lambdaEnergy, theta1, theta2);
    
    double F1F3F5 = exp(-muS*mT)*F1+F3-F5;
    double F2F4F5 = exp(-muS*mT)*F2+F4-F5;
    
    double xdLThinPart1 = -F1F3F5*(  (mEn*wn*cos(theta1)) / (2*mE*mE) );
    
    double onePlusExp = 1.0 - exp( - (mEjnu + mE/cos(theta1) )*mT );
    double onePlusmExp = 1.0 - exp( - (mEjnu - mE/cos(theta1) )*mT );
    double oneMinusFrac = 1.0 - mE/(mEjnu*cos(theta1));
    double onePlusFrac = 1.0 + mE/(mEjnu*cos(theta1));
    
    double bigBraketPart1 = -F1*mT*exp(-muS*mT);
    double bigBraketPart2 = -1.0*( (exp(-muS*mT)/mEjnu)/oneMinusFrac )*onePlusmExp;
    double bigBraketPart3 = ( (1.0/mEjnu)/onePlusFrac )*onePlusExp;
    double bigBraketPart4 = mT*E1Fun(mEjnu*mT)*exp(-(mE*mT)/cos(theta1));
    double frontOfBigBraket = ((wn*mEn)/(2.0*mE)) ;
    double bigBraketTerm = frontOfBigBraket*(bigBraketPart1 + bigBraketPart2 + bigBraketPart3 + bigBraketPart4);
    
    double smallBraket1 = F2*exp(-muS*mT);
    double smallBraket2 = -1.0*E1Fun(mEjnu*mT)*exp(-(mE*mT)/cos(theta1));
    double frontSmallBraket = -1.0*mT*( (mEn*wn*cos(theta2))/(2*mEil*cos(theta1)) );
    double smallBraketTerm = frontSmallBraket*(smallBraket1 + smallBraket2);
    
    double xdLThinValue = xdLThinPart1 + bigBraketTerm + smallBraketTerm;
    
    return xdLThinValue;
    
}

/*
 * C^IE for thin samples
 */

double CinIEThin(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
    
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
    
            
            double xLv = xLCurly(compound, energy, lambda->getEnergy(), (elementSet[i].getEnergyLines())[j].getEnergy(), theta1, theta2, mT );
            double xdLv = xdLThin(nElement, energy, lambda->getEnergy() , (elementSet[i].getEnergyLines())[j].getEnergy(), compound, theta1, theta2, mT);
                    
            double part1 = atomKroneckerDelta( &elementSet[i], nElement) -  XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  muStar(compound, energy, lambda->getEnergy(), theta1, theta2)  );
            double part2 = xdLv/xLv;
            double S = xSjnuThin(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2, mT) ;
            
            if(mxIsNaN(part1)==1) {continue;}
            if(mxIsNaN(part2)==1) {continue;}
            
            temp = temp +  S*(part1 + part2);
        }
    }
    
    // convenction with old dividing factor
    //double cieReturn = temp/(  XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) *dividingFactor(iElement, lambda, energy, compound, theta1, theta2)  );     
    double cieReturn = temp/( dividingFactorThin(iElement, lambda, energy, compound, theta1, theta2, mT) );     
 
    if (mxIsNaN(cieReturn)==1){
        return 0.0;
    } else {
        return cieReturn;
    }
    
}


/*
 * C^IE for thin samples H one plus
 */

double CinIEThinHOnePlus(Element *iElement, Element *nElement, double energy, ElementEnergyLine *lambda, Compound *compound, double theta1, double theta2, double mT){
    
    int numberOfElements = (compound->getElements()).size();
    std::vector<Element> elementSet = (compound->getElements());
    double temp = 0.0;
    
    for(int i = 0; i < numberOfElements; i++){
        int numberOfEnergyLines = (elementSet[i].getEnergyLines()).size();
        for(int j = 0; j < numberOfEnergyLines; j++){
    
            
            double xLv = xLCurly(compound, energy, lambda->getEnergy(), (elementSet[i].getEnergyLines())[j].getEnergy(), theta1, theta2, mT );
            double xdLv = xdLThin(nElement, energy, lambda->getEnergy() , (elementSet[i].getEnergyLines())[j].getEnergy(), compound, theta1, theta2, mT);
            
            double part1 = atomKroneckerDelta( &elementSet[i], nElement) -  XRLFunctions::xrlCSTotal(nElement->getZ(),energy)*(nElement->getW())/(cos(theta1)*  muStar(compound, energy, lambda->getEnergy(), theta1, theta2)  );
            double part2 = xdLv/xLv;
            double S = xSjnuThin(iElement, &elementSet[i], compound, energy, lambda, &(elementSet[i].getEnergyLines())[j], theta1, theta2, mT) ;
            
            if(mxIsNaN(part1)==1) {continue;}
            if(mxIsNaN(part2)==1) {continue;}
            
            temp = temp +  S*(part1 + part2);
        }
    }
    
    // convenction with old dividing factor
    //double cieReturn = temp/(  XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy) *dividingFactor(iElement, lambda, energy, compound, theta1, theta2)  );     
    double muS = muStar(compound, energy, lambda->getEnergy(), theta1, theta2);
    double tau = XRLFunctions::xrlCSFluorLineKissel(iElement->getZ(), lambda->getXrlint(), energy);
    double df = dividingFactorThinOnePlus(iElement, lambda, energy, compound, theta1, theta2, mT);
    double cieReturn = temp/( tau*(1.0 - exp(-1.0*muS*mT) )*df );
    
    if (mxIsNaN(cieReturn)==1){
        return 0.0;
    } else {
        return cieReturn;
    }
    
}
