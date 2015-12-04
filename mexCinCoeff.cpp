#include "mex.h"
#include "windows.h"
#include "Element.h"
#include "XRLFunctions.h"
#include "Compound.h"
#include "IncidentEnergyObj.h"
#include "CinFunctioncs.h"
#include <cmath>

typedef float (*xraylibPointer)(int,float);


std::vector<Element> makeElementVector(double *zArray, double* wArray, mwSize m);

/* The computational routine */
void thisCin(double *inZiMatrix, double *inZnMatrix, double inEnergy, int xrlintLambda, double *zArray, double *wArray, double *theta1, double theta2, double *out, mwSize n, mwSize m, mwSize nt, mwSize mt)
{
    
    IncidentEnergyObj *incidentEnergy = new IncidentEnergyObj(inEnergy);
    ElementEnergyLine* lambda = new ElementEnergyLine((int)inZiMatrix[0],xrlintLambda);
    
    std::vector<Element> elementsForCompound;
    elementsForCompound = makeElementVector(zArray, wArray, m);
  
    Compound *compound = new Compound( elementsForCompound );

    Element *iElement = new Element((int)inZiMatrix[0], inZiMatrix[1]);
    Element *nElement = new Element((int)inZnMatrix[0], inZnMatrix[1]);
    
    mwSize ii;
    mwSize jj;
    
    for(ii = 0; ii < nt; ii++){
        for(jj = 0; jj < mt; jj++){
            out[ii+jj*nt] = Cin(iElement, nElement ,  inEnergy, lambda, compound, theta1[ii+jj*nt], theta2);
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *outMatrix;
    
    double incidentEnergy;
    int lambdaXrlint;
    
    double *inZiMatrix;
    double *inZnMatrix;
    
    double *inZMatrix;
    double *inWMatrix;
    
    double *theta1Matrix;
    double theta2;
    
    size_t theta1MatrixRows;
    size_t theta1MatrixCol;
    
    size_t inZMatrixRows;
    size_t inZMatrixCol;
    
    size_t inWMatrixRows;
    size_t inWMatrixCol;
    
    inZiMatrix = mxGetPr(prhs[0]);
    inZnMatrix = mxGetPr(prhs[1]);
    
    incidentEnergy = mxGetScalar(prhs[2]);
    lambdaXrlint = mxGetScalar(prhs[3]);
    
    inZMatrix = mxGetPr(prhs[4]);
    inZMatrixCol = mxGetN(prhs[4]);
    inZMatrixRows = mxGetM(prhs[4]);
    
    inWMatrix = mxGetPr(prhs[5]);
    inWMatrixCol = mxGetN(prhs[5]);
    inWMatrixRows = mxGetM(prhs[5]);
    
    if (inZMatrixRows != 1 || inWMatrixRows != 1 || (inWMatrixCol != inZMatrixCol)) {
        mexErrMsgTxt("Error code (info): \n Somethings wrond with the input. \n Pleas check the compound array.");
    }
   
    theta1Matrix = mxGetPr(prhs[6]);
    theta1MatrixRows = mxGetN(prhs[6]);
    theta1MatrixCol = mxGetM(prhs[6]);
    
    theta2 = mxGetScalar(prhs[7]);

    
    plhs[0] = mxCreateDoubleMatrix((mwSize)theta1MatrixRows,(mwSize)theta1MatrixCol,mxREAL);
    outMatrix = mxGetPr(plhs[0]);


    thisCin(inZiMatrix, inZnMatrix, incidentEnergy, lambdaXrlint, inZMatrix, inWMatrix, theta1Matrix, theta2, outMatrix, (mwSize)inZMatrixRows, (mwSize)inZMatrixCol,(mwSize)theta1MatrixRows, (mwSize)theta1MatrixCol);

}

