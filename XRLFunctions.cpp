#include "XRLFunctions.h"
#include "windows.h"
// on linux use dlfcn.h instead of windows.h

/*
 * we are using the xraylib library dll file
 * see T. Schoonjans, A. Brunetti, B. Golosio, M. Sanchez del Rio, 
 * V. A. Solé, C. Ferrero and L. Vincze, "The xraylib library for X-ray--matter interactions. 
 * Recent developments", Spectrochimica Acta B 66 (2011) 776-784
 * and https://github.com/tschoonj/xraylib
 * for more details
 */

// The xraylib copyright/license is below.

//Copyright (c) 2009, Bruno Golosio, Antonio Brunetti, Manuel Sanchez del Rio, Tom Schoonjans and Teemu Ikonen
//All rights reserved.

//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//    * The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.

//THIS SOFTWARE IS PROVIDED BY Bruno Golosio, Antonio Brunetti, Manuel Sanchez del Rio, Tom Schoonjans and Teemu Ikonen ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Bruno Golosio, Antonio Brunetti, Manuel Sanchez del Rio, Tom Schoonjans and Teemu Ikonen BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


double XRLFunctions::xrlCSTotal(int Z,double energy) {
    
    typedef float (*xraylibPointer)(int,float);
    HINSTANCE xraylibDLL = LoadLibrary("libxrl.dll");
    xraylibPointer xraylibFunction;
    
    xraylibFunction = (xraylibPointer)GetProcAddress(xraylibDLL,"CS_Total");
    
    return (double)xraylibFunction(Z, (float)energy);
}


double XRLFunctions::xrlCSFluorLineKissel(int Z, int xrlint, double energy) {
    
    typedef float (*xraylibPointer)(int,int,float);
    HINSTANCE xraylibDLL = LoadLibrary("libxrl.dll");
    xraylibPointer xraylibFunction;
    
    xraylibFunction = (xraylibPointer)GetProcAddress(xraylibDLL,"CS_FluorLine_Kissel");
    
    return (double)xraylibFunction(Z, xrlint,(float)energy);
    
}


double XRLFunctions::xrlLineEnergy(int Z, int xrlint) {
    
    typedef float (*xraylibPointer)(int,int);
    HINSTANCE xraylibDLL = LoadLibrary("libxrl.dll");
    xraylibPointer xraylibFunction;
    
    xraylibFunction = (xraylibPointer)GetProcAddress(xraylibDLL,"LineEnergy");
    
    return (double)xraylibFunction(Z, xrlint);
}


