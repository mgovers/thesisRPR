/*!
 * \file TCGLNDiagnostic.cpp
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 
 * \copyright
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

/*!
 * \class TCGLNDiagnostic
 *
 * Class to calculate and internals of TMatrixElementCGLN
 * to evaluate the decomposition in the CGLN basis.
 * 
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */
 
#include "TCGLNDiagnostic.h"
#include <complex>

using std::complex;

TCGLNDiagnostic::TCGLNDiagnostic(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr) 
 : TMatrixElementCGLN( w, k, costhkcm, pK, particlesPtr, observPtr)
{}
GammaStructure TCGLNDiagnostic::CalcM(const int index, const FourVector< std::complex< double > >& epsilon) 
{
    return TCalculateCGLNMatrix::CalcM(index,epsilon,fk4vect,fp4vect,fpY4vect);
}

std::complex< double > TCGLNDiagnostic::CalcSandwichedM(const int Lp, const int Ly, const int index, const FourVector< std::complex< double > >& epsilon) 
{
    determineSpinors(Lp,Ly);
    return fHyperonSpinors[Ly]*CalcM(index,epsilon)*fNucleonSpinors[Lp];
}
GammaStructure TCGLNDiagnostic::CalcMwithkvect(const int index) 
{
  FourVector <complex< double > > k4vect(fk4vect[0],fk4vect[1],fk4vect[2],fk4vect[3]);
  return CalcM(index,k4vect);
}
