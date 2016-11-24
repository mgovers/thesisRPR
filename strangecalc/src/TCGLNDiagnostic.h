/*!
 * \file TCGLNDiagnostic.h
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

#ifndef TCGLNDIAGNOSTIC_H
#define TCGLNDIAGNOSTIC_H

#include <TMatrixElementCGLN.h>
#include <complex>

class TCGLNDiagnostic : public TMatrixElementCGLN
{

public:
    TCGLNDiagnostic(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr);
    TCGLNDiagnostic(const TCGLNDiagnostic& other);
    GammaStructure CalcM(const int index, const FourVector< std::complex< double > >& epsilon);
    GammaStructure CalcMwithkvect(const int index);
    std::complex< double > CalcSandwichedM(const int Lp, const int Ly, const int index, const FourVector< std::complex< double > >& epsilon);

};

#endif // TCGLNDIAGNOSTIC_H
