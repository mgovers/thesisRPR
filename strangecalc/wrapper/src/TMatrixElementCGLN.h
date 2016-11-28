/*!
 * \file TMatrixElementCGLN.h
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 
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

#ifndef TMATRIXELEMENTCGLN_H
#define TMATRIXELEMENTCGLN_H

#include "TMatrixElement.h"
#include "TCalculateCGLNCoeff.h"
#include "TCurrent.h"

class TMatrixElementCGLN: public TMatrixElement
{
public:
    static TMatrixElementCGLN* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< TMatrixElementCGLN Singleton factory
    static TMatrixElementCGLN* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label=-1); //!< TMatrixElementCGLN Singleton factory,  allowing off-shell kinematics
    virtual ~TMatrixElementCGLN(); //!< Destructor
    virtual void ReInitialise( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< Initialises the object as if it were newly constructed.
    virtual std::complex<double> calculateM ( int,int,int );//<! Calculates the complex value of the matrix element
    virtual const FourVector<GammaStructure>& GetCurrent(); //!< returns the hadronic current J^\mu
    virtual bool isGaugeInvariantForSpins (int Lp, int Ly); //!< Tests whether the amplitude is gauge invariant
    virtual void SetDataLabel(int label) {fDataLabel=label;} //!< Set the data label; used for caching the CGLN matrices.
    
protected:
    TMatrixElementCGLN(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< Protected constructor
    TMatrixElementCGLN(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label=-1); //!< Protected constructor, allowing off-shell kinematics
    TMatrixElementCGLN(const TMatrixElementCGLN& other); // Copy constructor
    TMatrixElementCGLN& operator=(const TMatrixElementCGLN& right);  // Assignment operator
    
    // Singleton implementation.
    static bool fgInstanceFlag;
    static TMatrixElementCGLN * fgMatrixElement;
    
    // Static consts
    static const std::vector< FourVector<std::complex< double > > > kUnitEpsilonFV;

    TCalculateCGLNCoeff* fCalculatorPtr; //!< the Coefficient Calculator, owned by this matrix element.
    TCurrent fCurrent; //!< The TCurrent object which holds the coefficients A_i in its state and calculates the rest jit
    FourVector<GammaStructure> fBareCurrent;
    int fDataLabel;
    bool fBareCurrentIsCalculated;
    
};

#endif  // #ifdef TMATRIXELEMENTCGLN_H
