/*!
 * \file TMatrixElementOld.h
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 
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

#ifndef TMATRIXELEMENTOLD_H
#define TMATRIXELEMENTOLD_H

#include "FourVector.h"
#include "GammaStructure.h"
#include "Structures.h"
#include "TMatrixElement.h"

class TMatrixElementOld : public TMatrixElement
{
public:
    static TMatrixElementOld* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr); //!< TMatrixElementOld singleton factory,  
    static TMatrixElementOld* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY); //!< TMatrixElementOld singleton factory,  allowing off-shell kinematics
    virtual ~TMatrixElementOld(); //!< Destructor
    virtual void ReInitialise( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< Initialises the object as if it were newly constructed.
    
    virtual std::complex<double> calculateM(int L, int Lp, int Ly); //!< Calculates the complex value of the matrix element

    /*! Retrieve the amputated current */
    const FourVector<GammaStructure>& GetCurrent() { return current;}
protected:
    TMatrixElementOld(double w, double k, double costhkcm, double pk,Class* particlesPtr, Observable* observPtr);

    /*! Calculates all components of the amputated current*/
    void determineCurrent();
    
    // Singleton implementation.
    static bool fgInstanceFlag;
    static TMatrixElementOld * fgMatrixElement;
    
    /*! The amputated current: a 4vector whose components are complex 4x4 matrices */
    FourVector<GammaStructure> current; 

};

#endif  // #ifdef TMATRIXELEMENTOLD_H
