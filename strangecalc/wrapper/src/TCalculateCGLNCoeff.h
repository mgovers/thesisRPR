/*!
 * \file TCalculateCGLNCoeff.h
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

#ifndef TCALCULATECGLNCOEFF_H
#define TCALCULATECGLNCOEFF_H
#include <complex>
#include <vector>
#include "Structures.h"
#include "TCachemap.h"
template <class T>
class FourVector;

class TCalculateCGLNCoeff
{

public:
    virtual ~TCalculateCGLNCoeff();

    virtual std::complex< double >  GetA ( int index,
                                           const char* cacheKey,
                                           int classindex,
                                           const Properties& particle,
                                           const Observable& observ,
                                           const double nucleon_charge,
                                           const double hyperon_charge,
                                           const double mN, const double mY, const double mK,
                                           const double w, const double k, const double costhkcm );
    
    static TCalculateCGLNCoeff* GetCGLNCoeffCalculator(std::string modeltype="consistent");
    virtual std::complex< double >  CalcA ( int index, int classindex,
					    const Properties& particle,
					    const Observable& observ,
					    const double nucleon_charge,
					    const double hyperon_charge,
					    const double mN, const double mY, const double mK,
					    const double w, const double k, const double costhkcm );
    

protected:
//     virtual GammaStructure GetAM ( int index, int classindex,
//                                    const Properties& particle,
//                                    const FourVector< std::complex<double> >& epsilon,
//                                   double su, double t ) const
    TCalculateCGLNCoeff();
    
    virtual std::complex< double >  CalcA1 ( int classindex,
            const Properties& particle,
            const Observable& observ) const = 0;
    virtual std::complex< double >  CalcA2 ( int classindex,
            const Properties& particle,
            const Observable& observ) const = 0;
    virtual std::complex< double >  CalcA3 ( int classindex,
            const Properties& particle,
            const Observable& observ) const = 0;
    virtual std::complex< double >  CalcA4 ( int classindex,
            const Properties& particle,
            const Observable& observ) const=0;
    virtual std::complex< double >  CalcA5 ( int classindex,
            const Properties& particle,
            const Observable& observ) const = 0;
    virtual std::complex< double >  CalcA6 ( int classindex,
            const Properties& particle,
            const Observable& observ) const = 0;

    // mandelstam variables
    double fS;
    double fT;
    double fU;

    // kinematic variables
    double fMass; //!< mass of exchanged particle
    double fkk; //!< (k4vect*k4vect)
    double fkp; //!< (k4vect*p4vect)
    double fppY; //!< (p4vect*pY4vect)
    double fkpY; //!< (k4vect*pY4vect)
    
    //propagators
    std::complex<double> fDenominator_s; //!< s-channel propagator
    std::complex<double> fDenominator_u; //!< u-channel propagator
    std::complex<double> fDenominator_t; //!< t-channel propagator
    
    //N,K,Y data
    double fNucleon_charge;
    double fHyperon_charge;
    double fmN;
    double fmY;
    double fmK;
    
    double fWidth; //!< width, in the case of a variable width (Lorentz formfactor)
    double* fWidthModifier;

private:
    //static std::vector< std::vector < TCachemap<std::complex< double > > > > fCachemapAMVector;
    static std::vector < TCachemap < std::complex< double > > > fCachemapAVector; //!< CLASSMAX x 6 matrix of cache containers for A^n_i*M_i's.


//std::vector< std::vector < TCachemap<GammaStructure> > > TCalculateCGLNCoeff::fCachemapAMVector ( CLASSMAX, std::vector < TCachemap < GammaStructure > > ( 6 ) );

// define some const members
};
#endif
