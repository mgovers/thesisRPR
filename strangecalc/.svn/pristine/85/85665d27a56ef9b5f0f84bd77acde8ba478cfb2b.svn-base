/*!
 * \file TCalculateConsistentCoeff.h
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

#ifndef TCALCULATECONSISTENTCOEFF_H
#define TCALCULATECONSISTENTCOEFF_H

#include "TCalculateCGLNCoeff.h"

class TCalculateConsistentCoeff : public TCalculateCGLNCoeff
{
public:
    static TCalculateConsistentCoeff* GetConsistentCalculator();
    static TCalculateConsistentCoeff* GetVarcutoff1Calculator();
    static TCalculateConsistentCoeff* GetVarcutoff2Calculator();
    static TCalculateConsistentCoeff* GetLorentzCalculator();
    virtual ~TCalculateConsistentCoeff();
    
protected:
    TCalculateConsistentCoeff(double exponent=0.0);
    std::complex< double >  CalcA1 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const;
    std::complex< double >  CalcA2 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const ;
    std::complex< double >  CalcA3 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const ;
    std::complex< double >  CalcA4 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const ;
    std::complex< double >  CalcA5 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const ;
    std::complex< double >  CalcA6 ( int classindex,
                                     const Properties& particle,
                                     const Observable& observ ) const ;
    double fExponent;
private:
    static bool fgInstanceFlag; //!< Singleton implementation.
    static TCalculateConsistentCoeff* fgCalculator;

};

#endif
