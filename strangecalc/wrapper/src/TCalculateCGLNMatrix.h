/*!
 * \file TCalculateCGLNMatrix.h
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

#ifndef TCALCULATECGLNMATRIX_H
#define TCALCULATECGLNMATRIX_H

#include <Structures.h>
#include <TCachemap.h>
#include <FourVector.h>
#include <GammaStructure.h>
#include <vector>
#include <complex>

class TCalculateCGLNMatrix
{
public:
//     TCalculateCGLNMatrix();
//     virtual ~TCalculateCGLNMatrix();

    static GammaStructure GetM(const int index,
			const char* cacheKey,
                        const FourVector< std::complex<double> >& epsilon,
                        const FourVector<double>& k4vect,
                        const FourVector<double>& p4vect,
                        const FourVector<double>& pY4vect);

    static GammaStructure CalcM(const int index,
                                const FourVector< std::complex<double> >& epsilon,
                                const FourVector<double>& k4vect,
                                const FourVector<double>& p4vect,
                                const FourVector<double>& pY4vect);
    
    static std::complex<double> GetSandwichedM ( const int index, 
					   const char* cacheKey,
					   const Matrix<4,1>& nucleonSpinor,
					   const Matrix<1,4>& hyperonSpinor,
					   const FourVector< std::complex< double > >& epsilon, 
					   const FourVector< double >& k4vect, 
					   const FourVector< double >& p4vect, 
					   const FourVector< double >& pY4vect );

    static std::complex<double> GetSandwichedM ( const int index, 
						 const int cachelabel, 
						 const int L, const int Lp, const int Ly,
						 const Matrix<4,1>& nucleonSpinor,
						 const Matrix<1,4>& hyperonSpinor,
						 const FourVector< std::complex< double > >& epsilon, 
						 const FourVector< double >& k4vect, 
						 const FourVector< double >& p4vect, 
						 const FourVector< double >& pY4vect );
					   
    static std::complex< double > CalcSandwichedM(const int index,
						  const Matrix< 4, 1 >& nucleonSpinor, 
						  const Matrix< 1, 4 >& hyperonSpinor, 
						  const FourVector< std::complex< double > >& epsilon, 
						  const FourVector< double >& k4vect, 
						  const FourVector< double >& p4vect, 
						  const FourVector< double >& pY4vect);
					   
    static void SetDataSize(int size);
    
    
protected:
    static std::complex< double >* fgSandwichedMCache;
    static bool* fgSandwichedMIsCalculated;
    static int fgDataSize;
    static std::vector<TCachemap<GammaStructure> > fgCachemapMVector; //!< vector of cache containers for M_i's.
    //static std::vector<TCachemap<std::complex< double > > > fgCachemapSandwichedMVector ;
    
    //! This is static because the M_i's do not depend on particle specifics, only kinematics, so can be used program-wide.
    static inline int GetCacheIndex(const int index, const int label, const int L, const int Lp, const int Ly) 
      { return label*72 + L*24+Lp*12+Ly*6+index; }
    
};


#endif
