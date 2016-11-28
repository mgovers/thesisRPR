/*!
 * \file TCalculateCGLNMatrix.cpp
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
 * \class TCalculateCGLNMatrix
 *
 * A class to calculate and cache the CGLN amplitude basis.
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */

#include "TCalculateCGLNMatrix.h"
#include <iostream>
#include <complex>
using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

int TCalculateCGLNMatrix::fgDataSize = 0;
complex< double > * TCalculateCGLNMatrix::fgSandwichedMCache= NULL;
bool* TCalculateCGLNMatrix::fgSandwichedMIsCalculated = NULL;
std::vector<TCachemap<GammaStructure> > TCalculateCGLNMatrix::fgCachemapMVector(6) ;
//std::vector<TCachemap<complex< double > > > TCalculateCGLNMatrix::fgCachemapSandwichedMVector(6) ;

//________________________________________________________________________
/*! \brief Memoizing function to retreive M_index for certain kinematics
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * \param index index of M; 0...3 for photoproduction, 0...5 for electroproduction
 * \param cachekey key to retreive the memoized result
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
GammaStructure
TCalculateCGLNMatrix::GetM(const int index,
			   const char* cacheKey,
                           const FourVector< std::complex<double> >& epsilon,
                           const FourVector<double>& k4vect,
                           const FourVector<double>& p4vect,
                           const FourVector<double>& pY4vect)
{
    // Create a temp object to hold the result
    GammaStructure result;

    // Look up cachekey in cacheMap and return its associated value if found.
    if ( fgCachemapMVector[index].find(cacheKey,result) ) return result;

    // if cacheKey was not found in the map, we evaluate the function
    result = CalcM(index, epsilon,k4vect, p4vect, pY4vect);

    // and store its result in cacheMap under cacheKey.
    fgCachemapMVector[index].insert(cacheKey,result);
    return result;
}

/*! \brief Calculate the GammaStructure M_i
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * \param index index of M; 0...3 for photoproduction, 0...5 for electroproduction
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
GammaStructure
TCalculateCGLNMatrix::CalcM(const int index,
                            const FourVector< std::complex<double> >& epsilon,
                            const FourVector<double>& k4vect,
                            const FourVector<double>& p4vect,
                            const FourVector<double>& pY4vect)
{
    switch (index)
    {
    case 0:
        return GammaStructure(0.0,1.0) * ((GMU*k4vect)*(GMU*epsilon));
    case 1:
        return 2.0 * GammaStructure(0.0,1.0) * ((pY4vect*k4vect)*(p4vect*epsilon) - (k4vect*p4vect)*(pY4vect*epsilon));
    case 2: //This term changes sign when channel is crossed.
        return GammaStructure(0.0,1.0) * ((GMU*epsilon) * (k4vect*p4vect) - (GMU*k4vect) * (p4vect*epsilon));
    case 3: //This term changes sign when channel is crossed.
        return GammaStructure(0.0,1.0) * ((GMU*epsilon) * (k4vect*pY4vect) - (GMU*k4vect) * (pY4vect*epsilon));
    case 4:
        return GammaStructure(0.0,1.0) * (pY4vect*epsilon);
    case 5: 
        return GammaStructure(0.0,1.0) * (GMU*epsilon);
    default:
        cerr << "Error in TCalculateCGLNMatrix::CalcM(int index, ... ): index " << index << " not in {0,...,5}."<< endl;
        exit(1);
        return GammaStructure(0.0,1.0);
    }
}

// /*! \brief Memoizing function to retreive the sandwiched <Ly|M_i|Lp> for certain kinematics
//  * Note that the kaon momentum fourvector is fixed by momentum conservation.
//  * \param index index of M; 0...3 for photoproduction, 0...5 for electroproduction
//  * \param cachekey key to retreive the memoized result
//  * \param nucleonSpinor Nucleon spinor (column matrix)
//  * \param hyperonSpinor Hyperon spinor (row matrix)
//  * \param epsilon photon polarization vector in CM frame
//  * \param k4vect photon momentum fourvector in CM frame
//  * \param p4vect proton momentum fourvector in CM frame
//  * \param pY4vect hyperon momentum fourvector in CM frame
//  * \return the (scalar,complex) sandwiched M_index
//  */
// complex< double > TCalculateCGLNMatrix::GetSandwichedM(const int index, 
// 						    const char* cacheKey, 
// 						    const Matrix< 4, 1 >& nucleonSpinor, 
// 						    const Matrix< 1, 4 >& hyperonSpinor, 
// 						    const FourVector< std::complex< double > >& epsilon, 
// 						    const FourVector< double >& k4vect, 
// 						    const FourVector< double >& p4vect, 
// 						    const FourVector< double >& pY4vect)
// {
//     complex<double> result;
//     // Look up cachekey in cacheMap and return its associated value if found.
//     if ( fgCachemapSandwichedMVector[index].find(cacheKey,result) ) return result;
// 
//     // if cacheKey was not found in the map, we evaluate the function
//     result = CalcSandwichedM(index,nucleonSpinor,hyperonSpinor,epsilon,k4vect, p4vect, pY4vect);
// 
//     // and store its result in cacheMap under cacheKey.
//     fgCachemapSandwichedMVector[index].insert(cacheKey,result);
//     return result;
// }


/*! \brief Memoizing function to retreive the sandwiched <Ly|M_i|Lp> for certain kinematics
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * \param index index of M; 0...3 for photoproduction, 0...5 for electroproduction
 * \param cachelabel integer to retreive the memoized result, corresponds to a (sub-)datapoint
 * \param nucleonSpinor Nucleon spinor (column matrix)
 * \param hyperonSpinor Hyperon spinor (row matrix)
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 * \return the (scalar,complex) sandwiched M_index
 */
complex< double > TCalculateCGLNMatrix::GetSandwichedM(const int index, 
							     const int cachelabel, 
							     const int L, const int Lp, const int Ly,
							     const Matrix< 4, 1 >& nucleonSpinor, 
							     const Matrix< 1, 4 >& hyperonSpinor, 
							     const FourVector< std::complex< double > >& epsilon, 
							     const FourVector< double >& k4vect, 
							     const FourVector< double >& p4vect, 
							     const FourVector< double >& pY4vect)
{
  if (cachelabel<0) return CalcSandwichedM(index,nucleonSpinor,hyperonSpinor,epsilon,k4vect, p4vect, pY4vect);
  int cache_index = GetCacheIndex(index,cachelabel,L,Lp,Ly);
  if (!fgSandwichedMIsCalculated[cache_index])
  {
    // if the cache_index was not yet calculated, we evaluate the function
    fgSandwichedMCache[cache_index] = TCalculateCGLNMatrix::CalcSandwichedM(index,nucleonSpinor,hyperonSpinor,epsilon,k4vect, p4vect, pY4vect);
    fgSandwichedMIsCalculated[cache_index]=true;
  }
  return fgSandwichedMCache[cache_index];
}

/*! \brief Calculate the sandwiched CGLN basis element <Ly|M_i|Lp>
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * \param index index of M; 0...3 for photoproduction, 0...5 for electroproduction
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 * \return the (scalar,complex) sandwiched M_index
 */
complex< double > TCalculateCGLNMatrix::CalcSandwichedM(const int index, 
							     const Matrix< 4, 1 >& nucleonSpinor, 
							     const Matrix< 1, 4 >& hyperonSpinor, 
							     const FourVector< std::complex< double > >& epsilon, 
							     const FourVector< double >& k4vect, 
							     const FourVector< double >& p4vect, 
							     const FourVector< double >& pY4vect)

{
    return hyperonSpinor*CalcM(index,epsilon,k4vect,p4vect,pY4vect).value()*nucleonSpinor;
}

void TCalculateCGLNMatrix::SetDataSize(const int size)
{
  fgDataSize = size;
  delete [] fgSandwichedMCache;
  fgSandwichedMCache = NULL;
  delete [] fgSandwichedMIsCalculated;
  fgSandwichedMIsCalculated = NULL;
  
  if (size>0)
  {
    fgSandwichedMCache= new complex< double > [size*72]; // FIXME -- this ends up as reachable  (is never deleted)
    fgSandwichedMIsCalculated = new bool[size*72]; // FIXME -- this ends up as reachable  (is never deleted)
    for (int i=0; i<size*72; i++)
    {
      fgSandwichedMIsCalculated[i]=false;
      fgSandwichedMCache[i]=0.0;
    }
  }
}
