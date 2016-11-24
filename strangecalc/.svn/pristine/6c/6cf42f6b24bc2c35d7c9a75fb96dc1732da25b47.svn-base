/*!
 * \file TCalculateCGLNCoeff.cpp
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
 * \class TCalculateCGLNCoeff
 *
 * Abstract base class to calculate and cache the coefficients of the
 * reaction-amplitude's decomposition in the CGLN basis.
 * 
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */
#include "TCalculateCGLNCoeff.h"
#include "TCalculateConsistentCoeff.h"
#include "FourVector.h"
#include "Structures.h"
#include <iostream>
#include "numtoa.h"
#include <vector>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::complex;
using std::vector;

vector < TCachemap < complex< double > > > TCalculateCGLNCoeff::fCachemapAVector ( 6 );

/*! \brief Private constructor for TCalculateCGLNCoeff
 */
TCalculateCGLNCoeff::TCalculateCGLNCoeff()
{   
  fWidthModifier = new double[CLASSMAX];  // FIXME -- this ends up as reachable  (is never deleted)
  for(int c=0; c < CLASSMAX; c++)
    fWidthModifier[c]= 1.0;
}

/*! \brief Destructor for TCalculateCGLNCoeff
 */
TCalculateCGLNCoeff::~TCalculateCGLNCoeff()
{
  delete fWidthModifier;
}

/*! \brief Create a new calculator calcname
* By default a new consistent calculator is returned
* \param calcname Name of the specific calculator impementation
*/
TCalculateCGLNCoeff* TCalculateCGLNCoeff::GetCGLNCoeffCalculator(std::string modeltype) 
{
  if (!modeltype.compare("consistent"))
    return TCalculateConsistentCoeff::GetConsistentCalculator();
  else if (!modeltype.compare("varcutoff1"))
    return TCalculateConsistentCoeff::GetVarcutoff1Calculator();
  else if (!modeltype.compare("varcutoff2"))
    return TCalculateConsistentCoeff::GetVarcutoff2Calculator();
  else if (!modeltype.compare("lorentz") || !modeltype.compare("lorentz_gaussian"))
    return TCalculateConsistentCoeff::GetLorentzCalculator();
  else
  {
    cerr<< "Error in TCalculateCGLNCoeff::GetCGLNCoeffCalculator(string modeltype): modeltype '"<< modeltype << "' was not recognised!"<< endl;
    exit(1);
  }
  
  return TCalculateConsistentCoeff::GetConsistentCalculator();

}


//________________________________________________________________________
/*! \brief Memoizing function to retreive A_{index,classindex} for certain particle & kinematics
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * \param index index of M and of A; 0...3 for photoproduction, 0...5 for electroproduction
 * \param cacheKey key for cache lookup; this is calculated beforehand as it is needed repeatedly for the same values
 * \param classindex classindex of A, determines the diagram type (e.g. Born, 1/2+ resonance, ...)
 * \param particle Properties struct containing the coupling constants, mass, width, etc. of the exchanged particle
 * \param observ Observ struct which contains information on the observable and model details
 * \param nucleon_charge Charge of the external nucleon involved
 * \param nucleon_charge Charge of the external hyperon involved
 * \param mN Mass of the external nucleon involved
 * \param mY Mass of the external hyperon involved
 * \param mK Mass of the external kaon involved
 * \param s mandelstam-s, actually four-momentum^2 going through the N*!!
 * \param t mandelstam-t
 * \param u mandelstam-u, actually four-momentum^2 going through the Y*!!
 */
complex< double >
TCalculateCGLNCoeff::GetA ( int index, 
			    const char* cacheKey,
			    int classindex,
                            const Properties& particle,
                            const Observable& observ,
                            const double nucleon_charge,
                            const double hyperon_charge,
                            const double mN, const double mY, const double mK,
                            const double w, const double k, const double costhkcm )
{
  // Create a temp object to hold the result
  complex< double > result;

  // Look up cachekey in cacheMap and return its associated value if found.
  if ( fCachemapAVector[index].find ( cacheKey,result ) ) return result;

  // if cacheKey was not found in the map, we evaluate the function
  result = CalcA ( index,classindex,particle,observ, nucleon_charge, hyperon_charge,mN,mY,mK,w,k,costhkcm );
  // and store its result in cacheMap under cacheKey.
  fCachemapAVector[index].insert ( cacheKey,result );
  return result;
}

//________________________________________________________________________
/*! \brief Function to calculate A_{index,classindex} for certain kinematics.
 * Note that the kaon momentum fourvector is fixed by momentum conservation.
 * This function simply redirects to CalcA${index}(classindex, particle, k4vect, p4vect, pY4vect)
 * which can be specified in the derived classes e.g. TCalculateConsistentCoeff
 * \param index index of M and of A; 1...4 for photoproduction, 1...6 for electroproduction
 * \param classindex classindex of A, determines the diagram type (e.g. Born, 1/2+ resonance, ...)
 * \param particle particle struct, contains coupling constants, Mass, width.
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
complex< double >
TCalculateCGLNCoeff::CalcA ( int index, int classindex,
                             const Properties& particle,
                             const Observable& observ,
                             const double nucleon_charge,
                             const double hyperon_charge,
                             const double mN, const double mY, const double mK,
                             const double w, const double k, const double costhkcm )
{
  // Assign the mandelstam variables
  double pK;

  fkk = w*w - k*k; // (k4vect*k4vect)
    
  fmN = mN; // Mass of the nucleon
  fmY = mY; // Mass of the hyperon
  
  fMass = particle.mass; // Mass of exchanged particle
  fWidth = fWidthModifier[classindex]*particle.width; // Width of the particle. 
  // In some cases, this has to be modified with a fWidthModifier factor
    
  if (!observ.kaoncapture)
    {
      fS = pow(w + sqrt(k*k + mN*mN),2.);
      pK = sqrt(pow(fS - mK*mK - mY*mY,2.) - 4.*pow(mK*mY,2.))/(2.*sqrt(fS));
      fT = pow(w - sqrt(pK*pK + mK*mK),2.) - k*k - pK*pK + 2.*k*pK*costhkcm;
      fU = mN*mN + mK*mK + mY*mY - fS - fT + fkk;
      fkp = 0.5* ( fS - mN*mN - fkk ); // (k4vect*p4vect)
      fkpY = -0.5* ( fU - mY*mY - fkk ); // (k4vect*pY4vect)
      fDenominator_s = fS - fMass*fMass 
	+ complex<double> ( 0.0,1.0 ) *fMass*fWidth; // s-channel propagator (resonant)
      fDenominator_u = fU - fMass*fMass; // u-channel propagator (non-resonant so width==0)
    }
  else
    {
      fU = pow(w + sqrt(k*k + mY*mY),2.);
      pK = sqrt(pow(fU - mK*mK - mY*mY,2.) - 4.*pow(mK*mY,2.))/(2.*sqrt(fU));
      fT = pow(w - sqrt(pK*pK + mK*mK),2.) - k*k - pK*pK + 2.*k*pK*costhkcm;
      fS = mN*mN + mK*mK + mY*mY - fU - fT + fkk;
      fkp = -0.5* ( fS - mN*mN - fkk ); // (k4vect*p4vect)
      fkpY = 0.5* ( fU - mY*mY - fkk ); // (k4vect*pY4vect)  
      fDenominator_s = fS - fMass*fMass; // s-channel propagator (non-resonant so width==0)
      fDenominator_u = fU - fMass*fMass 
	+ complex<double> ( 0.0,1.0 ) *fMass*fWidth; // u-channel propagator (resonant)
    
    }

  fppY = -0.5* ( fT - mN*mN - mY*mY ); // (p4vect*pY4vect)
  fDenominator_t = fT - fMass*fMass; //+ complex<double> ( 0.0,1.0 ) *fMass*particle.width; // t-channel propagator
    
  fNucleon_charge= nucleon_charge;
  fHyperon_charge= hyperon_charge;
  switch ( index )
    {
    case 0:
      return CalcA1 ( classindex, particle,observ );
      break;
    case 1:
      return CalcA2 ( classindex, particle,observ );
      break;
    case 2:
      return CalcA3 ( classindex, particle,observ );
      break;
    case 3:
      return CalcA4 ( classindex, particle,observ );
      break;
    case 4:
      return CalcA5 ( classindex, particle,observ );
      break;
    case 5:
      return CalcA6 ( classindex, particle,observ );
      break;
    default: //error
      cerr << "Error in TCalculateCGLNCoeff::CalcA(int index, ... ): index " << index << " not in {0,...,5}."<< endl;
      exit(1);
      return 0;
    }
}
