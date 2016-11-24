/*!
 * \file TCurrent.cpp
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
 * \class TCurrent
 *
 * A class to calculate the contracted current using a 
 * CGLN matrix decomposition and coefficient calculator.
 * 
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */

#include "TCurrent.h"
#include "TCalculateConsistentCoeff.h"
#include "FourVector.h"
#include "GammaStructure.h"
#include <complex>
#include <iostream>
#include "Matrix.h"
#include <numtoa.h>
using std::complex;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

typedef vector<int>::size_type size_type;

const GammaStructure TCurrent::kZeroGS = GammaStructure ( 0. );

/*! \brief Constructor with specific coefficient calculator.
 * \param coeffCalcPtr pointer to a TCalculateCGLNCoeff object
 * \param particlesPtr Pointer to array containing the Properties structs with particle properties
 * \param observPtr Observ struct which contains information on the observable and model details
 * \param s mandelstam-s, actually four-momentum^2 going through the N*!!
 * \param t mandelstam-t
 * \param u mandelstam-u, actually four-momentum^2 going through the Y*!!
 * 
 */
TCurrent::TCurrent (TCalculateCGLNCoeff* coeffCalcPtr, Class* particlesPtr, Observable* observPtr, const double w, const double k, const double costhkcm) 
: fCoeffCalcPtr(coeffCalcPtr)
, fParticlesPtr(particlesPtr)
, fObservPtr(observPtr)
, fOmega(w)
, fK(k)
, fCosthkcm(costhkcm)
, fCoefficients( vector<complex<double> > (6,complex<double>(0.0)))
{
  DetermineCoefficients();
}

/*! \brief Getter for the coefficient calculator pointer.
 */
TCalculateCGLNCoeff* TCurrent::GetCoeffCalculator() const
{
  return fCoeffCalcPtr;
}

/*! \brief Add TCurrents
 * Intelligently add two TCurrents simply by adding the coefficients.
 * \param other TCurrent on the right-hand side of the + operator
 */
TCurrent & TCurrent::operator += ( const TCurrent & other )
{
  for (size_type i=0;i<fCoefficients.size(); i++ )
    {
      fCoefficients[i] += other.GetCoefficient ( i );
    }
  return *this;
}

/*! \brief Return the coefficient A_index.
 * \param index Index of the matching M for this coefficient.
 */
complex<double> TCurrent::GetCoefficient ( int index ) const
{
  return fCoefficients[index];
}

/*! \brief Add the CGLN coefficients of a diagram to fCoefficients
 * \param classindex Classindex corresponding to the type of exchange diagram
 */
void TCurrent::AddDiagram ( const int classindex, const Properties& particle)

{
  for (size_type i=0;i<fCoefficients.size(); i++ )
    fCoefficients[i] += fCoeffCalcPtr->CalcA ( i,classindex,particle,*fObservPtr,
					       fParticlesPtr[0].partic[0].E, // nucleon_charge
					       fParticlesPtr[2].partic[0].E, // hyperon_charge
					       fParticlesPtr[0].partic[0].mass, // mN
					       fParticlesPtr[2].partic[0].mass, // mY
					       fParticlesPtr[1].partic[0].mass, // mK
					       fOmega,fK,fCosthkcm);
}

/*! \brief Clear the TCurrent coefficients and fCurrentGS
 *
 */
void TCurrent::Clear()
{
  for (size_type i=0;i<fCoefficients.size(); i++ )
    {
      fCoefficients[i] = 0.0;
    }
  fCurrentGS = kZeroGS;
}
/*! \brief (Re-)initialise the TCurrent.
 * Essentially does the same as the constructor but with an existing object.
 * \param coeffCalcPtr pointer to a TCalculateCGLNCoeff object
 * \param particlesPtr Pointer to array containing the Properties structs with particle properties
 * \param observPtr Observ struct which contains information on the observable and model details
 * \param s mandelstam-s, actually four-momentum^2 going through the N*!!
 * \param t mandelstam-t
 * \param u mandelstam-u, actually four-momentum^2 going through the Y*!!
 */
void TCurrent::Initialise(TCalculateCGLNCoeff* coeffCalcPtr, Class* particlesPtr, Observable* observPtr, const double w, const double k, const double costhkcm)
{
  fCoeffCalcPtr = coeffCalcPtr;
  fParticlesPtr = particlesPtr;
  fObservPtr = observPtr;
  fOmega = w;
  fK = k;
  fCosthkcm = costhkcm;
  DetermineCoefficients();
}
  
/*! \brief Calculate and return the GammaStructure that contains the current.
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
const GammaStructure& TCurrent::GetCurrentGS(const FourVector< std::complex<double> >& epsilon, // fEpsilon[L]  or // kUnitEpsilon[mu]
                                  const FourVector<double>& k4vect,
                                  const FourVector<double>& p4vect,
                                  const FourVector<double>& pY4vect)
{
  DetermineCurrentGS(epsilon,k4vect,p4vect,pY4vect);
  return fCurrentGS;
}

/*! \brief Determine the coefficients for given s/t and particle/fObservPtr info, using fCoeffCalculator.
 */
void TCurrent::DetermineCoefficients()
{
  // Clear the fCoefficients array
  Clear();
  // Loop over all types of diagrams (see Lagrangian.h for a list)
  for ( int diagram=0; diagram<CLASSMAX; diagram++ )
  {
    // Loop over all exchanged fParticlesPtr for a specific diagram
    for ( int particle=0; particle<fParticlesPtr[diagram].particount; particle++ )
    {
      AddDiagram ( diagram,fParticlesPtr[diagram].partic[particle]);
    } // particle loop
  } // diagram loop
}

/*! \brief Calculate the current GammaStructure using fCoefficients of the CGLN basis.
 * based on the CGLN coefficients and the CGLN matrices. The kinematicsKey is 
 * calculated here as the function GetM is called repeatedly for the same kinematics.
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
void TCurrent::DetermineCurrentGS ( const FourVector< std::complex<double> >& epsilon,
                                    const FourVector<double>& k4vect,
                                    const FourVector<double>& p4vect,
                                    const FourVector<double>& pY4vect)
{
  // set the currentGS to GS containing only zeroes
  fCurrentGS = kZeroGS;
    
  // make array for lookup key
  float kinematicsArray[fKinDF]= {(float)fOmega,(float)fK,(float)fCosthkcm,
  				  (float)epsilon[0].real(),(float)epsilon[1].real(),(float)epsilon[2].real(),(float)epsilon[3].real(),
  				  (float)epsilon[0].imag(),(float)epsilon[1].imag(),(float)epsilon[2].imag(),(float)epsilon[3].imag()};

  // Create unique key for cacheMap lookups
  char kinematicsKey[sizeof(float)*fKinDF*2 + 1];
  floatstochar(fKinDF, kinematicsArray, kinematicsKey);

  for ( int i=0;i< fCoefficients.size(); i++ )
  {
    fCurrentGS = fCurrentGS + TCalculateCGLNMatrix::GetM(i,kinematicsKey,epsilon,k4vect, p4vect, pY4vect) *fCoefficients[i];
  }
}


/*! \brief Determine the sandwiched current 
 * \param L Photon polarisation vector lambda (0=-1, 1=0, 2=+1)
 * \param Lp Nucleon spin (0=-1/2, 1=+1/2)
 * \param Ly Hyperon spin (0=-1/2, 1=+1/2)
 * \param nucleonSpinor Nucleon spinor (column matrix)
 * \param hyperonSpinor Hyperon spinor (row matrix)
 * \param epsilon photon polarization vector in CM frame
 * \param k4vect photon momentum fourvector in CM frame
 * \param p4vect proton momentum fourvector in CM frame
 * \param pY4vect hyperon momentum fourvector in CM frame
 */
complex< double > TCurrent::DetermineSandwichedCurrent(const int cachelabel,
						       const int L, const int Lp, const int Ly,
						       const Matrix< 4, 1 >& nucleonSpinor, 
						       const Matrix< 1, 4 >& hyperonSpinor, 
						       const FourVector< std::complex< double > >& epsilon, 
						       const FourVector< double >& k4vect, 
						       const FourVector< double >& p4vect, 
						       const FourVector< double >& pY4vect)
{
  complex< double > sandwichedCurrent = 0.0;
  
  for ( size_type i=0;i< fCoefficients.size(); i++ )
    sandwichedCurrent += TCalculateCGLNMatrix::GetSandwichedM(i,cachelabel,L,Lp,Ly,nucleonSpinor,hyperonSpinor,epsilon,k4vect, p4vect, pY4vect) *fCoefficients[i];
  return sandwichedCurrent;
}
