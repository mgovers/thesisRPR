/*!
 * \file TMatrixElementCGLN.cpp
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

#include "TMatrixElementCGLN.h"
#include <string>
#include <iostream>
using std::vector;
using std::complex;
using std::string;
using std::cout;
using std::endl;

/*!
 * \class TMatrixElementCGLN
 *
 * Class to calculate a matrix element for given kinematics,
 * photon polarisation, nucleon and hyperon spin.
 * The current contracted with epsilon is calculated
 * using its decomposition in the CGLN basis.
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 */

vector< FourVector < std::complex<  double > > > GetUnitEpsilonFV() {
    std::vector < FourVector < complex< double > > > unitEpsilon;
    unitEpsilon.push_back(FourVector< complex <double> >(1.0,0.0,0.0,0.0));
    unitEpsilon.push_back(FourVector< complex <double> >(0.0,-1.0,0.0,0.0));
    unitEpsilon.push_back(FourVector< complex <double> >(0.0,0.0,-1.0,0.0));
    unitEpsilon.push_back(FourVector< complex <double> >(0.0,0.0,0.0,-1.0));
    return unitEpsilon;
}

// initialise statics
const vector< FourVector<complex< double > > > TMatrixElementCGLN::kUnitEpsilonFV = GetUnitEpsilonFV();
bool TMatrixElementCGLN::fgInstanceFlag = false;
TMatrixElementCGLN * TMatrixElementCGLN::fgMatrixElement = NULL;

/*! \brief Constructor for TMatrixElementCGLN
* determines which TCalculateCGLNCoeff to use based on modelType in observPtr.
* \param w photon energy in the CMF
* \param k photon momentum in the CMF
* \param costhkcm cos(theta_K) in the CMF
* \param pK Kaon momentum in the CMF
* \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
* \param observPtr pointer to the Observ struct
* \param label integer label of the datapoint, corresponds to given kinematics and functions as a key for caching the CGLN matrices.
*/
TMatrixElementCGLN::TMatrixElementCGLN ( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label) 
: TMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr)
, fCalculatorPtr(TCalculateCGLNCoeff::GetCGLNCoeffCalculator(observPtr->modelType))
, fCurrent(TCurrent(fCalculatorPtr, particlesPtr,observPtr, w, k, costhkcm))
, fDataLabel(label)
, fBareCurrentIsCalculated(false)
{
}

/*! \brief Constructor for TMatrixElementCGLN, allowing off-shell kinematics
 *  which determines which TCalculateCGLNCoeff to use based on modelType in observPtr
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param kin_mN Nucleon mass (used only for off-shell kinematics!)
 * \param kin_mK Kaon mass (used only for off-shell kinematics!)
 * \param kin_mY Hyperon mass (used only for off-shell kinematics!)
 * \param label integer label of the datapoint, corresponds to given kinematics and functions as a key for caching the CGLN matrices.
 */
TMatrixElementCGLN::TMatrixElementCGLN ( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label) 
: TMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr)
, fCalculatorPtr(TCalculateCGLNCoeff::GetCGLNCoeffCalculator(observPtr->modelType))
, fCurrent(TCurrent(fCalculatorPtr, particlesPtr,observPtr, w, k, costhkcm))
, fDataLabel(label)
, fBareCurrentIsCalculated(false)
{
}

/*! \brief Destructor for TMatrixElementCGLN.
 * The CalculatorPtr is deleted.
 */
TMatrixElementCGLN::~TMatrixElementCGLN()
{
    delete fCalculatorPtr;
    fCalculatorPtr=NULL;
    fgInstanceFlag=false;
}

/*! \brief Singleton factory for TMatrixElementCGLN
 * Creates or gets the unique TMatrixElementCGLN instance and initialises it.
 * Leaves the object in a state as if it were newly constructed with the same arguments.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param label integer label of the datapoint, corresponds to given kinematics and functions as a key for caching the CGLN matrices.
 */
TMatrixElementCGLN* TMatrixElementCGLN::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)
{
  if(! fgInstanceFlag)
  {
    fgMatrixElement = new TMatrixElementCGLN(w,k,costhkcm,pK,particlesPtr,observPtr, label); // FIXME -- this ends up as reachable  (is never deleted)
    fgInstanceFlag= true;
    return fgMatrixElement;
  }
  else
  {
    fgMatrixElement->ReInitialiseWithMasses(w,k,costhkcm,pK,particlesPtr,observPtr, 
	particlesPtr[0].partic[0].mass, particlesPtr[1].partic[0].mass, particlesPtr[2].partic[0].mass,label);
    return fgMatrixElement;
  }
}
/*! \brief Singleton factory for TMatrixElementCGLN, allowing off-shell kinematics
 * Creates or gets the unique TMatrixElementCGLN instance and initialises it.
 * Leaves the object in a state as if it were newly constructed with the same arguments.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param kin_mN Nucleon mass (used only for off-shell kinematics!)
 * \param kin_mK Kaon mass (used only for off-shell kinematics!)
 * \param kin_mY Hyperon mass (used only for off-shell kinematics!)
 * \param label integer label of the datapoint, corresponds to given kinematics and functions as a key for caching the CGLN matrices.
 */
TMatrixElementCGLN* TMatrixElementCGLN::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label)
{
  if(! fgInstanceFlag)
  {
    fgMatrixElement = new TMatrixElementCGLN(w,k,costhkcm,pK,particlesPtr,observPtr,kin_mN, kin_mK, kin_mY, label); // FIXME -- this ends up as reachable  (is never deleted)
    fgInstanceFlag= true;
    return fgMatrixElement;
  }
  else
  {
    fgMatrixElement->ReInitialiseWithMasses(w,k,costhkcm,pK,particlesPtr,observPtr,kin_mN, kin_mK, kin_mY, label);
    return fgMatrixElement;
  }
}



/*! \brief (Re-)initialises the TMatrixElementCGLN as if it were newly constructed.
 * Essentially does the same as the constructor but with an existing object.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param label integer label of the datapoint, corresponds to given kinematics and functions as a key for caching the CGLN matrices.
 */
void TMatrixElementCGLN::ReInitialise( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)
{
  // initialise the TMatrixElement (kinematics, spinors, etc.)
  TMatrixElement::ReInitialise(w,k,costhkcm,pK,particlesPtr,observPtr);
  // set the datalabel
  SetDataLabel(label);  
  
  // set the correct calculator
  fCalculatorPtr = TCalculateCGLNCoeff::GetCGLNCoeffCalculator(observPtr->modelType);
  
  // initialise the current
  fCurrent.Initialise(fCalculatorPtr,fParticlesPtr,fObservPtr,w,k,costhkcm);
  
  // initialise the bare current
  fBareCurrentIsCalculated=false;
  for ( int mu=0; mu < 4; mu++)
    fBareCurrent[mu] = kZeroGS;
}

/*! \brief Calculates the actual complex value of the matrix element
 * for fixed photon and baryon polarisations L,Lp and Ly
 * This is done by calculating the sandwiched, contracted current with the L-polarisation vector fEpsilon[L]
 * and the spinors fNucleonSpinors[Lp] and fHyperonSpinors[Ly] 
 * \param L Photon polarisation
 * \param Lp Nucleon polarisation
 * \param Ly Hyperon polarisation
 */
complex<double> TMatrixElementCGLN::calculateM(int L, int Lp, int Ly)
{
    // if this matrixelement hasn't been determined before, proceed
    if (!fMEisCalculated[L][Lp][Ly])
    {
        // never calculate it again
        fMEisCalculated[L][Lp][Ly] = true;

        // Calculate the baryon spinors when necessary
        determineSpinors(Lp,Ly);
        /* We are now ready to calculate the matrix element
         * =(hyperonSpinor[Ly] * ContractedCurrent[L] * nucleonSpinor[Lp] )
         * We'll not use the 4vector product here, because that way all components
         * of the right-side 4vector would be evaluated. Since many components of
         * epsilon are zero, there would be a large amount of redundant calculations
         */

	fMatrixElement[L][Lp][Ly] = fCurrent.DetermineSandwichedCurrent(
	  fDataLabel,L,Lp,Ly,
	  fNucleonSpinors[Lp],fHyperonSpinors[Ly],
	  fEpsilon[L],fk4vect,fp4vect,fpY4vect);
    }
    return fMatrixElement[L][Lp][Ly];
}

/*! \brief Tests whether the amplitude is gauge invariant for a given nucleon and hyperon polarization
 * This implementation immediately replaces the epsilon of the M_i CGLN amplitudes by k4vect.
 * \param Lp Nucleon polarization
 * \param Ly Hyperon polarization
 */
bool TMatrixElementCGLN::isGaugeInvariantForSpins(int Lp, int Ly)
{
    // Calculate the baryon spinors when necessary
    determineSpinors(Lp,Ly);

    // The contraction of fk4vect with J is calculated by replacing epsilon by k4vect in the calculation 
    // of the contracted current \epsilon_\mu J^\mu. This should be the case per construction.
    FourVector< complex< double > >  k4vect(fk4vect[0],fk4vect[1],fk4vect[2],fk4vect[3]);    
    TCurrent kdotCurrent(fCalculatorPtr, fParticlesPtr, fObservPtr, fOmega, fK, fCosthkcm);
    GammaStructure kdotCurrentGS = kdotCurrent.GetCurrentGS(k4vect,fk4vect,fp4vect,fpY4vect);

    complex< double > test = fHyperonSpinors[Ly]*kdotCurrentGS.value()*fNucleonSpinors[Lp];

    // Return true when both real and imaginary part of test are <10^-10
    return ( (fabs(real(test))<=1e-10) && (fabs(imag(test))<=1e-10) );

}

/*! \brief Gets the hadronic four-current J^\mu in its GammaStructure form.
 * This is done by calculating the contracted current with an photon polarisation vector epsilon
 * containing 1.0 in the mu'th position and zeros elsewhere, for mu=0..3
  */
const FourVector<GammaStructure>& TMatrixElementCGLN::GetCurrent()
{
  if (!fBareCurrentIsCalculated)
    fBareCurrentIsCalculated=true;
    for (int mu=0;mu<4;mu++)
    {
      fBareCurrent[mu] = fCurrent.GetCurrentGS(kUnitEpsilonFV[mu],fk4vect,fp4vect,fpY4vect);
    }
  return fBareCurrent;
}
