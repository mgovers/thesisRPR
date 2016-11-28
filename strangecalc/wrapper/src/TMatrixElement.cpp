/*!
 * \file TMatrixElement.cpp
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

#include "TMatrixElement.h"
#include "TMatrixElementCGLN.h"
#include "TMatrixElementOld.h"
using std::vector;
using std::complex;
using std::cout;
using std::endl;

/*!
 * \class TMatrixElement
 *
 * Abstract base class to calculate a matrix element for given kinematics,
 * photon polarisation, nucleon and hyperon spin.
 *
 * This class represents a general matrix element M^{Lp,Ly}_{L}
 * with
 *\verbatim
  Lp = nucleon polarization
  Ly = hyperon polarization
  L  = photon polarization
 \endverbatim
 *
 * A general matrix element can be decomposed as:
 * - photon polarization vector (one for each polarization)
 * - 2 spinors (incoming nucleon and outgoing hyperon, one for each polarization)
 * - An amputated current
 *
 * The main feature of this class is the function calculateM():
 * Given specific kinematics and Lp,Ly and L this function returns
 * the complex value of the matrix element
 * (Note that for electroproduction the returned value is NOT the
 * full amplitude, but only the hadronic part)
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */

const Matrix<4,1> TMatrixElement::kZero41Matrix; //!< 4x1 matrix with all elements zero
const Matrix<1,4> TMatrixElement::kZero14Matrix; //!< 1x4 matrix with all elements zero
const GammaStructure TMatrixElement::kZeroGS(0.0); //!< zero-valued GammaStructure


/*! \brief Constructor for the base class TMatrixElement
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 */
TMatrixElement::TMatrixElement ( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr )
 : mN(particlesPtr[0].partic[0].mass) // this is the proton or neutron mass
 , mK(particlesPtr[1].partic[0].mass) // this is the mass of K+ or K0
 , mY(particlesPtr[2].partic[0].mass) // this is the mass of S+- or S0
 , fOmega(w)
 , fK(k)
 , fCosthkcm(costhkcm)
 , fPK(pK)
 , fParticlesPtr(particlesPtr)
 , fObservPtr(observPtr)
 , fMEisCalculated(vector<vector<vector<bool> > > (3,vector<vector<bool> >(2, vector<bool>(2,false))))
 , fNucleonSpinorIsCalculated(vector< bool> (2,false))
 , fHyperonSpinorIsCalculated(vector< bool> (2,false))
 , fEpsilon(vector<FourVector<complex< double > > > (3,FourVector<complex< double> >()))
 , fNucleonSpinors(vector<Matrix<4,1>  > (2,kZero41Matrix))
 , fHyperonSpinors(vector<Matrix<1,4>  > (2,kZero14Matrix))
 , fMatrixElement(vector< vector< vector< complex< double > > > >  (3,vector< vector< complex< double> > >( 2,vector< complex < double> >(2 ,complex<double>(0.0) ) ) ))
{

    CalculateKinematics();
    determineEpsilon();
}

/*! \brief Constructor for the base class TMatrixElement that supports off-shell kinematics
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param kin_mN Nucleon mass (used only for off-shell kinematics!)
 * \param kin_mK Kaon mass (used only for off-shell kinematics!)
 * \param kin_mY Hyperon mass (used only for off-shell kinematics!)
 */
TMatrixElement::TMatrixElement ( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY)
: mN(kin_mN)
, mK(kin_mK)
, mY(kin_mY)
, fOmega(w)
, fK(k)
, fCosthkcm(costhkcm)
, fPK(pK)
, fParticlesPtr(particlesPtr)
, fObservPtr(observPtr)
, fMEisCalculated(vector<vector<vector<bool> > > (3,vector<vector<bool> >(2, vector<bool>(2,false))))
, fNucleonSpinorIsCalculated(vector< bool> (2,false))
, fHyperonSpinorIsCalculated(vector< bool> (2,false))
, fEpsilon(vector<FourVector<complex< double > > > (3,FourVector<complex< double> >()))
, fNucleonSpinors(vector<Matrix<4,1>  > (2,kZero41Matrix))
, fHyperonSpinors(vector<Matrix<1,4>  > (2,kZero14Matrix))
, fMatrixElement(vector< vector< vector< complex< double > > > >  (3,vector< vector< complex< double> > >( 2,vector< complex < double> >(2 ,complex<double>(0.0) ) ) ))
{
  CalculateKinematics();
  determineEpsilon();
}

/*! \brief Singleton factory for the base class TMatrixElement
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param label (optional) data label for caching of <|M_i|> cgln 
 */
TMatrixElement* TMatrixElement::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label) 
{
//  cout << "TMatrixElement* TMatrixElement::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)\n"; //DEBUG 
  if(observPtr->cgln)
  {
    return TMatrixElementCGLN::GetMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr,label);
  }
  else
  {
    return TMatrixElementOld::GetMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr);
  }
}

/*! \brief Singleton factory for the base class TMatrixElement that supports off-shell kinematics
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 * \param kin_mN Nucleon mass (used only for off-shell kinematics!)
 * \param kin_mK Kaon mass (used only for off-shell kinematics!)
 * \param kin_mY Hyperon mass (used only for off-shell kinematics!)
 * \param label (optional) data label for caching of <|M_i|> cgln 
 */
TMatrixElement* TMatrixElement::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr,double kin_mN, double kin_mK, double kin_mY, int label) 
{
  //  cout << "TMatrixElement* TMatrixElement::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)\n"; //DEBUG 
  if(observPtr->cgln)
  {
    return TMatrixElementCGLN::GetMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr,kin_mN, kin_mK, kin_mY, label);
  }
  else
  {
    return TMatrixElementOld::GetMatrixElement(w,k,costhkcm,pK,particlesPtr,observPtr,kin_mN, kin_mK, kin_mY);
  }
}

/*! \brief Calculates all kinematical quantities and initialises vectors and matrices.
 * - Initialises all bool vectors to false;
 * - Initialises the matrix elements etc. to zero
 * - Initialises kinematics and assigns particles and observ pointers.
*/
void TMatrixElement::CalculateKinematics()
{
    /* Initialize kinematics */
    if ( !fObservPtr->kaoncapture )
    {
        fk4vect = FourVector<double> ( fOmega,0,0,fK );             // photon 4vector
        fp4vect = FourVector<double> ( sqrt ( fK*fK+mN*mN ),0,0,-fK ); // nucleon 4vector
        fpK4vect = FourVector<double> ( sqrt ( fPK*fPK+mK*mK ),
                                        fPK*sqrt ( 1-fCosthkcm*fCosthkcm ),
                                        0,fPK*fCosthkcm );      // kaon 4vector
        fpY4vect = FourVector<double> ( sqrt ( fPK*fPK+mY*mY ),
                                        -1*fPK*sqrt ( 1-fCosthkcm*fCosthkcm ),
                                        0,-1*fPK*fCosthkcm );   // hyperon 4vector
    }
    else /* Use crossing symmetry: kaon and proton along z-axis */
    {
        fk4vect = FourVector<double> ( fOmega,fK*sqrt ( 1-fCosthkcm*fCosthkcm ),
                                       0,fK*fCosthkcm );           // photon 4vector
        fp4vect = FourVector<double> ( sqrt ( fPK*fPK+mN*mN ),0,0,-fPK );   // nucleon 4vector
        fpK4vect = FourVector<double> ( sqrt ( fPK*fPK+mK*mK ),0,0,fPK ); // kaon 4vector
        fpY4vect = FourVector<double> ( sqrt ( fK*fK+mY*mY ),
                                        -1*fK*sqrt ( 1-fCosthkcm*fCosthkcm ),
                                        0,-1*fK*fCosthkcm );   // hyperon 4vector
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
 */
void TMatrixElement::ReInitialise( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)
{
  fOmega = w;
  fK = k;
  fCosthkcm = costhkcm;
  fPK = pK;
  fParticlesPtr = particlesPtr;
  fObservPtr = observPtr;
  // Initialise fMatrixElement
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<2;j++)
    {
      for (int l=0;l<2;l++)
      {
	fMEisCalculated[i][j][l]=false;
	fMatrixElement[i][j][l]=0.0;
      }
    }
  }
  
  //Initialise spinors
  for(int l=0;l<2;l++)
  {
    fNucleonSpinorIsCalculated[l]=false;
    fHyperonSpinorIsCalculated[l]=false;
    fNucleonSpinors[l] = kZero41Matrix;
    fHyperonSpinors[l] = kZero14Matrix;
    
  }
  //Initialise kinematical variables
  CalculateKinematics();
  
  //Initialise epsilon_L
  determineEpsilon();
}

/*! \brief (Re-)initialises the TMatrixElementCGLN as if it were newly constructed, allowing off-shell kinematics.
 * Essentially does the same as the constructor but with an existing object.
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
void TMatrixElement::ReInitialiseWithMasses( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr,double kin_mN, double kin_mK, double kin_mY, int label)
{
  mN=kin_mN;
  mK=kin_mK;
  mY=kin_mY;
  ReInitialise(w,k,costhkcm,pK,particlesPtr,observPtr,label); // this will be the ReInitialise of the derived class!
}


/*! \brief Destructor.
*
*/
TMatrixElement::~TMatrixElement() {}


/*! \brief Tests whether the amplitude is gauge invariant for a given nucleon and hyperon polarization.
 * \param Lp Nucleon polarization
 * \param Ly Hyperon polarization
 */
bool TMatrixElement::isGaugeInvariantForSpins(int Lp, int Ly)
{
    // Calculate the baryon spinors when necessary
    determineSpinors(Lp,Ly);

    /* Now we will calculate k^(mu)*M^(Lp,Ly)_(mu) and store it in test.
     * For the reason stated in calculateM() we do not use the 4vector product
     */

    complex<double> test(0.0,0.0); // The contraction of fk4vect with M

    test += ( fk4vect[0]==complex<double>(0,0) ? 0 :
              fk4vect[0] * ( fHyperonSpinors[Ly] *
                             GetCurrent()[0].value() *
                             fNucleonSpinors[Lp] ) );

    for (int mu=1; mu<=3; mu++)
    {
        test -= ( fk4vect[mu]==complex<double>(0,0) ? 0 :
                  fk4vect[mu] * ( fHyperonSpinors[Ly] *
                                  GetCurrent()[mu].value() *
                                  fNucleonSpinors[Lp] ) );
    }


    /* Return true when both real and imaginary part of test
     * are <10^-10
     */
  
    return ( (fabs(real(test))<=1e-10) && (fabs(imag(test))<=1e-10) );
}

/*! \brief Tests whether the amplitude is gauge invariant for all nucleon and hyperon polarizations.
 */
bool TMatrixElement::isGaugeInvariant()
{
    for (int Lp=0; Lp<2; ++Lp)
        for (int Ly=0; Ly<2; ++Ly)
            if (!isGaugeInvariantForSpins(Lp,Ly)) 
            {
	      return false;
	    }

    return true;
}


/*! \brief Sets all calculated-flags to false.
 * In case settings in the spin specification part of observ change, all
 * spin-dependent quantities need to be re-evaluated.
 * In practice we:
 * - redetermine the photon polarization vectors
 * - set fNucleonSpinorIsCalculated, fHyperonSpinorIsCalculated and fMEisCalculated to false.
 */
void TMatrixElement::updateSpinDependencies()
{
    // redetermine the photon polarizations
    determineEpsilon();

    // set all calc*'s to false
    for (int i=0; i<2; i++)
    {
        fNucleonSpinorIsCalculated[i] = false;
        fHyperonSpinorIsCalculated[i] = false;
    }
    for (int i=0; i<3; i++)
        for (int j=0; j<2; j++)
            for (int k=0; k<2; k++)
                fMEisCalculated[i][j][k] = false;

}

/*! \brief Calculates the spinor wavefunction of a particle.
 * We use the definition of Peskin and Schroeder p.803
 * Spinors are normalised to 2m.
 * \param fp4vect energy-momentum 4vector:
 * \param theta zenital angle pinning down the quantization axis in the CMF
 * \param phi azimuthal angle pinning down the quantization axis in the CMF
 * \param L polarization of the particle (0 corresponds to -1/2, 1 to +1/2)
 */
Matrix<4,1> TMatrixElement::CalculateSpinor(const FourVector<double>& fp4vect,
        double theta, double phi, int L)
{
    Matrix<4,1> spinor; // Initialize a 4x1 matrix containing zeros


    /* Define the 3 Pauli matrices
     * ************************/
    static Matrix<2,2> sigmaPauli1(0.0,1.0,1.0,0.0);
    static Matrix<2,2> sigmaPauli2(0.0,
                                   complex<double>(0,-1.0),
                                   complex<double>(0,1.0),
                                   0.0);
    static Matrix<2,2> sigmaPauli3(1.0,0.0,0.0,-1.0);

    /* Determine the 2-component spinor which
     * depends upon the polarization L and the
     * quantization axis defined by theta and phi.
     * along z-axis (defined by the incoming photon) this
     * spinor is either (1,0) or (0,1) depending on L.
     * For other polarization axis, the spinor
     * is obtained through combinations
     * with Wigner D-matrices
     * See my(=Pieter's) notes
     * ********************************************/
    Matrix<2,1> spinQuantSpinor;

    // polarization = -1/2
    if (L==0)
    {
        spinQuantSpinor(0,0) = -sin(theta/2.0)
                               *exp(complex<double>(0,-1.0*phi));
        spinQuantSpinor(1,0) = cos(theta/2.0);
    }
    // polarization = +1/2
    else if (L==1)
    {
        spinQuantSpinor(0,0) = cos(theta/2.0);
        spinQuantSpinor(1,0) = sin(theta/2.0)
                               *exp(complex<double>(0,phi));
    }

    /* Determine the baryon spinor using the
     * definition of Bjorken and Drell or Gross
     * ********************************************/
    Matrix<2,2> pSigma
    = fp4vect[1] * sigmaPauli1
      + fp4vect[2] * sigmaPauli2
      + fp4vect[3] * sigmaPauli3;
    double norm = sqrt(fp4vect[0]+sqrt(fp4vect*fp4vect)); // on-mass shell 4vectors

    // the first two components of the spinor
    spinor(0,0) = norm * spinQuantSpinor(0,0);
    spinor(1,0) = norm * spinQuantSpinor(1,0);

    // the two last components of the spinor
    spinQuantSpinor = pSigma * spinQuantSpinor;
    spinor(2,0) = spinQuantSpinor(0,0) / norm;
    spinor(3,0) = spinQuantSpinor(1,0) / norm;

    return spinor;
}


/*! \brief Calculates the polarization vectors epsilon_L.
 * Assign the correct polarization vectors to each component of fEpsilon[L]
 * L is the photon polarization + 1
 *
 */
void TMatrixElement::determineEpsilon()
{
    /* Assign L=0 polarization vector
     * *****************************
     * The photon only has a longitudinal component
     * in electroproduction (when Q^2!=0) */
    double sqrtQ2 = sqrt(fK*fK-fOmega*fOmega);
    if (fObservPtr->electroprod && sqrtQ2!=0.0)
        fEpsilon[1] = FourVector< complex<double> >(fK/sqrtQ2,0,0,fOmega/sqrtQ2);
    else
        fEpsilon[1] = FourVector< complex<double> >(0,0,0,0);


    /* Assign transverse polarization vectors
     * **************************************
     * circular polarization for electroproduction
     *                           circ. double pol. photoprod
     * linear polarization in all other cases
     * Attention: lin. double pol. photoprod has a different
     *            definition for the polarization vectors
     *
     * radiative kaon capture: the kaon travels along the z-axis,
     * but the photon does not! Make sure that the spatial part of the
     * polarization vectors are orthogonal to the spatial part of the
     * photon momentum \vec{k} = w * (sintheta_k, 0, fCosthkcm).
     *
     */

    // Circular polarization
    if ( (fObservPtr->photoprod && (fObservPtr->photo.pol.dpol
                                    && fObservPtr->photo.pol.doubpol.circbeam))
            || fObservPtr->electroprod)
    {
        // L=+1 component
        fEpsilon[2] = FourVector< complex<double> >(0,-1/sqrt(2.0),
                      complex<double>(0,-1/sqrt(2.0)),
                      0);
        // L=-1 component
        fEpsilon[0] = FourVector< complex<double> >(0,1/sqrt(2.0),
                      complex<double>(0,-1/sqrt(2.0)),
                      0);
    }
    // Linear polarization (first kind)
    else if (fObservPtr->photoprod && (fObservPtr->photo.pol.dpol
                                       && fObservPtr->photo.pol.doubpol.linbeam))
    {
        // L=+1 component
        fEpsilon[2] = FourVector< complex<double> >(0,1/sqrt(2.0),1/sqrt(2.0),0);
        // L=-1 component
        fEpsilon[0] = FourVector< complex<double> >(0,1/sqrt(2.0),-1/sqrt(2.0),0);
    }
    else if (fObservPtr->kaoncapture)
    {
        double sintheta_k = sqrt(1.0-fCosthkcm*fCosthkcm);
        /*
        // This part is for circularly polarized photons
        // L=+1 component
        fEpsilon[2] = FourVector< complex<double> >(0,-fCosthkcm/sqrt(2.0),
          				 complex<double>(0,-sqrt(0.5)),
          				 sintheta_k/sqrt(2.0));
        // L=-1 component
        fEpsilon[0] = FourVector< complex<double> >(0,fCosthkcm/sqrt(2.0),
          				 complex<double>(0,-sqrt(0.5)),
          				 -sintheta_k/sqrt(2.0));
        */
        // This part is for linearly polarized photons
        // L='x' component (actually in xz-plane)
        fEpsilon[2] = FourVector< complex<double> >(0,fCosthkcm,0,-sintheta_k);
        // L=y component
        fEpsilon[0] = FourVector< complex<double> >(0,0,1,0);
    }
    // Linear polarization (second kind)
    else
    {
        // L=+1 component
        fEpsilon[2] = FourVector< complex<double> >(0,1,0,0);
        // L=-1 component
        fEpsilon[0] = FourVector< complex<double> >(0,0,1,0);
    }
}

/*! \brief Calculates the spinors for the nucleon and the hyperon
 * given specific polarizations Lp and Ly.
 * When the spinor has previously been calculated, the process is not repeated.
 */
void TMatrixElement::determineSpinors(int Lp, int Ly)
{
    double theta, phi; /* these angles pin down the axis along which
		      * the baryons are polarized
		      * (to be determined later) */

    /* if the baryon spinor with polarization=Lp does not exist yet,
     * it has to be determined */
    if (!fNucleonSpinorIsCalculated[Lp])
    {
        fNucleonSpinorIsCalculated[Lp] = true; // Never calculate it again

        /* We need to determine the axis along which the nucleon is
         * polarized ( x, y or z ).
         * Most of the time we need the y-axis. First we select all
         * cases with a different axis (x or z) and then take
         * theta=90,phi=90 (=y) in all other cases.
         */
        theta = PI/2.0;
        phi = PI/2.0;

        /* in electroproduction with a polarized nucleon
         * the x- or z-axis are possible */
        if (fObservPtr->electroprod && fObservPtr->elec.baryon_pol
                && fObservPtr->elec.bar_pol.tarpol)
        {
            // if polarized along -x
            if (fObservPtr->elec.bar_pol.x_pol)
            {
                theta = PI/2.0;
                phi = 0;
            }
            // if polarized along -z
            else if (fObservPtr->elec.bar_pol.z_pol)
            {
                theta = 0.0;
                phi = 0.0;
            }
        }
        /* In double polarized photoproduction the nucleon is always
         * polarized along the (-x)- or (-z)-axis */
        else if (fObservPtr->photoprod && fObservPtr->photo.pol.dpol)
        {
            // if polarized along -x
            if (fObservPtr->photo.pol.doubpol.x_barvec)
            {
                theta = PI/2.0;
                phi = PI;
            }
            // if polarized along -z
            else if (fObservPtr->photo.pol.doubpol.z_barvec)
            {
                theta = PI;
                phi = 0.0;
            }
        }

        // Finally we calculate the spinor with the chosen quantization axis
        fNucleonSpinors[Lp] = CalculateSpinor(fp4vect,theta,phi,Lp);
    }

    /* if the hyperon spinor with polarization=Ly does not exist yet,
     * it has to be determined */
    if (!fHyperonSpinorIsCalculated[Ly])
    {
        fHyperonSpinorIsCalculated[Ly] = true; // Never calculate it again

        /* Choosing the quantization axis of the hyperon is completely analog to
         * the nucleon case (since simultanious polarization of both baryons is not
         * implemented at the present time). The code is copy-pasted
         * (with omission of comments).
         * One should note however that we choose the hyperons spin axis in the
         * primed reference frame, thus (theta,phi) change. In a first step we
         * take z' along the kaon momentum (Drechsel convention).
         */
        theta = PI/2.0;
        phi = PI/2.0;

        if (fObservPtr->electroprod && fObservPtr->elec.baryon_pol
                && fObservPtr->elec.bar_pol.recpol)
        {
            if (fObservPtr->elec.bar_pol.x_pol)
            {
                theta = PI/2.0 + acos(fCosthkcm);
                if (theta<=PI)
                    phi = 0.0;
                else
                {
                    theta *= -1.0;
                    phi = PI;
                }
            }
            else if (fObservPtr->elec.bar_pol.z_pol)
            {
                theta = acos(fCosthkcm);
                phi = 0.0;
            }
        }
        else if (fObservPtr->photoprod && fObservPtr->photo.pol.dpol)
        {
            if (fObservPtr->photo.pol.doubpol.x_barvec)
            {
                theta = PI/2.0 + acos(fCosthkcm);
                if (theta<=PI)
                    phi = 0.0;
                else
                {
                    theta *= -1.0;
                    phi = PI;
                }
            }
            else if (fObservPtr->photo.pol.doubpol.z_barvec)
            {
                theta = acos(fCosthkcm);
                phi = 0.0;
            }
        }


        /* The primed reference frame can be defined in two different ways.
         * Supra we have taken the Drechsel convention which chooses z'
         * along the kaon momentum.
         * The user is given the option to choose the z' axis along the
         * hyperon momentum (original convention of Stijn's thesis).
         * In this case there's a rotation of PI around the y=y' axis
         */
        if (fObservPtr->quant_axis_thesis)
        {
            theta = PI - theta;
            phi = PI - phi;
        }

        /* Finally we calculate the spinor with the chosen quantization axis.
         * We want the conjugate spinor = hermitian conjugate of the spinor
         * times gamma^0.
         */
        fHyperonSpinors[Ly] = (CalculateSpinor(fpY4vect,theta,phi,Ly)).H()*G0;
    }
}

/*! \brief  Parses a matrix with a fixed underflow (hard-coded to 1e-14) */
template<int K, int L>
Matrix<K,L> cleanMatrix(const Matrix<K,L>& x)
{
  const double underflow = 1.e-14;

  Matrix<K,L> result;
  for(int i=0; i<K; ++i)
    for(int j=0; j<L; ++j)
      result(i,j) = complex<double>( std::fabs(x(i,j).real()) < underflow ? 0. : x(i,j).real(),
				     std::fabs(x(i,j).imag()) < underflow ? 0. : x(i,j).imag());
  return result;
}

/*! \brief  Prints information about the object to the screen */
void TMatrixElement::Print()
{
    cout << "*******************" << endl
         << "Matrix Element info" << endl
         << "*******************" << endl
         << "Photon polarization vectors:" << endl
         << " L=-1:\t" << fEpsilon[0] << endl
         << " L= 0:\t" << fEpsilon[1] << endl
         << " L=+1:\t" << fEpsilon[2] << endl << endl
         << "Nucleon spinor:" << endl
         << " L=-1/2:\t"
         << fNucleonSpinors[0] << endl
         << " L=+1/2:\t"
         << fNucleonSpinors[1] << endl
         << endl
         << "Hyperon spinor:" << endl
         << " L=-1/2:\t"
         << fHyperonSpinors[0] << endl
         << " L=+1/2:\t"
         << fHyperonSpinors[1] << endl
         << endl
         << "Amputated current:" << endl
         << "[0]: " << cleanMatrix(GetCurrent()[0].value()) << endl
         << "[1]: " << cleanMatrix(GetCurrent()[1].value()) << endl
         << "[2]: " << cleanMatrix(GetCurrent()[2].value()) << endl
         << "[3]: " << cleanMatrix(GetCurrent()[3].value()) << endl
         << endl
         << "Current" << endl
         << "Lp=-1/2, Ly=-1/2: " << endl
         << "[0]: " << fHyperonSpinors[0]*GetCurrent()[0].value()*fNucleonSpinors[0] << endl
         << "[1]: " << fHyperonSpinors[0]*GetCurrent()[1].value()*fNucleonSpinors[0] << endl
         << "[2]: " << fHyperonSpinors[0]*GetCurrent()[2].value()*fNucleonSpinors[0] << endl
         << "[3]: " << fHyperonSpinors[0]*GetCurrent()[3].value()*fNucleonSpinors[0] << endl
         << endl
         << "Lp=-1/2, Ly=+1/2: " << endl
         << "[0]: " << fHyperonSpinors[1]*GetCurrent()[0].value()*fNucleonSpinors[0] << endl
         << "[1]: " << fHyperonSpinors[1]*GetCurrent()[1].value()*fNucleonSpinors[0] << endl
         << "[2]: " << fHyperonSpinors[1]*GetCurrent()[2].value()*fNucleonSpinors[0] << endl
         << "[3]: " << fHyperonSpinors[1]*GetCurrent()[3].value()*fNucleonSpinors[0] << endl
         << endl
         << "Lp=+1/2, Ly=-1/2: " << endl
         << "[0]: " << fHyperonSpinors[0]*GetCurrent()[0].value()*fNucleonSpinors[1] << endl
         << "[1]: " << fHyperonSpinors[0]*GetCurrent()[1].value()*fNucleonSpinors[1] << endl
         << "[2]: " << fHyperonSpinors[0]*GetCurrent()[2].value()*fNucleonSpinors[1] << endl
         << "[3]: " << fHyperonSpinors[0]*GetCurrent()[3].value()*fNucleonSpinors[1] << endl
         << endl
         << "Lp=+1/2, Ly=+1/2: " << endl
         << "[0]: " << fHyperonSpinors[1]*GetCurrent()[0].value()*fNucleonSpinors[1] << endl
         << "[1]: " << fHyperonSpinors[1]*GetCurrent()[1].value()*fNucleonSpinors[1] << endl
         << "[2]: " << fHyperonSpinors[1]*GetCurrent()[2].value()*fNucleonSpinors[1] << endl
         << "[3]: " << fHyperonSpinors[1]*GetCurrent()[3].value()*fNucleonSpinors[1] << endl
         << endl
         << "M^(Lp=-1/2 , Ly=-1/2) "
         << ( isGaugeInvariantForSpins(0,0) ? "is " : "is not " ) << "gauge invariant" << endl
         << "M^(Lp=+1/2 , Ly=-1/2) "
         << ( isGaugeInvariantForSpins(1,0) ? "is " : "is not " ) << "gauge invariant" << endl
         << "M^(Lp=-1/2 , Ly=+1/2) "
         << ( isGaugeInvariantForSpins(0,1) ? "is " : "is not " ) << "gauge invariant" << endl
         << "M^(Lp=+1/2 , Ly=+1/2) "
         << ( isGaugeInvariantForSpins(1,1) ? "is " : "is not " ) << "gauge invariant" << endl
         << endl;
}
