/*!
 * \file TMatrixElementOld.cpp
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

#include "TMatrixElementOld.h"
#include <iostream>
#include "Lagrangian.h"

using std::complex;
using std::cout;
using std::cerr;
using std::endl;

/*!
 * \class TMatrixElementOld
 *
 * Class to calculate a matrix element for given kinematics,
 * photon polarisation, nucleon and hyperon spin.
 * The current is calculated by explictly adding the 
 * individual diagrams (EMvertex*propagator*Svertex).
 * 
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 *
 */

// initialise statics
bool TMatrixElementOld::fgInstanceFlag = false;
TMatrixElementOld * TMatrixElementOld::fgMatrixElement = NULL;


/*! \brief The constructor, initialises kinematics, spinors, amputated current.
 * - the kinematics of the process are specified
 * - the 4 components of the amputated current are calculated.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 */
TMatrixElementOld::TMatrixElementOld(double w, double k, double costhkcm, double pk, Class* particlesPtr, Observable* observPtr)
  : TMatrixElement(w,k,costhkcm,pk,particlesPtr,observPtr)
  , current(FourVector<GammaStructure>(kZeroGS,kZeroGS,kZeroGS,kZeroGS))
{
    determineCurrent();
}

/*! \brief Destructor
*/
TMatrixElementOld::~TMatrixElementOld()
{
  fgInstanceFlag=false;
}

/*! \brief Singleton factory for TMatrixElementOld
 * Creates or get the unique TMatrixElementOld instance and initialises it.
 * Leaves the object in a state as if it were newly constructed with the same arguments.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 */
TMatrixElementOld* TMatrixElementOld::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr)
{
  if(! fgInstanceFlag)
  {
    fgMatrixElement = new TMatrixElementOld(w,k,costhkcm,pK,particlesPtr,observPtr);
    fgInstanceFlag= true;
    return fgMatrixElement;
  }
  else
  {
    fgMatrixElement->ReInitialiseWithMasses(w,k,costhkcm,pK,particlesPtr,observPtr,particlesPtr[0].partic[0].mass, particlesPtr[1].partic[0].mass, particlesPtr[2].partic[0].mass);
    return fgMatrixElement;
  }
}

/*! \brief Singleton factory for TMatrixElementOld, allowing off-shell kinematics
 * Creates or gets the unique TMatrixElementOld instance and initialises it, 
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
 */
TMatrixElementOld* TMatrixElementOld::GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY)
{
  if(! fgInstanceFlag)
  {
    fgMatrixElement = new TMatrixElementOld(w,k,costhkcm,pK,particlesPtr,observPtr);
    fgInstanceFlag= true;
    return fgMatrixElement;
  }
  else
  {
    fgMatrixElement->ReInitialiseWithMasses(w,k,costhkcm,pK,particlesPtr,observPtr,kin_mN, kin_mK, kin_mY);
    return fgMatrixElement;
  }
}

//--------------------------------------------------------------------------

/*! \brief (Re-)initialises the TMatrixElementOld as if it were newly constructed.
 * Essentially does the same as the constructor but with an existing object.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param particlesPtr pointer to the array of Class structs which contain the particles' parameters
 * \param observPtr pointer to the Observ struct
 */
void TMatrixElementOld::ReInitialise (double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label)
{
  TMatrixElement::ReInitialise(w,k,costhkcm,pK,particlesPtr,observPtr);
  for(int mu=0;mu<4;mu++)
    current[mu] = kZeroGS;
  
  determineCurrent();
}


/*! \brief Calculates the matrix element
 * This is the main function of the MatrixElement class
 * It calculates the numerical value of the matrix element
 * for fixed photon and baryon polarizations L,Lp and Ly
 */
complex<double> TMatrixElementOld::calculateM(int L, int Lp, int Ly)
{
    // if this matrixelement hasn't been determined before, proceed
    if (!fMEisCalculated[L][Lp][Ly])
    {
        // never calculate it again
        fMEisCalculated[L][Lp][Ly] = true;

        // Calculate the baryon spinors when necessary
        determineSpinors(Lp,Ly);

        /* We are now ready to calculate the matrix element
         * = epsilon[L] * ( hyperonSpinor[Ly] * current * nucleonSpinor[Lp] )
         * We'll not use the 4vector product here, because that way all components
         * of the right-side 4vector would be evaluated. Since many components of
         * epsilon are zero, there would be a large amount of redundant calculations
         */

        fMatrixElement[L][Lp][Ly] = 0.0; // set initial value of the matrix element to zero

        fMatrixElement[L][Lp][Ly] += ( (fEpsilon[L])[0]==complex<double>(0,0) ? 0 :
                                       (fEpsilon[L])[0] * ( fHyperonSpinors[Ly] *
                                                            current[0].value() *
                                                            fNucleonSpinors[Lp] ) );

        for (int mu=1; mu<=3; mu++)
        {
            fMatrixElement[L][Lp][Ly] -= ( (fEpsilon[L])[mu]==complex<double>(0,0) ? 0 :
                                           (fEpsilon[L])[mu] * ( fHyperonSpinors[Ly] *
                                                                 current[mu].value() *
                                                                 fNucleonSpinors[Lp] ) );
        }
    }

    return fMatrixElement[L][Lp][Ly];

}


//--------------------------------------------------------------------------


/*! \brief Calculates the components of the amputated current 
 * This function calculates all components of the amputated current
 * We add all contributions from all possible diagrams.
 * In case of Regge trajectory exchange we multiply with the
 * Regge propagator. In t-channel Reggeization a gauge-invariance
 * restoration term is added for certain particles
 */
void TMatrixElementOld::determineCurrent()
{
    /* tempCurrent[] will be used to add different contributions to the
     * amputated current.
     * It is of type FourVector<GammaStructure> to make use of the
     * matrix multiplication optimalization of the GammaStructure class.
     * It is only at the very end, that the GammaStructure-components
     * are convertex to complex 4x4 matrices. See isGaugeInvariant(int,int)
     * for example. */

    // diagramContribution temporarely holds the contribution of one diagram
    FourVector<GammaStructure> diagramContribution;


    // Loop over all types of diagrams (see Lagrangian.h for a list)
    for (int diagram=0; diagram<CLASSMAX; diagram++)
    {
        // Loop over all exchanged particles for a specific diagram
        for (int particle=0; particle<fParticlesPtr[diagram].particount; particle++)
        {
            if (!fObservPtr->kaoncapture)
                diagramContribution =
                    calculateDiagram[diagram](fParticlesPtr[diagram].partic[particle],
                                              fk4vect,fp4vect,fpK4vect,fpY4vect,fObservPtr->kaoncapture);
            else /* Use crossing symmetry k -> -k; pk -> -pk */
                diagramContribution =
                    calculateDiagram[diagram](fParticlesPtr[diagram].partic[particle],
                                              -fk4vect,fp4vect,-fpK4vect,fpY4vect,fObservPtr->kaoncapture);

            /* In case of Regge trajectory exchange
             * which only occurs in diagrams S,T,U,A,B,C */
            if ( fObservPtr->regge && diagram<=5 )
            {
                /* If we're Reggeizing a t-channel diagram
                 * the electric part of the s-channel or u-channel diagram
                 * needs to be added to restore gauge invariance
                 *
                 * When the charge of the exchanged t-channel particle
                 * is zero, this gauge invariance restauration procedure
                 * can obviously be ommited.
                 *
                 * This procedure remains unaffected in the case of radiative
                 * kaon capture.
                 */
                if ( diagram==1 && fParticlesPtr[diagram].partic[particle].E!=0.0 )
                {
                    // Add electric part of s-channel
                    // when there's a proton in the initial state
                    if ( fParticlesPtr[0].partic[0].E!=0.0 ) {

                        /* We copy the s-channel exchange particle */
                        Properties sChannelParticle = fParticlesPtr[0].partic[0];

                        /* Set the magnetic explicitely to zero
                         * and copy the strong coupling constant g_KYN
                         * from the t-channel particle */
                        sChannelParticle.G = 0;
                        sChannelParticle.H = fParticlesPtr[diagram].partic[particle].H;

                        /* And use the same formfactor as the original particle  */
                        sChannelParticle.formfactorE =
                            fParticlesPtr[diagram].partic[particle].formfactorE;

                        /* Now we can add this diagram to the amputated current */
                        if (!fObservPtr->kaoncapture)
                            diagramContribution +=
                                calculateDiagram[0](sChannelParticle,
                                                    fk4vect,fp4vect,fpK4vect,fpY4vect,fObservPtr->kaoncapture);
                        else /* Use crossing symmetry k -> -k; pk -> -pk */
                            diagramContribution +=
                                calculateDiagram[0](sChannelParticle,
                                                    -fk4vect,fp4vect,-fpK4vect,fpY4vect,fObservPtr->kaoncapture);
                    } // end add s-channel electric

                    // Add electric part of u-channel
                    // when there's a charged hyperon in the final state
                    else if ( fParticlesPtr[2].partic[0].E!=0.0 ) {

                        /* We copy the u-channel exchange particle */
                        Properties uChannelParticle = fParticlesPtr[2].partic[0];

                        /* Set the magnetic explicitely to zero
                         * and copy the strong coupling constant g_KYN
                         * from the t-channel particle */
                        uChannelParticle.G = 0;
                        uChannelParticle.H = fParticlesPtr[diagram].partic[particle].H;

                        /* And use the same formfactor as the original particle  */
                        uChannelParticle.formfactorE =
                            fParticlesPtr[diagram].partic[particle].formfactorE;

                        /* Now we can add this diagram to the amputated current */
                        if (!fObservPtr->kaoncapture)
                            diagramContribution +=
                                calculateDiagram[2](uChannelParticle,
                                                    fk4vect,fp4vect,fpK4vect,fpY4vect,fObservPtr->kaoncapture);
                        else /* Use crossing symmetry k -> -k; pk -> -pk */
                            diagramContribution +=
                                calculateDiagram[2](uChannelParticle,
                                                    -fk4vect,fp4vect,-fpK4vect,fpY4vect,fObservPtr->kaoncapture);
                    } // end add u-channel electric

                    else {
                        cerr << "ERROR in TMatrixElementOld::determineCurrent(): "
                             << "Problem restoring gauge invariance for "
                             << fParticlesPtr[diagram].partic[particle].nickname
                             << " exchange.\n";
                        exit(1);
                    }

                } // end restore gauge-invariance

                // We multiply with the correct Regge propagator
		double sign = (!fObservPtr->kaoncapture ? 1. : -1.);
		double 
		  s = (fp4vect + sign*fk4vect)*(fp4vect + sign*fk4vect),
		  t = (fk4vect - fpK4vect)*(fk4vect - fpK4vect),
		  u = (fk4vect - sign*fpY4vect)*(fk4vect - sign*fpY4vect),
		  Q2 = -1.*(fk4vect*fk4vect); 

                if (!fObservPtr->kaoncapture)
                    diagramContribution = 
		      diagramContribution
		      *propagatorRegge(fParticlesPtr[diagram].partic[particle],
				       s, u, t, Q2, 0., fObservPtr);
                else /* Use crossing symmetry s <-> u */
                    diagramContribution = 
		      diagramContribution
		      *propagatorRegge(fParticlesPtr[diagram].partic[particle],
				       u, s, t, Q2, 0., fObservPtr);
            } // end of Regge case

            // Add the contribution of this diagram to the current
            current += diagramContribution;

        } // particle loop
    } // diagram loop
}
