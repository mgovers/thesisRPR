/*!
 * \file TMatrixElement.h
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

#ifndef TMATRIXELEMENT_H
#define TMATRIXELEMENT_H

#include <complex>
#include <vector>
#include "FourVector.h"
#include "GammaStructure.h"
#include "Matrix.h"
#include "Structures.h"

class TMatrixElement
{
public:
    static TMatrixElement* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< MatrixElement singleton factory
    static TMatrixElement* GetMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label=-1);//!< MatrixElement singleton factory,  allowing off-shell kinematics
    virtual ~TMatrixElement();

    virtual void ReInitialise( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, int label=-1); //!< Initialises the object as if it were newly constructed.
    virtual void ReInitialiseWithMasses( double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY, int label=-1);//!< Initialises the object as if it were newly constructed, allowing off-shell kinematics.
    virtual const FourVector<GammaStructure>& GetCurrent() =0; //!< returns the hadronic current J^\mu
    virtual std::complex<double> calculateM(int L, int Lp, int Ly)=0;//!< Calculates the complex value of the matrix element
    virtual bool isGaugeInvariant(); //!< Tests whether the amplitude is gauge invariant
    virtual bool isGaugeInvariantForSpins (int Lp, int Ly); //!< Tests whether the amplitude is gauge invariant
    virtual void updateSpinDependencies();
    virtual void SetDataLabel(int label){};
    virtual void CalculateKinematics();
    virtual void Print();
    


protected:
    TMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr); //!< protected constructor
    TMatrixElement(double w, double k, double costhkcm, double pK, Class* particlesPtr, Observable* observPtr, double kin_mN, double kin_mK, double kin_mY); //!< protected constructor,  allowing off-shell kinematics
    TMatrixElement ( const TMatrixElement& other);    
    TMatrixElement& operator=(const TMatrixElement& right); 
    
    // Masses of the incoming and outgoing mesons and baryons
    double mN; //!< nucleon mass
    double mK; //!< kaon mass
    double mY; //!< hyperon mass

    /* // The Mandelstam variables */
    /* double fS; */
    /* double fT; */
    /* double fU; */

    // Kinematic input
    double fOmega;
    double fK;
    double fCosthkcm;
    double fPK;

    /* Energy-momentum 4vector of the incoming and outgoing particles
    * in the CM frame */
    FourVector<double> fk4vect;  //!< photon energy-momentum 4vector in the CM frame
    FourVector<double> fp4vect;  //!< nucleon energy-momentum 4vector in the CM frame
    FourVector<double> fpK4vect; //!< kaon energy-momentum 4vector in the CM frame
    FourVector<double> fpY4vect; //!< hyperon energy-momentum 4vector in the CM frame

    // Particles and Observ pointers
    Class* fParticlesPtr;
    Observable* fObservPtr;

    // bools & vector<bools> to check whether something has been calculated
    std::vector< std::vector< std::vector< bool> > > fMEisCalculated;
    std::vector< bool> fNucleonSpinorIsCalculated;
    std::vector< bool> fHyperonSpinorIsCalculated;

    // Spinors and fourvectors
    std::vector<FourVector<std::complex<double> > > fEpsilon;
    static const Matrix<4,1> kZero41Matrix;
    static const Matrix<1,4> kZero14Matrix;
    static const GammaStructure kZeroGS;

    /*! 2-component array containing for each Lp-value
     * the spinor of the incoming nucleon
     *\verbatim
     index -> Lp: 0->-1/2
                  1->+1/2 
     \endverbatim
    */
    std::vector< Matrix<4,1> > fNucleonSpinors; 

    /*! 2-component array containing for each Ly-value
     * the conjugate spinor of the outgoing hyperon
     *\verbatim
     index -> Ly: 0->-1/2
                  1->+1/2 
     \endverbatim
    */
    std::vector< Matrix<1,4> > fHyperonSpinors;

    std::vector<std::vector< std::vector <std::complex< double > > > > fMatrixElement;  //!< the actual matrix element

    /*! Assigns a correct polarization vector to epsilon for each photon polarization */
    void determineEpsilon(); 

    /*! Calculates nucleon and hyperon spinors */
    void determineSpinors(int Lp, int Ly);
    

 private:
    /*! Calculates the spinor wavefunction of a particle
     * given its energy-momentum 4vector, angles of its
     * quantization axis in the CMF and polarization.*/
    Matrix<4,1> CalculateSpinor(const FourVector<double>&,double,double,int);



};

#endif // TMATRIXELEMENT_H
