/*!
 * \file TCurrent.h
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

#ifndef TCURRENT_H
#define TCURRENT_H

#include "TCalculateCGLNMatrix.h"
class TCalculateCGLNCoeff;

class TCurrent
{

public:
    TCurrent(TCalculateCGLNCoeff* coeffCalculator, Class* particles, Observable* observ, const double w, const double k, const double costhkcm); // Constructor for TCurrent object
    virtual ~TCurrent() {}; // Destructor for TCurrent object
    TCurrent& operator+= ( const TCurrent& ); // Add two TCurrents by adding their respective coefficients.
    void Clear(); // Clear the coefficients vector and the fCurrentGS.
    void Initialise(TCalculateCGLNCoeff* coeffCalcPtr, Class* particlesPtr, Observable* observPtr, const double w, const double k, const double costhkcm);// Re-initialise TCurrent object resetting pointers, kinematics, etc. and determining the coefficients again
    void AddDiagram ( int classlabel, const Properties & particle); // Add a diagram of class classlabel with properties particle
    
    const GammaStructure& GetCurrentGS( const FourVector<std::complex<double> >& epsilon,
                            const FourVector<double>& k4vect,
                            const FourVector<double>& p4vect,
                            const FourVector<double>& pY4vect);// Get the GammaStructure containing the current
        
    void DetermineCurrentGS(const FourVector< std::complex< double > >& epsilon,
                            const FourVector< double >& k4vect,
                            const FourVector< double >& p4vect,
                            const FourVector< double >& pY4vect); // Determine the GammaStructure containing the current
    
    std::complex<double> GetCoefficient ( int index ) const; // Get the index'th coefficient (of M_{index})
    std::complex<double> DetermineSandwichedCurrent( const int cachelabel,
				    const int L, const int Lp, const int Ly,
				    const Matrix< 4, 1 >& nucleonSpinor, 
				    const Matrix< 1, 4 >& hyperonSpinor, 
				    const FourVector<std::complex<double> >& epsilon,
				    const FourVector<double>& k4vect,
				    const FourVector<double>& p4vect,
				    const FourVector<double>& pY4vect); // Determine <Ly|J_mu.epsilon_L^mu|Lp>
    TCalculateCGLNCoeff* GetCoeffCalculator() const; // Get the calculator ( i.e. the theoretical model) for the CGLN coefficients


private:
  
  //! \remark Pointers do NOT belong to TCurrent!
  //! CGLN coefficient calculator pointer
  TCalculateCGLNCoeff* fCoeffCalcPtr;
  
  // Particles and Observ pointers
  Class* fParticlesPtr;
  Observable* fObservPtr;
  
  // kinematics variables
  double fOmega;
  double fK;
  double fCosthkcm;
  
  //! the calculated CGLN coefficients
  std::vector<std::complex<double> > fCoefficients;

  GammaStructure fCurrentGS;
//  const static int fVarDF = 8; // s,t, 3 coupling constants, classindex, mass, width, TODO: extend with 3 ratios
  static const int fKinDF = 11;// s and t or u and t, and epsilon
  static const GammaStructure kZeroGS; //!< 4x4 matrix containing zeros
  void DetermineCoefficients();

};

#endif
