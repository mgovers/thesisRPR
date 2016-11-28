/*!
 * \file TCalculateConsistentCoeff.cpp
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

/*!
 * \class TCalculateConsistentCoeff
 *
 * Calculates the coefficients of the reaction-amplitude's 
 * decomposition in the CGLN basis.
 * This is an implementation for the version of the RPR model 
 * which features consistent, gauge-invariant couplings for
 * spin 3/2 and 5/2 particles.
 * 
 * \author Tom Vrancx <Tom.Vrancx@UGent.be>
 * \author Lesley De Cruz <Lesley.DeCruz@UGent.be>
 */ 

#include "TCalculateConsistentCoeff.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
#include <complex>
using std::complex;
#include <cmath>
#include <cstring>
#include "Lagrangian.h"
#include "FormFactorParametrization.h"

bool TCalculateConsistentCoeff::fgInstanceFlag = false;
TCalculateConsistentCoeff* TCalculateConsistentCoeff::fgCalculator = NULL;

// local coupling constants
double e_; //!< electrical charge divided by |e_electron|
double kappa_; //!< magnetic moment for diagrams S, T, U, A, B, C, D, E, F and G // first EM coupling constant for diagrams H, I, J, L, M, N, O and Q
double g_; //!< strong coupling constant divided by sqrt(4PI) for diagrams S, T, U, A, B, C, D, E, F and G // second EM coupling constant for diagrams H, I, J, L, M, N, O and Q
double f_; //!< strong coupling constant for diagrams H, I, J, L, M, N, O and Q
double gv_; //!< strong vector coupling constant
double gt_; //!< strong tensor coupling constant
complex<double> ReggePropagator; //!< Regge propagator
complex<double> coefficient; //!< variable used by diagram T to restore gauge invariance
complex<double> Fpi; //!< variable used by diagram R
complex<double> Fs; //!< variable used by diagram R
double mandel; //!< variable used by diagram R


/*! \brief Singleton factory for the consistent model, features no spin-dependence in the cutoff.
 * The exponent of the correction factor for the strong cutoff is set to 0.
 */
TCalculateConsistentCoeff* TCalculateConsistentCoeff::GetConsistentCalculator()
{
  if(! fgInstanceFlag)
    {
      fgCalculator= new TCalculateConsistentCoeff(0.0); // FIXME -- this ends up as reachable  (is never deleted)
      fgInstanceFlag=true;
    }
    else
    {
      fgCalculator->fExponent=0.0;
      for(int c=0; c < CLASSMAX; c++)
	fgCalculator->fWidthModifier[c]= 1.0;
    }
    
    return fgCalculator;
}

/*! \brief Singleton factory for model with variable cutoff, exponent 1.
 * Exponent n= 1 for the correction factor (J'/J)^n for the strong cutoff.
 */
TCalculateConsistentCoeff* TCalculateConsistentCoeff::GetVarcutoff1Calculator()
{
  if(! fgInstanceFlag)
  {
    fgCalculator= new TCalculateConsistentCoeff(1.0); // FIXME -- this ends up as reachable  (is never deleted)
    fgInstanceFlag=true;
  }
  else
  {
    fgCalculator->fExponent=1.0;
    for(int c=0; c < CLASSMAX; c++)
      fgCalculator->fWidthModifier[c]= 1.0;
  }

  return fgCalculator;
}

/*! \brief Singleton factory for model with variable cutoff, exponent 1/2.
 * Exponent n= 1/2 for the correction factor (J'/J)^n for the strong cutoff.
 */
TCalculateConsistentCoeff* TCalculateConsistentCoeff::GetVarcutoff2Calculator()
{
  if(! fgInstanceFlag)
  {
    fgCalculator= new TCalculateConsistentCoeff(0.5); // FIXME -- this ends up as reachable  (is never deleted)
    fgInstanceFlag=true;
  }
  else
  {
    fgCalculator->fExponent=0.5;
    for(int c=0; c < CLASSMAX; c++)
      fgCalculator->fWidthModifier[c]= 1.0;
  }
  return fgCalculator;
}

/*! \brief Singleton factory for model with Lorentz-type dependence of the width.
 * No correction factor for the cutoff, but for the width.
 * Assuming that strong_formfactor is set to Lorentz! (2)
 */
TCalculateConsistentCoeff* TCalculateConsistentCoeff::GetLorentzCalculator()
{  
  if(! fgInstanceFlag)
  {
    fgCalculator= new TCalculateConsistentCoeff(0.0); // FIXME -- this ends up as reachable  (is never deleted)
    fgInstanceFlag=true;
  }
  else
  {
    fgCalculator->fExponent= 0.0;
  }
  for(int c=10; c < 14; c++)
    fgCalculator->fWidthModifier[c]= 1.0/sqrt(pow(2.,1./3.)-1.);
  for(int c=14; c < 18; c++)
    fgCalculator->fWidthModifier[c]= 1.0/sqrt(pow(2.,1./5.)-1.);
  return fgCalculator;
}


/*! \brief Constructor for model with variable cutoff.
 * \param exponent Exponent n for the correction factor (J'/J)^n for the strong cutoff.
 */
TCalculateConsistentCoeff::TCalculateConsistentCoeff(double exponent)
: TCalculateCGLNCoeff(), fExponent(exponent)
 {
   fgInstanceFlag=true;
 }

/*! \brief Destructor
 */
TCalculateConsistentCoeff::~TCalculateConsistentCoeff()
{
  fgInstanceFlag=false;
}


complex< double >
TCalculateConsistentCoeff::CalcA1 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*g_*(2.0*e_*M_P+kappa_*(fMass+fmN)))/(2.0*M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1: {
      double nucleoncharge_ = fNucleon_charge;
      double hyperoncharge_ = fHyperon_charge;
      if (particle.formfactorE!=NULL)
      {
	nucleoncharge_ = (*particle.formfactorE).value(fNucleon_charge,-1.0*fkk);
	hyperoncharge_ = (*particle.formfactorE).value(fHyperon_charge,-1.0*fkk);
      }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fT);
      coefficient = 0.0;
      if(observ.regge)
      	{
      	  ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      	  if (particle.E!=0.0)
      	    {
      	      // Add electric part of s-channel
      	      // when there's a proton in the initial state
      	      if (fNucleon_charge!=0.0)
		coefficient += (ELEC*g_*nucleoncharge_)/(fS-fmN*fmN);

      	      // Add electric part of u-channel
      	      // when there's a charged hyperon in the final state
      	      else if (fHyperon_charge!=0.0)
		coefficient += (ELEC*g_*hyperoncharge_)/(fU-fmY*fmY);

      	      else
      		{
      		  cerr << "ERROR in TCalculateConsistentCoeff::CalcA1(...): "
      		       << "Problem restoring gauge invariance for "
      		       << particle.nickname
      		       << " exchange.\n";
      		  exit(1);
      		}
      	    }
        }
      else
	ReggePropagator = 1;
      return ReggePropagator*coefficient;
    }

      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*g_*(2.0*e_*M_P+kappa_*(fMass+fmY)))/(2.0*M_P*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram A
    case 3:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE!=NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*g_*kappa_*(fMass+fmY))/(2.0*M_P*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram B
    case 4:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(-1.0)*ELEC*kappa_*(gt_*fT+gv_*(M_L+M_P)*(fmY+fmN)) / (1000.0*(M_L+M_P)*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      return (ELEC*g_*kappa_*(fMass+fmN))/(2.0*M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      return (ELEC*g_*kappa_*(fMass-fmN))/(2.0*M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      return -1.0*((ELEC*kappa_*g_*(fMass+fmY))/(2.0*M_P*(-1.0*fDenominator_u)));

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      return (ELEC*kappa_*g_*(fmY-fMass))/(2.0*M_P*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(-1.0*fkp*(fmY*(-6.0*kappa_*fMass*M_P-8.0*kappa_*fmN*M_P+g_*(fmN*fmN))+fppY*(g_*fMass-12.0*kappa_*M_P)+g_*fMass*fkpY)+fkk*(2.0*kappa_*M_P*(fkpY+fMass*fmY+2.0*fmN*fmY+4.0*fppY)-g_*fmY*fkp)+2.0*kappa_*fmN*M_P*((2.0*fMass-fmN)*fkpY+2.0*(fMass+fmN)*(fmN*fmY+fppY))-2.0*g_*fmY*(fkp*fkp));

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(fkp*(fmY*(6.0*kappa_*fMass*M_P-8.0*kappa_*fmN*M_P+g_*(fmN*fmN))-fppY*(12.0*kappa_*M_P+g_*fMass)-g_*fMass*fkpY)+fkk*(2.0*kappa_*M_P*(-1.0*fkpY+fMass*fmY-2.0*fmN*fmY-4.0*fppY)+g_*fmY*fkp)+2.0*kappa_*fmN*M_P*((2.0*fMass+fmN)*fkpY+2.0*(fMass-fmN)*(fmN*fmY+fppY))+2.0*g_*fmY*(fkp*fkp));

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(fmN*fkpY*(-6.0*kappa_*fMass*M_P-8.0*kappa_*fmY*M_P-2.0*g_*fkpY+g_*(fmY*fmY))+fkk*(-2.0*kappa_*M_P*fkp+2.0*kappa_*M_P*(fMass*fmN+2.0*fmN*fmY+4.0*fppY)+g_*fmN*fkpY)-fkp*(2.0*kappa_*fmY*M_P*(2.0*fMass-fmY)+g_*fMass*fkpY)+fppY*(fkpY*(g_*fMass-12.0*kappa_*M_P)+4.0*kappa_*fmY*M_P*(fMass+fmY))+4.0*kappa_*fmN*(fmY*fmY)*M_P*(fMass+fmY));

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(fmN*fkpY*(-6.0*kappa_*fMass*M_P+8.0*kappa_*fmY*M_P+2.0*g_*fkpY-g_*(fmY*fmY))+fkk*(2.0*kappa_*M_P*(fkp-4.0*fppY)+2.0*kappa_*fmN*M_P*(fMass-2.0*fmY)-g_*fmN*fkpY)-fkp*(2.0*kappa_*fmY*M_P*(2.0*fMass+fmY)+g_*fMass*fkpY)+fppY*(fkpY*(12.0*kappa_*M_P+g_*fMass)+4.0*kappa_*fmY*M_P*(fMass-fmY))+4.0*kappa_*fmN*(fmY*fmY)*M_P*(fMass-fmY));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*((fkp*fkp)*(-1.0*fppY*(fmY*(8.0*kappa_*fMass*M_P-12.0*kappa_*fmN*M_P+g_*(fmN*fmN))+24.0*kappa_*M_P*fkpY)+fkpY*(fmY*(3.0*fmN*(g_*fmN-4.0*kappa_*M_P)+8.0*kappa_*fMass*M_P)-g_*fMass*fkpY)+(fppY*fppY)*(16.0*kappa_*M_P+g_*fMass)+2.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P)+(fkk*fkk)*(fppY*(2.0*kappa_*M_P*(fkpY-fMass*fmY+3.0*fmN*fmY+6.0*fppY)-g_*fmY*fkp)-2.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P)+fkk*(fkp*(fppY*(fkpY*(g_*fMass-12.0*kappa_*M_P)+2.0*kappa_*fmY*M_P*(9.0*fmN-5.0*fMass)-g_*(fmN*fmN)*fmY)+fkpY*(fmY*(2.0*kappa_*M_P*(fMass-3.0*fmN)+g_*(fmN*fmN))-2.0*kappa_*M_P*fkpY)+(fppY*fppY)*(30.0*kappa_*M_P+g_*fMass)-4.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P)+fmY*(fkp*fkp)*(2.0*kappa_*fmY*M_P+g_*fkpY-3.0*g_*fppY)+2.0*kappa_*fmN*M_P*(fkpY*(fmN*fmY*(fMass-3.0*fmN)-3.0*(fMass+2.0*fmN)*fppY)-(fmN*fmY+fppY)*((3.0*fMass-4.0*fmN)*fppY+(fmN*fmN)*fmY)))+fmN*fkp*(fmN*fkpY*(fmY*(2.0*kappa_*M_P*(7.0*fMass-9.0*fmN)+g_*(fmN*fmN))-fppY*(28.0*kappa_*M_P+g_*fMass))+(fkpY*fkpY)*(6.0*kappa_*M_P*(fMass+fmN)-g_*fMass*fmN)-6.0*kappa_*M_P*(fMass-fmN)*fppY*(fmN*fmY+fppY))+2.0*fmY*pow(fkp,3.0)*(2.0*kappa_*fmY*M_P+g_*fkpY-g_*fppY)+2.0*kappa_*pow(fmN,3.0)*M_P*fkpY*((3.0*fMass+2.0*fmN)*fkpY+3.0*(fMass-fmN)*(fmN*fmY+fppY)));

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*((fkp*fkp)*(fkpY*(fmY*(8.0*kappa_*fMass*M_P+12.0*kappa_*fmN*M_P-3.0*g_*(fmN*fmN))+24.0*kappa_*M_P*fppY)+fppY*(fppY*(g_*fMass-16.0*kappa_*M_P)-4.0*kappa_*fmY*M_P*(2.0*fMass+3.0*fmN)+g_*(fmN*fmN)*fmY)-2.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P-g_*fMass*(fkpY*fkpY))+(fkk*fkk)*(2.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P-fppY*(2.0*kappa_*M_P*(fkpY+fMass*fmY+3.0*fmN*fmY+6.0*fppY)-g_*fmY*fkp))+fkk*(fkp*(fkpY*(fppY*(12.0*kappa_*M_P+g_*fMass)+2.0*kappa_*fmY*M_P*(fMass+3.0*fmN)-g_*(fmN*fmN)*fmY)+fppY*(fppY*(g_*fMass-30.0*kappa_*M_P)-2.0*kappa_*fmY*M_P*(5.0*fMass+9.0*fmN)+g_*(fmN*fmN)*fmY)+2.0*kappa_*M_P*(fkpY*fkpY)+4.0*kappa_*(fmN*fmN)*(fmY*fmY)*M_P)-fmY*(fkp*fkp)*(2.0*kappa_*fmY*M_P+g_*fkpY-3.0*g_*fppY)+2.0*kappa_*fmN*M_P*(fkpY*(fmN*fmY*(fMass+3.0*fmN)-3.0*(fMass-2.0*fmN)*fppY)+(fmN*fmY+fppY)*((fmN*fmN)*fmY-(3.0*fMass+4.0*fmN)*fppY)))-fmN*fkp*(fmN*fkpY*(fppY*(g_*fMass-28.0*kappa_*M_P)-2.0*kappa_*fmY*M_P*(7.0*fMass+9.0*fmN)+g_*(fmN*fmN)*fmY)+(fkpY*fkpY)*(6.0*kappa_*M_P*(fmN-fMass)+g_*fMass*fmN)+6.0*kappa_*M_P*(fMass+fmN)*fppY*(fmN*fmY+fppY))-2.0*fmY*pow(fkp,3.0)*(2.0*kappa_*fmY*M_P+g_*fkpY-g_*fppY)+2.0*kappa_*pow(fmN,3.0)*M_P*fkpY*((3.0*fMass-2.0*fmN)*fkpY+3.0*(fMass+fmN)*(fmN*fmY+fppY)));

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fkk,2.0)*(fppY*(2.0*kappa_*M_P*(fkp+fmN*(fMass-3.0*fmY))-g_*fmN*fkpY)+2.0*kappa_*pow(fmN,2.0)*pow(fmY,2.0)*M_P-12.0*kappa_*M_P*pow(fppY,2.0))+fkk*(fkpY*(fppY*(fppY*(30.0*kappa_*M_P+g_*fMass)+2.0*kappa_*fmN*M_P*(9.0*fmY-5.0*fMass)-g_*fmN*pow(fmY,2.0))-4.0*kappa_*pow(fmN,2.0)*pow(fmY,2.0)*M_P)+fkp*(-fmN*(pow(fmY,2.0)-fkpY)*(g_*fkpY-2.0*kappa_*M_P*(fMass-3.0*fmY))-fppY*(fkpY*(g_*fMass-12.0*kappa_*M_P)+6.0*kappa_*fmY*M_P*(fMass+2.0*fmY)))+fmN*pow(fkpY,2.0)*(3.0*g_*fppY-2.0*kappa_*fmN*M_P)-2.0*kappa_*M_P*pow(fkp,2.0)*fkpY+2.0*kappa_*fmY*M_P*(fmN*fmY+fppY)*((3.0*fMass-4.0*fmY)*fppY+fmN*pow(fmY,2.0)))+fkpY*(fkpY*(fmN*fppY*(8.0*kappa_*fMass*M_P-12.0*kappa_*fmY*M_P+g_*pow(fmY,2.0))-pow(fppY,2.0)*(16.0*kappa_*M_P+g_*fMass)-2.0*kappa_*pow(fmN,2.0)*pow(fmY,2.0)*M_P)+2.0*fmN*pow(fkpY,2.0)*(2.0*kappa_*fmN*M_P-g_*fppY)-6.0*kappa_*fmY*M_P*(fMass-fmY)*fppY*(fmN*fmY+fppY))+fkp*(fmN*(pow(fmY,2.0)-fkpY)*(fkpY*(4.0*kappa_*M_P*(3.0*fmY-2.0*fMass)+2.0*g_*fkpY-g_*pow(fmY,2.0))+6.0*kappa_*pow(fmY,2.0)*M_P*(fMass-fmY))+fppY*(pow(fmY,2.0)*fkpY*(28.0*kappa_*M_P+g_*fMass)-24.0*kappa_*M_P*pow(fkpY,2.0)+6.0*kappa_*pow(fmY,3.0)*M_P*(fMass-fmY)))+pow(fkp,2.0)*(fmY*fkpY*(6.0*kappa_*M_P*(fMass+fmY)-g_*fMass*fmY)-2.0*kappa_*pow(fmY,3.0)*M_P*(3.0*fMass+2.0*fmY)+g_*fMass*pow(fkpY,2.0)));

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fkk,2.0)*(2.0*kappa_*pow(fmN,2.0)*pow(fmY,2.0)*M_P-fppY*(-2.0*kappa_*M_P*fkp+2.0*kappa_*M_P*(fMass*fmN+3.0*fmN*fmY+6.0*fppY)+g_*fmN*fkpY))+fkpY*(fkpY*(fppY*(fppY*(g_*fMass-16.0*kappa_*M_P)-4.0*kappa_*fmN*M_P*(2.0*fMass+3.0*fmY)+g_*fmN*pow(fmY,2.0))-2.0*kappa_*pow(fmN,2.0)*pow(fmY,2.0)*M_P)+2.0*fmN*pow(fkpY,2.0)*(2.0*kappa_*fmN*M_P-g_*fppY)+6.0*kappa_*fmY*M_P*(fMass+fmY)*fppY*(fmN*fmY+fppY))+fkk*(fmN*fppY*(fkpY*(2.0*kappa_*M_P*(5.0*fMass+9.0*fmY)+3.0*g_*fkpY-g_*pow(fmY,2.0))-6.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY))+fkp*(fppY*(fkpY*(12.0*kappa_*M_P+g_*fMass)+6.0*kappa_*fmY*M_P*(fMass-2.0*fmY))-fmN*(pow(fmY,2.0)-fkpY)*(2.0*kappa_*M_P*(fMass+3.0*fmY)+g_*fkpY))-pow(fppY,2.0)*(fkpY*(g_*fMass-30.0*kappa_*M_P)+2.0*kappa_*fmY*M_P*(3.0*fMass+4.0*fmY))+2.0*kappa_*pow(fmN,2.0)*M_P*(pow(fmY,4.0)-fkpY*(fkpY+2.0*pow(fmY,2.0)))-2.0*kappa_*M_P*pow(fkp,2.0)*fkpY)+fkp*(-fmN*(pow(fmY,2.0)-fkpY)*(fkpY*(-4.0*kappa_*M_P*(2.0*fMass+3.0*fmY)-2.0*g_*fkpY+g_*pow(fmY,2.0))+6.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY))-fppY*(pow(fmY,2.0)*fkpY*(g_*fMass-28.0*kappa_*M_P)+24.0*kappa_*M_P*pow(fkpY,2.0)+6.0*kappa_*pow(fmY,3.0)*M_P*(fMass+fmY)))+pow(fkp,2.0)*(fmY*fkpY*(6.0*kappa_*M_P*(fmY-fMass)+g_*fMass*fmY)+2.0*kappa_*pow(fmY,3.0)*M_P*(3.0*fMass-2.0*fmY)-g_*fMass*pow(fkpY,2.0)));

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      if(particle.formfactorE == NULL)
	Fpi = particle.E;
      else
	Fpi = (*particle.formfactorE).value(particle.E, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);

      Fs = protonEMTFF(-fkk, particle.E > 0 ? fS : fU, particle.Y)*
	propagatorRegge(particle, fS, fU, fT, -fkk, particle.Z, &observ);

      f_ = particle.I;

      if(particle.E < 0)
	return -1.0*((sqrt(2.0)*ELEC*f_*Fs)/(2.0*fkpY-fkk+pow(fmN,2.0)-pow(fmY,2.0)));

      return (ELEC*sqrt(2.0)*Fs*f_)/(2.0*fkp+fkk);
	  

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return sqrt(2.0)*ELEC*f_*Fpi*(fmN*(kappa_*fmN+fmN+fmY)-kappa_*fppY)/(fmN*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      return 0.0;

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA1(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


complex< double >
TCalculateConsistentCoeff::CalcA2 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(e_*ELEC*g_)/(fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1: {
      double nucleoncharge_ = fNucleon_charge;
      double hyperoncharge_ = fHyperon_charge;
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
      {
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
	nucleoncharge_ = (*particle.formfactorE).value(fNucleon_charge,-1.0*fkk);
	hyperoncharge_ = (*particle.formfactorE).value(fHyperon_charge,-1.0*fkk);
      }
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fT);
      coefficient = (e_*ELEC*g_)/(fkpY*fDenominator_t);
      if(observ.regge)
	{
	  ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
	  if (particle.E!=0.0)
	    {
	      // Add electric part of s-channel
	      // when there's a proton in the initial state
	      if (fNucleon_charge!=0.0)
		coefficient += (nucleoncharge_*ELEC*g_)/(fkpY*(fS-fmN*fmN));

	      else
		{
		  cerr << "ERROR in TCalculateConsistentCoeff::CalcA1(...): "
		       << "Problem restoring gauge invariance for "
		       << particle.nickname
		       << " exchange.\n";
		  exit(1);
		}
	    }
	}
      else
	ReggePropagator = 1;
      return ReggePropagator*coefficient;
    }
      
      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      return 0.0;

      //Diagram A
    case 3:
      return 0.0;

      //Diagram B
    case 4:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL || !observ.electroprod)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(-1.0)*ELEC*kappa_*gt_ / (1000.0*(M_L+M_P)*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*((complex<double>(0.0,1.0)*ELEC*kappa_*(fkp-fkpY)*(gt_*(fMass*fMass)-gv_*(M_L+M_P)*(fmN+fmY)))/(2000.0*(fMass*fMass)*(M_L+M_P)*fkpY*fDenominator_t));

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*(fkk*(fkpY*(8.0*kappa_*M_P+g_*fmN)+fMass*fmY*(2.0*kappa_*M_P+g_*fmN)+fppY*(-4.0*kappa_*M_P+3.0*g_*fMass-2.0*g_*fmN))-3.0*fkpY*(2.0*fkp+(fmN*fmN))*(g_*(fMass-fmN)-4.0*kappa_*M_P)))/(48.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*(3.0*fkpY*(2.0*fkp+(fmN*fmN))*(4.0*kappa_*M_P+g_*(fMass+fmN))-fkk*(-1.0*fkpY*(8.0*kappa_*M_P+g_*fmN)+fMass*fmY*(2.0*kappa_*M_P+g_*fmN)+fppY*(4.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmN))))/(48.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+(fmY*fmY))*(fkpY*(4.0*kappa_*M_P+g_*(fmY-fMass))-2.0*kappa_*fkk*M_P))/(16.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+(fmY*fmY))*(fkpY*(4.0*kappa_*M_P+g_*(fMass+fmY))-2.0*kappa_*fkk*M_P))/(16.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*(-1.0*fDenominator_u)));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*((fkk*fkk)*(fmY*fkp*(10.0*kappa_*fmY*M_P+g_*fmY*(fMass+fmN)-g_*fppY)+2.0*fppY*(fkpY*(6.0*kappa_*M_P+g_*fmN)+fmY*(kappa_*M_P*(fmN-3.0*fMass)-g_*fMass*fmN))+(fmN*fmN)*(fmY*fmY)*(6.0*kappa_*M_P+g_*(fMass+fmN))-(fppY*fppY)*(6.0*kappa_*M_P+5.0*g_*fMass+3.0*g_*fmN)+2.0*kappa_*M_P*fkpY*(fmN*fmY-fkpY))+fkk*(fkp*(fkpY*(2.0*fppY*(24.0*kappa_*M_P+7.0*g_*fMass+6.0*g_*fmN)+4.0*kappa_*fmY*M_P*(2.0*fMass+fmN)+g_*fmN*fmY*(3.0*fMass+fmN))-(fkpY*fkpY)*(24.0*kappa_*M_P+g_*(fMass+4.0*fmN))+(fmN*fmN)*(fmY*fmY)*(14.0*kappa_*M_P+3.0*g_*(fMass+fmN))+fppY*(-1.0*fppY*(8.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmN)+4.0*kappa_*fmY*M_P*(fmN-2.0*fMass)-g_*fmN*fmY*(2.0*fMass+fmN)))+fmN*(2.0*fkpY*(fmN*fmY*(kappa_*M_P*(2.0*fMass+fmN)+g_*fMass*fmN)+fppY*(g_*fmN*(4.0*fMass+3.0*fmN)-2.0*kappa_*M_P*(fMass-5.0*fmN)))+(fkpY*fkpY)*(-1.0*(2.0*kappa_*M_P*(fMass+9.0*fmN)+g_*fmN*(fMass+3.0*fmN)))+(fmN*fmY+fppY)*((fmN*fmN)*fmY*(4.0*kappa_*M_P+g_*(fMass+fmN))-(fMass+fmN)*fppY*(2.0*kappa_*M_P+g_*fmN)))+2.0*fmY*(fkp*fkp)*(6.0*kappa_*fmY*M_P+g_*fkpY+g_*fmY*(fMass+fmN)-g_*fppY))-fkpY*(2.0*fkp+(fmN*fmN))*(5.0*(fmN*fmN)*fkpY*(4.0*kappa_*M_P+g_*(fMass+fmN))+fkp*(2.0*(3.0*fkpY-2.0*fppY)*(4.0*kappa_*M_P+g_*(fMass+fmN))-4.0*kappa_*fMass*fmY*M_P-g_*fmN*fmY*(fMass+fmN))-2.0*g_*fmY*(fkp*fkp))+2.0*kappa_*pow(fkk,3.0)*(fmY*fmY)*M_P))/(320.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*((fkk*fkk)*(fmY*fkp*(fmY*(-10.0*kappa_*M_P+g_*fMass-g_*fmN)+g_*fppY)-2.0*fppY*(fkpY*(6.0*kappa_*M_P+g_*fmN)+kappa_*fmY*M_P*(3.0*fMass+fmN)+g_*fMass*fmN*fmY)+(fmN*fmN)*(fmY*fmY)*(g_*(fMass-fmN)-6.0*kappa_*M_P)+(fppY*fppY)*(6.0*kappa_*M_P-5.0*g_*fMass+3.0*g_*fmN)+2.0*kappa_*M_P*fkpY*(fkpY-fmN*fmY))+fkk*(fkp*(fkpY*(fmY*(8.0*kappa_*fMass*M_P-4.0*kappa_*fmN*M_P+3.0*g_*fMass*fmN-g_*(fmN*fmN))+2.0*fppY*(-24.0*kappa_*M_P+7.0*g_*fMass-6.0*g_*fmN))+(fkpY*fkpY)*(24.0*kappa_*M_P-g_*fMass+4.0*g_*fmN)+(fmN*fmN)*(fmY*fmY)*(3.0*g_*(fMass-fmN)-14.0*kappa_*M_P)+fmY*fppY*(g_*fmN*(fmN-2.0*fMass)-4.0*kappa_*M_P*(2.0*fMass+fmN))+(fppY*fppY)*(8.0*kappa_*M_P-5.0*g_*fMass+4.0*g_*fmN))-fmN*(fkpY*(fppY*(4.0*kappa_*M_P*(fMass+5.0*fmN)+2.0*g_*fmN*(3.0*fmN-4.0*fMass))-2.0*fmN*fmY*(kappa_*M_P*(2.0*fMass-fmN)+g_*fMass*fmN))+(fkpY*fkpY)*(2.0*kappa_*M_P*(fMass-9.0*fmN)+g_*fmN*(fMass-3.0*fmN))+(fmN*fmY+fppY)*((fmN*fmN)*fmY*(4.0*kappa_*M_P+g_*(fmN-fMass))+(fMass-fmN)*fppY*(2.0*kappa_*M_P+g_*fmN)))+2.0*fmY*(fkp*fkp)*(fmY*(-6.0*kappa_*M_P+g_*fMass-g_*fmN)-g_*fkpY+g_*fppY))+fkpY*(2.0*fkp+(fmN*fmN))*(5.0*(fmN*fmN)*fkpY*(4.0*kappa_*M_P+g_*(fmN-fMass))+fkp*(-2.0*(3.0*fkpY-2.0*fppY)*(g_*(fMass-fmN)-4.0*kappa_*M_P)+4.0*kappa_*fMass*fmY*M_P+g_*fmN*fmY*(fMass-fmN))-2.0*g_*fmY*(fkp*fkp))-2.0*kappa_*pow(fkk,3.0)*(fmY*fmY)*M_P))/(320.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+pow(fmY,2.0))*(fkpY*(fkpY*(4.0*fppY*(4.0*kappa_*M_P+g_*(fMass+fmY))+4.0*kappa_*fMass*fmN*M_P-2.0*g_*fmN*fkpY+g_*fmN*fmY*(fMass+fmY))-fkp*(5.0*pow(fmY,2.0)-6.0*fkpY)*(4.0*kappa_*M_P+g_*(fMass+fmY)))-fkk*(fkpY*(fppY*(28.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmY)+2.0*kappa_*fmN*M_P*(fMass+2.0*fmY)+g_*fMass*fmN*fmY)+fkp*(fkpY*(12.0*kappa_*M_P+g_*fmY)-2.0*kappa_*fmY*M_P*(fMass+4.0*fmY))+2.0*kappa_*fmY*M_P*(fMass-fmY)*(fmN*fmY+fppY)-g_*fmN*pow(fkpY,2.0))+2.0*kappa_*pow(fkk,2.0)*M_P*(fmN*fmY+5.0*fppY)))/(320.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+pow(fmY,2.0))*(fkpY*(fkpY*(4.0*fppY*(4.0*kappa_*M_P+g_*(fmY-fMass))-fmN*(4.0*kappa_*fMass*M_P+2.0*g_*fkpY+g_*fmY*(fMass-fmY)))+fkp*(5.0*pow(fmY,2.0)-6.0*fkpY)*(g_*(fMass-fmY)-4.0*kappa_*M_P))+fkk*(fkpY*(fppY*(-28.0*kappa_*M_P+5.0*g_*fMass-4.0*g_*fmY)+2.0*kappa_*fmN*M_P*(fMass-2.0*fmY)+g_*fMass*fmN*fmY)+fkp*(-1.0*fkpY*(12.0*kappa_*M_P+g_*fmY)-2.0*kappa_*fmY*M_P*(fMass-4.0*fmY))+2.0*kappa_*fmY*M_P*(fMass+fmY)*(fmN*fmY+fppY)+g_*fmN*pow(fkpY,2.0))+2.0*kappa_*pow(fkk,2.0)*M_P*(fmN*fmY+5.0*fppY)))/(320.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      if(particle.formfactorE == NULL)
	Fpi = particle.E;
      else
	Fpi = (*particle.formfactorE).value(particle.E, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);

      Fs = protonEMTFF(-fkk, particle.E > 0 ? fS : fU, particle.Y)*
	propagatorRegge(particle, fS, fU, fT, -fkk, particle.Z, &observ);
      
      f_ = particle.I;

      if(particle.E < 0)
	return (sqrt(2.0)*ELEC*f_*Fpi)/(fkpY*fDenominator_t);

      return ELEC*(sqrt(2.0)*f_*(2.0*Fpi*fkp+Fpi*fkk+Fs*(fDenominator_t)))/((2.0*fkp+fkk)*fkpY*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return ELEC*f_*Fpi*kappa_/(sqrt(2.0)*fmN*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return -ELEC*f_*Fpi*kappa_/(sqrt(2.0)*fmN*fDenominator_t);

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA2(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


complex< double >
TCalculateConsistentCoeff::CalcA3 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*kappa_*g_)/(M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram A
    case 3:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram B
    case 4:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*ELEC*kappa_*(gt_*(fmY-fmN)-gv_*(M_L+M_P)) / (1000.0*(M_L+M_P)*(fDenominator_t));

      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);

      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      return (ELEC*kappa_*g_)/(M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      return -1.0*((ELEC*kappa_*g_)/(M_P*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp*(6.0*kappa_*fmY*M_P+2.0*g_*fkpY-g_*fMass*fmY+g_*fmN*fmY-g_*fppY)+fmN*(fkpY*(-4.0*kappa_*M_P+g_*fMass+2.0*g_*fmN)+fppY*(-4.0*kappa_*M_P+g_*fMass-g_*fmN)+2.0*kappa_*fmY*M_P*(fMass+3.0*fmN)+g_*fmN*fmY*(fmN-fMass))+fkk*(6.0*kappa_*fmY*M_P+g_*fmN*fmY-3.0*g_*fppY));

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return (1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp*(6.0*kappa_*fmY*M_P+2.0*g_*fkpY+g_*fmY*(fMass+fmN)-g_*fppY)+fmN*(-1.0*fkpY*(4.0*kappa_*M_P+g_*fMass-2.0*g_*fmN)-fppY*(4.0*kappa_*M_P+g_*(fMass+fmN))-2.0*kappa_*fmY*M_P*(fMass-3.0*fmN)+g_*fmN*fmY*(fMass+fmN))+fkk*(6.0*kappa_*fmY*M_P+g_*fmN*fmY-3.0*g_*fppY));

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+(fmY*fmY))*(2.0*kappa_*M_P*(fMass+fmY)+g_*fkpY))/(8.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+(fmY*fmY))*(2.0*kappa_*M_P*(fmY-fMass)+g_*fkpY))/(8.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fmN,3.0)*fkpY*(fkpY*(8.0*kappa_*M_P+2.0*g_*fMass-3.0*g_*fmN)+2.0*fppY*(4.0*kappa_*M_P+g_*(fMass+fmN))+2.0*kappa_*fmY*M_P*(fMass-5.0*fmN)-2.0*g_*fmN*fmY*(fMass+fmN))+(fkk*fkk)*((fmY*fmY)*(2.0*kappa_*M_P*(fMass-fmN)+g_*fkp+g_*(fmN*fmN))+2.0*fmY*fppY*(5.0*kappa_*M_P+g_*fmN)-5.0*g_*(fppY*fppY))+fkk*(fkp*(fppY*(fmY*(28.0*kappa_*M_P+3.0*g_*fMass+5.0*g_*fmN)+10.0*g_*fkpY)-fkpY*(3.0*fmY*(4.0*kappa_*M_P+g_*fmN)+g_*fkpY)+3.0*(fmY*fmY)*(2.0*kappa_*M_P*(fMass-fmN)+g_*(fmN*fmN))-9.0*g_*(fppY*fppY))+fkpY*(-1.0*fkpY*(2.0*kappa_*M_P*(fMass-2.0*fmN)+g_*(fmN*fmN))-2.0*fmN*fmY*(kappa_*M_P*(fMass+5.0*fmN)+g_*(fmN*fmN)))+2.0*fppY*(fmN*fmY*(kappa_*M_P*(5.0*fmN-2.0*fMass)+g_*fmN*(fMass+fmN))-fkpY*(2.0*kappa_*fMass*M_P+g_*fmN*(fMass-3.0*fmN)))+(fmN*fmN)*(fmY*fmY)*(2.0*kappa_*M_P*(fMass-fmN)+g_*(fmN*fmN))-(fppY*fppY)*(2.0*kappa_*M_P*(fMass+2.0*fmN)+g_*fmN*(2.0*fMass+3.0*fmN))+2.0*g_*(fmY*fmY)*(fkp*fkp))+fmN*fkp*(fppY*(2.0*fkpY*(4.0*kappa_*M_P+g_*fMass+5.0*g_*fmN)-4.0*kappa_*fmY*M_P*(fMass-2.0*fmN)+g_*fmN*fmY*(fMass+fmN))+fkpY*(3.0*fkpY*(4.0*kappa_*M_P+g_*fMass-3.0*g_*fmN)-fmN*fmY*(32.0*kappa_*M_P+6.0*g_*fMass+7.0*g_*fmN))-(fppY*fppY)*(4.0*kappa_*M_P+g_*(fMass+fmN))+2.0*kappa_*fmN*(fmY*fmY)*M_P*(fMass-fmN))+2.0*(fkp*fkp)*(fppY*(8.0*kappa_*fmY*M_P+6.0*g_*fkpY+g_*fmY*(fMass+fmN))-fkpY*(12.0*kappa_*fmY*M_P+3.0*g_*fkpY+2.0*g_*fMass*fmY+3.0*g_*fmN*fmY)+2.0*kappa_*(fmY*fmY)*M_P*(fMass-fmN)-g_*(fppY*fppY)));

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fmN,3.0)*fkpY*(fkpY*(-8.0*kappa_*M_P+2.0*g_*fMass+3.0*g_*fmN)+2.0*fmY*(kappa_*M_P*(fMass+5.0*fmN)+g_*fmN*(fmN-fMass))+2.0*fppY*(-4.0*kappa_*M_P+g_*fMass-g_*fmN))-(fkk*fkk)*((fmY*fmY)*(-2.0*kappa_*M_P*(fMass+fmN)+g_*fkp+g_*(fmN*fmN))+2.0*fmY*fppY*(5.0*kappa_*M_P+g_*fmN)-5.0*g_*(fppY*fppY))+fkk*(fkp*(fkpY*(12.0*kappa_*fmY*M_P+3.0*g_*fmN*fmY-10.0*g_*fppY)+3.0*(fmY*fmY)*(2.0*kappa_*M_P*(fMass+fmN)-g_*(fmN*fmN))+fmY*fppY*(-28.0*kappa_*M_P+3.0*g_*fMass-5.0*g_*fmN)+g_*(fkpY*fkpY)+9.0*g_*(fppY*fppY))+fkpY*(fkpY*(g_*(fmN*fmN)-2.0*kappa_*M_P*(fMass+2.0*fmN))+2.0*fmN*fmY*(-kappa_*fMass*M_P+5.0*kappa_*fmN*M_P+g_*(fmN*fmN)))-2.0*fppY*(fkpY*(2.0*kappa_*fMass*M_P+g_*fmN*(fMass+3.0*fmN))+fmN*fmY*(kappa_*M_P*(2.0*fMass+5.0*fmN)+g_*fmN*(fmN-fMass)))+(fmN*fmN)*(fmY*fmY)*(2.0*kappa_*M_P*(fMass+fmN)-g_*(fmN*fmN))+(fppY*fppY)*(g_*fmN*(3.0*fmN-2.0*fMass)-2.0*kappa_*M_P*(fMass-2.0*fmN))-2.0*g_*(fmY*fmY)*(fkp*fkp))+fmN*fkp*(fppY*(2.0*fkpY*(-4.0*kappa_*M_P+g_*fMass-5.0*g_*fmN)-4.0*kappa_*fmY*M_P*(fMass+2.0*fmN)+g_*fmN*fmY*(fMass-fmN))+fkpY*(3.0*fkpY*(-4.0*kappa_*M_P+g_*fMass+3.0*g_*fmN)+fmN*fmY*(32.0*kappa_*M_P-6.0*g_*fMass+7.0*g_*fmN))+(fppY*fppY)*(4.0*kappa_*M_P-g_*fMass+g_*fmN)+2.0*kappa_*fmN*(fmY*fmY)*M_P*(fMass+fmN))+2.0*(fkp*fkp)*(fppY*(-8.0*kappa_*fmY*M_P-6.0*g_*fkpY+g_*fMass*fmY-g_*fmN*fmY)+fkpY*(12.0*kappa_*fmY*M_P+3.0*g_*fkpY-2.0*g_*fMass*fmY+3.0*g_*fmN*fmY)+2.0*kappa_*(fmY*fmY)*M_P*(fMass+fmN)+g_*(fppY*fppY)));

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+pow(fmY,2.0))*(fkpY*(-fmN*fkpY*(4.0*kappa_*M_P+g_*(fMass+2.0*fmY))-fppY*(4.0*kappa_*M_P*(3.0*fMass-2.0*fmY)-6.0*g_*fkpY+g_*fmY*(fMass+fmY))+fmN*fmY*(fMass+fmY)*(2.0*kappa_*M_P+g_*fmY))+fkk*(fkpY*(2.0*kappa_*fmN*M_P+g_*fmN*fmY-5.0*g_*fppY)-2.0*kappa_*fmY*M_P*fkp+2.0*kappa_*M_P*((5.0*fMass-4.0*fmY)*fppY-fMass*fmN*fmY))+fkp*(fkpY*(4.0*kappa_*M_P*(3.0*fmY-2.0*fMass)+4.0*g_*fkpY+g_*fmY*(fMass-4.0*fmY))+10.0*kappa_*pow(fmY,2.0)*M_P*(fMass-fmY)));

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fkpY+fkk+pow(fmY,2.0))*(fkpY*(fmN*fkpY*(g_*(fMass-2.0*fmY)-4.0*kappa_*M_P)+fppY*(4.0*kappa_*M_P*(3.0*fMass+2.0*fmY)+6.0*g_*fkpY+g_*fmY*(fMass-fmY))-fmN*fmY*(fMass-fmY)*(2.0*kappa_*M_P+g_*fmY))+fkk*(fkpY*(2.0*kappa_*fmN*M_P+g_*fmN*fmY-5.0*g_*fppY)-2.0*kappa_*fmY*M_P*fkp+2.0*kappa_*M_P*(fMass*fmN*fmY-(5.0*fMass+4.0*fmY)*fppY))+fkp*(fkpY*(4.0*kappa_*M_P*(2.0*fMass+3.0*fmY)+4.0*g_*fkpY-g_*fmY*(fMass+4.0*fmY))-10.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY)));

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return sqrt(2.0)*ELEC*f_*Fpi*(kappa_*fmN-kappa_*fmY+fmN)/(fmN*fDenominator_t);

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return sqrt(2.0)*ELEC*f_*Fpi*(1.0 + kappa_)/fDenominator_t;

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA3(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


complex< double >
TCalculateConsistentCoeff::CalcA4 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*kappa_*g_)/(M_P*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram A
    case 3:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE!=NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*kappa_*g_)/(M_P*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram B
    case 4:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(-1.0*ELEC*kappa_*(gt_*(fmY-fmN)+gv_*(M_L+M_P)))/(1000.0*(M_L+M_P)*(fDenominator_t));

      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      return -1.0*((ELEC*kappa_*g_)/(M_P*(-1.0*fDenominator_u)));

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      return (ELEC*kappa_*g_)/(M_P*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp+fkk+(fmN*fmN))*(2.0*kappa_*M_P*(fMass+fmN)-g_*fkp))/(8.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp+fkk+(fmN*fmN))*(2.0*kappa_*M_P*(fMass-fmN)+g_*fkp))/(8.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(2.0*fmN*fkpY*(6.0*kappa_*M_P-g_*fMass+g_*fmY)+fkp*(fmY*(-4.0*kappa_*M_P+g_*fMass+2.0*g_*fmY)-4.0*g_*fkpY)+fppY*(fmY*(4.0*kappa_*M_P-g_*fMass+g_*fmY)-2.0*g_*fkpY)-fkk*(6.0*kappa_*fmN*M_P+g_*fmN*fmY-3.0*g_*fppY)-fmN*fmY*(2.0*kappa_*M_P*(fMass+3.0*fmY)+g_*fmY*(fmY-fMass)));

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (1.0/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(-2.0*fmN*fkpY*(6.0*kappa_*M_P+g_*(fMass+fmY))+fkp*(fmY*(4.0*kappa_*M_P+g_*fMass-2.0*g_*fmY)+4.0*g_*fkpY)-fppY*(4.0*kappa_*fmY*M_P-2.0*g_*fkpY+g_*fmY*(fMass+fmY))+fkk*(6.0*kappa_*fmN*M_P+g_*fmN*fmY-3.0*g_*fppY)+fmN*fmY*(g_*fmY*(fMass+fmY)-2.0*kappa_*M_P*(fMass-3.0*fmY)));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp+fkk+(fmN*fmN))*((fkp*fkp)*(4.0*kappa_*fmY*M_P+4.0*g_*fkpY+g_*fMass*fmY+2.0*g_*fmN*fmY-6.0*g_*fppY)+fkp*(fkpY*(4.0*kappa_*M_P*(2.0*fMass-3.0*fmN)-g_*fmN*(fMass-4.0*fmN))+fmN*fmY*(fMass+fmN)*(2.0*kappa_*M_P+g_*fmN)-fppY*(4.0*kappa_*M_P*(3.0*fMass-2.0*fmN)+g_*fmN*(fMass+fmN)))+fkk*(fkp*(2.0*kappa_*fmY*M_P+g_*fmN*fmY-5.0*g_*fppY)+2.0*kappa_*M_P*(-fmN*fkpY+fMass*fmN*fmY+(4.0*fmN-5.0*fMass)*fppY))+10.0*kappa_*(fmN*fmN)*M_P*(fMass-fmN)*fkpY);

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      return (1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*(2.0*fkp+fkk+(fmN*fmN))*((fkp*fkp)*(4.0*kappa_*fmY*M_P+4.0*g_*fkpY-g_*fMass*fmY+2.0*g_*fmN*fmY-6.0*g_*fppY)+fkp*(fkpY*(g_*fmN*(fMass+4.0*fmN)-4.0*kappa_*M_P*(2.0*fMass+3.0*fmN))+fmN*fmY*(-1.0*(fMass-fmN))*(2.0*kappa_*M_P+g_*fmN)+fppY*(4.0*kappa_*M_P*(3.0*fMass+2.0*fmN)+g_*fmN*(fMass-fmN)))+fkk*(fkp*(2.0*kappa_*fmY*M_P+g_*fmN*fmY-5.0*g_*fppY)+2.0*kappa_*M_P*((5.0*fMass+4.0*fmN)*fppY-fmN*(fkpY+fMass*fmY)))-10.0*kappa_*(fmN*fmN)*M_P*(fMass+fmN)*fkpY);

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fkk,2.0)*(pow(fmN,2.0)*(2.0*kappa_*M_P*(fMass-fmY)-g_*fkpY+g_*pow(fmY,2.0))+2.0*fmN*fppY*(5.0*kappa_*M_P+g_*fmY)-5.0*g_*pow(fppY,2.0))+fkpY*(fmN*fppY*(2.0*fkpY*(8.0*kappa_*M_P+g_*(fMass+fmY))-fmY*(g_*fmY*(fMass+fmY)-4.0*kappa_*M_P*(fMass-2.0*fmY)))+pow(fppY,2.0)*(4.0*kappa_*fmY*M_P-2.0*g_*fkpY+g_*fmY*(fMass+fmY))-2.0*kappa_*pow(fmN,2.0)*M_P*(fMass-fmY)*(pow(fmY,2.0)-2.0*fkpY))+fkk*(pow(fmN,2.0)*(fkpY*(6.0*kappa_*M_P*(fmY-fMass)+2.0*g_*fkpY-3.0*g_*pow(fmY,2.0))+2.0*kappa_*pow(fmY,2.0)*M_P*(fMass-fmY)+g_*pow(fmY,4.0))+fkp*(2.0*fppY*(2.0*kappa_*fMass*M_P+5.0*g_*fkpY+g_*fmY*(fMass-3.0*fmY))-3.0*fmN*fkpY*(4.0*kappa_*M_P+g_*fmY)+2.0*fmN*fmY*(kappa_*M_P*(fMass+5.0*fmY)+g_*pow(fmY,2.0)))+fmN*fppY*(2.0*fmY*(kappa_*M_P*(5.0*fmY-2.0*fMass)+g_*fmY*(fMass+fmY))-fkpY*(28.0*kappa_*M_P+3.0*g_*fMass+5.0*g_*fmY))+pow(fkp,2.0)*(-2.0*kappa_*fMass*M_P+4.0*kappa_*fmY*M_P+g_*fkpY-g_*pow(fmY,2.0))-pow(fppY,2.0)*(2.0*kappa_*M_P*(fMass+2.0*fmY)-9.0*g_*fkpY+g_*fmY*(2.0*fMass+3.0*fmY)))+fkp*(fmY*fkpY*(2.0*fppY*(4.0*kappa_*M_P+g_*fMass+5.0*g_*fmY)-fmN*fmY*(32.0*kappa_*M_P+6.0*g_*fMass+7.0*g_*fmY))+2.0*pow(fkpY,2.0)*(12.0*kappa_*fmN*M_P+2.0*g_*fMass*fmN+3.0*g_*fmN*fmY-6.0*g_*fppY)+2.0*pow(fmY,3.0)*(-1.0*fppY*(4.0*kappa_*M_P+g_*(fMass+fmY))-kappa_*fmN*M_P*(fMass-5.0*fmY)+g_*fmN*fmY*(fMass+fmY)))+pow(fkp,2.0)*(-3.0*fmY*fkpY*(4.0*kappa_*M_P+g_*(fMass-3.0*fmY))+pow(fmY,3.0)*(8.0*kappa_*M_P+2.0*g_*fMass-3.0*g_*fmY)-6.0*g_*pow(fkpY,2.0)));

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      return (1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*(pow(fkk,2.0)*(pow(fmN,2.0)*(-2.0*kappa_*M_P*(fMass+fmY)-g_*fkpY+g_*pow(fmY,2.0))+2.0*fmN*fppY*(5.0*kappa_*M_P+g_*fmY)-5.0*g_*pow(fppY,2.0))+fkpY*(fmN*fppY*(2.0*fkpY*(8.0*kappa_*M_P-g_*fMass+g_*fmY)+fmY*(g_*fmY*(fMass-fmY)-4.0*kappa_*M_P*(fMass+2.0*fmY)))+pow(fppY,2.0)*(fmY*(4.0*kappa_*M_P-g_*fMass+g_*fmY)-2.0*g_*fkpY)+2.0*kappa_*pow(fmN,2.0)*M_P*(fMass+fmY)*(pow(fmY,2.0)-2.0*fkpY))+fkk*(pow(fmN,2.0)*(fkpY*(6.0*kappa_*M_P*(fMass+fmY)+2.0*g_*fkpY-3.0*g_*pow(fmY,2.0))-2.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY)+g_*pow(fmY,4.0))+fkp*(-2.0*fppY*(2.0*kappa_*fMass*M_P-5.0*g_*fkpY+g_*fmY*(fMass+3.0*fmY))-3.0*fmN*fkpY*(4.0*kappa_*M_P+g_*fmY)+2.0*fmN*fmY*(g_*pow(fmY,2.0)-kappa_*M_P*(fMass-5.0*fmY)))+fmN*fppY*(fkpY*(-28.0*kappa_*M_P+3.0*g_*fMass-5.0*g_*fmY)+2.0*fmY*(kappa_*M_P*(2.0*fMass+5.0*fmY)+g_*fmY*(fmY-fMass)))+pow(fkp,2.0)*(2.0*kappa_*M_P*(fMass+2.0*fmY)+g_*fkpY-g_*pow(fmY,2.0))+pow(fppY,2.0)*(2.0*kappa_*M_P*(fMass-2.0*fmY)+9.0*g_*fkpY+g_*fmY*(2.0*fMass-3.0*fmY)))+fkp*(fmY*fkpY*(fmN*fmY*(-32.0*kappa_*M_P+6.0*g_*fMass-7.0*g_*fmY)+2.0*fppY*(4.0*kappa_*M_P-g_*fMass+5.0*g_*fmY))+pow(fkpY,2.0)*(24.0*kappa_*fmN*M_P-4.0*g_*fMass*fmN+6.0*g_*fmN*fmY-12.0*g_*fppY)+2.0*pow(fmY,3.0)*(fppY*(-4.0*kappa_*M_P+g_*fMass-g_*fmY)+kappa_*fmN*M_P*(fMass+5.0*fmY)+g_*fmN*fmY*(fmY-fMass)))+pow(fkp,2.0)*(3.0*fmY*fkpY*(-4.0*kappa_*M_P+g_*fMass+3.0*g_*fmY)+pow(fmY,3.0)*(8.0*kappa_*M_P-2.0*g_*fMass-3.0*g_*fmY)-6.0*g_*pow(fkpY,2.0)));

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return sqrt(2.0)*ELEC*f_*Fpi/fDenominator_t;

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      if(particle.formfactorG == NULL)
	Fpi = particle.G;
      else
	Fpi = (*particle.formfactorG).value(particle.G, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ)/1.e3;

      kappa_ = particle.H;
      f_ = particle.I;

      return -1.0*sqrt(2.0)*ELEC*f_*Fpi*(1.0 + kappa_)/fDenominator_t;

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA4(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


complex< double >
TCalculateConsistentCoeff::CalcA5 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(2.0*e_*ELEC*g_*fkp)/(fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1: {
      double nucleoncharge_ = fNucleon_charge;
      double hyperoncharge_ = fHyperon_charge;
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
      {
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
	nucleoncharge_ = (*particle.formfactorE).value(fNucleon_charge,-1.0*fkk);
	hyperoncharge_ = (*particle.formfactorE).value(fHyperon_charge,-1.0*fkk);
      }
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fT);
      coefficient = (-2.0*e_*ELEC*g_*(fkpY-fkp))/(fkpY*fDenominator_t);
      if(observ.regge)
	{
	  ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
	  if (particle.E!=0.0)
	    {
	      // Add electric part of s-channel
	      // when there's a proton in the initial state
	      if (fNucleon_charge!=0.0)
		coefficient += (2.0*nucleoncharge_*ELEC*g_*fkp)/(fkpY*(fS-fmN*fmN));

	      // Add electric part of u-channel
	      // when there's a charged hyperon in the final state
	      else if (fHyperon_charge!=0.0)		
		coefficient += (2.0*hyperoncharge_*ELEC*g_)/(fU-fmY*fmY);
		
	      else
		{
		  cerr << "ERROR in TCalculateConsistentCoeff::CalcA5(...): "
		       << "Problem restoring gauge invariance for "
		       << particle.nickname
		       << " exchange.\n";
		  exit(1);
		}
	    }
	}
      else
	ReggePropagator = 1;
      return ReggePropagator*coefficient;
    }
      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(2.0*e_*ELEC*g_)/(fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram A
    case 3:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram B
    case 4:
      return 0.0;
 
      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*((complex<double>(0.0,1.0)*ELEC*kappa_*(fkp-fkpY)*(fkp*(gt_*(fMass*fMass)-gv_*(M_L+M_P)*(fmN+fmY))+fkpY*(gt_*(fMass*fMass)+gv_*(M_L+M_P)*(fmN+fmY))))/(1000.0*(fMass*fMass)*(M_L+M_P)*fkpY*fDenominator_t));

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkp*(fkpY*(16.0*kappa_*M_P-3.0*g_*fMass+2.0*g_*fmN)-fMass*fmY*(2.0*kappa_*M_P+g_*fmN)+fppY*(4.0*kappa_*M_P-3.0*g_*fMass+2.0*g_*fmN))+6.0*kappa_*M_P*(fkk+(fmN*fmN))*fkpY))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkp*(fkpY*(16.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmN)+fMass*fmY*(2.0*kappa_*M_P+g_*fmN)+fppY*(4.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmN))+6.0*kappa_*M_P*(fkk+(fmN*fmN))*fkpY))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkpY*(fppY*(4.0*kappa_*M_P-3.0*g_*fMass+2.0*g_*fmY)-fMass*fmN*(2.0*kappa_*M_P+g_*fmY))+fkp*(fkpY*(-16.0*kappa_*M_P+3.0*g_*fMass-2.0*g_*fmY)+6.0*kappa_*M_P*(fkk+(fmY*fmY)))))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*(-1.0*fDenominator_u)));

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkpY*(fMass*fmN*(2.0*kappa_*M_P+g_*fmY)+fppY*(4.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmY))+fkp*(6.0*kappa_*M_P*(fkk+(fmY*fmY))-fkpY*(16.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmY))))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
      	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*((fkp*fkp)*(fkk*fmY*(10.0*kappa_*fmY*M_P-g_*(fkpY+fppY)+g_*fmY*(fMass+fmN))+(fkpY*fkpY)*(24.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmN)-fkpY*(fmN*fmY*(4.0*kappa_*M_P+g_*fmN)+24.0*kappa_*M_P*fppY)+(fmN*fmN)*(fmY*fmY)*(14.0*kappa_*M_P+3.0*g_*(fMass+fmN))+fppY*(-1.0*fppY*(8.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmN)+4.0*kappa_*fmY*M_P*(fmN-2.0*fMass)-g_*fmN*fmY*(2.0*fMass+fmN)))+fkp*(fkk*(-1.0*fkpY*(fppY*(36.0*kappa_*M_P+5.0*g_*fMass+2.0*g_*fmN)+2.0*kappa_*fmY*M_P*(fMass+3.0*fmN)+g_*fMass*fmN*fmY)+(fkpY*fkpY)*(10.0*kappa_*M_P+g_*fmN)+(fmN*fmN)*(fmY*fmY)*(6.0*kappa_*M_P+g_*(fMass+fmN))-2.0*fmY*fppY*(kappa_*M_P*(3.0*fMass-fmN)+g_*fMass*fmN)-(fppY*fppY)*(6.0*kappa_*M_P+5.0*g_*fMass+3.0*g_*fmN))+fmN*(fmN*fkpY*(fppY*(-12.0*kappa_*M_P+3.0*g_*fMass+2.0*g_*fmN)+6.0*kappa_*fmY*M_P*(fMass-fmN)+g_*fMass*fmN*fmY)+(fkpY*fkpY)*(2.0*kappa_*M_P*(fMass+15.0*fmN)+g_*fmN*(4.0*fMass+3.0*fmN))+(fmN*fmY+fppY)*((fmN*fmN)*fmY*(4.0*kappa_*M_P+g_*(fMass+fmN))-(fMass+fmN)*fppY*(2.0*kappa_*M_P+g_*fmN)))+2.0*kappa_*(fkk*fkk)*(fmY*fmY)*M_P)+2.0*fmY*pow(fkp,3.0)*(6.0*kappa_*fmY*M_P-g_*(fkpY+fppY)+g_*fmY*(fMass+fmN))+2.0*kappa_*M_P*(fkk+(fmN*fmN))*fkpY*(fmN*((fMass+4.0*fmN)*fkpY+(fMass-fmN)*(fmN*fmY+fppY))-fkk*(fmN*fmY+5.0*fppY))))/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*((fkp*fkp)*(fkk*fmY*(fmY*(g_*(fMass-fmN)-10.0*kappa_*M_P)+g_*(fkpY+fppY))+(fkpY*fkpY)*(-24.0*kappa_*M_P+5.0*g_*fMass-4.0*g_*fmN)+fkpY*(fmN*fmY*(4.0*kappa_*M_P+g_*fmN)+24.0*kappa_*M_P*fppY)+(fmN*fmN)*(fmY*fmY)*(3.0*g_*(fMass-fmN)-14.0*kappa_*M_P)+fmY*fppY*(g_*fmN*(fmN-2.0*fMass)-4.0*kappa_*M_P*(2.0*fMass+fmN))+(fppY*fppY)*(8.0*kappa_*M_P-5.0*g_*fMass+4.0*g_*fmN))+fkp*(fkk*(-1.0*fppY*(fkpY*(-36.0*kappa_*M_P+5.0*g_*fMass-2.0*g_*fmN)+2.0*fmY*(kappa_*M_P*(3.0*fMass+fmN)+g_*fMass*fmN))-fkpY*(fkpY*(10.0*kappa_*M_P+g_*fmN)+2.0*kappa_*fmY*M_P*(fMass-3.0*fmN)+g_*fMass*fmN*fmY)+(fmN*fmN)*(fmY*fmY)*(g_*(fMass-fmN)-6.0*kappa_*M_P)+(fppY*fppY)*(6.0*kappa_*M_P-5.0*g_*fMass+3.0*g_*fmN))+fmN*(fmN*fkpY*(fppY*(12.0*kappa_*M_P+3.0*g_*fMass-2.0*g_*fmN)+6.0*kappa_*fmY*M_P*(fMass+fmN)+g_*fMass*fmN*fmY)+(fkpY*fkpY)*(2.0*kappa_*M_P*(fMass-15.0*fmN)+g_*fmN*(4.0*fMass-3.0*fmN))+(fmN*fmY+fppY)*((fmN*fmN)*fmY*(g_*(fMass-fmN)-4.0*kappa_*M_P)+(fmN-fMass)*fppY*(2.0*kappa_*M_P+g_*fmN)))-2.0*kappa_*(fkk*fkk)*(fmY*fmY)*M_P)+2.0*fmY*pow(fkp,3.0)*(fmY*(g_*(fMass-fmN)-6.0*kappa_*M_P)+g_*(fkpY+fppY))+2.0*kappa_*M_P*(fkk+(fmN*fmN))*fkpY*(fmN*((fMass-4.0*fmN)*fkpY+(fMass+fmN)*(fmN*fmY+fppY))+fkk*(fmN*fmY+5.0*fppY))))/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*(complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkpY*(pow(fmN,2.0)*(-2.0*fkpY+fkk+pow(fmY,2.0))*(-1.0*fkpY*(6.0*kappa_*M_P+g_*(fMass+fmY))+pow(fmY,2.0)*(4.0*kappa_*M_P+g_*(fMass+fmY))+2.0*kappa_*fkk*M_P)+fmN*fppY*(fkpY*(-4.0*kappa_*M_P*(fmY-2.0*fMass)-2.0*g_*fkpY+g_*fmY*(2.0*fMass+fmY))+fkk*(2.0*kappa_*M_P*(fmY-3.0*fMass)+g_*fkpY-2.0*g_*fMass*fmY)+2.0*kappa_*pow(fmY,2.0)*M_P*(fmY-fMass))-pow(fppY,2.0)*(-1.0*fkpY*(8.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmY)+fkk*(6.0*kappa_*M_P+5.0*g_*fMass+3.0*g_*fmY)+fmY*(fMass+fmY)*(2.0*kappa_*M_P+g_*fmY)))+fkp*(fmN*(fkpY*(fkk*(2.0*kappa_*M_P*(fMass+3.0*fmY)+g_*fMass*fmY)-pow(fmY,2.0)*(6.0*kappa_*M_P*(fMass-fmY)+g_*fMass*fmY))-pow(fkpY,2.0)*(fmY*(4.0*kappa_*M_P+g_*fmY)+g_*fkk)-2.0*kappa_*fmY*M_P*(fkk+pow(fmY,2.0))*(fkk+fmY*(fmY-fMass))+2.0*g_*pow(fkpY,3.0))+fppY*(fkpY*(fkk*(36.0*kappa_*M_P+5.0*g_*fMass+2.0*g_*fmY)+pow(fmY,2.0)*(12.0*kappa_*M_P-3.0*g_*fMass-2.0*g_*fmY))+2.0*kappa_*M_P*(fkk+pow(fmY,2.0))*(fmY*(fMass-fmY)-5.0*fkk)-24.0*kappa_*M_P*pow(fkpY,2.0)))+pow(fkp,2.0)*(pow(fkpY,2.0)*(-1.0*(24.0*kappa_*M_P+5.0*g_*fMass+4.0*g_*fmY))+fkpY*(fkk*(10.0*kappa_*M_P+g_*fmY)+fmY*(2.0*kappa_*M_P*(fMass+15.0*fmY)+g_*fmY*(4.0*fMass+3.0*fmY)))-2.0*kappa_*fmY*M_P*(fkk+pow(fmY,2.0))*(fMass+4.0*fmY))))/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkpY*(pow(fmN,2.0)*(-2.0*fkpY+fkk+pow(fmY,2.0))*(fkpY*(6.0*kappa_*M_P-g_*fMass+g_*fmY)+pow(fmY,2.0)*(-4.0*kappa_*M_P+g_*fMass-g_*fmY)-2.0*kappa_*fkk*M_P)-fmN*fppY*(fkpY*(-4.0*kappa_*M_P*(2.0*fMass+fmY)+g_*fkk+g_*fmY*(fmY-2.0*fMass))+2.0*fkk*(kappa_*M_P*(3.0*fMass+fmY)+g_*fMass*fmY)+2.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY)-2.0*g_*pow(fkpY,2.0))+pow(fppY,2.0)*(fkpY*(-8.0*kappa_*M_P+5.0*g_*fMass-4.0*g_*fmY)+fkk*(6.0*kappa_*M_P-5.0*g_*fMass+3.0*g_*fmY)+fmY*(-1.0*(fMass-fmY))*(2.0*kappa_*M_P+g_*fmY)))+fkp*(fmN*(-1.0*fkpY*(pow(fmY,2.0)*(6.0*kappa_*M_P*(fMass+fmY)+g_*fMass*fmY)-fkk*(2.0*kappa_*M_P*(fMass-3.0*fmY)+g_*fMass*fmY))+pow(fkpY,2.0)*(fmY*(4.0*kappa_*M_P+g_*fmY)+g_*fkk)+2.0*kappa_*fmY*M_P*(fkk+pow(fmY,2.0))*(fkk+fmY*(fMass+fmY))-2.0*g_*pow(fkpY,3.0))+fppY*(-1.0*fkpY*(fkk*(36.0*kappa_*M_P-5.0*g_*fMass+2.0*g_*fmY)+pow(fmY,2.0)*(12.0*kappa_*M_P+3.0*g_*fMass-2.0*g_*fmY))+2.0*kappa_*M_P*(fkk+pow(fmY,2.0))*(5.0*fkk+fmY*(fMass+fmY))+24.0*kappa_*M_P*pow(fkpY,2.0)))+pow(fkp,2.0)*(pow(fkpY,2.0)*(24.0*kappa_*M_P-5.0*g_*fMass+4.0*g_*fmY)+fkpY*(fmY*(2.0*kappa_*M_P*(fMass-15.0*fmY)+g_*fmY*(4.0*fMass-3.0*fmY))-fkk*(10.0*kappa_*M_P+g_*fmY))-2.0*kappa_*fmY*M_P*(fkk+pow(fmY,2.0))*(fMass-4.0*fmY))))/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fkpY*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      if(particle.formfactorE == NULL)
	Fpi = particle.E;
      else
	Fpi = (*particle.formfactorE).value(particle.E, -1.0*fkk);

      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);

      Fs = protonEMTFF(-fkk, particle.E > 0 ? fS : fU, particle.Y)*
	propagatorRegge(particle, fS, fU, fT, -fkk, particle.Z, &observ);

      f_ = particle.I;

      if(particle.E < 0)
	return (2.0*sqrt(2.0)*ELEC*f_*(fkpY*(-2.0*Fpi*fkpY+Fpi*fkk+Fpi*(pow(fmY,2.0)-pow(fmN,2.0))+Fs*(pow(fMass,2.0)-pow(fmN,2.0)-pow(fmY,2.0))+2.0*Fs*fppY)+Fpi*fkp*(2.0*fkpY-fkk+pow(fmN,2.0)-pow(fmY,2.0))))/(fkpY*(2.0*fkpY-fkk+pow(fmN,2.0)-pow(fmY,2.0))*(fDenominator_t));

      return ELEC*(2.0*sqrt(2.0)*f_*(fkp*(2.0*Fpi*fkpY-Fpi*fkk+Fs*((-1.0*fDenominator_t)))-2.0*Fpi*pow(fkp,2.0)+Fpi*fkk*fkpY))/((2.0*fkp+fkk)*fkpY*((-1.0*fDenominator_t)));

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      return 0.0;

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA5(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


complex< double >
TCalculateConsistentCoeff::CalcA6 ( int classindex,
                                    const Properties& particle,
                                    const Observable& observ) const
{
  switch ( classindex )
    {
      //----------------------------------------------------------------------------------------
      //Diagram S
    case 0:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(ELEC*g_*(2.0*e_*M_P*(fMass-fmN)+kappa_*fkk))/(2.0*M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram T
    case 1:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram U
    case 2:
      if (particle.formfactorE==NULL)
	e_ = particle.E;
      else
	e_ = (*particle.formfactorE).value(particle.E,-1.0*fkk);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(-1.0*((ELEC*g_*(2.0*e_*M_P*(fMass-fmY)+kappa_*fkk))/(2.0*M_P*fDenominator_u)));

      //----------------------------------------------------------------------------------------
      //Diagram A
    case 3:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE!=NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      if(!observ.electroprod)
	return 0.0;
      return ReggePropagator*(-1.0*(ELEC*g_*kappa_*fkk)/(2.0*M_P*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram B
    case 4:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram C
    case 5:
      if (particle.formfactorH==NULL)
	gv_ = particle.H ;
      else
	gv_ = (*particle.formfactorH).value(particle.H,fT);
      if (particle.formfactorI==NULL)
	gt_ = particle.I ;
      else
	gt_ = (*particle.formfactorI).value(particle.I,fT);
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (!observ.regge)
	ReggePropagator = 1;
      else
	ReggePropagator = propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      return ReggePropagator*(complex<double>(0.0,1.0)*ELEC*kappa_*(fkp-fkpY)*(gt_*(fmY-fmN)+gv_*(M_L+M_P)))/(1000.0*(M_L+M_P)*(fDenominator_t));

      //----------------------------------------------------------------------------------------
      //Diagram D, N*(1/2+)
    case 6:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if(!observ.electroprod)
	return 0.0;
      return (ELEC*g_*kappa_*fkk)/(2.0*M_P*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram E, N*(1/2-)
    case 7:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fS,particle.spin,fExponent,fWidth);
      if(!observ.electroprod)
	return 0.0;
      return -1.0*((ELEC*g_*kappa_*fkk)/(2.0*M_P*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram F, Y*(1/2+)
    case 8:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if(!observ.electroprod)
	return 0.0;
      return -1.0*((ELEC*kappa_*g_*fkk)/(2.0*M_P*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram G, Y*(1/2-)
    case 9:
      /* When no longitudinal coupling is desired the the
       * formfactor should be set to NULL */
      if (particle.formfactorE != NULL && !particle.long_coupling)
        {
	  cout << "Form factor electric term  for " << particle.nickname
	       << "-exchange without longitudinal coupling is non-zero." << endl << endl;
	  exit(1);
        }
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H ;
      else
	g_ = (*particle.formfactorH).value(particle.H,fU);
      if(!observ.electroprod)
	return 0.0;
      return (ELEC*kappa_*g_*fkk)/(2.0*M_P*fDenominator_u);

      //----------------------------------------------------------------------------------------
      //Diagram H, N*(3/2+)
    case 10:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*((complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkp*(4.0*kappa_*fmY*M_P+2.0*g_*(fkpY+fppY)-g_*fMass*fmY)+2.0*kappa_*M_P*(fMass+fmN)*(-2.0*fkpY+fmN*fmY-2.0*fppY)+2.0*kappa_*fkk*fmY*M_P))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s));

      //----------------------------------------------------------------------------------------
      //Diagram I, N*(3/2-)
    case 11:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkp*(fmY*(4.0*kappa_*M_P+g_*fMass)+2.0*g_*(fkpY+fppY))+2.0*kappa_*M_P*(fMass-fmN)*(2.0*fkpY-fmN*fmY+2.0*fppY)+2.0*kappa_*fkk*fmY*M_P))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*fDenominator_s);

      //----------------------------------------------------------------------------------------
      //Diagram J, Y*(3/2+)
    case 12:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fmN*fkpY*(g_*fMass-4.0*kappa_*M_P)+2.0*fkp*(2.0*kappa_*M_P*(fMass+fmY)+g_*fkpY)-2.0*fppY*(2.0*kappa_*M_P*(fMass+fmY)+g_*fkpY)+2.0*kappa_*fkk*fmN*M_P+2.0*kappa_*fmN*fmY*M_P*(fMass+fmY)))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram L, Y*(3/2-)
    case 13:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return (complex<double>(0.0,1.0)*ELEC*f_*fkk*(fmN*fkpY*(4.0*kappa_*M_P+g_*fMass)+fkp*(4.0*kappa_*M_P*(fMass-fmY)-2.0*g_*fkpY)+2.0*fppY*(2.0*kappa_*M_P*(fmY-fMass)+g_*fkpY)-2.0*kappa_*fkk*fmN*M_P+2.0*kappa_*fmN*fmY*M_P*(fMass-fmY)))/(24.0*(M_KP*M_KP)*pow(M_P,3.0)*(-1.0*fDenominator_u));

      //----------------------------------------------------------------------------------------
      //Diagram M, N*(5/2+)
    case 14:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
      	return 0.0;
      return (1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*fkk*((fkp*fkp)*(2.0*fkpY*(4.0*kappa_*fmY*M_P+g_*fkpY+g_*fmY*(fMass+fmN))+(fmY*fmY)*(4.0*kappa_*M_P*(fMass-fmN)+3.0*g_*(fmN*fmN))+2.0*g_*fppY*(fmN*fmY-2.0*fkpY)-6.0*g_*(fppY*fppY))+fkp*(fkpY*((fmN*fmN)*fmY*(8.0*kappa_*M_P+2.0*g_*fMass+g_*fmN)-(fMass-2.0*fmN)*fkpY*(g_*fmN-4.0*kappa_*M_P))+fppY*(fmN*fmY*(4.0*kappa_*fMass*M_P+g_*fmN*(fMass+fmN))-2.0*fMass*fkpY*(4.0*kappa_*M_P+g_*fmN))+(fmN*fmN)*(fmY*fmY)*(6.0*kappa_*M_P*(fMass-fmN)+g_*(fmN*fmN))-(fppY*fppY)*(4.0*kappa_*M_P*(3.0*fMass-2.0*fmN)+g_*fmN*(fMass+2.0*fmN)))+fkk*(fkp*(fkpY*(4.0*kappa_*fmY*M_P+g_*fmN*fmY-4.0*g_*fppY)+(fmY*fmY)*(2.0*kappa_*M_P*(fMass-fmN)+g_*(fmN*fmN))-fppY*(fmY*(4.0*kappa_*M_P+g_*fMass-g_*fmN)+4.0*g_*fppY))+2.0*kappa_*M_P*(fmN*(fmY*(fMass+fmN)*fkpY-(fkpY*fkpY)+fmN*(fmY*fmY)*(fMass-fmN))+(2.0*fMass-fmN)*fppY*(fmN*fmY-2.0*fkpY)+(3.0*fmN-4.0*fMass)*(fppY*fppY))+g_*(fmY*fmY)*(fkp*fkp))+2.0*kappa_*(fmN*fmN)*M_P*(fMass-fmN)*(fkpY*(2.0*fppY-fmN*fmY)+3.0*(fkpY*fkpY)+(fmN*fmN)*(fmY*fmY)-(fppY*fppY))-2.0*kappa_*(fkk*fkk)*fmY*M_P*fppY+2.0*g_*(fmY*fmY)*pow(fkp,3.0));

      //----------------------------------------------------------------------------------------
      //Diagram N, N*(5/2-)
    case 15:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fS,particle.spin,fExponent,fWidth) ;
      if(!observ.electroprod)
	return 0.0;
      return (1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*fDenominator_s))*complex<double>(0.0,1.0)*ELEC*f_*fkk*((fkp*fkp)*(-2.0*fkpY*(fmY*(4.0*kappa_*M_P+g_*(fmN-fMass))+g_*fkpY)+(fmY*fmY)*(4.0*kappa_*M_P*(fMass+fmN)-3.0*g_*(fmN*fmN))-2.0*g_*fppY*(fmN*fmY-2.0*fkpY)+6.0*g_*(fppY*fppY))-fkp*(fkpY*((fmN*fmN)*fmY*(8.0*kappa_*M_P-2.0*g_*fMass+g_*fmN)+2.0*fMass*fppY*(4.0*kappa_*M_P+g_*fmN))+(fMass+2.0*fmN)*(fkpY*fkpY)*(g_*fmN-4.0*kappa_*M_P)+(fmN*fmN)*(fmY*fmY)*(g_*(fmN*fmN)-6.0*kappa_*M_P*(fMass+fmN))+fmN*fmY*fppY*(g_*fmN*(fmN-fMass)-4.0*kappa_*fMass*M_P)+(fppY*fppY)*(4.0*kappa_*M_P*(3.0*fMass+2.0*fmN)+g_*fmN*(fMass-2.0*fmN)))+fkk*(-1.0*fkp*(fkpY*(4.0*kappa_*fmY*M_P+g_*fmN*fmY-4.0*g_*fppY)+(fmY*fmY)*(g_*(fmN*fmN)-2.0*kappa_*M_P*(fMass+fmN))+fppY*(-4.0*kappa_*fmY*M_P+g_*fmY*(fMass+fmN)-4.0*g_*fppY))+2.0*kappa_*M_P*((2.0*fMass+fmN)*fppY*(fmN*fmY-2.0*fkpY)+fmN*fkpY*(fkpY+fmY*(fMass-fmN))+(fmN*fmN)*(fmY*fmY)*(fMass+fmN)-(4.0*fMass+3.0*fmN)*(fppY*fppY))-g_*(fmY*fmY)*(fkp*fkp))+2.0*kappa_*(fmN*fmN)*M_P*(fMass+fmN)*(fkpY*(2.0*fppY-fmN*fmY)+3.0*(fkpY*fkpY)+(fmN*fmN)*(fmY*fmY)-(fppY*fppY))+2.0*kappa_*(fkk*fkk)*fmY*M_P*fppY-2.0*g_*(fmY*fmY)*pow(fkp,3.0));

      //----------------------------------------------------------------------------------------
      //Diagram O, Y*(5/2+)
    case 16:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkk*(pow(fmN,2.0)*(pow(fmY,2.0)-fkpY)*(2.0*kappa_*M_P*(fmY-fMass)+g_*fkpY)+fmN*fppY*(fkpY*(-4.0*kappa_*M_P-g_*fMass+g_*fmY)+2.0*kappa_*fmY*M_P*(fmY-2.0*fMass))+fkp*(4.0*fppY*(kappa_*M_P*(fmY-2.0*fMass)+g_*fkpY)-fmN*fkpY*(4.0*kappa_*M_P+g_*fmY)+2.0*kappa_*fmN*fmY*M_P*(fMass+fmY))+2.0*pow(fppY,2.0)*(4.0*kappa_*fMass*M_P-3.0*kappa_*fmY*M_P-2.0*g_*fkpY)+2.0*kappa_*fmY*M_P*pow(fkp,2.0))+pow(fmN,2.0)*(-3.0*pow(fmY,2.0)*fkpY+2.0*pow(fkpY,2.0)+pow(fmY,4.0))*(-1.0*(2.0*kappa_*M_P*(fMass-fmY)-g_*fkpY))+fkp*(fkpY*(2.0*fMass*fppY*(4.0*kappa_*M_P+g_*fmY)-fmN*pow(fmY,2.0)*(8.0*kappa_*M_P+2.0*g_*fMass+g_*fmY))+2.0*pow(fkpY,2.0)*(4.0*kappa_*fmN*M_P+g_*fmN*(fMass+fmY)-2.0*g_*fppY)+2.0*kappa_*pow(fmY,2.0)*M_P*(fmY-fMass)*(fmN*fmY-2.0*fppY))+fmN*fmY*fkpY*fppY*(4.0*kappa_*fMass*M_P-2.0*g_*fkpY+g_*fmY*(fMass+fmY))+pow(fppY,2.0)*(fkpY*(4.0*kappa_*M_P*(2.0*fmY-3.0*fMass)+6.0*g_*fkpY-g_*fmY*(fMass+2.0*fmY))+2.0*kappa_*pow(fmY,2.0)*M_P*(fMass-fmY))+pow(fkp,2.0)*((fMass-2.0*fmY)*fkpY*(4.0*kappa_*M_P-g_*fmY)+6.0*kappa_*pow(fmY,2.0)*M_P*(fmY-fMass)-2.0*g_*pow(fkpY,2.0))+2.0*kappa_*pow(fkk,2.0)*fmN*M_P*fppY);

      //----------------------------------------------------------------------------------------
      //Diagram Q, Y*(5/2-)
    case 17:
      if (particle.formfactorG==NULL)
	kappa_ = particle.G;
      else
	kappa_ = (*particle.formfactorG).value(particle.G,-1.0*fkk);
      if (particle.formfactorH==NULL)
	g_ = particle.H;
      else
	g_ = (*particle.formfactorH).value(particle.H,-1.0*fkk);
      if(particle.formfactorI==NULL)
	f_ = particle.I ;
      else
	f_ = (*particle.formfactorI).value(particle.I,fU) ;
      if(!observ.electroprod)
	return 0.0;
      return -1.0*(1.0/(160.0*pow(M_KP,4.0)*pow(M_P,5.0)*(-1.0*fDenominator_u)))*complex<double>(0.0,1.0)*ELEC*f_*fkk*(fkk*(pow(fmN,2.0)*(pow(fmY,2.0)-fkpY)*(2.0*kappa_*M_P*(fMass+fmY)+g_*fkpY)+fmN*fppY*(fkpY*(g_*(fMass+fmY)-4.0*kappa_*M_P)+2.0*kappa_*fmY*M_P*(2.0*fMass+fmY))+fkp*(4.0*fppY*(kappa_*M_P*(2.0*fMass+fmY)+g_*fkpY)-fmN*fkpY*(4.0*kappa_*M_P+g_*fmY)+2.0*kappa_*fmN*fmY*M_P*(fmY-fMass))-2.0*pow(fppY,2.0)*(4.0*kappa_*fMass*M_P+3.0*kappa_*fmY*M_P+2.0*g_*fkpY)+2.0*kappa_*fmY*M_P*pow(fkp,2.0))+pow(fmN,2.0)*(-3.0*pow(fmY,2.0)*fkpY+2.0*pow(fkpY,2.0)+pow(fmY,4.0))*(2.0*kappa_*M_P*(fMass+fmY)+g_*fkpY)+fkp*(fkpY*(-fmN*pow(fmY,2.0)*(8.0*kappa_*M_P-2.0*g_*fMass+g_*fmY)-2.0*fMass*fppY*(4.0*kappa_*M_P+g_*fmY))+pow(fkpY,2.0)*(8.0*kappa_*fmN*M_P-2.0*g_*fMass*fmN+2.0*g_*fmN*fmY-4.0*g_*fppY)+2.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY)*(fmN*fmY-2.0*fppY))+fmN*fmY*fkpY*fppY*(-4.0*kappa_*fMass*M_P-2.0*g_*fkpY+g_*fmY*(fmY-fMass))+pow(fppY,2.0)*(fkpY*(4.0*kappa_*M_P*(3.0*fMass+2.0*fmY)+6.0*g_*fkpY+g_*fmY*(fMass-2.0*fmY))-2.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY))+pow(fkp,2.0)*((fMass+2.0*fmY)*fkpY*(g_*fmY-4.0*kappa_*M_P)+6.0*kappa_*pow(fmY,2.0)*M_P*(fMass+fmY)-2.0*g_*pow(fkpY,2.0))+2.0*kappa_*pow(fkk,2.0)*fmN*M_P*fppY);

      //----------------------------------------------------------------------------------------
      //Diagram R
    case 18:
      if(particle.E > 0)
	return 0.0;
      
      if(particle.formfactorE == NULL)
	Fpi = particle.E;
      else
	Fpi = (*particle.formfactorE).value(particle.E, -1.0*fkk);
      
      Fpi *= propagatorRegge(particle, fS, fU, fT, 0., 0., &observ);
      
      Fs = protonEMTFF(-fkk, particle.E > 0 ? fS : fU, particle.Y)*
	propagatorRegge(particle, fS, fU, fT, -fkk, particle.Z, &observ);

      f_ = particle.I;

      return (sqrt(2.0)*ELEC*f_*Fs*(fmN-fmY))/(2.0*fkpY-fkk+pow(fmN,2.0)-pow(fmY,2.0));

      //----------------------------------------------------------------------------------------
      //Diagram V
    case 19:
      return 0.0;

      //----------------------------------------------------------------------------------------
      //Diagram W
    case 20:
      return 0.0;

    default: //error
      cerr << "Error in TCalculateConsistentCoeff::CalcA6(index, ...): index out of range";
      exit(1);
    }
  return 0;
}


//aconst calcAM TCalculateConsistentCoeff::fCalcMVector[CLASSMAX][6] = ;
