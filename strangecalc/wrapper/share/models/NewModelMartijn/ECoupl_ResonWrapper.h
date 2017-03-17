/*!
 * \file ECoupl_ResonWrapper.h
 * \ingroup wrapper
 *
 * \author Martijn Govers <martinus.govers@ugent.be> //:::ADDED BY MARTIJN:::
 
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

#ifndef ECOUPLRESONWRAPPER_H
#define ECOUPLRESONWRAPPER_H

#include "version.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include "strange_func.h"
#include <FormFactor.h>
#include "io_specific.h"
#include </home/mgovers/Software/strangecalcWorkDir/wrapper/src/FormFactorParametrization.h>

using std::string;
using std::strcmp;
using std::strncmp;
using std::strcpy;
using std::fixed;
using std::scientific;
using std::cerr;
using std::cout;
using std::endl;

namespace external{ // build a safe zone around the fortran function
  extern "C" {
    /* call Fortran code for including measured helicity amplitudes from https://userweb.jlab.org/~isupov/couplings */
    float ecoupl_(int &i, int &j, float &Q2);
  }
}

double ecoupling(int i,int j, double Q2)
{
  float Q2f = (float) (Q2/1000000.); // unit of energy is GeV instead of MeV in implementation of Fortran code
  return ((double) external::ecoupl_(i,j,Q2f))*sqrt(1000); // helicity amplitudes in Fortran code are in units of GeV^(-1/2)
}

/* Spin-1/2 Nucleon Resonance Parametrizations
 * ------------------------------------------- */

/*! Parametrization using measured helicity amplitudes for S11(1535) 1/2^- resonance
 * Dirac form factor
 *
 * = 0 in p(e,eK)Y, because transitioning (can be calculated if this is not the case)
 *
 */

double measured_S11_1535_dirac(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) {return 0.;} // not implemented

/*! Parametrization using measured helicity amplitudes for S11(1535) 1/2^- resonance
 * cfr. Van Cauteren et al. Eur. Phys. J. A (2004) 20: 283. doi:10.1140/epja/i2003-10158-3
 * Dirac form factor
 *
 * \verbatim
 
           1        / Q^2 |\  /|0               Mres - Mp |\  /|+              \
   ________________ | ___ | \/ |        (Q^2) + _________ | \/ |         (Q^2) |
          2         |     |    |1/2,1/2                   |    |1/2,-1/2       |
   Qm(Q^2)  Qp(Q^2) \ |p|                           2                          /
 = _____________________________________________________________________________

                              Mres - Mp     |\  /|+ 
                          _________________ | \/ |          (0)
                                   2        |    |1/2,-1/2
                          2 ( Qm(0)  Qp(0) )

 where
          ____________________
   Qp = \/ Q^2 + (Mres + Mp)^2  ,
          ____________________
   Qm = \/ Q^2 + (Mres - Mp)^2  , 

        / / Mp^2+Mres^2+Q^2 \2         \(1/2)
  |p| = | | ______________  |   - Mp^2 |
        \ \      2Mres      /          /

 \endverbatim
 *
 * -> no parameters
 */

double measured_S11_1535_pauli(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0)
{
  int j=18;
  double Mp = 938.272;
  double Mres = 1535.;
  double Qplus = sqrt(Q2+(Mres+Mp)*(Mres+Mp));
  double Qmin2 = Q2+(Mres-Mp)*(Mres-Mp);
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);

  double result = ((2.*Q2*ecoupling(3,j,Q2)) + ((Mres-Mp)*absp*ecoupling(1,j,Q2)))/(Qmin2*Qplus*absp);
  return (result*(Mres-Mp)*(Mres-Mp)*(Mres+Mp)/((Mres-Mp)*ecoupling(1,j,0.)));
}

/************************************************************************/

/*! Parametrization using measured helicity amplitudes for S11(1650) 1/2^- resonance
 * Dirac form factor
 *
 * = 0 in p(e,eK)Y, because transitioning (can be calculated if this is not the case)
 *
 */

double measured_S11_1650_dirac(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) {return 0.;} // not implemented

/*! Parametrization using measured helicity amplitudes for S11(1650) 1/2^- resonance
 * cfr. Van Cauteren et al. Eur. Phys. J. A (2004) 20: 283. doi:10.1140/epja/i2003-10158-3
 * Dirac form factor
 *
 * \verbatim
 
           1        / Q^2 |\  /|0               Mres - Mp |\  /|+              \
   ________________ | ___ | \/ |        (Q^2) + _________ | \/ |         (Q^2) |
          2         |     |    |1/2,1/2                   |    |1/2,-1/2       |
   Qm(Q^2)  Qp(Q^2) \ |p|                           2                          /
 = _____________________________________________________________________________

                              Mres - Mp     |\  /|+ 
                          _________________ | \/ |          (0)
                                   2        |    |1/2,-1/2
                          2 ( Qm(0)  Qp(0) )

 where
          ____________________
   Qp = \/ Q^2 + (Mres + Mp)^2  ,
          ____________________
   Qm = \/ Q^2 + (Mres - Mp)^2  , 

        / / Mp^2+Mres^2+Q^2 \2         \(1/2)
  |p| = | | ______________  |   - Mp^2 |
        \ \      2Mres      /          /

 \endverbatim
 *
 * -> no parameters
 */

double measured_S11_1650_pauli(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0)
{
  int j=19;
  double Mp = 938.272;  
  double Mres = 1650.; 
  double Qplus = sqrt(Q2+(Mres+Mp)*(Mres+Mp));
  double Qmin2 = Q2+(Mres-Mp)*(Mres-Mp);
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);

  double result = ((2.*Q2*ecoupling(3,j,Q2)) + ((Mres-Mp)*absp*ecoupling(1,j,Q2)))/(Qmin2*Qplus*absp);
  return (result*(Mres-Mp)*(Mres-Mp)*(Mres+Mp)/((Mres-Mp)*ecoupling(1,j,0.)));
}

/************************************************************************/

/*! Parametrization using measured helicity amplitudes for S11(1535) 1/2^- resonance
 * Dirac form factor
 *
 * = 0 in p(e,eK)Y, because transitioning (can be calculated if this is not the case)
 *
 */

double measured_P11_1440_dirac(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) {return 0.;} // not implemented

/*! Parametrization using measured helicity amplitudes for S11(1535) 1/2^- resonance
 * cfr. Van Cauteren et al. Eur. Phys. J. A (2004) 20: 283. doi:10.1140/epja/i2003-10158-3
 * Dirac form factor
 *
 * \verbatim
 
          -1        / Q^2 |\  /|0               Mres + Mp |\  /|+              \
   ________________ | ___ | \/ |        (Q^2) - _________ | \/ |         (Q^2) |
          2         |     |    |1/2,1/2                   |    |1/2,-1/2       |
   Qp(Q^2)  Qm(Q^2) \ |p|                           2                          /
 = _____________________________________________________________________________

                              Mres + Mp     |\  /|+ 
                          _________________ | \/ |          (0)
                                   2        |    |1/2,-1/2
                          2 ( Qp(0)  Qm(0) )

 where
          ____________________
   Qp = \/ Q^2 + (Mres + Mp)^2  ,
          ____________________
   Qm = \/ Q^2 + (Mres - Mp)^2  , 

        / / Mp^2+Mres^2+Q^2 \2         \(1/2)
  |p| = | | ______________  |   - Mp^2 |
        \ \      2Mres      /          /

 \endverbatim
 *
 * -> no parameters
 */

double measured_P11_1440_pauli(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0)
{
  int j=20;
  double Mp = 938.272;
  double Mres = 1440.0;
  double Qplus2 = Q2+(Mres+Mp)*(Mres+Mp);
  double Qmin = sqrt(Q2+(Mres-Mp)*(Mres-Mp));
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);

  double result = -1.*((2.*Q2*ecoupling(3,j,Q2)) - ((Mres+Mp)*absp*ecoupling(1,j,Q2)))/(Qplus2*Qmin*absp);
  return (result*(Mres+Mp)*(Mres+Mp)*(Mres-Mp)/((Mres+Mp)*ecoupling(1,j,0.)));
}

/* Spin-3/2 Nucleon Resonance Parametrizations
 * ------------------------------------------- */

/*! Parametrization using measured helicity amplitudes for P11(1720) 3/2^+ resonance
 * kappa_1 form factor
 *
 * \verbatim

           1        /   ___ |\  /|+		   |\  /|+	       \
   ________________ | \/ 3  | \/ |         (Q^2) + | \/ |	 (Q^2) |
   |p(Q^2)| Qp(Q^2) \       |    |1/2,-1/2	   |    |3/2,1/2       /
 = _____________________________________________________________________

              1      /   ___ |\  /|+              |\  /|+          \
        ____________ | \/ 3  | \/ |         (0) + | \/ |       (0) |
                     |       |    |1/2,-1/2       |    |3/2,1/2    |
        |p(0)| Qp(0) \				                   /
   
 where
          ____________________
   Qp = \/ Q^2 + (Mres + Mp)^2  ,

        / / Mp^2+Mres^2+Q^2 \2        \(1/2)
  |p| = | | _______________ |  - Mp^2 |
        \ \      2Mres      /         /

 \endverbatim
 *
 */

double measured_P13_1720_1(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) 
{ 
  int j=17;
  double Mp = 938.272;
  double Mres = 1720.0;
  double Qplus = sqrt(Q2+(Mres+Mp)*(Mres+Mp));
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);

  double result = (sqrt(3.)*ecoupling(1,j,Q2) + ecoupling(2,j,Q2))/(Qplus*absp);
  return (result*(Mres+Mp)*(Mres*Mres-Mp*Mp)) / (2.*Mres*(sqrt(3.)*ecoupling(1,j,0.) + ecoupling(2,j,0.)));
}

/*! Parametrization using measured helicity amplitudes for P13(1720) 3/2^+ resonance
 * kappa_2 form factor
 *
 * \verbatim

           1          /       |\  /|0		                       |\  /|+	           
   __________________ | 6 Q^2 | \/ |       (Q^2) + 3 (Mp-Mres)|p(Q^2)| | \/ |	    (Q^2) 
   |p(Q^2)|^3 Qm(Q^2) \       |    |1/2,1/2                            |    |1/2,-1/2      

                                          ___ / Q^2 + Mp^2 - MpMres \          |\  /|+            \
                                      - \/ 3  | ___________________ | |p(Q^2)| | \/ |       (Q^2) |
                                              \         Mres        /          |    |3/2,1/2      /
 = ________________________________________________________________________________________________

              1       /             |\  /|+                ___  Mp            |\  /|+          \
       ______________ | 3 (Mp-Mres) | \/ |         (0) - \/ 3  ____ (Mp-Mres) | \/ |       (0) |
       |p(0)|^2 Qm(0) \             |    |1/2,-1/2             Mres           |    |3/2,1/2    /
   
 where
          ____________________
   Qm = \/ Q^2 + (Mres - Mp)^2  ,

        / / Mp^2+Mres^2+Q^2 \2        \(1/2)
  |p| = | | _______________ |  - Mp^2 |
        \ \      2Mres      /         /

 \endverbatim
 *
 */

double measured_P13_1720_2(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) 
{ 
  int j=17;
  double Mp = 938.272;
  double Mres = 1720.0;
  double Qmin = sqrt(Q2+(Mres-Mp)*(Mres-Mp));
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);
  double absp0 = (Mres*Mres-Mp*Mp)/(2.*Mres);

  double result = (6.*Q2*ecoupling(3,j,Q2) + 3.*(Mp-Mres)*absp*ecoupling(1,j,Q2)-sqrt(3.)*(Q2+Mp*Mp-Mp*Mres)*absp*ecoupling(2,j,Q2)/Mres) / (absp*absp*absp*Qmin);
  return (result * absp0*absp0*(Mres-Mp) / (3.*(Mp-Mres)*ecoupling(1,j,0.)-sqrt(3.)*Mp*(Mp-Mres)*ecoupling(2,j,0.)/Mres) );
}

/************************************************************************/

/*! Parametrization using measured helicity amplitudes for D13(1520) 3/2^- resonance
 * kappa_1 form factor
 *
 * \verbatim

           1        /   ___ |\  /|+		   |\  /|+	       \
   ________________ | \/ 3  | \/ |         (Q^2) - | \/ |	 (Q^2) |
   |p(Q^2)| Qm(Q^2) \       |    |1/2,-1/2	   |    |3/2,1/2       /
 = _____________________________________________________________________

              1      /   ___ |\  /|+              |\  /|+          \
        ____________ | \/ 3  | \/ |         (0) - | \/ |       (0) |
                     |       |    |1/2,-1/2       |    |3/2,1/2    |
        |p(0)| Qm(0) \				                   /
   
 where
          ____________________
   Qm = \/ Q^2 + (Mres - Mp)^2  ,

        / / Mp^2+Mres^2+Q^2 \2        \(1/2)
  |p| = | | _______________ |  - Mp^2 |
        \ \      2Mres      /         /

 \endverbatim
 *
 */

double measured_D13_1520_1(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) 
{ 
  int j=15;
  double Mp = 938.272;
  double Mres = 1520.0;
  double Qmin = sqrt(Q2+(Mres-Mp)*(Mres-Mp));
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);

  double result = (sqrt(3.)*ecoupling(1,j,Q2) - ecoupling(2,j,Q2))/(Qmin*absp);
  return (result*(Mres-Mp)*(Mres*Mres-Mp*Mp)) / (2.*Mres*(sqrt(3.)*ecoupling(1,j,0.) - ecoupling(2,j,0.)));
}


/*! Parametrization using measured helicity amplitudes for D13(1520) 3/2^- resonance
 * kappa_2 form factor
 *
 * \verbatim

           1          /       |\  /|0		                       |\  /|+	           
   __________________ | 6 Q^2 | \/ |       (Q^2) + 3 (Mp+Mres)|p(Q^2)| | \/ |	    (Q^2) 
   |p(Q^2)|^3 Qp(Q^2) \       |    |1/2,1/2                            |    |1/2,-1/2      

                                          ___ / Q^2 + Mp^2 + MpMres \          |\  /|+            \
                                      - \/ 3  | ___________________ | |p(Q^2)| | \/ |       (Q^2) |
                                              \         Mres        /          |    |3/2,1/2      /
 = ________________________________________________________________________________________________

              1       /             |\  /|+                ___  Mp            |\  /|+          \
       ______________ | 3 (Mp+Mres) | \/ |         (0) - \/ 3  ____ (Mp+Mres) | \/ |       (0) |
       |p(0)|^2 Qp(0) \             |    |1/2,-1/2             Mres           |    |3/2,1/2    /
   
 where
          ____________________
   Qp = \/ Q^2 + (Mres + Mp)^2  ,

        / / Mp^2+Mres^2+Q^2 \2        \(1/2)
  |p| = | | _______________ |  - Mp^2 |
        \ \      2Mres      /         /

 \endverbatim
 *
 */

double measured_D13_1520_2(double Q2,double p1=0,double p2=0,double p3=0,
		       double p4=0,double p5=0,double p6=0,double p7=0,
		       double p8=0,double p9=0,double p10=0,double p11=0,
		       double p12=0,double p13=0,double p14=0,double p15=0) 
{ 
  int j=15;
  double Mp = 938.272;
  double Mres = 1520.0;
  double Qplus = sqrt(Q2+(Mres+Mp)*(Mres+Mp));
  double absp = sqrt(((Mp*Mp+Mres*Mres+Q2)*(Mp*Mp+Mres*Mres+Q2)/(4.*Mres*Mres))-Mp*Mp);
  double absp0 = (Mres*Mres-Mp*Mp)/(2.*Mres);

  double result = (6.*Q2*ecoupling(3,j,Q2) + 3.*(Mp+Mres)*absp*ecoupling(1,j,Q2)-sqrt(3.)*(Q2+Mp*Mp+Mp*Mres)*absp*ecoupling(2,j,Q2)/Mres) / (absp*absp*absp*Qplus);
  return (result * absp0*absp0*(Mres+Mp) / (3.*(Mp+Mres)*ecoupling(1,j,0.)-sqrt(3.)*Mp*(Mp+Mres)*ecoupling(2,j,0.)/Mres) );
}

/*! Specify the formfactor obtained from the measured helicity amplitudes for a certain resonance state.
 * Returns NULL if and only if it should do nothing (for instance, if the parametrization is not measured)
 */
FormFactor* specify_external_ff(int resonanceValue, int formfactorType)
{
  ffCalculator ff_parametrization;
  switch (resonanceValue)
    {
      case 18: // S11(1535)
        if (formfactorType==0) // dirac
          ff_parametrization = measured_S11_1535_dirac;
        else if (formfactorType==1) // pauli
          ff_parametrization = measured_S11_1535_pauli;
        else // do nothing
          return NULL;
	break;
      case 19: // S11(1650)
        if (formfactorType==0) // dirac
          ff_parametrization = measured_S11_1650_dirac;
        else if (formfactorType==1) // pauli
          ff_parametrization = measured_S11_1650_pauli;
        else // do nothing
          return NULL;
	break;
      case 20: // P11(1440)
        if (formfactorType==0) // dirac
          ff_parametrization = measured_P11_1440_dirac;
        else if (formfactorType==1) // pauli
          ff_parametrization = measured_P11_1440_pauli;
        else // do nothing
          return NULL;
	break;
      case 15: // D13(1520)
        if (formfactorType==1) // kappa1
	  ff_parametrization = measured_D13_1520_1;
	else if (formfactorType==2) // kappa2
	  ff_parametrization = measured_D13_1520_2;
	else // do nothing
          return NULL;
	break;
      case 17: // P13(1720)
	if (formfactorType==1) // kappa1
	  ff_parametrization = measured_P13_1720_1;
	else if (formfactorType==2) // kappa2
	  ff_parametrization = measured_P13_1720_2;
	else // do nothing
          return NULL;
	break;
      case 6: // F15(1685)
	cerr << "Measured EM Form Factor of F15(1685) exists but has not been implemented yet.\n\n";
	exit(1);
      case 21: // S31(1620)
	cerr << "Measured EM Form Factor of S31(1620) exists but has not been implemented yet.\n\n";
	exit(1);
      case 25: // P33(1232)
	cerr << "Measured EM Form Factor of P33(1232) exists but has not been implemented yet.\n\n";
	exit(1);
      case 13: // D33(1700)
	cerr << "Measured EM Form Factor of D33(1700) exists but has not been implemented yet.\n\n";
	exit(1);
      case 14: // P13 missing (1720)
	cerr << "Measured EM Form Factor of P13 missing (1720) exists but has not been implemented yet.\n\n";
	exit(1);
      default:
	// the helicity amplitudes have not been measured yet. Do nothing
	return NULL;
    }
  if (ff_parametrization==NULL)
  {
    cout << "Unexpected Error in recognizing parametrization of measured form factors.";
    exit(1);
  }

  // Allocate some dynamic memory and create the correct FormFactor object
  FormFactor* newformfactor = new FormFactor(ff_parametrization, 0., 0., 0.,
						    0., 0., 0., 0., 0., 0.,
						    0., 0., 0., 0., 0., 0.);
    // Notify when something went wrong. If function returned NULL then nothing happened.
  if (newformfactor==NULL)
    {
      cout << "Error in allocation of formfactor!" << endl << endl;
      exit(1);
    }

  // change the passed pointer of the FormFactor to the correct FormFactor
  return newformfactor;
}




/*!\brief Method determining the input parameter value for the Fortran code ecoupl
 * \returns value if existing else -1 
 */
int determineResonanceValue(char* nickname){
/*
Class  NickNam    Name            J^pi     Mass_(MeV)    Width_(MeV) 	Value
-----  -------    ----            ----     ---------     ----------- 	-----
							  
 D        N1      N(1440)        1/2^+     1440.0        300.0      	20
 I        N2      N(1520)        3/2^-     1520.0        115.0      	15
 E        N3      N(1535)        1/2^-     1535.0        150.0       	18
 E        N4      N(1650)        1/2^-     1650.0        150.0      	19
 I        N5      N(1700)        3/2^-     1700.0        100.0      	 
 D        N6      N(1710)        1/2^+     1710.0        100.0      	
 H        N7      N(1720)        3/2^+     1720.0        150.0      	17
 I        N8      m-D13(1895)    3/2^-     1895.0        200.0		 
 H        N9      P13(1900)      3/2^+     1900.0        500.0		 
 E        N10     m-S11(1895)	 1/2^-     1895.0        200.0		 
 D        N11     m-P11(1895)  	 1/2^+     1895.0	 200.0		 
 E        N12     S11_Capstick   1/2^-     1945.0	 595.0		
 D        N13     P11_Capstick   1/2^+     1975.0         45.0      	
 I        N14     D13_Capstick   3/2^-     1950.0        140.0		 
 H	  N15	  P13_Capstick	 3/2^+     1960.0        535.0		 
 I	  N19	  N(2050)	 3/2^-	   2050.0	 300.0		
 H	  N20	  N(2050)	 3/2^+	   2050.0	 300.0		 
 E	  N21	  N(2050)	 1/2^-	   2050.0	 300.0		 
 D	  N22	  N(2050)	 1/2^+	   2050.0	 300.0		
 M	  N23	  N(1680)	 5/2^+	   1685.0	 130.0		6
 N	  N24	  N(1675)	 5/2^-	   1675.0	 150.0		
 M	  N25	  N(2000)	 5/2^+	   2000.0	 140.0		 
 N	  N26	  N(2200)	 5/2^-	   2200.0	 250.0		
 E        N27     N(1690)_2010   1/2^-     1690.0         80.0           
 H        N28     N(1920)_2010   3/2^+     1920.0        440.0	         
 I        N29     N(2100)_2010   3/2^-     2100.0        200.0          
 E        D1      Delta(1620)    1/2^-     1620.0        150.0    	21
 E        D2      Delta(1900)    1/2^-     1900.0        200.0     	 
 D        D3      Delta(1910)    1/2^+     1910.0        250.0     	
 H        D4      Delta(1232)    3/2^+     1232.0        120.0      	25
 H        D5      Delta(1600)    3/2^+     1600.0        350.0      	 
 I        D6      Delta(1700)    3/2^-     1700.0        300.0      	13
 H        D7      Delta(1920)    3/2^+     1920.0        200.0		 
 I        D8      D33(1940)      3/2^-     1940.0        286.1		 
 E        D9      D13(2150)      1/2^-     2150.0	 200.0		 
 E        D10     S13(2150)-II   1/2^-     2150.0	 300.0		
 M	  D11	  Delta(1905)	 5/2^+	   1890.0	 330.0		
 N	  D12	  Delta(1930)	 5/2^-	   1960.0	 360.0		
 M	  D13	  Delta(2000)	 5/2^+	   2000.0	 195.0		 
*/
  if (!strcmp(nickname,"N1")) return 20;
  else if (!strcmp(nickname,"N2")) return 15;
  else if (!strcmp(nickname,"N3")) return 18;
  else if (!strcmp(nickname,"N4")) return 19;
  else if (!strcmp(nickname,"N7")) return 17;
  else if (!strcmp(nickname,"N23")) return 6;
  else if (!strcmp(nickname,"D1")) return 21;
  else if (!strcmp(nickname,"D4")) return 54;
  else if (!strcmp(nickname,"D6")) return 13;
  else return -1;
}

#endif
