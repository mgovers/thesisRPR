/*!
 * \file calcMatrixElement.cpp
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
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

#include <iostream>
#include <complex>
#include "calcMatrixElement.h"
#include "TMatrixElement.h"
#include <TCalculateCGLNMatrix.h>
using namespace std;

/*!
 * Determine the squared matrix element for photopolarization
 * given the specific kinematics and all info about the process
 * and the observables.
 */
double determine_M2_photo(double w, double k,double costheta_k, double pk, 
			  Class particles[], Observable* observ, int label)
{
  /* Construct the matrixelement
   * ***************************/
  TMatrixElement* me = TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles, observ, label);

  /* Calculate the squared matrix element
   * ************************************/
  double m2 = determine_M2_photo_withME(me,observ); 
  
  return m2;
}


/**********************************************************************************/


/*!
 * Determine the squared matrix element for photopolarization
 * given a matrix element and a specification of the observables.
 */
double determine_M2_photo_withME( TMatrixElement* me, Observable* observ)
{
  double m2 = 0; // this variable will hold the final result

  // Structure containing all info about polarizations of the process
  Polar polarization = observ->photo.pol;

  if(polarization.nopol)
    {
      /* Determine unpolarized matrixelement
       * ***********************************/
      // Loop over proton polarization
      for( int Lp=0; Lp<=1; Lp++)
	{
	  // Loop over hyperon polarization
	  for( int Ly=0; Ly<=1; Ly++)
	    {
	      // Loop over photon polarization
	      for( int L=0; L<=2; L++)
		{
		  // add squared matrixelement contribution for
		  // specific polarizations
		  m2 += norm(me->calculateM(L,Lp,Ly));

		} // L loop
	    } // Ly loop
	} // Lp loop
    }
  else if(polarization.spol)
    {
      if(polarization.sinpol.phopol)
	{
	  /* Determine matrixelement <g_5> for polarized photon
	   * **************************************************/
	  // Loop over proton polarization
	  for( int Lp=0; Lp<=1; Lp++)
	    {
	      // Loop over hyperon polarization
	      for( int Ly=0; Ly<=1; Ly++)
		{
		  // add squared matrixelement contribution for L=-1 (y-axis)
		  m2 += norm(me->calculateM(0,Lp,Ly));

		  // subtract squared matrixelement contrib for L=+1 (x-axis)
		  m2 -= norm(me->calculateM(2,Lp,Ly));

		} // Ly loop
	    } // Lp loop
	}
      else if(polarization.sinpol.tarpol)
	{
	  /* Determine matrixelement <g_5> for polarized target
	   * **************************************************/
	  // Loop over hyperon polarization
	  for( int Ly=0; Ly<=1; Ly++)
	    {
	      // Loop over photon polarization
	      for( int L=0; L<=2; L++)
		{ 
		  // add squared matrixelement contribution for Lp=+1/2
		  m2 += norm(me->calculateM(L,1,Ly));
		  
		  // subtract squared matrixelement contrib for Lp=-1/2
		  m2 -= norm(me->calculateM(L,0,Ly));

		} // L loop
	    } // Ly loop
	}
      else if(polarization.sinpol.recpol)
	{
	  /* Determine matrixelement <g_5> for polarized recoil
	   * **************************************************/
	  // Loop over target polarization
	  for( int Lp=0; Lp<=1; Lp++)
	    {
	      // Loop over photon polarization
	      for( int L=0; L<=2; L++)
		{ 
		  // add squared matrixelement contribution for Ly=+1/2
		  m2 += norm(me->calculateM(L,Lp,1));
		  
		  // subtract squared matrixelement contrib for Ly=-1/2
		  m2 -= norm(me->calculateM(L,Lp,0));

		} // L loop
	    } // Lp loop
	} 
    } // single polarization
  else if(polarization.dpol)
    {
      if(polarization.doubpol.beamtar)
	{
	  /* Determine matrixelement <g_5^+> for polarized beam/target
	   * *********************************************************/
	  // Loop over hyperon polarization
	  for( int Ly=0; Ly<=1; Ly++)
	    {
	      // add squared matrixelement contribution for L=+1 Lp=+1/2
	      m2 += norm(me->calculateM(2,1,Ly));
	      
	      // subtract squared matrixelement contrib for L=+1 Lp=-1/2
	      m2 -= norm(me->calculateM(2,0,Ly));
	    }
	  // if the beam helicity is -1, there's an additional minus-sign
	  m2 *= polarization.doubpol.beamhel;
	}
      else if(polarization.doubpol.beamrec)
	{
	  /* Determine matrixelement <g_5^+> for polarized beam/recoil
	   * *********************************************************/
	  // Loop over target polarization
	  for( int Lp=0; Lp<=1; Lp++)
	    {
	      // add squared matrixelement contribution for L=+1 Ly=+1/2
	      m2 += norm(me->calculateM(2,Lp,1));
	      
	      // subtract squared matrixelement contrib for L=+1 Ly=-1/2
	      m2 -= norm(me->calculateM(2,Lp,0));
	    }
	  // if the beam helicity is -1, there's an additional minus-sign
	  m2 *= polarization.doubpol.beamhel;
	}
    } // double polarization
  
  //me->print();

  return m2;

}


/**********************************************************************************/

/*!
 * Calculates the virtual photon cross-sections 
 * dsigma_{..}/dOmega_K* for the electroproduction,
 * in units of mubarn.
 *
 * When calc_response_func() is called in electrodiffcross(),
 * all variables (w, k, ...) etc should be given in MeV.
 * 
 * We contruct a MatrixElement object me.
 * Using this object we can calculate the value of
 * H_{lambda_e, lambda_ec} for specific kinematics.
 * 
 * IMPORTANT note on symmetry of the hadronic tensor under 
 * lambda_c <-> lambda:
 * 
 * Since H_{l,l'} = Re H_{l,l'} + i Im H_{l,l'} 
 * = sum_{l1,l2} M_l^{l1,l2} (M_l'^{l1,l2})^+ 
 * (82 stijn notes), we have
 *
 * H_{l',l} = (H_{l,l'})* = Re H_{l,l'} - i Im H_{l,l'}.
 *
 * Response functions R should be REAL:
 *       -> R, c^R = real lin. comb. of H_{l,l'},
 *          containing both H_{l,l'} and H_{l',l}, with a relative + sign 
 *          => Im parts cancel, real parts are equal
 *          so R, c^R = real lin. comb. of Re H_{l,l'} 
 *       -> s^R = imaginary lin. comb. of H_{l,l'},
 *          containing both H_{l,l'} and H_{l',l}, with a relative - sign
 *          => Re parts cancel, Im parts are equal
 *          so s^R = i. lin. comb. of i Im H_{l,l'}
 * 
 * SO: in R, c^R: H_{l,l'} + H_{l',l} = 2 Re H_{l,l'}
 *        s^R: H_{l,l'} - H_{l',l} = 2 (i Im H_{l,l'})
 *
 * Other consequence: H_{lambda,lambda} is real.
 *
 * The above symmetry properties reduces the number of matrix elements
 * we need to calculate.
 *
 * In 3 cases (s_TT, s_TL, s_TL_pol) the response function is of the form
 * i * Im(...). The response functions need to be real and in (136) we notice
 * an additional i in the response function definitions. Therefor we add
 * an additional factor (-1.0) when we store those 3 response functions
 * in the Response_func structure.
 *
 * \param  respcalc tells calc_response_func() which response functions need to be calculated.
 */
int calc_response_func(Response_func* resp, Response_func_calc* respcalc,
		       double w, double k, double costheta_k,
		       double pk, double mandel_s, double mp,
		       Class particles[], Observable* observ, int label)
{
  double chiphotofac, equivphotonenergy;
  double convertcoeff;
  double invfourpisquared = 1./(16.*PI*PI);


  /* Construct the matrixelement
   * ***************************/
  TMatrixElement* me = TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles, observ, label);

  // Structure containing all info about polarizations of the process
  Polar polarization = observ->photo.pol;


  /*-----------------------------------
   * Construction of H_{ll'} lin comb's    (dimensionless)
   *-----------------------------------*/

  // Initialize all response functions to zero
  resp->L = 0;
  resp->T = 0;
  resp->c_TT = 0;
  resp->c_TL = 0;
  resp->s_TT = 0;
  resp->s_TL = 0;
  resp->TT_pol = 0;
  resp->s_TL_pol = 0;
  resp->c_TL_pol = 0;


  /* We need to sum over all nucleon and hyperon polarizations
   * Two situations occur for each type of polarization:
   * - when the nucleon/hyperon is not polarized the different
   *   polarization contributions are summed
   * - In case of polarization we take the difference
   *
   * Before starting the polarization loops we determine what
   * type of combinations (summations/subtractions) of hadronic
   * matrix elements we need and store this in the hadron_comp[4]
   * array.
   * Then it suffices to multiply the H_{ll'} lin comb's within
   * the loops with the corresponding array-component.
   */

  /* Define the hadron_comp[] array
   * 4 components are enough because there
   * are 2x2 nucleon/hyperon polarizations 
   *
   * To access the correct element: hadron_comp[Lp*2+Ly]
   * Lp Ly -> index
   * -- --    -----
   *  0  0        0
   *  0  1        1
   *  1  0        2
   *  1  1        3
   */
  double hadron_comp[4]; 

  // In case of nucleon polarization
  if(observ->elec.bar_pol.tarpol)
    {
      hadron_comp[0] = -1;
      hadron_comp[1] = -1;
      hadron_comp[2] = +1;
      hadron_comp[3] = +1;
    }

  // In case of hyperon polarization
  else if(observ->elec.bar_pol.recpol)
    {
      hadron_comp[0] = -1;
      hadron_comp[1] = +1;
      hadron_comp[2] = -1;
      hadron_comp[3] = +1;
    }

  /* Simultanious target/recoil polarization
   * is not yet implemented
   * So only the unpolarized case is left */
  else
    {
      hadron_comp[0] = +1;
      hadron_comp[1] = +1;
      hadron_comp[2] = +1;
      hadron_comp[3] = +1;
    }

  // Initialize hadronic currents
  complex<double> h[3*3];
  for(int L=0; L<3; ++L)
    for(int LL=0; LL<3; ++LL)
      h[L*3+LL] = 0.0;
  
  // Calculate hadronic currents
  for( int Lp=0; Lp<=1; Lp++)   // Loop over nucleon polarization
    for( int Ly=0; Ly<=1; Ly++) // Loop over hyperon polarization
      for(int L=0; L<3; ++L) {
	for(int LL=0; LL<3; ++LL) {
	  h[L*3+LL] += (hadron_comp[Lp*2+Ly] * me->calculateM(L,Lp,Ly) * conj(me->calculateM(LL,Lp,Ly)));
	}
      }
  

  /* Construction of the "regular" response functions,
   * appearing in the unpolarized cross-section */
	  
  /* R_L ~ H_{0,0} (real) */
  if(respcalc->L)
    resp->L = real( h[1*3+1] );
	    	  
  /* R_T ~ H_{1,1} + H_{-1,-1} (both real) */
  if(respcalc->T)
    resp->T = real( h[2*3+2] + h[0*3+0] );
	    	  
  /* c^R_TT ~ H_{1,-1} + H_{-1,1} = 2 Re H_{1,-1} */
  if(respcalc->c_TT)
    resp->c_TT = 2.0 * real( h[2*3+0] );
	  
  /* c^R_TL ~ H_{0,1} + H_{1,0} - H_{-1,0} - H{0,-1} 
   * = 2 Re [H_{0,1} - H_{-1,0}] */
  if(respcalc->c_TL)
    resp->c_TL = 2.0 * 
      real ( h[1*3+2] - h[0*3+1] );
	  

  /* The response functions which are accessible
   * by polarized baryons are calculated. */
  if(observ->elec.baryon_pol)
    {
      /* s^R_{TT} ~ H_{1,-1} - H_{-1,1} 
       *          = 2 i Im H_{1,-1} 
       * -> extra factor (-1.0) because complex !*/
      if(respcalc->s_TT)
	resp->s_TT = -2.0 * imag( h[2*3+0] );
	      	      
      /* s^R{TL} ~ H_{1,0} - H_{0,1} + H_{-1,0} - H_{0,-1}
       *         = 2 i Im [H_{1,0} + H_{-1,0}] 
       * -> extra factor (-1.0) because complex ! */
      if(respcalc->s_TL)
	resp->s_TL = -2.0 * imag( h[2*3+1] + h[0*3+1] );
    }
	  	  
	  	  
  /* The response functions which are accessible
   * by polarized electrons are calculated. */
  if(observ->elec.elec_pol)
    {
      /* R_TT' ~ H_{1,1} - H{-1,-1} (both real) */
      if(respcalc->TT_pol)
	resp->TT_pol = real( h[2*3+2] - h[0*3+0] );
	      	      
      /* s^R_TL' ~ H_{1,0} - H_{0,1} - H_{-1,0} + H_{0,-1}
       * = 2 i Im [H_{1,0} - H_{-1,0}] 
       * -> extra factor (-1.0) because complex ! */
      if(respcalc->s_TL_pol)
	resp->s_TL_pol = -2.0 * imag( h[2*3+1] - h[0*3+1] );
	      	      
	      
      /* The response functions which are accessible by polarized
       * electrons and polarized baryons are calculated. */
      if(observ->elec.baryon_pol)
	{
	  /* c^R_TL' ~ H_{1,0} + H_{0,1} + H_{-1,0} + H_{0,-1}
	   * = 2 Re [H_{1,0} + H_{-1,0}] */     
	  if(respcalc->c_TL_pol)
	    resp->c_TL_pol = 2.0 * real( h[2*3+1] + h[0*3+1] );
	}      
    }
  
  /*-------------------------------------- 
   * Some useful variables are constructed 
   *--------------------------------------*/

  /* K_H (stijn notes 96) (MeV) */  
  equivphotonenergy = (mandel_s - mp*mp) / (2 * mp);

  /* chi (stijn notes 95), but with 1/4 instead of 1/16 (so not including
   * spinaverage) (MeV^{-2}) */
  chiphotofac = pk / (4. * sqrt(mandel_s) * mp * equivphotonenergy);


  /*---------------------
   * Conversion to mubarn
   *---------------------*/

  /* Factor for conversion from MeV^-2 
   * [ from H_{ll'}(dimensionless).\chi(MeV^-2).cste ] to mu barn (10^-34 m^2) 
   * == (hbar*c)^2 * 10^4 */
  
  convertcoeff = 197.32 * 197.32 * 1e4;  

   

  /*------------------------------------------
   * Virtual photon cross-sec dsigma/dOmega_K*
   *------------------------------------------*/

  /*
   * The virtual photon cross sections d(sigma)/d(omega_K*) are obtained 
   * after multiplying with the correct front factors: stijn notes (136).
   * 
   */

  resp->L *= 2 * chiphotofac * convertcoeff *invfourpisquared;
  resp->T *= chiphotofac * convertcoeff *invfourpisquared; 
  resp->c_TT *= - chiphotofac * convertcoeff*invfourpisquared;
  resp->c_TL *= - chiphotofac *  convertcoeff*invfourpisquared;
  

  if(observ->elec.baryon_pol)
    {
      resp->s_TT *= - chiphotofac * convertcoeff*invfourpisquared;
      resp->s_TL *= - chiphotofac * convertcoeff*invfourpisquared;
    }

  if(observ->elec.elec_pol)
    {
      resp->TT_pol *= chiphotofac * convertcoeff*invfourpisquared;
      resp->s_TL_pol *= - chiphotofac * convertcoeff*invfourpisquared;
      
      if(observ->elec.baryon_pol)
	{
	  resp->c_TL_pol *= - chiphotofac * convertcoeff*invfourpisquared;
        }
    }
  
  return 0;
}


/**********************************************************************************/


/*!
 * This function initializes an object of the MatrixElement class given the
 * kinematics and the specifications (particles[] and observ) it is given
 * as input parameters.
 * A pointer to the new object is returned.
 *
 * \warning since the memory for the object is dynamically allocated, one
 *          has to make sure the memory is released when no longer needed.
 *          This can be done using release_matrixelement() defined infra.
 */
TMatrixElement* construct_matrixelement(double w, double k,double costheta_k, double pk, 
					Class particles[], Observable* observ, int label)
{
  return TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles, observ, label);
}

/**********************************************************************************/

/*! Wrapper function: gives user in C-code access to the transversity amplitudes.
 */
void get_transversity_amplitude(TMatrixElement* matrixelement, double* amplitude, int index)
{
  complex<double> a;

  switch(index)
    {
    case 1:
      a = matrixelement->calculateM(0,1,1);
      break;
    case 2:
      a = matrixelement->calculateM(0,0,0);
      break;
    case 3:
      a = matrixelement->calculateM(2,0,1);
      break;
    case 4:
      a = matrixelement->calculateM(2,1,0);
      break;
    default:
      cerr << "ERROR: invalid transversity amplitude index (1, 2, 3, or 4)\n";
      exit(1);
    }
  
  amplitude[0] = real(a);
  amplitude[1] = imag(a);
}


/*! Wrapper function: gives user in C-code access to the helicity amplitudes.
 */
void get_helicity_amplitude(TMatrixElement* matrixelement, double* amplitude, int index)
{
  complex<double> a;

  switch(index)
    {
    case 1:
      a = matrixelement->calculateM(2,0,1);
      break;
    case 2:
      a = matrixelement->calculateM(2,1,1);
      break;
    case 3:
      a = matrixelement->calculateM(2,0,0);
      break;
    case 4:
      a = matrixelement->calculateM(2,1,0);
      break;
    default:
      cerr << "ERROR: invalid helicity amplitude index (1, 2, 3, or 4)\n";
      exit(1);
    }
  
  amplitude[0] = real(a);
  amplitude[1] = imag(a);
}


/**********************************************************************************/


/*!
 * The memory that's allocated by construct_matrixelement() should be released
 * at some point in time. This function takes care of that.
 */
int release_matrixelement( TMatrixElement* matrixelement)
{

  // Release dynamically allocated memory
  //delete matrixelement;

  return 0;
}


/**********************************************************************************/


/*! Wrapper function:
 * Makes it possible to access member function updateSpinDependencies() of the
 * MatrixElement class
 * Note that when the spin dependencies (== epsilon, among others) are updated, 
 * the data (cache) index should be incremented!
 */
int updateSpinDependencies_matrixelement( TMatrixElement* matrixelement, int* label)
{
  if (*label>=0)(*label)++; 
  matrixelement->updateSpinDependencies();
  matrixelement->SetDataLabel(*label);
  return 0;
}


/*!
 * Determine the squared matrix element for photopolarization
 * given the specific kinematics and all info about the process
 * and the observables.
 */
void set_datasize(int size)
{
TCalculateCGLNMatrix::SetDataSize(size);	 
}
