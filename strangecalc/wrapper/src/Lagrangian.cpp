/*!
 * \file Lagrangian.cpp
 * \ingroup wrapper
 *
 * Defines all types of Lagrangians and propagators
 * necessary to construct the amputated current.
 *
 * We define a number of functions (calculateDiagram*()) that put
 * together the necessary lagrangians and propagators to determine 
 * the contribution of a specific diagram to the amputated current.
 *
 * To improve speed these calculateDiagram*() functions are wrapped
 * in a TCalculateDiagram class that can cache the results if
 * memoizing is turned on (on a per-diagram basis, hard-coded).
 *
 * All references are to Stijn's notes
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 *
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

#include <complex>
#include <iostream>
#include <cstring>
#include "BetaIncomplete.h"
#include "Matrix.h"
#include "Tensor.h"
#include "FourVectorHandle.h"
#include "Lagrangian.h"
#include <numtoa.h>
#include <string.h>
#include <map>
#include <stdarg.h>
#include <TCachemap.h>
#include "FormFactorParametrization.h"
using std::cout; using std::cerr; using std::endl;
using std::complex;

#define pi 4.*atan(1.) //!< Local macro definition of pi

/*!
 * \return the natural logarithm of the gamma function 
 * from Numerical recipes in C 
 */
float gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);

}

/*! \return the gamma function */
double gamma_function(double z)
{
  const static double gamma_function_pi = 3.141592654;
 
  if(z > 0.)
    return exp(gammln(z));
  
  if(z < 0.)
    return gamma_function_pi / ( sin(gamma_function_pi * z) * exp(gammln(1. - z)) );   
  
  return 1.0;
}


/* ************** *
 * EM Lagrangians *
 * ************** */

/*! Photon (with 4vector k4vect) couples to a 1/2^+ particle (particle)
 * Lagrangian = (160) on p.29
 */
FourVector<GammaStructure> EMvertex_12p_12p(const Properties& particle,
					    const FourVector<double>& k4vect)
{
  double e; // electrical charge divided by |e_electron|
  if(particle.formfactorE==NULL)
    e = particle.E;
  else 
    e = (*particle.formfactorE).value(particle.E,-1.0*k4vect*k4vect);


  double kappa; // magnetic moment
  if(particle.formfactorG==NULL)
    kappa = particle.G;
  else
    kappa = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);


  return complex<double>(0.0,-1.0) * ELEC * 
    ( e * GMU
      + (-1.0*kappa)/(2.0*M_P) * ( GMU * (GMU*k4vect) 
				   - GammaStructure(1.0)*k4vect));
}

/* ******************************************************************************* */

/*! Photon (with 4vector k4vect) couples to a 1/2^+ particle 
 * and becomes a 1/2^- particle (particle)
 * Lagrangian = (168) on p.30
 */
FourVector<GammaStructure> EMvertex_12p_12m(const Properties& particle,
					    const FourVector<double>& k4vect)
{
  double e; // electrical charge divided by |e_electron|
  if(particle.formfactorE==NULL)
    e = particle.E;
  else 
    e = (*particle.formfactorE).value(particle.E,-1.0*k4vect*k4vect);


  double kappa; // magnetic moment
  if(particle.formfactorG==NULL)
    kappa = particle.G;
  else
    kappa = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);


  return complex<double>(0.0,-1.0) * ELEC * 
    ( e * GMU
      + (-1.0*kappa)/(2.0*M_P) * ( GMU * (GMU*k4vect)  
        - GammaStructure(1.0)*k4vect)) 
    * GammaStructure(0.0,1.0);
}

/* ******************************************************************************* */

/*! Photon (with 4vector k4vect) couples to a 0^- particle (particle) with
 * 4vector pK4vect.
 * Lagrangian = (163) on p.29
 */
FourVector< complex<double> > EMvertex_0m_0m(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& pK4vect)
{
  double e; // electrical charge divided by |e_electron|
  if(particle.formfactorE==NULL)
    e = particle.E;
  else 
    e = (*particle.formfactorE).value(particle.E,-1.0*k4vect*k4vect);


  FourVector<double> realVertex = ELEC * e * (2.0 * pK4vect - k4vect);

  return FourVector< complex<double> >(complex<double>(0.0,-1.0)*realVertex[0],
				       complex<double>(0.0,-1.0)*realVertex[1],
				       complex<double>(0.0,-1.0)*realVertex[2],
				       complex<double>(0.0,-1.0)*realVertex[3]);
}

/* ******************************************************************************* */

/*! Photon (with 4vector k) couples to a 1^- particle (particle) with
 * energy-momentum 4vector q.
 * The outgoing particle is a 0^-.
 */
TensorRank2 EMvertex_1m_0m(const Properties& particle,
			   const FourVector<double>& k,
			   const FourVector<double>& q)
{
  double g; // EM coupling constant
  if(particle.formfactorG==NULL)
    g = particle.G;
  else
    g = (*particle.formfactorG).value(particle.G,-1.0*k*k);


  double m = 1e3; // Normalization = 1 GeV by convention

  // Lagrangian = (165) on p.30
  // vertex = (-ieg/M)*epsilon^(abcd)*k4vect_b*q_d
  TensorRank2 propagator(        0.0         , k[3]*q[2]-k[2]*q[3],
			 k[1]*q[3]-k[3]*q[1] , k[2]*q[1]-k[1]*q[2],
			 k[2]*q[3]-k[3]*q[2] ,         0.0,
			 k[0]*q[3]-k[3]*q[0] , k[2]*q[0]-k[0]*q[2],
			 k[3]*q[1]-k[1]*q[3] , k[3]*q[0]-k[0]*q[3],
				 0.0         , k[0]*q[1]-k[1]*q[0],
			 k[1]*q[2]-k[2]*q[1] , k[0]*q[2]-k[2]*q[0],
			 k[1]*q[0]-k[0]*q[1] ,         0.0        );
  
  return (complex<double>(0.0,-1.0) * ELEC * g / m ) * propagator;
}

/* ******************************************************************************* */

/*! Photon (with 4vector k) couples to a 1^+ particle (particle) with
 * energy-momentum 4vector q.
 * The outgoing particle is a 0^-.
 */
TensorRank2 EMvertex_1p_0m(const Properties& particle,
			   const FourVector<double>& k,
			   const FourVector<double>& q)
{
  double g; // EM coupling constant
  if(particle.formfactorG==NULL)
    g = particle.G;
  else
    g = (*particle.formfactorG).value(particle.G,-1.0*k*k);


  double m = 1.0e3; // Normalization = 1 GeV by convention

  TensorRank2 metric;
  TensorRank2 qk(q,k); // direct product q^{mu}k_{nu}

  // Lagrangian = (166) on p.30
  // vertex = (ieg/M)*( (k.q) g^{mu,ksi}g_{ksi,nu} - q^{mu}k_{nu} )
  //   

  return complex<double>(ELEC*g/m) * ( (k*q)*metric*metric - qk );
}

/* ******************************************************************************* */

/*! Photon (with 4vector k4vect) excites a 1/2^+ particle (with 4vector p4vect)
 * to a 3/2^+ state (particle)
 *
 * When the argument 'timeReversed' is true (default is false), we consider the vertex
 * the other way around: a 3/2^+ particle decaying into a 1/2^+ by radiating a photon.
 *
 * A small remark about the 4vectors 'k4vect' and 'p4vect'. These are the momenta of
 * the photon and 1/2^+ particle respectively. It's the callers responsability to
 * make sure they point in the correct direction. Out of the vertex when 
 * timeReversed=true and vice versa.
 */
TensorRank2 EMvertex_12p_32p(const Properties& particle,
			     const FourVector<double>& k4vect,
			     const FourVector<double>& p4vect,
			     const FourVector<double>& q,
			     bool timeReversed)
{
  double kappa1; // first EM coupling constant
  if(particle.formfactorG==NULL)
    kappa1 = particle.G;
  else
    kappa1 = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);

  double kappa2; // second EM coupling constant
  if(particle.formfactorH==NULL)
    kappa2 = particle.H;
  else
    kappa2 = (*particle.formfactorH).value(particle.H,-1.0*k4vect*k4vect);


  const static TensorRank2 metric = 
    TensorRank2()*GammaStructure(1.0); // = the metric
  const static TensorRank2 GmuGmu = TensorRank2(GMU,GMU);

  if(!particle.gic){
    TensorRank2 offshellX = metric - (particle.X+0.5) * GmuGmu;

    // if(timeReversed) {
    //   // We have an extra minus sign because we need to change the direction
    //   // of the photon 4vector.
    //   return - ELEC / (2.0*M_P) * GammaStructure(0.0,1.0) *
    // 	( kappa1 * ( TensorRank2(GMU,k4vect)
    // 		     - (particle.Y+.5) * TensorRank2(GMU,(GMU*k4vect)*GMU)
    // 		     - metric * (GMU*k4vect)
    // 		     + (particle.Y+.5) * TensorRank2((GMU*k4vect)*GMU,GMU) )
    // 	  -
    // 	  kappa2/ (2.0*M_P) * ( TensorRank2(p4vect,k4vect) * offshellX
    // 				- (p4vect*k4vect) * offshellX )
    // 	  );
    // }
  
    // formula (171) and (172) on p.30
    return ELEC / (2.0*M_P) *
      ( kappa1 * ( TensorRank2(k4vect,GMU)
		   - (particle.Y+.5) * TensorRank2(GMU*(GMU*k4vect),GMU)
		   - metric * (GMU*k4vect)
		   + (particle.Y+.5) * TensorRank2(GMU,GMU*(GMU*k4vect)) )
	-
	kappa2/ (2.0*M_P) * ( offshellX * TensorRank2(k4vect,p4vect)
			      - offshellX * (p4vect*k4vect) )
	) *GammaStructure(0.0,1.0);
  }

  TensorRank2 T1 = (GMU*k4vect)*metric - TensorRank2(k4vect,GMU);
  TensorRank2 T2 = (p4vect*k4vect)*metric - TensorRank2(k4vect,p4vect)*GammaStructure(1.0); 

  return complex<double>(0.0,1.0) * ELEC / (4.0*M_P*M_P) * (GMU*q) * GammaStructure(0.0,1.0) * (kappa1 * T1 + kappa2 / (2.0*M_P) * T2);
}

/* ******************************************************************************* */

/*! Photon (with 4vector k4vect) excites a 1/2^+ particle (with 4vector p4vect)
 * to a 3/2^- state (particle)
 *
 * When the argument 'timeReversed' is true (default is false), we consider the vertex
 * the other way around: a 3/2^- particle decaying into a 1/2^+ by radiating a photon.
 *
 * A small remark about the 4vectors 'k4vect' and 'p4vect'. These are the momenta of
 * the photon and 1/2^+ particle respectively. It's the callers responsability to
 * make sure they point in the correct direction. Out of the vertex when 
 * timeReversed=true and vice versa.
 */
TensorRank2 EMvertex_12p_32m(const Properties& particle,
			     const FourVector<double>& k4vect,
			     const FourVector<double>& p4vect,
			     const FourVector<double>& q,
			     bool timeReversed)
{
  double kappa1; // first EM coupling constant
  if(particle.formfactorG==NULL)
    kappa1 = particle.G;
  else
    kappa1 = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);

  double kappa2; // second EM coupling constant
  if(particle.formfactorH==NULL)
    kappa2 = particle.H;
  else
    kappa2 = (*particle.formfactorH).value(particle.H,-1.0*k4vect*k4vect);


  const static TensorRank2 metric = 
    TensorRank2()*GammaStructure(1.0); // = the metric
  const static TensorRank2 GmuGmu = TensorRank2(GMU,GMU);

  if(!particle.gic){
    TensorRank2 offshellX = metric - (particle.X+0.5) * GmuGmu;
    
    // if(timeReversed) {
    //   // We have an extra minus sign because we need to change the direction
    //   // of the photon 4vector.
    //   return - ELEC / (2.0*M_P) *
    // 	( kappa1 * ( TensorRank2(GMU,k4vect)
    // 		     - (particle.Y+.5) * TensorRank2(GMU,(GMU*k4vect)*GMU)
    // 		     - metric * (GMU*k4vect)
    // 		     + (particle.Y+.5) * TensorRank2((GMU*k4vect)*GMU,GMU) )
    // 	  -
    // 	  kappa2/ (2.0*M_P) * ( TensorRank2(p4vect,k4vect) * offshellX
    // 				- (p4vect*k4vect) * offshellX )
    // 	  );
    // }
		// formula (171) and (172) on p.30
      return ELEC / 2.0 / M_P *
	( kappa1 * ( TensorRank2(k4vect,GMU)
		    - (particle.Y+.5) * TensorRank2(GMU*(GMU*k4vect),GMU)
		    - metric * (GMU*k4vect)
		    + (particle.Y+.5) * TensorRank2(GMU,GMU*(GMU*k4vect)) )
	  -
	  kappa2/ 2.0 / M_P * ( offshellX * TensorRank2(k4vect,p4vect)
				- offshellX * (p4vect*k4vect) )
	  );
      
    }
 
  TensorRank2 T1 = (GMU*k4vect)*metric - TensorRank2(k4vect,GMU);
  TensorRank2 T2 = (p4vect*k4vect)*metric - TensorRank2(k4vect,p4vect)*GammaStructure(1.0);

  return complex<double>(0.0,1.0) * ELEC / (4.0*M_P*M_P) 
    * (GMU*q) * (kappa1 * T1 + kappa2 / (2.0*M_P) * T2);
}

/* ******************************************************************************* */

/* Photon (with 4vector k4vect) excites a 1/2^+ particle (with 4vector p4vect)
 * to a 5/2^+ state (particle)
 */
TensorRank3 EMvertex_12p_52p(const Properties& particle,
			     const FourVector<double>& k4vect,
			     const FourVector<double>& p4vect,
			     const FourVector<double>& q)
{
  double kappa1; // first EM coupling constant
  if(particle.formfactorG==NULL)
    kappa1 = particle.G;
  else
    kappa1 = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);

  double kappa2; // second EM coupling constant
  if(particle.formfactorH==NULL)
    kappa2 = particle.H;
  else
    kappa2 = (*particle.formfactorH).value(particle.H,-1.0*k4vect*k4vect);


  const static TensorRank2 metric = 
    TensorRank2()*GammaStructure(1.0); // = the metric

  // always use gic
  TensorRank3 T1 = TensorRank3(p4vect,(GMU*k4vect)*metric - TensorRank2(k4vect,GMU));
  TensorRank3 T2 = TensorRank3(p4vect,(p4vect*k4vect)*metric - TensorRank2(k4vect,p4vect)*GammaStructure(1.0));

  return complex<double>(0.0,1.0) * ELEC / (16.0*M_P*M_P*M_P*M_P) * (q*q) * (kappa1 * T1 + kappa2 / (2.0*M_P) * T2);
}

/* ******************************************************************************* */

/*! Photon (with 4vector k4vect) excites a 1/2^+ particle (with 4vector p4vect)
 * to a 5/2^- state (particle)
 */
TensorRank3 EMvertex_12p_52m(const Properties& particle,
			     const FourVector<double>& k4vect,
			     const FourVector<double>& p4vect,
			     const FourVector<double>& q)
{
  double kappa1; // first EM coupling constant
  if(particle.formfactorG==NULL)
    kappa1 = particle.G;
  else
    kappa1 = (*particle.formfactorG).value(particle.G,-1.0*k4vect*k4vect);

  double kappa2; // second EM coupling constant
  if(particle.formfactorH==NULL)
    kappa2 = particle.H;
  else
    kappa2 = (*particle.formfactorH).value(particle.H,-1.0*k4vect*k4vect);


  const static TensorRank2 metric = 
    TensorRank2()*GammaStructure(1.0); // = the metric

  // always use gic for spin 5/2 particles
  TensorRank3 T1 = TensorRank3(p4vect,(GMU*k4vect)*metric - TensorRank2(k4vect,GMU));
  TensorRank3 T2 = TensorRank3(p4vect,(p4vect*k4vect)*metric - TensorRank2(k4vect,p4vect)*GammaStructure(1.0));

  return complex<double>(0.0,1.0) * ELEC / (16.0*M_P*M_P*M_P*M_P)*
    (q*q) * GammaStructure(0.0,1.0) * (kappa1 * T1 + kappa2 / (2.0*M_P) * T2);
}


/* ******************************************************************************* */
/* ******************************************************************************* */


/* ****************** *
 * Strong Lagrangians *
 * ****************** */

/*! Born s,t,u-channel:
 * A 1/2^+ particle(particle) decays into a 0^- and 1/2^+ particle
 * Lagrangian = (164) on p.29
 */
GammaStructure StrongVertex_12p_0m_12p(const Properties& particle,
				       double mandel)
{
  double g; // strong coupling constant divided by sqrt(4PI)
  if(particle.formfactorH==NULL)
    g = particle.H;
  else
    g = (*particle.formfactorH).value(particle.H,mandel);


  return g*GammaStructure(0.0,1.0);
}

/* ******************************************************************************* */

/*! A 1/2^- particle(particle) decays into a 0^- and 1/2^+ particle
 * Lagrangian = (169) on p.29 -> pseudoscalar coupling
 */
GammaStructure StrongVertex_12p_0m_12m(const Properties& particle,
				       double mandel)
{
  double g; // strong coupling constant divided by sqrt(4PI)
  if(particle.formfactorH==NULL)
    g = particle.H ;
  else
    g = (*particle.formfactorH).value(particle.H,mandel);


  return g*GammaStructure(1.0);
}

/* ******************************************************************************* */

/*! A 1^- particle(particle) with energy-momentum 4vector q4vect
 * decays into a 1/2^+ (with mass mN) and 1/2^+ particle (with mass mY)
 */
FourVector<GammaStructure> StrongVertex_12p_1m_12p(const Properties& particle,
						   const FourVector<double>& q4vect,
						   double mandel)
{
  double gv; // strong vector coupling constant
  if(particle.formfactorH==NULL)
    gv = particle.H ;
  else
    gv = (*particle.formfactorH).value(particle.H,mandel);

  double gt; // strong tensor coupling constant
  if(particle.formfactorI==NULL)
    gt = particle.I ;
  else
    gt = (*particle.formfactorI).value(particle.I,mandel);


  return ( gv * GMU +
	   gt/(M_P+M_L) * ( (GMU*q4vect)*GMU - q4vect*GammaStructure(1.0) ) )
    * complex<double>(0.0,1.0);
}

/* ******************************************************************************* */

/*! A 1^+ particle(particle) with energy-momentum 4vector q4vect
 * decays into a 1/2^+ (with mass mN) and 1/2^+ particle (with mass mY)
 */
FourVector<GammaStructure> StrongVertex_12p_1p_12p(const Properties& particle,
						   const FourVector<double>& q4vect,
						   double mandel)
{
  /* Lagrangian = (167) on p.30
   * it is the same as for coupling to a positive parity vector meson
   * except for the gamma^5 on the right */

  return StrongVertex_12p_1m_12p(particle,q4vect,mandel) * GammaStructure(0,1.0);
}

/* ******************************************************************************* */

/*! A 3/2^+ particle (particle) decays into a 1/2^+ particle and a
 * 0^- particle with energy-momentum 4vector pK4vect
 */
FourVector<GammaStructure> StrongVertex_32p_0m_12p(const Properties& particle,
						   const FourVector<double>& pK4vect,
						   const FourVector<double>& q,
						   double mandel)
{
  double f; /* strong coupling constant (divided by square root 4pi) */
  if(particle.formfactorI==NULL)
    f = particle.I ;
  else
    f = (*particle.formfactorI).value(particle.I,mandel) ;

  if(!particle.gic){
    // formula (173) p.30
  return complex<double>(0.0,-1.0) * f/M_KP * 
      ( pK4vect*GammaStructure(1.0) - (particle.Z+0.5) * (GMU*pK4vect) *GMU);
  }

  return complex<double>(0.0,1.0) * f / M_KP / M_KP * (GMU*q) * pK4vect;
}

/* ******************************************************************************* */

/*! A 3/2^- particle (particle) decays into a 1/2^+ particle and a
 * 0^- particle with energy-momentum 4vector pK4vect
 */
FourVector<GammaStructure> StrongVertex_32m_0m_12p(const Properties& particle,
						   const FourVector<double>& pK4vect,
						   const FourVector<double>& q,
						   double mandel)
{
  double f; /* strong coupling constant (divided by square root 4pi) */
  if(particle.formfactorI==NULL)
    f = particle.I ;
  else
    f = (*particle.formfactorI).value(particle.I,mandel) ;

  if(!particle.gic){
    // formula (173) p.30
  return complex<double>(0.0,1.0) * f/M_KP * GammaStructure(0.0,1.0) *
      ( pK4vect*GammaStructure(1.0) - (particle.Z+0.5) * (GMU*pK4vect) *GMU);
  }

  return complex<double>(0.0,1.0) * f / M_KP / M_KP * GammaStructure(0.0,1.0) * (GMU*q) * pK4vect;
}

/* ******************************************************************************* */

/*! A 5/2^+ particle (particle) decays into a 1/2^+ particle and a
 * 0^- particle with energy-momentum 4vector pK4vect
 */
TensorRank2 StrongVertex_52p_0m_12p(const Properties& particle,
				    const FourVector<double>& pK4vect,
				    const FourVector<double>& q,
				    double mandel)
{
  double f; /* strong coupling constant (divided by square root 4pi) */
  if(particle.formfactorI==NULL)
    f = particle.I ;
  else
    f = (*particle.formfactorI).value(particle.I,mandel) ;

  // always use gic couplings for spin 5/2  
  return complex<double>(0.0,1.0) * f / (M_KP*M_KP*M_KP*M_KP) * GammaStructure(0.0,1.0) * (q*q) * TensorRank2(pK4vect,pK4vect);
}


/* ******************************************************************************* */

/*! A 5/2^- particle (particle) decays into a 1/2^+ particle and a
 * 0^- particle with energy-momentum 4vector pK4vect
 */
TensorRank2 StrongVertex_52m_0m_12p(const Properties& particle,
				    const FourVector<double>& pK4vect,
				    const FourVector<double>& q,				    
	                            double mandel)
{
  double f; /* strong coupling constant (divided by square root 4pi) */
  if(particle.formfactorI==NULL)
    f = particle.I ;
  else
    f = (*particle.formfactorI).value(particle.I,mandel) ;
  
  // always use gic couplings for spin 5/2
  return complex<double>(0.0,1.0) * f / (M_KP*M_KP*M_KP*M_KP) * (q*q) * TensorRank2(pK4vect,pK4vect);

}

/* ******************************************************************************* */
/* ******************************************************************************* */


/* *********** *
 * PROPAGATORS *
 * *********** */

/*! \return the progagator of a spin-0 particle (particle) with
 * energy-momentum 4vector p4vect.
 */
complex<double> propagatorSpin0(const Properties& particle,
				const FourVector<double>& p4vect,
				const bool nonresonant)
{
  // formula is (193) of Stijns notes p.33
  complex<double> propagator = p4vect*p4vect - particle.mass * particle.mass
  + (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width));
  return (complex<double>(0.0,1.0)/propagator);
}

/* ******************************************************************************* */

/*! \return the progagator of a spin-1/2 particle (particle) with
 * energy-momentum 4vector p4vect.
 */
GammaStructure propagatorSpin12(const Properties& particle,
				const FourVector<double>& p4vect,
				const bool nonresonant)
{
  // formula is (190) of Stijns notes p.32
  GammaStructure propagator = GMU*p4vect + GammaStructure(1.0)*particle.mass;
  return propagator*(complex<double>(0.0,1.0)/
	(p4vect*p4vect - particle.mass * particle.mass
	+ (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width))));
}

/* ******************************************************************************* */

/*! \return the propagator of a spin-1 particle (particle) with
 * energy-momentum 4vect p4vect.
 * The return type is a TensorRank2, meaning it is the
 * direct product of two 4vectors
 */
TensorRank2 propagatorSpin1(const Properties& particle,
			    const FourVector<double>& p4vect,
			    const bool nonresonant)
{
  const static TensorRank2 metric;
  TensorRank2 pp(p4vect,p4vect); // direct product of p4vect with itself

  // formula is (194) on p.33
  return 
    ( metric - (1.0/(particle.mass*particle.mass)) * pp )
    * ( complex<double>(0.0,-1.0) / 
      ( p4vect*p4vect - particle.mass*particle.mass 
	+ (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width))) );
}

/* ******************************************************************************* */

/*! \return the propagator of a spin-3/2 particle (particle) with
 * energy-momentum 4vector q.
 * We use the Rarita-Schwinger propagator.
 */
TensorRank2 propagatorSpin32(const Properties& particle,
			     const FourVector<double>& q,
			     const bool nonresonant)
{
  // formula is (191) on p.32
  const static TensorRank2 metric = 
    TensorRank2() * GammaStructure(1.0); // the metric
  const static TensorRank2 GmuGmu(GMU,GMU); // gamma^{mu} gamma^{nu}
  TensorRank2 qq(q,q); // q^{mu} q^{nu}
  TensorRank2 gammaq(GMU,q); // gamma^{mu} q^{nu}
  TensorRank2 qgamma(q,GMU); // q^{mu} gamma^{nu}
  double q2 = q*q;

  if(!particle.gic)
  {
    return complex<double>(0.0,1.0) * (GMU*q + particle.mass*GammaStructure(1.0)) *
    ( 1.0 / 3.0 / ( q*q -particle.mass*particle.mass 
		    + (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width))  )
    )
    * ( 3.0 * metric
	- GmuGmu
	- 2.0/particle.mass/particle.mass * qq * GammaStructure(1.0) 
	- 1.0/particle.mass * ( gammaq - qgamma )
	);
  }
  //We will only use the spin-3/2 part of the spin-3/2 propagator
  //in describing consistent couplings

  return complex<double>(0.0,1.0) * (GMU*q + particle.mass*GammaStructure(1.0)) *
    ( 1.0 / 3.0 / ( q*q -particle.mass*particle.mass 
		       + (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width))  )  )
    * ( 3.0*metric - 2.0/q2*qq*GammaStructure(1.0) 
	- 1.0/q2*(GMU*q)*(gammaq - qgamma) - GmuGmu);
}

/* ******************************************************************************* */

/*! \return the propagator of a spin-5/2 particle (particle) with
 * energy-momentum 4vector q.
 * We use the 5/2-propagator.
 */
TensorRank4 propagatorSpin52(const Properties& particle,
			     const FourVector<double>& q,
			     const bool nonresonant)
{
  const static TensorRank2 metric = 
    TensorRank2() * GammaStructure(1.0); // the metric
  double q2 = q*q;
  const static FourVector<double> t = FourVector<double>(1.0,0.0,0.0,0.0); // t unit vector
  const static FourVector<double> x = FourVector<double>(0.0,1.0,0.0,0.0); // x unit vector
  const static FourVector<double> y = FourVector<double>(0.0,0.0,1.0,0.0); // y unit vector
  const static FourVector<double> z = FourVector<double>(0.0,0.0,0.0,1.0); // z unit vector
    
  TensorRank2 P = -1.0 * metric + 1.0 / q2 * TensorRank2(q,q) * GammaStructure(1.0); // spin-1 projection operator
  FourVector<GammaStructure> gP = GMU * P;

  TensorRank4 T;
  TensorRank2 rank2;


  //T
  rank2 = TensorRank2(t,x);
  rank2 += TensorRank2(x,t);
  T = TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(t,y);
  rank2 += TensorRank2(y,t);
  T += TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(t,z);
  rank2 += TensorRank2(z,t);
  T += TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(t,q);
  rank2 += TensorRank2(q,t);
  T += (1.0/q2) * TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(t,gP);
  rank2 += TensorRank2(gP,t);
  T += 1.0/5.0 * TensorRank4(rank2,rank2);
  
  T *= -1.0; 
  

  //X
  rank2 = TensorRank2(x,y);
  rank2 += TensorRank2(y,x);
  T += TensorRank4(rank2,rank2);

  rank2 = TensorRank2(x,z);
  rank2 += TensorRank2(z,x);
  T += TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(x,q);
  rank2 += TensorRank2(q,x);
  T += (1.0/q2) * TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(x,gP);
  rank2 += TensorRank2(gP,x);
  T += 1.0/5.0 * TensorRank4(rank2,rank2);
  
  
  //Y
  rank2 = TensorRank2(y,z);
  rank2 += TensorRank2(z,y);
  T += TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(y,q);
  rank2 += TensorRank2(q,y);
  T += (1.0/q2) * TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(y,gP);
  rank2 += TensorRank2(gP,y);
  T += 1.0/5.0 * TensorRank4(rank2,rank2);
  
  
  //Z
  rank2 = TensorRank2(z,q);
  rank2 += TensorRank2(q,z);
  T += (1.0/q2) * TensorRank4(rank2,rank2);
  
  rank2 = TensorRank2(z,gP);
  rank2 += TensorRank2(gP,z);
  T += 1.0/5.0 * TensorRank4(rank2,rank2);
  
  
  //Q
  rank2 = TensorRank2(q,gP);
  rank2 += TensorRank2(gP,q);
  T += 1.0/(5.0*q2) * TensorRank4(rank2,rank2);
  
  TensorRank4 S = TensorRank4(t,t,t,t) + TensorRank4(x,x,x,x) + TensorRank4(y,y,y,y)
    + TensorRank4(z,z,z,z) + 1.0/q2/q2 * TensorRank4(q,q,q,q);

  return complex<double>(0.0,1.0) * (particle.mass*GammaStructure(1.0) + GMU*q) *
    ( 1.0 / 10.0 / ( q2 - particle.mass*particle.mass 
                     + (nonresonant ? 0.: complex<double>(0.0,particle.mass*particle.width)) ) )
    * (10.0 * S + 5.0 * T - 2.0 * TensorRank4(P,P) );
}

/* ******************************************************************************* */

/*! \return the Regge propagator, where:
 * - particle is the first materialization of the exchanged trajectory
 * - s and u (MeV^2) are the well known mandelstam variables
 * - mandel = t or u (MeV^2) depending on which channel is being reggeized
 * 
 * The expression for the propagator is (409-419) on p.67-68
 */
complex<double> propagatorRegge(const Properties& particle, double s, double u,
				double mandel, double Q2, double a, const Observable* observ)
{
  /* First we need to select the correct slope and intercept
   * of the trajectory: alpha(t) = a0 + a1 * (t-mass^2)
   * and the spin of the exchanged particles.
   */
  double a0, a1, s0; // intercept, slope, scale factor

  a0 = particle.spin;
  s0 = 1.e6;

  // K(494) trajectory
  if(!strcmp(particle.nickname,"K+") || !strcmp(particle.nickname,"K0"))
    {
      a1 = 0.7e-6/(1. + a*Q2/s);
     }
  // K*(892) trajectory
  else if(!strcmp(particle.nickname,"B1") || !strcmp(particle.nickname,"B2"))
    {
      a1 = 0.85e-6;
    }
  // K1(1270) trajectory
  else if(!strcmp(particle.nickname,"K1") || !strcmp(particle.nickname,"K10"))
    {
      a1 = 0.7e-6;
    }
  // K1(1400) trajectory
  else if(!strcmp(particle.nickname,"C1") || !strcmp(particle.nickname,"C2"))
    {
      /* This trajectory is not on solid ground.
       * We assume the members are K1(1400), K2(1820) and K4(2500). 
       * The latter, however, is not well known, so we define a 
       * trajectory through the 1st two points.
       */
      a1 = 0.75e-6;
    }
  // K*(1410) trajectory
  else if(!strcmp(particle.nickname,"B3") || !strcmp(particle.nickname,"B4"))
    {
      /* This trajectory is not well known.
       * We take the first materialization at the PDG-mass
       * of the K*(1410).
       * The slope is calculated with the Bonn-model
       * see Phys. Rev. C 75, 045204 (2007) for more info.
       */
      // a1 = 0.85e-6;
      a1 = 0.83e-6;
    }

  // pi(140)/b1(1235) trajectory
  else if(!strcmp(particle.nickname,"P+") || !strcmp(particle.nickname,"P-"))
    {
      a1 = 0.74e-6/(1. + a*Q2/s);
      s0 = 1./a1;
    }
  // rho(770)/a2(1320) trajectory
  else if(!strcmp(particle.nickname,"R+") || !strcmp(particle.nickname,"R-"))
    {
      a1 = 0.85e-6;
      a0 = 0.53 + a1*particle.mass*particle.mass;
      s0 = 1./a1;
    }
  // a1(1260) trajectory
  else if(!strcmp(particle.nickname,"A1+") || !strcmp(particle.nickname,"A1-"))
    {
      a1 = 0.85e-6;
      a0 = 0.53 - 1. + a1*particle.mass*particle.mass;
      s0 = 1./a1;
    }
  else
    {
      cout << "Regge trajectory exchange for "
	   << particle.nickname
	   << " is not implemented!" << endl << endl;
      exit(1);
    }

  // calculate alpha(mandel)
  double alpha = a0 + a1*(mandel - particle.mass*particle.mass);

  /* There're 2 ways to deal with spin in the Regge propagator.
   * The user is given the option in the input.
   * Stijn:  all occurences of alpha(mandel) are replaced by
   *         alpha(mandel) - spin of the exchanged particle
   * Guidal: alpha(mandel) is only replaced by alpha(mandel)-spin
   *         in the gamma-function and in the power of s
   *
   * Tamara has studied the implications of the different conventions.
   * She concludes that:
   * - for non-degenerate trajectories Stijns recipe doesn't produce
   *   the correct poles. Meaning in that case the Guidal recipe is the
   *   only correct option. For the time being, non-degenerate trajectories
   *   are not implemented however.
   * - for degenerate trajectories with a rotating phase both recipes
   *   are identical
   * - for degenerate trajectories with constant phases, both recipes
   *   differ by a factor (-1.0). 
   */

  /* calculate alpha(mandel) - spin of the exchange particles
   * to be used in gamma function and as power for s  */
  double alphaMinSpin = alpha - particle.spin;

  /* When we use the Janssen recipe for spin, all occurences of alpha(mandel)
   * are replaced by alpha(mandel) - spin
   */
  if(!observ->reg.spinshift_guidal)
    alpha = alphaMinSpin;


  /* When the trajectory is degenerate, one can choose between
   * a constant or a rotating phase.
   * When the trajector is non-degenerate, it has either 
   * positive or negative signature. */
  complex<double> phase;
  short regge_phase = (observ->iso.isospin == observ->iso.iso_base ?
		       particle.regge_phase : particle.regge_phase_nonbase);

  if(regge_phase==0) // constant
    phase = 1.0;
  else if(regge_phase==1) // rotating
    phase = exp(complex<double>(0.0,-1.0*PI*alpha));
  else if(regge_phase==2) // non-degenerate
    phase = 0.5 + 0.5*exp(complex<double>(0.0,-1.0*PI*alphaMinSpin));
  else
    cerr << "Error in propagatorRegge(..): "
	 << "unknown phase for Regge trajectory" << endl;
  
  /* It is possible to change the asymptotic behaviour
   * of the Regge propagator */
  if(observ->reg.s_modif)
    s = (s-u)/2.0;


  /* Calculate and return the Regge propagator correction factor
   * = 1/Feynman propagator * Regge propagator
   */

  return ( (mandel - particle.mass*particle.mass) * pow(s/s0 , alpha) * 
  	   PI * a1 / sin(PI * alpha) * phase / gamma_function(1.0 + alpha) );
}


/* ******************************************************************************* */
/* ******************************************************************************* */


/* ******************************* *
 * Calculate Diagram contributions *
 * ******************************* */

/*! Calculates the contribution of a Born s-channel diagram with exchanged
 * particle (particle)
 */

FourVector<GammaStructure> calculateDiagramS(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return ( StrongVertex_12p_0m_12p(particle,(k4vect+p4vect)*(k4vect+p4vect)) *
	   propagatorSpin12(particle,k4vect+p4vect,true) *
	   EMvertex_12p_12p(particle,k4vect) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of a Born t-channel diagram with exchanged
 * particle (particle)
 */
FourVector<GammaStructure> calculateDiagramT(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  FourVector< complex<double> > EMvertex = EMvertex_0m_0m(particle,k4vect,pK4vect);
  
  FourVector<GammaStructure> EMvertex2(GammaStructure(1.0)*EMvertex[0],
				       GammaStructure(1.0)*EMvertex[1],
				       GammaStructure(1.0)*EMvertex[2],
				       GammaStructure(1.0)*EMvertex[3]);

  return ( EMvertex2 *
	   propagatorSpin0(particle,pK4vect+(-1.0*k4vect),true) * 
	   StrongVertex_12p_0m_12p(particle,(k4vect-pK4vect)*(k4vect-pK4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of a Born u-channel diagram with exchanged
 * particle (particle)
 */
FourVector<GammaStructure> calculateDiagramU(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return ( EMvertex_12p_12p(particle,k4vect) *
	   propagatorSpin12(particle,pY4vect+(-1.0*k4vect),true) *
	   StrongVertex_12p_0m_12p(particle,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of extended Born u-channel diagram with exchanged
 * particle (particle)
 */
FourVector<GammaStructure> calculateDiagramA(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  /* since the baryons are not the same, there's no electric term in the
   * EM vertex.
   * -> We check whether the relevant coupling constant is zero 
   * When no longitudinal coupling is desired the the
   * formfactor should be set to NULL */
  if(particle.E != 0.0)
    {
      cout << "Coupling constant E for " << particle.nickname
	   << " in mass.in should be set to zero." << endl << endl;
      exit(1);
    }
  if(particle.formfactorE != NULL && !particle.long_coupling)
    {
      cout << "Form factor electric term  for " << particle.nickname
	   << "-exchange without longitudinal coupling is non-zero." << endl << endl;
      exit(1);
    }


  return ( EMvertex_12p_12p(particle,k4vect) *
	   propagatorSpin12(particle,pY4vect+(-1.0*k4vect),true) *
	   StrongVertex_12p_0m_12p(particle,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of vector meson(negative parity) exchange (B-diagram)
 */
FourVector<GammaStructure> calculateDiagramB(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{

  FourVector<double> q4vect = k4vect - pK4vect; // 4momentum of the exchanged particle
  return (EMvertex_1m_0m(particle,k4vect,q4vect) * 
	  propagatorSpin1(particle,q4vect,true) ) * 
	  StrongVertex_12p_1m_12p(particle,q4vect,(k4vect-pK4vect)*(k4vect-pK4vect));
	
}

/* ******************************************************************************* */

/*! Calculates the contribution of vector meson(positive parity) exchange (C-diagram)
 */
FourVector<GammaStructure> calculateDiagramC(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  FourVector<double> q4vect = k4vect - pK4vect; // 4momentum of the exchanged particle

  return ( ( EMvertex_1p_0m(particle,k4vect,q4vect) *
	     propagatorSpin1(particle,q4vect,true)) *
	   StrongVertex_12p_1p_12p(particle,q4vect,
				   (k4vect-pK4vect)*(k4vect-pK4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of nucleon resonance (positive parity) exchange
 * (D-diagram)
 */
FourVector<GammaStructure> calculateDiagramD(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  /* for resonance exchange, there's no electric term in the
   * EM vertex.
   * -> We check whether the relevant coupling constant is zero 
   * When no longitudinal coupling is desired the the
   * formfactor should be set to NULL */
  if(particle.E != 0.0)
    {
      cout << "Coupling constant E for " << particle.nickname
	   << " in mass.in should be set to zero." << endl << endl;
      exit(1);
    }
  if(particle.formfactorE != NULL && !particle.long_coupling)
    {
      cout << "Form factor electric term  for " << particle.nickname
	   << "-exchange without longitudinal coupling is non-zero." << endl << endl;
      exit(1);
    }

  return ( StrongVertex_12p_0m_12p(particle,(k4vect+p4vect)*(k4vect+p4vect)) *
	    // In the case of kaon capture, N* is non-resonant
	   propagatorSpin12(particle,k4vect+p4vect,kaoncapture) * 
	   EMvertex_12p_12p(particle,k4vect) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of nucleon resonance (negative parity) exchange
 * (E-diagram)
 */
FourVector<GammaStructure> calculateDiagramE(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  /* for resonance exchange, there's no electric term in the
   * EM vertex.
   * -> We check whether the relevant coupling constant is zero 
   * When no longitudinal coupling is desired the the
   * formfactor should be set to NULL */
  if(particle.E != 0.0)
    {
      cout << "Coupling constant E for " << particle.nickname
	   << " in mass.in should be set to zero." << endl << endl;
      exit(1);
    }
  if(particle.formfactorE != NULL && !particle.long_coupling)
    {
      cout << "Form factor electric term  for " << particle.nickname
	   << "-exchange without longitudinal coupling is non-zero." << endl << endl;
      exit(1);
    }

  return ( StrongVertex_12p_0m_12m(particle,(k4vect+p4vect)*(k4vect+p4vect)) *
	   // In the case of kaon capture, N* is non-resonant
	   propagatorSpin12(particle,k4vect+p4vect,kaoncapture) *
	   EMvertex_12p_12m(particle,k4vect) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of hyperon resonance (positive parity) exchange
 * (F-diagram)
 */
FourVector<GammaStructure> calculateDiagramF(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  /* since the baryons are not the same, there's no electric term in the
   * EM vertex.
   * -> We check whether the relevant coupling constant is zero 
   * When no longitudinal coupling is desired the the
   * formfactor should be set to NULL */
  if(particle.E != 0.0)
    {
      cout << "Coupling constant E for " << particle.nickname
	   << " in mass.in should be set to zero." << endl << endl;
      exit(1);
    }
  if(particle.formfactorE != NULL && !particle.long_coupling)
    {
      cout << "Form factor electric term  for " << particle.nickname
	   << "-exchange without longitudinal coupling is non-zero." << endl << endl;
      exit(1);
    }
  // In the case of kaon capture, Y* is resonant, 
  // so kaoncapture=true =>  nonresonant=false => nonresonant = !kaoncapture
  return (GMU[0]*((Fourvector)calculateDiagramD(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];

  // return ( EMvertex_12p_12p(particle,k4vect) *
  // 	  // In the case of kaon capture, Y* is resonant
  // 	   propagatorSpin12(particle,pY4vect+(-1.0*k4vect),!kaoncapture) *
  // 	   StrongVertex_12p_0m_12p(particle,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of hyperon resonance (negative parity) exchange
 * (G-diagram)
 */
FourVector<GammaStructure> calculateDiagramG(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  /* since the baryons are not the same, there's no electric term in the
   * EM vertex.
   * -> We check whether the relevant coupling constant is zero 
   * When no longitudinal coupling is desired the the
   * formfactor should be set to NULL */
  if(particle.E != 0.0)
    {
      cout << "Coupling constant E for " << particle.nickname
	   << " in mass.in should be set to zero." << endl << endl;
      exit(1);
    }
  if(particle.formfactorE != NULL && !particle.long_coupling)
    {
      cout << "Form factor electric term  for " << particle.nickname
	   << "-exchange without longitudinal coupling is non-zero." << endl << endl;
      exit(1);
    }
  // In the case of kaon capture, Y* is resonant, 
  // so kaoncapture=true =>  nonresonant=false => nonresonant = !kaoncapture
  return (GMU[0]*((Fourvector)calculateDiagramE(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];

  // return ( EMvertex_12p_12m(particle,k4vect) *
  // 	   // In the case of kaon capture, Y* is resonant
  // 	   propagatorSpin12(particle,pY4vect+(-1.0*k4vect),!kaoncapture) *
  // 	   StrongVertex_12p_0m_12m(particle,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-3/2 nucleon resonance (positive parity)
 * exchange (H-diagram)
 */
FourVector<GammaStructure> calculateDiagramH(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return ( StrongVertex_32p_0m_12p(particle,pK4vect,k4vect+p4vect,(k4vect+p4vect)*(k4vect+p4vect)) *
	   // In the case of kaon capture, N* is non-resonant
	   propagatorSpin32(particle,k4vect+p4vect,kaoncapture) * 
	   EMvertex_12p_32p(particle,k4vect,p4vect,k4vect+p4vect) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-3/2 nucleon resonance (negative parity)
 * exchange (I-diagram)
 */
FourVector<GammaStructure> calculateDiagramI(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return ( StrongVertex_32m_0m_12p(particle,pK4vect,k4vect+p4vect,(k4vect+p4vect)*(k4vect+p4vect)) *
	   // In the case of kaon capture, N* is non-resonant
	   propagatorSpin32(particle,k4vect+p4vect,kaoncapture) * 
	   EMvertex_12p_32m(particle,k4vect,p4vect,k4vect+p4vect) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-3/2 hyperon resonance (positive parity)
 * exchange (J-diagram)
 */
FourVector<GammaStructure> calculateDiagramJ(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  // In the case of kaon capture, Y* is resonant, 
  // so kaoncapture=true =>  nonresonant=false => nonresonant = !kaoncapture
  return (GMU[0]*((Fourvector)calculateDiagramH(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];

  // return ( EMvertex_12p_32p(particle,-k4vect,pY4vect,k4vect-pY4vect,true) *
  // 	   // In the case of kaon capture, Y* is resonant
  // 	   propagatorSpin32(particle,pY4vect-k4vect,!kaoncapture) *
  // 	   StrongVertex_32p_0m_12p(particle,pK4vect,k4vect-pY4vect,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-3/2 hyperon resonance (negative parity)
 * exchange (L-diagram)
 */
FourVector<GammaStructure> calculateDiagramL(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  // In the case of kaon capture, Y* is resonant, 
  // so kaoncapture=true =>  nonresonant=false => nonresonant = !kaoncapture
  return (GMU[0]*((Fourvector)calculateDiagramI(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];

  // return ( EMvertex_12p_32m(particle,-k4vect,pY4vect,k4vect-pY4vect,true) *
  // 	   // In the case of kaon capture, Y* is resonant
  // 	   propagatorSpin32(particle,pY4vect-k4vect,!kaoncapture) *
  // 	   StrongVertex_32m_0m_12p(particle,pK4vect,k4vect-pY4vect,(k4vect-pY4vect)*(k4vect-pY4vect)) );
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-5/2 nucleon resonance (positive parity)
 * exchange (M-diagram)
 */
// Memoized wrapper function
FourVector<GammaStructure> calculateDiagramM(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  TensorRank2 rank2 = StrongVertex_52p_0m_12p(particle,pK4vect,k4vect+p4vect,(k4vect+p4vect)*(k4vect+p4vect)) % propagatorSpin52(particle,k4vect+p4vect,kaoncapture);

  return rank2 % EMvertex_12p_52p(particle,k4vect,p4vect,k4vect+p4vect);
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-5/2 nucleon resonance (negative parity)
 * exchange (N-diagram)
 */
// Memoized wrapper function
FourVector<GammaStructure> calculateDiagramN(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  TensorRank2 rank2 = StrongVertex_52m_0m_12p(particle,pK4vect,k4vect+p4vect,(k4vect+p4vect)*(k4vect+p4vect)) % propagatorSpin52(particle,k4vect+p4vect,kaoncapture);

  return rank2 % EMvertex_12p_52m(particle,k4vect,p4vect,k4vect+p4vect);
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-5/2 hyperon resonance (positive parity)
 * exchange (O-diagram)
 */
// Memoized wrapper function
FourVector<GammaStructure> calculateDiagramO(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return (GMU[0]*((Fourvector)calculateDiagramM(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];
}

/* ******************************************************************************* */

/*! Calculates the contribution of spin-5/2 hyperon resonance (negative parity)
 * exchange (Q-diagram)
 */
// Memoized wrapper function
FourVector<GammaStructure> calculateDiagramQ(const Properties& particle,
					     const FourVector<double>& k4vect,
					     const FourVector<double>& p4vect,
					     const FourVector<double>& pK4vect,
					     const FourVector<double>& pY4vect,
					     const bool kaoncapture)
{
  return (GMU[0]*((Fourvector)calculateDiagramN(particle,-k4vect,pY4vect,-pK4vect,p4vect,!kaoncapture)).HermitianConjugate())*GMU[0];
}
