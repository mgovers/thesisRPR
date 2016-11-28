/*!
 * \file FormFactorParametrization.cpp
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
#include <cmath>
#include <cstdlib>
#include <Structures.h>
#include "FormFactorParametrization.h"
#include <complex>

using std::complex;

using namespace std;


/*! MONOPOLE formfactor
 * has one parameter = cutoff
 */
double monopole(double Q2,double cutoff,double p2,double p3,
		double p4,double p5,double p6,double p7,
		double p8,double p9,double p10,double p11,
		double p12,double p13,double p14,double p15)
{
  return 1.0/ ( 1.0 + Q2/cutoff/cutoff );
}


/*! MONOPOLE formfactor which is zero at Q^2=0
 * has two parameters = cutoff and cutoff2
 */
double monopole0(double Q2,double cutoff,double cutoff2,double p3,
		double p4,double p5,double p6,double p7,
		double p8,double p9,double p10,double p11,
		double p12,double p13,double p14,double p15)
{
  return monopole(Q2,cutoff) * Q2 / ( Q2 + cutoff2*cutoff2 );
}

/************************************************************************/

/*! DIPOLE formfactor
 * has one parameter = cutoff
 */
double dipole(double Q2,double cutoff,double p2,double p3,
	      double p4,double p5,double p6,double p7,
	      double p8,double p9,double p10,double p11,
	      double p12,double p13,double p14,double p15)
{
  return 1.0/ pow( 1.0 + Q2/cutoff/cutoff, 2 );
}


/*! DIPOLE formfactor which is zero at Q^2=0
 * has two parameters = cutoff and cutoff2
 */
double dipole0(double Q2,double cutoff,double cutoff2,double p3,
	       double p4,double p5,double p6,double p7,
	       double p8,double p9,double p10,double p11,
	       double p12,double p13,double p14,double p15)
{
  return dipole(Q2,cutoff) * Q2 / ( Q2 + cutoff2*cutoff2 );
}

/************************************************************************/

/*! CONSTANT formfactor
 * has one parameter = constant
 * This parametrization is used when the vertex does not have a formfactor
 */
double constant(double Q2,double constant,double p2,double p3,
	      double p4,double p5,double p6,double p7,
	      double p8,double p9,double p10,double p11,
	      double p12,double p13,double p14,double p15)
{
  return constant;
}


/*! Function returning message: "Not defined"
 */

double not_defined(double Q2,double p1,double p2,double p3,
		   double p4,double p5,double p6,double p7,
		   double p8,double p9,double p10,double p11,
		   double p12,double p13,double p14,double p15)
{
  cout << "Chosen formfactor parametrization does not exist for coupling constant" << endl;
  exit(1);

  return 1.0;
}

/************************************************************************/

/*! proton Dirac formfactor parameterized by Gari and Kruempelmann
 * Phys.Lett.B,173:10 and Phys.Rev.D,45:1817
 *
 * All formula references are to p.49-50 of Stijn's notes
 *
 * 13 parameters
 */
double gariKruempelmann_proton_dirac(double Q2,double k_iv,double k_is,double grfr,
				     double kr,double gwfw,double kw, double gfff,
				     double kf,double muf,double l1rw,double l1D,
				     double l2,double lQCD,double p14,double p15)
{
  // Q~^2, (316)
  double Qt2 = Q2*log( (l2*l2+Q2)/lQCD/lQCD ) / log( l2*l2/lQCD/lQCD );

  // (312)
  double f1 = l1rw*l1rw / (l1rw*l1rw+Qt2) * l2*l2 / (l2*l2+Qt2);

  // (312)
  double f1D = l1D*l1D / (l1D*l1D+Qt2) * l2*l2 / (l2*l2+Qt2);

  // (314)
  double f1Phi = f1 * pow(Q2/(l1rw*l1rw+Q2),1.5);


  // (310)
  double f1_is = 
    gwfw * M_OMEGA*M_OMEGA / (M_OMEGA*M_OMEGA+Q2) * f1 +
    gfff * M_PHI*M_PHI / (M_PHI*M_PHI+Q2) * f1Phi +
    (1.0 - gwfw) * f1D;

  // (308)
  double f1_iv = 
    grfr * M_RHO*M_RHO / (M_RHO*M_RHO+Q2) * f1 +
    (1.0 - grfr) * f1D;

  // (306)
  return (f1_is + f1_iv) / 2.0;  
}


/*! Proton Pauli formfactor parameterized by Gari and Kruempelmann
 * Phys.Lett.B,173:10 and Phys.Rev.D,45:1817
 *
 * All formula references are to p.49-50 of Stijn's notes
 *
 * 13 parameters
 */
double gariKruempelmann_proton_pauli(double Q2,double k_iv,double k_is,double grfr,
				     double kr,double gwfw,double kw, double gfff,
				     double kf,double muf,double l1rw,double l1D,
				     double l2,double lQCD,double p14,double p15)
{
  double kappa = (k_is + k_iv ) / 2.0; // proton magnetic moment

  // Q~^2, (316)
  double Qt2 = Q2*log( (l2*l2+Q2)/lQCD/lQCD ) / log( l2*l2/lQCD/lQCD );

  // (313)
  double f2 = pow(l1rw*l1rw / (l1rw*l1rw+Qt2),2) * l2*l2 / (l2*l2+Qt2);

  // (313)
  double f2D = pow(l1D*l1D / (l1D*l1D+Qt2),2) * l2*l2 / (l2*l2+Qt2);

  // (315)
  double f2Phi = f2 * pow(l1rw*l1rw/M_PHI/M_PHI * (M_PHI*M_PHI+Q2)/(l1rw*l1rw+Q2) ,1.5);


  // (311)
  double f2_is = 
    kw*gwfw * M_OMEGA*M_OMEGA/(M_OMEGA*M_OMEGA+Q2) * f2 +
    kf*gfff * M_PHI*M_PHI/(M_PHI*M_PHI+Q2) * f2Phi +
    (k_is-kw*gwfw-kf*gfff) * f2D;

  // (309)
  double f2_iv = 
    kr*grfr * M_RHO*M_RHO/(M_RHO*M_RHO+Q2) * f2 +
    (k_iv-kr*grfr) * f2D;

  // (307)
  return (f2_is+f2_iv) / 2.0 / kappa;
}

/************************************************************************/

/*! Neutron Dirac formfactor parameterized by Gari and Kruempelmann
 * Phys.Lett.B,173:10 and Phys.Rev.D,45:1817
 *
 * All formula references are to p.49-50 of Stijn's notes
 *
 * 13 parameters
 */
double gariKruempelmann_neutron_dirac(double Q2,double k_iv,double k_is,double grfr,
				     double kr,double gwfw,double kw, double gfff,
				     double kf,double muf,double l1rw,double l1D,
				     double l2,double lQCD,double p14,double p15)
{
  // Q~^2, (316)
  double Qt2 = Q2*log( (l2*l2+Q2)/lQCD/lQCD ) / log( l2*l2/lQCD/lQCD );

  // (312)
  double f1 = l1rw*l1rw / (l1rw*l1rw+Qt2) * l2*l2 / (l2*l2+Qt2);

  // (312)
  double f1D = l1D*l1D / (l1D*l1D+Qt2) * l2*l2 / (l2*l2+Qt2);

  // (314)
  double f1Phi = f1 * pow(Q2/(l1rw*l1rw+Q2),1.5);


  // (310)
  double f1_is = 
    gwfw * M_OMEGA*M_OMEGA / (M_OMEGA*M_OMEGA+Q2) * f1 +
    gfff * M_PHI*M_PHI / (M_PHI*M_PHI+Q2) * f1Phi +
    (1.0 - gwfw) * f1D;

  // (308)
  double f1_iv = 
    grfr * M_RHO*M_RHO / (M_RHO*M_RHO+Q2) * f1 +
    (1.0 - grfr) * f1D;

  // (306)
  return (f1_is - f1_iv) / 2.0;  
}


/*! Neutron Pauli formfactor parameterized by Gari and Kruempelmann
 * Phys.Lett.B,173:10 and Phys.Rev.D,45:1817
 *
 * All formula references are to p.49-50 of Stijn's notes
 *
 * 13 parameters
 */
double gariKruempelmann_neutron_pauli(double Q2,double k_iv,double k_is,double grfr,
				     double kr,double gwfw,double kw, double gfff,
				     double kf,double muf,double l1rw,double l1D,
				     double l2,double lQCD,double p14,double p15)
{
  double kappa = (k_is - k_iv ) / 2.0; // neutron magnetic moment

  // Q~^2, (316)
  double Qt2 = Q2*log( (l2*l2+Q2)/lQCD/lQCD ) / log( l2*l2/lQCD/lQCD );

  // (313)
  double f2 = pow(l1rw*l1rw / (l1rw*l1rw+Qt2),2) * l2*l2 / (l2*l2+Qt2);

  // (313)
  double f2D = pow(l1D*l1D / (l1D*l1D+Qt2),2) * l2*l2 / (l2*l2+Qt2);

  // (315)
  double f2Phi = f2 * pow(l1rw*l1rw/M_PHI/M_PHI * (M_PHI*M_PHI+Q2)/(l1rw*l1rw+Q2) ,1.5);


  // (311)
  double f2_is = 
    kw*gwfw * M_OMEGA*M_OMEGA/(M_OMEGA*M_OMEGA+Q2) * f2 +
    kf*gfff * M_PHI*M_PHI/(M_PHI*M_PHI+Q2) * f2Phi +
    (k_is-kw*gwfw-kf*gfff) * f2D;

  // (309)
  double f2_iv = 
    kr*grfr * M_RHO*M_RHO/(M_RHO*M_RHO+Q2) * f2 +
    (k_iv-kr*grfr) * f2D;

  // (307)
  return (f2_is-f2_iv) / 2.0 / kappa;
}

/************************************************************************/

/*! Bonn model parametrization of Lambda elastic Dirac formfactor
 * (no paper yet - see Tim for more info)
 *
 * G_E = sum of 3 dipoles
 * G_M = dipole
 *
 * -> 7 parameters
 */
double bonn_lambda_dirac(double Q2,double cc1,double cc2,double cc3,
			 double le1,double le2,double le3,double lm,
			 double p8,double p9,double p10,double p11,
			 double p12,double p13,double p14,double p15)
{
  // kinematic factor
  double tau = Q2 / (4.0 * M_P * M_P); 

  // electric Sachs formfactor
  double g_e = 
    cc1 * dipole(Q2,le1) +
    cc2 * dipole(Q2,le2) +
    cc3 * dipole(Q2,le3);

  // magnetic Sachs formfactor
  double g_m = KAPPA_L * dipole(Q2,lm);


  // Dirac formfactor
  return 1.0 / (1.0 + tau) * (tau * g_m + g_e);
}


/*! Bonn model parametrization of Lambda elastic Pauli formfactor
 * (no paper yet - see Tim for more info)
 *
 * G_E = sum of 3 dipoles
 * G_M = dipole
 *
 * -> 7 parameters
 */
double bonn_lambda_pauli(double Q2,double cc1,double cc2,double cc3,
			 double le1,double le2,double le3,double lm,
			 double p8,double p9,double p10,double p11,
			 double p12,double p13,double p14,double p15)
{
  // kinematic factor
  double tau = Q2 / (4.0 * M_P * M_P); 

  // electric Sachs formfactor
  double g_e = 
    cc1 * dipole(Q2,le1) +
    cc2 * dipole(Q2,le2) +
    cc3 * dipole(Q2,le3);

  // magnetic Sachs formfactor
  double g_m = KAPPA_L * dipole(Q2,lm);


  // Pauli formfactor
  return 1.0 / (KAPPA_L * (1.0 + tau)) * (g_m - g_e);
}

/************************************************************************/

/*! Bonn model parametrization of Sigma^0 elastic Dirac formfactor
 * (no paper yet - see Tim for more info)
 *
 * G_E = sum of 2 dipoles
 * G_M = dipole
 *
 * -> 5 parameters
 */
double bonn_sigma0_dirac(double Q2,double cc1,double cc2,double le1,
			 double le2,double lm,double p6,double p7,
			 double p8,double p9,double p10,double p11,
			 double p12,double p13,double p14,double p15)
{
  // kinematic factor
  double tau = Q2 / (4.0 * M_P * M_P); 

  // electric Sachs formfactor
  double g_e = 
    cc1 * dipole(Q2,le1) +
    cc2 * dipole(Q2,le2);

  // magnetic Sachs formfactor
  double g_m = KAPPA_S0 * dipole(Q2,lm);


  // Dirac formfactor
  return 1.0 / (1.0 + tau) * (tau * g_m + g_e);
}


/*! Bonn model parametrization of Sigma^0 elastic Pauli formfactor
 * (no paper yet - see Tim for more info)
 *
 * G_E = sum of 2 dipoles
 * G_M = dipole
 *
 * -> 5 parameters
 */
double bonn_sigma0_pauli(double Q2,double cc1,double cc2,double le1,
			 double le2,double lm,double p6,double p7,
			 double p8,double p9,double p10,double p11,
			 double p12,double p13,double p14,double p15)
{
  // kinematic factor
  double tau = Q2 / (4.0 * M_P * M_P); 

  // electric Sachs formfactor
  double g_e = 
    cc1 * dipole(Q2,le1) +
    cc2 * dipole(Q2,le2);

  // magnetic Sachs formfactor
  double g_m = KAPPA_S0 * dipole(Q2,lm);


  // Pauli formfactor
  return 1.0 / (KAPPA_S0 * (1.0 + tau)) * (g_m - g_e);
}

/************************************************************************/

/*! Bonn model parametrization of Lambda-Sigma^0 transition Dirac formfactor
 * (no paper yet - see Tim for more info)
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_lambda_sigma0_dirac(double Q2,double p1,double p2,double p3,
				double p4,double p5,double p6,double p7,
				double p8,double p9,double p10,double p11,
				double p12,double p13,double p14,double p15)
{
  // convert Q^2 to GeV^2
  Q2 /= 1e6;

  // return the formfactor
  return 
    Q2 / ( Q2 + 3.10095 ) +  // makes sure that =0 @ Q^2=0
    0.400414 * dipole(Q2,0.944691) +
    2.25112  * dipole(Q2,0.842582) -
    0.205185 * dipole(Q2,1.772250);
}


/*! Bonn model parametrization of Lambda-Sigma^0 transition Pauli formfactor
 * (no paper yet - see Tim for more info)
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_lambda_sigma0_pauli(double Q2,double p1,double p2,double p3,
				double p4,double p5,double p6,double p7,
				double p8,double p9,double p10,double p11,
				double p12,double p13,double p14,double p15)
{
  // convert Q^2 to GeV^2
  Q2 /= 1e6;

  // return the formfactor
  return 
    0.334385  * dipole(Q2,0.859923) -
    0.0816798 * dipole(Q2,2.010060) +
    0.747295  * dipole(Q2,0.731168);
}

/************************************************************************/

/*! Kaon^+ formfactor parametrization by David et al. (Phys.Rev.C,53:2613)
 * based on a light-front CQM
 *
 * 3 parameters
 */
double david_kaonPlus(double Q2,double a,double l1,double l2,
		      double p4,double p5,double p6,double p7,
		      double p8,double p9,double p10,double p11,
		      double p12,double p13,double p14,double p15)
{
  // formula = (320) on p.51 of Stijn's notes
  return a * monopole(Q2,l1) + (1.0-a) * dipole(Q2,l2);
}

/************************************************************************/

/*! Kaon^0 formfactor parametrization by Bennhold et al. based on VMD
 * (proceedings of the 7th International Conference on the Structure of Baryons,
 *  Baryons '95* 323-326)
 *
 * no parameters
 */
double bennhold_kaon0(double Q2,double p1,double p2,double p3,
		      double p4,double p5,double p6,double p7,
		      double p8,double p9,double p10,double p11,
		      double p12,double p13,double p14,double p15)
{
  // formula = (322) on p.51 of Stijn's notes
  return -1.0/3.0 * monopole(Q2,M_OMEGA) + 1.0/3.0 * monopole(Q2,M_PHI);
}

/************************************************************************/

/*! Kaon^*+ - K^+ transition formfactor calculated by Muenz et al.
 * (Phys.Rev.C,52:2110)
 *
 * There's no parametrization.
 * We apply linear interpolation between their calculated values
 */
double muenz_kaonStarPlus(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  /* Array containing the calculated formfactor values: {Q^2 in GeV^2, F(Q^2)}
   * The array is declared static. This way it's only initialized once */
  
  static double ff[35][2] = {{0.0, 1},         
			     {0.187721, 0.782656},  
			     {0.294446, 0.576122},    
			     {0.437009, 0.398894},    
			     {0.651903, 0.270631},  
			     {0.938997, 0.176278},  
			     {1.15389, 0.119597},  
			     {1.40449, 0.0747761}, 
			     {1.72716, 0.0439776}, 
			     {1.9414, 0.0269385}, 
			     {2.15564, 0.0165012}, 
			     {2.33351, 0.00990273},
			     {2.43893, 0.00594213},
			     {2.54381, 0.00328569},
			     {2.57676, 0.00197134},
			     {2.68112, 0.00100448},
			     {2.70000, 0.00000000},
			     {2.78995, 0.00102541},
			     {2.82868, 0.00151203},
			     {2.86793, 0.00241951},
			     {2.94341, 0.00387187},
			     {3.12719, 0.00582861},
			     {3.34655, 0.00792242},
			     {3.71031, 0.0099255}, 
			     {4.14574, 0.0110014}, 
			     {4.61701, 0.0114693}, 
			     {5.05166, 0.0112455}, 
			     {5.70356, 0.0108069}, 
			     {6.46389, 0.00997122},
			     {7.07917, 0.0090119}, 
			     {7.73068, 0.00814537},
			     {8.27349, 0.00736082},
			     {8.74385, 0.00665104},
			     {9.25056, 0.00613415},
			     {9.68468, 0.00554231}};

  // convert Q^2 from MeV^2 to GeV^2
  Q2 /= 1e6;

  if(Q2 >= 9.68468)
    {
      cout << "(K*^+)-(K^+) transition form factor calculation by Muenz et al. "
	   << "does not exist for Q^2=" << Q2 << " GeV^2" << endl;
      exit(1);
    }


  /* Determine form factor value for given Q^2 */
  bool match = false; // becomes true when correct position in ff data table is found
  double formfactor = 0; // this variable will hold the formfactor value for Q^2
  
  for(int i=0; i<35 && !match; i++)
    {
      if(i > 0 && Q2 < ff[i][0] && Q2 >= ff[i-1][0])
	{
	  match = true; // correct position has been found
	  
	  /* Linear interpolation between points 
	   * (x1,y1) = (ff[i-1][0], ff{i-1][1])
	   * and (x2,y2) = (ff[i][0], ff[i][1]), with x1 <= Q^2 <= x2:
	   * y = y1 + (y2-y1)/(x2-x1) * (x-x1) */
	  
	  formfactor =  ff[i-1][1] + (ff[i][1] - ff[i-1][1]) * 
	    (Q2 - ff[i-1][0]) / (ff[i][0] - ff[i-1][0]);
	}
    }
  
  return formfactor;
}

/************************************************************************/

/*! Kaon^*0 - K^0 transition formfactor calculated by Muenz et al.
 * (Phys.Rev.C,52:2110)
 *
 * There's no parametrization.
 * We apply linear interpolation between their calculated values
 */
double muenz_kaonStar0(double Q2,double p1,double p2,double p3,
			double p4,double p5,double p6,double p7,
			double p8,double p9,double p10,double p11,
			double p12,double p13,double p14,double p15)
{
  /* Array containing the calculated formfactor values: {Q^2 in GeV^2, F(Q^2)}
   * The array is declared static. This way it's only initialized once */
  
  static double ff[25][2] = {{0.0, 1},         
			     {0.260577, 0.832243},  
			     {0.367696, 0.65136},    
			     {0.583772, 0.531156},    
			     {0.835948, 0.424398},  
			     {1.12422,  0.332258},  
			     {1.41263,  0.265493},  
			     {1.70104,  0.212144}, 
			     {2.09827,  0.173047}, 
			     {2.38668,  0.138274}, 
			     {2.78404, 0.11512}, 
			     {3.21738, 0.0920094},
			     {3.68707, 0.0750614},
			     {4.15703, 0.0637901},
			     {4.62699, 0.0542112},
			     {5.09695, 0.0460707},
			     {5.56703, 0.039961},
			     {6.14569, 0.0339666},
			     {6.61604, 0.0306913},
			     {7.15859, 0.0266243},
			     {7.80984, 0.0231005},
			     {8.42472, 0.0196365},
			     {9.0039,  0.0181126},
			     {9.43815, 0.0167029}, 
			     {9.79994, 0.0154011}};

  // convert Q^2 from MeV^2 to GeV^2
  Q2 /= 1e6;


  if(Q2 >= 9.79994)
    {
      cout << "(K*^0)-(K^0) transition form factor calculation by Muenz et al. "
	   << "does not exist for Q^2=" << Q2 << " GeV^2" << endl;
      exit(1);
    }

  /* Determine form factor value for given Q^2 */
  bool match = false; // becomes true when correct position in ff data table is found
  double formfactor = 0; // this variable will hold the formfactor value for Q^2
  
  for(int i=0; i<25 && !match; i++)
    {
      if(i > 0 && Q2 < ff[i][0] && Q2 >= ff[i-1][0])
	{
	  match = true; // correct position has been found
	  
	  /* Linear interpolation between points 
	   * (x1,y1) = (ff[i-1][0], ff{i-1][1])
	   * and (x2,y2) = (ff[i][0], ff[i][1]), with x1 <= Q^2 <= x2:
	   * y = y1 + (y2-y1)/(x2-x1) * (x-x1) */
	  
	  formfactor = ff[i-1][1] + (ff[i][1] - ff[i-1][1]) * 
	    (Q2 - ff[i-1][0]) / (ff[i][0] - ff[i-1][0]);
	}
    }
  
  return formfactor;
}

/************************************************************************/

/*! Bonn Parametrization for P_11(1710) 1/2^+ resonance
 * Dirac form factor
 *
 * = sum of 3 dipoles times monopole going to zero
 *
 * -> no parameters
 */
double bonn_P11_1710_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return Q2 / ( Q2 + 2.33818e6 ) / (-0.335529) *
    ( 0.64583  * dipole(Q2,967.961) -
      0.800948 * dipole(Q2,1341.01) +
      0.680674 * dipole(Q2,1168.22) );
}


/*! Bonn Parametrization for P_11(1710) 1/2^+ resonance
 * Pauli form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_P11_1710_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return ( -0.380619 * dipole(Q2,2907.89) -
	    1.883170 * dipole(Q2,2777.56) +
	    3.263790 * dipole(Q2,2278.12) );
}

/*! Bonn Parametrization for P_11(1710) 1/2^+ resonance
 * Dirac form factor off the neutron
 *\verbatim
              3                     
       2    _____         a         
      Q     \              i        
 =  _____ *  \     ______________   
     2       /             2        
    Q + C   /             Q   2     
            _____    (1 + ___)      
             i=1            2       
                          /\        
                            i       
 \endverbatim
 *
 * -> no parameters
 */
double bonn_P11_1710_neutron_dirac(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return Q2 / ( Q2 + 1.33691e6 ) *
    ( 0.640258   * dipole(Q2,1212.1881)
      -0.884962  * dipole(Q2,1159.2972)
      -0.0715926 * dipole(Q2,1087.6994) );
}


/*! Bonn Parametrization for P_11(1710) 1/2^+ resonance
 * Pauli form factor off the neutron
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_P11_1710_neutron_pauli(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return 
    3.45861   * dipole(Q2,1146.9438)
    -0.607387 * dipole(Q2,1524.4015)
    -1.85122  * dipole(Q2,989.02983);
}

/************************************************************************/

/*! Bonn Parametrization for P_31(1910) 1/2^+ resonance
 * Dirac form factor
 *
 * = sum of 4 dipoles
 *
 * -> no parameters
 */
double bonn_P31_1910_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return ( -7.66308  * dipole(Q2,1002.287) +
	    5.01613  * dipole(Q2,1077.418) +
	    3.15644  * dipole(Q2, 898.497) -
	    0.509131 * dipole(Q2,0855.260) ) / 0.944;
}


/* Bonn Parametrization for P_31(1910) 1/2^+ resonance
 * Pauli form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_P31_1910_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return ( -6.00467 * dipole(Q2,1676.732) +
	    2.00972 * dipole(Q2,2087.499) +
	    4.99548 * dipole(Q2,1332.333) );
}

/************************************************************************/

/*! Bonn Parametrization for S_11(1535) 1/2^- resonance
 * Dirac form factor
 *
 * = sum of 3 dipoles times monopole going to zero
 *
 * -> no parameters
 */
double bonn_S11_1535_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return Q2 / ( Q2 + 1.59678e6 ) / (-0.869181) *
    ( 0.0308628 * dipole(Q2,3227.30) +
      0.0259375 * dipole(Q2,848.205) -
      1.3463000 * dipole(Q2,784.119) );
}


/*! Bonn Parametrization for S_11(1535) 1/2^- resonance
 * Pauli form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_S11_1535_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return ( 1.1591100 * dipole(Q2,1052.19) -
	   0.2321680 * dipole(Q2,1859.03) +
	   0.0730628 * dipole(Q2,435.696) );
}


/*! Bonn Parametrization for S_11(1535) 1/2^- resonance
 * Dirac form factor off the neutron
 *\verbatim
      4
    _____     
    \        
     \             j
     /      a  *  Q
    /        j 
    _____     
     j=1      
 = ___________________________  
      
      7
    _____     
    \        
     \             i
     /      a  *  Q
    /        i 
    _____     
     i=0 
 \endverbatim
 *
 * -> no parameters
 */
double bonn_S11_1535_neutron_dirac(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  double Q=sqrt(Q2)/1e3; // in GeV

  return 
    (  2.70155 * Q + 3.27058 * Q*Q 
       - 0.25473 * Q*Q*Q - 0.105993 * Q*Q*Q*Q )
    /
    ( 3.13816 + 1.81194 * Q + 9.06668 * Q*Q + 16.5776 * Q*Q*Q
      + 13.7903 * Q*Q*Q*Q + 5.48482 * Q*Q*Q*Q*Q
      + 6.752 * Q*Q*Q*Q*Q*Q + 18.5677 * Q*Q*Q*Q*Q*Q*Q);
}


/*! Bonn Parametrization for S_11(1535) 1/2^- resonance
 * Pauli form factor off the neutron
 *\verbatim
       2
                            2
     _____    a * exp(-e * Q )
     \         i        i
      \     _________________
 =    /             2
     /             Q   2
     _____    (1 + ___)
      i=1            2
                   /\
                     i
\endverbatim
 *
 * -> no parameters
 */
double bonn_S11_1535_neutron_pauli(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return 
    ( -6.37776 * exp(-19.2722e-6 * Q2) / pow( 1.-Q2/-1.19205e-6 , 2) )
    +( 0.996767  * exp(-0.882475e-6 * Q2) / pow( 1.+Q2/0.866405e6  , 2) );
}

/************************************************************************/

/*! Bonn Parametrization for S_11(1650) 1/2^- resonance
 * Dirac form factor
 *
 * = sum of 3 dipoles times monopole going to zero
 *
 * -> no parameters
 */
double bonn_S11_1650_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return Q2 / ( Q2 + 3.08025e6 ) / (0.0132487) *
    ( 0.160509 * dipole(Q2,1585.47) -
      0.479094 * dipole(Q2,1431.48) +
      0.185153 * dipole(Q2,1406.48) );
}


/*! Bonn Parametrization for S_11(1650) 1/2^- resonance
 * Pauli form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_S11_1650_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return ( -3.872870 * dipole(Q2,1760.74) +
	    0.205359 * dipole(Q2,6791.44) +
	    4.667510 * dipole(Q2,1121.89) );
}


/*! Bonn Parametrization for S_11(1650) 1/2^- resonance
 * Dirac form factor off the neutron
 *\verbatim
      4
    _____     
    \        
     \             j
     /      a  *  Q
    /        j 
    _____     
     j=1      
 = ___________________________  
      
      7
    _____     
    \        
     \             i
     /      a  *  Q
    /        i 
    _____     
     i=0 
 \endverbatim
 *
 * -> no parameters
 */
double bonn_S11_1650_neutron_dirac(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  double Q=sqrt(Q2)/1e3; // in GeV

  return 
    ( 0.0366997 * Q + 0.117971 * Q*Q 
      -0.00533585 * Q*Q*Q + 0.00660971 * Q*Q*Q*Q )
    /
    ( 2.05416 + 3.55349 * Q + 6.07053 * Q*Q + 8.59447 * Q*Q*Q
      + 10.6793 * Q*Q*Q*Q + 12.2934 * Q*Q*Q*Q*Q
      + 12.4344 * Q*Q*Q*Q*Q*Q + 8.74409 * Q*Q*Q*Q*Q*Q*Q);
}


/*! Bonn Parametrization for S_11(1650) 1/2^- resonance
 * Pauli form factor off the neutron
 *\verbatim
       2
                            2
     _____    a * exp(-e * Q )
     \         i        i
      \     _________________
 =    /             2
     /             Q   2
     _____    (1 + ___)
      i=1            2
                   /\
                     i
 \endverbatim
 *
 * -> no parameters
 */
double bonn_S11_1650_neutron_pauli(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return 
    ( -0.027604 * exp(-19.8939e-6 * Q2) / pow( 1.-Q2/1.52548e-1 , 2) )
    +( 1.00088  * exp(-0.654771e-6 * Q2) / pow( 1.+Q2/3.69578e6  , 2) );
}

/************************************************************************/

/*! Bonn Parametrization for S_31(1900) 1/2^- resonance
 * Dirac form factor
 *
 * = sum of 3 dipoles times monopole going to zero
 *
 * -> no parameters
 */
double bonn_S31_1900_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of S_31(1900)"
       << " is not to be trusted!" << endl << endl;

  return Q2 / ( Q2 + 1.32098e6 ) / (0.0165148) *
    ( 0.995750 * dipole(Q2,866.548) -
      0.731828 * dipole(Q2,785.985) -
      0.220265 * dipole(Q2,775.505) );
}


/*! Bonn Parametrization for S_31(1900) 1/2^- resonance
 * Pauli form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_S31_1900_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of S_31(1900)"
       << " is not to be trusted!" << endl << endl;

  return ( 11.08070 * dipole(Q2,1357.65) +
	   0.698568 * dipole(Q2,5218.43) +
	   -10.7793 * dipole(Q2,1749.51) );
}

/************************************************************************/

/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P13_1720_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -1.37364 * dipole(Q2,1123.886) * exp(-1.514540e-6 * Q2) +
      2.37382 * dipole(Q2,1028.718) * exp(-0.143659e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P13_1720_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( 0.580233 * dipole(Q2,766.260) * exp(-0.358284e-6 * Q2) +
      0.419813 * dipole(Q2,991.895) * exp(-1.515110e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_1 form factor off the neutron
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P13_1720_neutron_1(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    ( -5.10196 * dipole(Q2,1153.) * exp(-0.270413e-6 * Q2)
      +6.10146 * dipole(Q2,1254.) * exp(-0.268665e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_2 form factor off the neutron
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P13_1720_neutron_2(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    ( 0.584837 * dipole(Q2,774.3) * exp(-0.213376e-6 * Q2) +
      0.415224 * dipole(Q2,911.5) * exp(-0.228213e-6 * Q2) );
}

/************************************************************************/

/*! Bonn Parametrization for P_13(1900) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P13_1900_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of P_13(1900)"
       << " is not to be trusted!" << endl << endl;

  return 
    ( -3.82382 * dipole(Q2,7246.123) * exp(-1.467270e-6 * Q2) +
      4.82588 * dipole(Q2,1875.857) * exp(-0.238605e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1900) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles  and 1 dipole with exponential
 *
 * -> no parameters
 */
double bonn_P13_1900_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of P_13(1900)"
       << " is not to be trusted!" << endl << endl;

  return 
    ( 4.68599 * dipole(Q2,2322.815) * exp(-0.158133e-6 * Q2) +
      5.00753 * dipole(Q2,1033.601) -
      8.69102 * dipole(Q2,1328.119) );
}

/************************************************************************/

/*! Bonn Parametrization for P_33(1920) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_P33_1920_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  /* exc3 uit Bonn-model, 1866 MeV */

  return ( -5.38627  * dipole(Q2,1644.600) +
	   5.49611  * dipole(Q2,1504.729) +
	   0.892168* dipole(Q2,1504.939) );
}



/*! Bonn Parametrization for P_33(1920) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_P33_1920_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  /* exc4 uit Bonn-model, 1951 MeV */
  /*
    return 
    ( 4.02690 * dipole(Q2,1068.836) * exp(-1.327680e-6 * Q2) -
      3.02677 * dipole(Q2,1110.086) * exp(-0.209678e-6 * Q2) );
  */

  /* exc3 uit Bonn-model, 1871 MeV */

  return ( -5.23423 * dipole(Q2,1444.479) +
	   1.80138 * dipole(Q2,1344.288) +
	   4.43416 * dipole(Q2,1346.945) );
}


/************************************************************************/

/*! Bonn Parametrization for D_13(1700) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1700_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -4.38057 * dipole(Q2,1134.976) +
      5.38332 * dipole(Q2,1037.642) );
}


/*! Bonn Parametrization for D_13(1700) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1700_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -8.48840 * dipole(Q2,998.747) +
      9.48879 * dipole(Q2,952.542) );
}

/************************************************************************/

/*! Bonn Parametrization for D_33(1700) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_D33_1700_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -3.98846 * dipole(Q2,1056.650) * exp(-0.271889e-6 * Q2) +
      4.98866 * dipole(Q2,1135.487) * exp(-0.325962e-6 * Q2) );

//     ( -4.38057 * dipole(Q2,1134.976) +
//        5.38332 * dipole(Q2,1037.642) );
}


/*! Bonn Parametrization for D_33(1700) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double bonn_D33_1700_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( //-2.84609e-11 * dipole(Q2,0.0) * exp(-19.5074e-6 * Q2) +
     1.03325     * dipole(Q2,24277.665) * exp(-1.09488e-6 * Q2) );
}

/************************************************************************/

/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1520_1(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
        		   double p8,double p9,double p10,double p11,
		           double p12,double p13,double p14,double p15)
{
  return 
    ( -4.44054 * dipole(Q2, 1163.53341) +
      5.44026 * dipole(Q2, 1087.99816) );
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1520_2(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
        		   double p8,double p9,double p10,double p11,
		           double p12,double p13,double p14,double p15)
{
  return 
    (  -4.93527 * dipole(Q2, 561.846954) +
       5.93212 * dipole(Q2, 588.288195) );
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_1 form factor off the neutron
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1520_neutron_1(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    -4.39709  * dipole(Q2,1155.0844)
    + 5.39697 * dipole(Q2,1062.7135);
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_2 form factor off the neutron
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_D13_1520_neutron_2(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    -2.22492  * dipole(Q2,947.69035)
    + 3.22323 * dipole(Q2,751.65085);
}

/************************************************************************/

/*! Bonn Parametrization for F_15(1680) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_F15_1680_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    5.83527 * dipole(Q2,1166.69)
		- 4.83527 * dipole(Q2,1234.39);
}


/*! Bonn Parametrization for F_15(1680) 5/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double bonn_F15_1680_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    13.609 * dipole(Q2,925.001)
		- 12.609 * dipole(Q2,946.076);
}

/************************************************************************/

/*! Bonn Parametrization for D_15(1675) 5/2^- resonance
 * kappa_1 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double bonn_D15_1675_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    1.01618 * dipole(Q2,1681.09) * exp(-0.731708e-6*Q2)
		-0.01618 * dipole(Q2,1463.54);
}


/*! Bonn Parametrization for D_15(1675) 5/2^- resonance
 * kappa_2 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double bonn_D15_1675_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    1.025 * dipole(Q2,900.054) * exp(-1.00011e-6*Q2)
		-0.025 * dipole(Q2,869.888);
}

/************************************************************************/

/*! Bonn Parametrization for F_15(2000) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_F15_2000_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for F_15(2000) 5/2^+ resonance
 * kappa_2 form factor
 *
 * dipole
 *
 * -> no parameters
 */
double bonn_F15_2000_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for D_15(2200) 5/2^- resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_D15_2200_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for D_15(2200) 5/2^- resonance
 * kappa_2 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_D15_2200_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for F_35(1905) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_F35_1905_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for F_35(1905) 5/2^+ resonance
 * kappa_2 form factor
 *
 * dipole
 *
 * -> no parameters
 */
double bonn_F35_1905_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for D_35(1930) 5/2^- resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_D35_1930_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for D_35(1930) 5/2^- resonance
 * kappa_2 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_D35_1930_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for F_35(2000) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double bonn_F35_2000_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for F_35(2000) 5/2^+ resonance
 * kappa_2 form factor
 *
 * dipole
 *
 * -> no parameters
 */
double bonn_F35_2000_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/************************************************************************/

/*! Proton electromagnetic transition form factor of the VR model
 * 
 * -> 1 parameter
 */
double protonEMTFF(double Q2, double mandel, double Lambda_infinity)
{
  if(mandel < M_P*M_P)  // n(e,e'pi^-)p (s = u)
    mandel = 2*M_P*M_P - mandel;

  return dipole(Q2, 840. + (Lambda_infinity - 840.)*(1 - M_P*M_P/mandel));
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/************************************************************************/

/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_P13_1720_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -1.37364 * dipole(Q2,1123.886) * exp(-1.514540e-6 * Q2) +
      2.37382 * dipole(Q2,1028.718) * exp(-0.143659e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_P13_1720_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    (2.46148 * dipole(Q2,726.758) * exp(-0.369751e-6 * Q2)
     - 1.4648 * dipole(Q2,775.538) * exp(-2.38661e-6 * Q2));
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_1 form factor off the neutron
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_P13_1720_neutron_1(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    ( -5.10196 * dipole(Q2,1153.) * exp(-0.270413e-6 * Q2)
      +6.10146 * dipole(Q2,1254.) * exp(-0.268665e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1720) 3/2^+ resonance
 * kappa_2 form factor off the neutron
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_P13_1720_neutron_2(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    ( 0.584837 * dipole(Q2,774.3) * exp(-0.213376e-6 * Q2) +
      0.415224 * dipole(Q2,911.5) * exp(-0.228213e-6 * Q2) );
}

/************************************************************************/

/*! Bonn Parametrization for P_13(1900) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_P13_1900_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of P_13(1900)"
       << " is not to be trusted!" << endl << endl;

  return 
    ( -3.82382 * dipole(Q2,7246.123) * exp(-1.467270e-6 * Q2) +
      4.82588 * dipole(Q2,1875.857) * exp(-0.238605e-6 * Q2) );
}


/*! Bonn Parametrization for P_13(1900) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles and 1 dipole with exponential
 *
 * -> no parameters
 */
double gic_bonn_P13_1900_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  cout << "The Bonn-parametrization of the EM form factor of P_13(1900)"
       << " is not to be trusted!" << endl << endl;

  return 
    ( 4.08501 * dipole(Q2,4612.74) * exp(-0.528292e-6 * Q2) +
      0.363665 * dipole(Q2,2890.6) -
      3.448675 * dipole(Q2,1165.11) );
}

/************************************************************************/

/*! Bonn Parametrization for P_33(1920) 3/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double gic_bonn_P33_1920_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  /* exc3 uit Bonn-model, 1866 MeV */

  return ( -5.38627  * dipole(Q2,1644.600) +
	   5.49611  * dipole(Q2,1504.729) +
	   0.892168* dipole(Q2,1504.939) );
}



/*! Bonn Parametrization for P_33(1920) 3/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_P33_1920_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  /* exc3 uit Bonn-model, 1871 MeV */

  return ( 1.1358 * dipole(Q2,1280.96) * exp(-0.894633e-6 * Q2)
	   - 0.1358 * dipole(Q2,1095.46));
}


/************************************************************************/

/*! Bonn Parametrization for D_13(1700) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_D13_1700_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -4.38057 * dipole(Q2,1134.976) +
      5.38332 * dipole(Q2,1037.642) );
}


/*! Bonn Parametrization for D_13(1700) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_D13_1700_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( 1.44604 * dipole(Q2,1330.09) * exp(-1.37536e-6 * Q2)
      - 0.44604 * dipole(Q2,736.894) );
}

/************************************************************************/

/*! Bonn Parametrization for D_33(1700) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_D33_1700_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( -3.98846 * dipole(Q2,1056.650) * exp(-0.271889e-6 * Q2) +
      4.98866 * dipole(Q2,1135.487) * exp(-0.325962e-6 * Q2) );

//     ( -4.38057 * dipole(Q2,1134.976) +
//        5.38332 * dipole(Q2,1037.642) );
}


/*! Bonn Parametrization for D_33(1700) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_D33_1700_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( 1.36855 * dipole(Q2,836.969) * exp(-0.753232e-6 * Q2)
      - 0.36855 * dipole(Q2,819.342) );
}

/************************************************************************/

/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_D13_1520_1(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
        		   double p8,double p9,double p10,double p11,
		           double p12,double p13,double p14,double p15)
{
  return 
    ( -4.44054 * dipole(Q2, 1163.53341) +
      5.44026 * dipole(Q2, 1087.99816) );
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_2 form factor
 *
 * = sum of a dipole and a dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_D13_1520_2(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
        		   double p8,double p9,double p10,double p11,
		           double p12,double p13,double p14,double p15)
{
  return 
    ( 1.45955 * dipole(Q2,1055.0) * exp(-0.850355e-6 * Q2)
      - 0.45955 * dipole(Q2,555.6) );
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_1 form factor off the neutron
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_D13_1520_neutron_1(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    -4.39709  * dipole(Q2,1155.0844)
    + 5.39697 * dipole(Q2,1062.7135);
}


/*! Bonn Parametrization for D_13(1520) 3/2^- resonance
 * kappa_2 form factor off the neutron
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_D13_1520_neutron_2(double Q2,double p1,double p2,double p3,
			       double p4,double p5,double p6,double p7,
			       double p8,double p9,double p10,double p11,
			       double p12,double p13,double p14,double p15)
{
  return 
    -2.22492  * dipole(Q2,947.69035)
    + 3.22323 * dipole(Q2,751.65085);
}

/************************************************************************/

/*! Bonn Parametrization for F_15(1680) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_F15_1680_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( dipole(Q2,1060.3) * exp(-0.463387e-6 * Q2) );
}


/*! Bonn Parametrization for F_15(1680) 5/2^+ resonance
 * kappa_2 form factor
 *
 * = sum of 2 dipoles with exponentials
 *
 * -> no parameters
 */
double gic_bonn_F15_1680_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( 6.95476 * dipole(Q2,958.306) * exp(-0.548062e-6 * Q2)
      - 5.95476 * dipole(Q2,901.666) * exp(-0.549094e-6 * Q2) );
}

/************************************************************************/

/*! Bonn Parametrization for D_15(1675) 5/2^- resonance
 * kappa_1 form factor
 *
 * = dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_D15_1675_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( dipole(Q2,1718.63) * exp(-0.769309e-6 * Q2) );
}


/*! Bonn Parametrization for D_15(1675) 5/2^- resonance
 * kappa_2 form factor
 *
 * = dipole with an exponential
 *
 * -> no parameters
 */
double gic_bonn_D15_1675_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    ( dipole(Q2,1592.03) * exp(-0.695637e-6 * Q2) );
}

/************************************************************************/

/*! Bonn Parametrization for F_15(2000) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_F15_2000_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    //exc2
    1.13741 * dipole(Q2,800.743) - 0.13741 * dipole(Q2,1739.44);

    //exc3
    //0.812835 * dipole(Q2,1259.89) * exp(-0.874599e-6 * Q2)
    //+ 0.187165 * dipole(Q2,576.272);

    //exc4
    //dipole(Q2,919.391) * exp(-0.131374e-6 * Q2);
    //was: dipole(Q2,840.);
}


/*! Bonn Parametrization for F_15(2000) 5/2^+ resonance
 * kappa_2 form factor
 *
 * sum of 2 dipoles
 *
 * -> no parameters
 */
double gic_bonn_F15_2000_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    //exc2
    1.48043 * dipole(Q2,1040.79) - 0.48043 * dipole(Q2,1613.09);

    //exc3
    //1.86768 * dipole(Q2,767.867) * exp(-0.476268e-6 * Q2)
    //- 0.86768 * dipole(Q2,4327.38) * exp(-0.904051e-6 * Q2);

    //exc4
    //dipole(Q2,1093.63) * exp(-0.333112e-6 * Q2);
    // was: dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for D_15(2200) 5/2^- resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_D15_2200_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for D_15(2200) 5/2^- resonance
 * kappa_2 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_D15_2200_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for F_35(1905) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_F35_1905_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for F_35(1905) 5/2^+ resonance
 * kappa_2 form factor
 *
 * dipole
 *
 * -> no parameters
 */
double gic_bonn_F35_1905_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for D_35(1930) 5/2^- resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_D35_1930_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for D_35(1930) 5/2^- resonance
 * kappa_2 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_D35_1930_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for F_35(2000) 5/2^+ resonance
 * kappa_1 form factor
 *
 * = dipole
 *
 * -> no parameters
 */
double gic_bonn_F35_2000_1(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}


/*! Bonn Parametrization for F_35(2000) 5/2^+ resonance
 * kappa_2 form factor
 *
 * dipole
 *
 * -> no parameters
 */
double gic_bonn_F35_2000_2(double Q2,double p1,double p2,double p3,
		       double p4,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  return 
    dipole(Q2,840.);
}

/************************************************************************/

/*! Bonn Parametrization for P_11(1440) 1/2^+ resonance
 * Dirac form factor
 *\verbatim
      4
    _____     
    \        
     \             j
     /      a  *  Q
    /        j 
    _____     
     j=1      
 = ___________________________  
      
      7
    _____     
    \        
     \             i
     /      a  *  Q
    /        i 
    _____     
     i=0 
 \endverbatim
 *
 * -> no parameters
 */
double bonn_P11_1440_dirac(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
double Q=sqrt(Q2)/1e3; // in GeV

  return 
    ( -0.345722 * Q + 0.190066 * Q*Q 
      + 0.15152 * Q*Q*Q - 0.00649536 * Q*Q*Q*Q )
    /
    ( 1.07772 + 1.49531 * Q - 0.919313 * Q*Q + 4.32276 * Q*Q*Q
      + 4.60165 * Q*Q*Q*Q - 0.426772 * Q*Q*Q*Q*Q
      - 3.76203 * Q*Q*Q*Q*Q*Q + 4.54486 * Q*Q*Q*Q*Q*Q*Q);
}


/*! Bonn Parametrization for P_11(1440) 1/2^+ resonance
 * Pauli form factor
 *\verbatim
       3
                            2
     _____    a * exp(-e * Q )
     \         i        i
      \     _________________
 =    /             2
     /             Q   2
     _____    (1 + ___)
      i=1            2
                   /\
                     i
 \endverbatim
 *
 * -> no parameters
 */
double bonn_P11_1440_pauli(double Q2,double p1,double p2,double p3,
			   double p4,double p5,double p6,double p7,
			   double p8,double p9,double p10,double p11,
			   double p12,double p13,double p14,double p15)
{
  return 
    -1.01714   * dipole(Q2,1237.655) * exp(-2.37857e-6 * Q2)
    + 2.0087   * dipole(Q2,747.1078) * exp(-0.55817e-6 * Q2)
    + 0.008442 * dipole(Q2,6385.648) * exp(0.396131e-6 * Q2);
}

/*! Bonn Parametrization for P_11(1440) 1/2^+ resonance
 * Dirac form factor off the neutron
 *\verbatim
              3                     
       2    _____         a         
      Q     \              i        
 =  _____ *  \     ______________   
     2       /             2        
    Q + C   /             Q   2     
            _____    (1 + ___)      
             i=1            2       
                          /\        
                            i       
 \endverbatim
 *
 * -> no parameters
 */
double bonn_P11_1440_neutron_dirac(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return Q2 / ( Q2 + 4.01161e6) *
    ( 0.19634    * dipole(Q2,1558.034)
      + 0.841408 * dipole(Q2,983.43378)
      - 1.10078  * dipole(Q2,1240.7054));
}


/*! Bonn Parametrization for P_11(1440) 1/2^+ resonance
 * Pauli form factor off the neutron
 *
 * = sum of 3 dipoles
 *
 * -> no parameters
 */
double bonn_P11_1440_neutron_pauli(double Q2,double p1,double p2,double p3,
				   double p4,double p5,double p6,double p7,
				   double p8,double p9,double p10,double p11,
				   double p12,double p13,double p14,double p15)
{
  return 
    10.277    * dipole(Q2,1293.6344)
    - 1.79556 * dipole(Q2,940.71515)
    - 7.48144 * dipole(Q2,1365.3608);
}

/************************************************************************/

/*! Dipole formfactor for strong interaction vertex
 *
 * formula is (17) of Phys.Rev.C,73:045107
 */
double strong_dipole(double mandel,double mass,double cutoff,double p3,
		     double p4,double p5,double p6,double p7,
		     double p8,double p9,double p10,double p11,
		     double p12,double p13,double p14,double p15)
{
  return pow(cutoff,4) / ( pow(cutoff,4) + pow(mandel - mass*mass,2) );
}

/************************************************************************/
/*! Gaussian formfactor for strong interaction vertex
 * 
 * formula is (18) of Phys.Rev.C,73:045107
 * Modification for high-spin particles: the strong cut-off is a function
 * of the spin of the exchanged resonance.
 *	 \Lambda_{J'} = \Lambda_{J} * (J'/J)^{n}
 * where e.g. n = 1 or 1/2 or 1/4....
 * 
 */
double strong_gaussian(double mandel,double mass,double cutoff,double spin,
		       double exponent,double p5,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{
  if (spin!=0.5) cutoff /= pow(2.0*spin,exponent);
  cutoff *=cutoff;
  cutoff *=cutoff;
  return exp( -1.0 * (mandel - mass*mass)* (mandel - mass*mass) / cutoff);
}

/************************************************************************/
/*! "Lorentz" formfactor for strong interaction vertex
 * 
 * The maximum is ensured to remain at the value s = m^2
 * 
 */
double strong_lorentz(double mandel,double mass,double cutoff,double spin,
		       double exponent,double width,double p6,double p7,
		       double p8,double p9,double p10,double p11,
		       double p12,double p13,double p14,double p15)
{  
  double mw2=mass*width; 
  mw2 *= mw2;
  return pow( mw2/((mandel - mass*mass)*(mandel - mass*mass) + mw2) ,(spin-0.5));
}

/************************************************************************/
/*! "Lorentz" x Gaussian formfactor for strong interaction vertex
 * 
 * The maximum is ensured to remain at the value s = m^2 for higher
 * spin particles, and spin=1/2  particles' heavy tails are suppressed 
 * by the Gaussian ffac.
 * 
 * 
 * formula is (18) of Phys.Rev.C,73:045107
 * Modification for high-spin particles: the strong cut-off is a function
 * of the spin of the exchanged resonance.
 *	 \Lambda_{J'} = \Lambda_{J} * (J'/J)^{n}
 * where e.g. n = 1 or 1/2 or 1/4....
 * 
 */
double strong_lorentz_gaussian(double mandel,double mass,double cutoff,double spin,
		      double exponent,double width,double p6,double p7,
		      double p8,double p9,double p10,double p11,
		      double p12,double p13,double p14,double p15)
{  
  
  if (spin!=0.5) cutoff /= pow(2.0*spin,exponent);
  cutoff *=cutoff;
  cutoff *=cutoff;
  double mw2=mass*width; 
  mw2 *= mw2;
  double mandelmm = (mandel - mass*mass)*(mandel - mass*mass);
  return pow( mw2/ (mandelmm + mw2) ,(spin-0.5)) * exp( -1.0 * mandelmm / cutoff);
  
  //return strong_lorentz( mandel, mass, cutoff, spin, exponent, width, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
    //  * strong_gaussian( mandel, mass, cutoff, spin, exponent, width, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15);
}
