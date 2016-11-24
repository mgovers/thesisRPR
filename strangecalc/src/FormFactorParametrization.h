/*!
 * \file FormFactorParametrization.h
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

#include <FormFactor.h>
#include <complex>


#ifndef FORMFACTORPARAMETRIZATION_H
#define FORMFACTORPARAMETRIZATION_H


/* *************************** *
 * FORMFACTOR PARAMETRIZATIONS *
 * *************************** */


/* Standard Parametrizations
 * ------------------------- */

// monopole formfactor
double monopole(double,double,double=0,double=0,
		double=0,double=0,double=0,double=0,
		double=0,double=0,double=0,double=0,
		double=0,double=0,double=0,double=0);

// monopole formfactor, zero for Q^2=0
double monopole0(double,double,double,double=0,
		 double=0,double=0,double=0,double=0,
		 double=0,double=0,double=0,double=0,
		 double=0,double=0,double=0,double=0);

// dipole formfactor
double dipole(double,double,double=0,double=0,
	      double=0,double=0,double=0,double=0,
	      double=0,double=0,double=0,double=0,
	      double=0,double=0,double=0,double=0);

// dipole formfactor, zero for Q^2=0
double dipole0(double,double,double,double=0,
	       double=0,double=0,double=0,double=0,
	       double=0,double=0,double=0,double=0,
	       double=0,double=0,double=0,double=0);

// constant formfactor
double constant(double,double,double=0,double=0,
		double=0,double=0,double=0,double=0,
		double=0,double=0,double=0,double=0,
		double=0,double=0,double=0,double=0);

// Function returns "Not Defined"
double not_defined(double=0,double=0,double=0,double=0,
	   double=0,double=0,double=0,double=0,
	   double=0,double=0,double=0,double=0,
	   double=0,double=0,double=0,double=0);



/* Gari-Kruempelmann Parametrization
 * --------------------------------- */

// Gari-Kruempelmann proton Dirac formfactor
double gariKruempelmann_proton_dirac(double,double,double,double,
				     double,double,double,double,
				     double,double,double,double,
				     double,double,double=0,double=0);

// Gari-Kruempelmann proton Pauli formfactor
double gariKruempelmann_proton_pauli(double,double,double,double,
				     double,double,double,double,
				     double,double,double,double,
				     double,double,double=0,double=0);

// Gari-Kruempelmann neutron Dirac formfactor
double gariKruempelmann_neutron_dirac(double,double,double,double,
				      double,double,double,double,
				      double,double,double,double,
				      double,double,double=0,double=0);

// Gari-Kruempelmann neutron Pauli formfactor
double gariKruempelmann_neutron_pauli(double,double,double,double,
				      double,double,double,double,
				      double,double,double,double,
				      double,double,double=0,double=0);

// Proton electromagnetic transition form factor of the VR model
double protonEMTFF(double, double, double);



/* Bonn Parametrizations
 * --------------------- */

// Bonn Lambda Dirac formfactor
double bonn_lambda_dirac(double,double,double,double,
			 double,double,double,double,
			 double=0,double=0,double=0,double=0,
			 double=0,double=0,double=0,double=0);

// Bonn Lambda Pauli formfactor
double bonn_lambda_pauli(double,double,double,double,
			 double,double,double,double,
			 double=0,double=0,double=0,double=0,
			 double=0,double=0,double=0,double=0);

// Bonn Sigma^0 Dirac formfactor
double bonn_sigma0_dirac(double,double,double,double,
			 double,double,double=0,double=0,
			 double=0,double=0,double=0,double=0,
			 double=0,double=0,double=0,double=0);

// Bonn Sigma^0 Pauli formfactor
double bonn_sigma0_pauli(double,double,double,double,
			 double,double,double=0,double=0,
			 double=0,double=0,double=0,double=0,
			 double=0,double=0,double=0,double=0);

// Bonn Lambda-Sigma^0 transition Dirac formfactor
double bonn_lambda_sigma0_dirac(double,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0);

// Bonn Lambda-Sigma^0 transition Pauli formfactor
double bonn_lambda_sigma0_pauli(double,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0,
				double=0,double=0,double=0,double=0);




/* Kaon formfactor Parametrizations
 * -------------------------------- */

// David Kaon^+ formfactor
double david_kaonPlus(double,double,double,double,
		      double=0,double=0,double=0,double=0,
		      double=0,double=0,double=0,double=0,
		      double=0,double=0,double=0,double=0);

// Bennhold Kaon^0 formfactor
double bennhold_kaon0(double,double=0,double=0,double=0,
		      double=0,double=0,double=0,double=0,
		      double=0,double=0,double=0,double=0,
		      double=0,double=0,double=0,double=0);

// Muenz Kaon^*+ formfactor
double muenz_kaonStarPlus(double,double=0,double=0,double=0,
			  double=0,double=0,double=0,double=0,
			  double=0,double=0,double=0,double=0,
			  double=0,double=0,double=0,double=0);

// Muenz Kaon^*0 formfactor
double muenz_kaonStar0(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);


/* Spin-1/2 Nucleon Resonance Parametrizations
 * ------------------------------------------- */

double bonn_P11_1710_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P11_1710_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P11_1710_neutron_dirac(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_P11_1710_neutron_pauli(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_P31_1910_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P31_1910_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S11_1535_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S11_1535_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S11_1535_neutron_dirac(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_S11_1535_neutron_pauli(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_S11_1650_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S11_1650_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S11_1650_neutron_dirac(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_S11_1650_neutron_pauli(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_S31_1900_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_S31_1900_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P11_1440_dirac(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P11_1440_pauli(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_P11_1440_neutron_dirac(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);

double bonn_P11_1440_neutron_pauli(double,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0,
				   double=0,double=0,double=0,double=0);


/* Spin-3/2 Nucleon Resonance Parametrizations
 * ------------------------------------------- */

double bonn_P13_1720_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_P13_1720_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_P13_1720_neutron_1(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double bonn_P13_1720_neutron_2(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double bonn_P13_1900_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_P13_1900_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_P33_1920_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_P33_1920_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_D13_1700_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_D13_1700_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_D33_1700_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_D33_1700_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double bonn_D13_1520_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D13_1520_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D13_1520_neutron_1(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double bonn_D13_1520_neutron_2(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1720_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1720_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1720_neutron_1(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1720_neutron_2(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1900_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_P13_1900_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_P33_1920_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_P33_1920_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_D13_1700_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_D13_1700_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_D33_1700_1(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_D33_1700_2(double,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

double gic_bonn_D13_1520_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D13_1520_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D13_1520_neutron_1(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);

double gic_bonn_D13_1520_neutron_2(double,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0,
			       double=0,double=0,double=0,double=0);


/* Spin-5/2 Nucleon Resonance Parametrizations
 * ------------------------------------------- */
 
 
double bonn_D15_1675_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D15_1675_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);
 
double bonn_F15_1680_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F15_1680_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F15_2000_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F15_2000_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D15_2200_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D15_2200_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F35_1905_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F35_1905_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D35_1930_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_D35_1930_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F35_2000_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double bonn_F35_2000_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F15_1680_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F15_1680_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D15_1675_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D15_1675_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F15_2000_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F15_2000_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D15_2200_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D15_2200_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F35_1905_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F35_1905_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D35_1930_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_D35_1930_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F35_2000_1(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);

double gic_bonn_F35_2000_2(double,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0,
			   double=0,double=0,double=0,double=0);



/* Strong formfactor Parametrizations
 * ---------------------------------- */

// Dipole formfactor for strong interaction vertex
double strong_dipole(double,double,double,double=0,
		     double=0,double=0,double=0,double=0,
		     double=0,double=0,double=0,double=0,
		     double=0,double=0,double=0,double=0);
// Gaussian formfactor for strong interaction vertex
double strong_gaussian(double,double,double,double=0.5,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);
// Lorentz formfactor for strong interaction vertex
double strong_lorentz (double,double,double,double=0.5,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);
// Lorentz*gaussian formfactor for strong interaction vertex
double strong_lorentz_gaussian (double,double,double,double=0.5,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0,
		       double=0,double=0,double=0,double=0);

/* ************************************* *
 * ARRAY OF FORM FACTOR PARAMETRIZATIONS *
 * ************************************* */

// Define an enumeration list CC for the type of coupling constants
const int NrOfCC = 7;
enum CC{ E = 0, G , H , I , X , Y , Z };

/* Define a 2-dim array of ffCalculator structures
 * The first coordinate is the parameterization type
 * the second coordinate is the coupling constant 
 *
 * If you add or remove parametrizations, make sure to make
 * the necessary changes to TStrangeModel::Streamer(TBuffer&,FormFactor*).
 */
const int NrOfParametrizations = 38;
const ffCalculator ff_parametrizations[NrOfParametrizations][NrOfCC] = 
  { { monopole, monopole, monopole, monopole,
      monopole, monopole, monopole },

    { monopole0, monopole, monopole, monopole,
      monopole, monopole, monopole },

    { dipole, dipole, dipole, dipole,
      dipole, dipole, dipole},

    { dipole0, dipole, dipole, dipole,
      dipole, dipole, dipole},

    { gariKruempelmann_proton_dirac, gariKruempelmann_proton_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { gariKruempelmann_neutron_dirac, gariKruempelmann_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { david_kaonPlus, not_defined, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bennhold_kaon0, not_defined, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, muenz_kaonStarPlus, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, muenz_kaonStar0, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_lambda_dirac, bonn_lambda_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_sigma0_dirac, bonn_sigma0_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_lambda_sigma0_dirac, bonn_lambda_sigma0_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_P11_1710_dirac, bonn_P11_1710_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_P31_1910_dirac, bonn_P31_1910_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_S11_1535_dirac, bonn_S11_1535_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_S11_1650_dirac, bonn_S11_1650_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_S31_1900_dirac, bonn_S31_1900_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_P13_1720_1, bonn_P13_1720_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_P13_1900_1, bonn_P13_1900_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_P33_1920_1, bonn_P33_1920_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_D13_1700_1, bonn_D13_1700_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_D33_1700_1, bonn_D33_1700_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_D13_1520_1, bonn_D13_1520_2,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_P11_1710_neutron_dirac, bonn_P11_1710_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_P13_1720_neutron_1, bonn_P13_1720_neutron_2,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_S11_1650_neutron_dirac, bonn_S11_1650_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_S11_1535_neutron_dirac, bonn_S11_1535_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_D13_1520_neutron_1, bonn_D13_1520_neutron_2,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_P11_1440_dirac, bonn_P11_1440_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { bonn_P11_1440_neutron_dirac, bonn_P11_1440_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_F15_1680_1, bonn_F15_1680_2,
      not_defined, not_defined, not_defined, not_defined },

    { not_defined, bonn_D15_1675_1, bonn_D15_1675_2,
      not_defined, not_defined, not_defined, not_defined},
    
    { not_defined, bonn_F15_2000_1, bonn_F15_2000_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, bonn_D15_2200_1, bonn_D15_2200_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, bonn_F35_1905_1, bonn_F35_1905_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, bonn_D35_1930_1, bonn_D35_1930_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, bonn_F35_2000_1, bonn_F35_2000_2,
      not_defined, not_defined, not_defined, not_defined}};


//! EM Formfactors for gauge-invariant interactions
const ffCalculator gic_ff_parametrizations[NrOfParametrizations][NrOfCC] = 
  { { monopole, monopole, monopole, monopole,
      monopole, monopole, monopole},

    { monopole0, monopole, monopole, monopole,
      monopole, monopole, monopole },

    { dipole, dipole, dipole, dipole,
      dipole, dipole, dipole},

    { dipole0, dipole, dipole, dipole,
      dipole, dipole, dipole},

    { gariKruempelmann_proton_dirac, gariKruempelmann_proton_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { gariKruempelmann_neutron_dirac, gariKruempelmann_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { david_kaonPlus, not_defined, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bennhold_kaon0, not_defined, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, muenz_kaonStarPlus, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, muenz_kaonStar0, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_lambda_dirac, bonn_lambda_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_sigma0_dirac, bonn_sigma0_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_lambda_sigma0_dirac, bonn_lambda_sigma0_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_P11_1710_dirac, bonn_P11_1710_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_P31_1910_dirac, bonn_P31_1910_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_S11_1535_dirac, bonn_S11_1535_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_S11_1650_dirac, bonn_S11_1650_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_S31_1900_dirac, bonn_S31_1900_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_P13_1720_1, gic_bonn_P13_1720_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_P13_1900_1, gic_bonn_P13_1900_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_P33_1920_1, gic_bonn_P33_1920_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D13_1700_1, gic_bonn_D13_1700_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D33_1700_1, gic_bonn_D33_1700_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D13_1520_1, gic_bonn_D13_1520_2,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_P11_1710_neutron_dirac, bonn_P11_1710_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_P13_1720_neutron_1, gic_bonn_P13_1720_neutron_2,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_S11_1650_neutron_dirac, bonn_S11_1650_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_S11_1535_neutron_dirac, bonn_S11_1535_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D13_1520_neutron_1, gic_bonn_D13_1520_neutron_2,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_P11_1440_dirac, bonn_P11_1440_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { bonn_P11_1440_neutron_dirac, bonn_P11_1440_neutron_pauli, not_defined,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_F15_1680_1, gic_bonn_F15_1680_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D15_1675_1, gic_bonn_D15_1675_2,
      not_defined, not_defined, not_defined, not_defined},
    
    { not_defined, gic_bonn_F15_2000_1, gic_bonn_F15_2000_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D15_2200_1, gic_bonn_D15_2200_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_F35_1905_1, gic_bonn_F35_1905_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_D35_1930_1, gic_bonn_D35_1930_2,
      not_defined, not_defined, not_defined, not_defined},

    { not_defined, gic_bonn_F35_2000_1, gic_bonn_F35_2000_2,
      not_defined, not_defined, not_defined, not_defined}};


#endif
