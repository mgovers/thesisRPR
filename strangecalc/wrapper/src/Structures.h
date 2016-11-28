/*!
 * \file Structures.h
 * \ingroup wrapper
 *  
 * This file declares all user-defined structures
 * and a number of physical and user-defined constants
 *
 * \author Tamara Corthals
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Stijn Janssen
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

#include "FormFactor.h"
#include <math.h>

#ifndef STRUCTURES_H
#define STRUCTURES_H



////////////////////////////////////
// A number of physical constants //
////////////////////////////////////

#define PI 3.14159265

/*! e= sqrt(4 pi /137) = sqrt(4 pi * alpha_e) */
#define ELEC  0.302861901

/*! Kappa_{K*0}->K^0} = KSP2KS0 * Kappa_{K*+-> K^+} 
 *
 * The sign is taken from Singer&Miller,
 * Phys. Rev. D 33, 141 - 158 (1986). 
 * The sign is confirmed by the BetheMeson
 * code by Tim Van Cauteren.
 */
#define KSP2KS0 -1.52

/*! Kappa_{K1(1270)^0}->K^0} = K1P2K10 * Kappa_{K1(1270)^+-> K^+}.
 *
 * This number was calculated with the BetheMeson code by Tim Van Cauteren. 
 */
#define K1P2K10 -0.97

/*! Kappa_{K1(1400)^0}->K^0} = K1400P2K14000 * Kappa_{K1(1400)^+-> K^+}.
 *
 * This number was calculated with the BetheMeson code by Pieter Vancraeyveld.
 */
#define K1400P2K14000 2.5277

/*! Kappa_{K_1410^0}->K^0} = B3P2B30 * Kappa_{K_1410^+-> K^+}.
 *
 * This number was calculated with the BetheMeson code by Tim Van Cauteren. 
 */
#define B4toB3 -1.1592

/* Masses */
#define M_P 938.272       //!< proton mass

#define M_N 939.565       //!< neutron mass

#define M_KP 493.677      //!< kaon^+ mass

#define M_K0 497.672      //!< kaon^0 mass

#define M_L 1115.683      //!< lambda mass

#define M_S0 1192.642     //!< sigma^0 mass

#define M_SP 1189.37      //!< sigma^+ mass

#define M_SM 1197.449     //!< sigma^- mass

#define M_RHO 775.5       //!< rho mass

#define M_OMEGA 782.65    //!< omega mass

#define M_PHI 1019.46     //!< phi mass

#define M_PIP 139.57      //!< pi^+ mass

#define M_PIM 139.57      //!< pi^- mass

#define M_PI0 134.977     //!< pi^0 mass

/* Anomalous magnetic moment */
#define KAPPA_P 1.793     //!< proton anomalous magnetic moment

#define KAPPA_N -1.913    //!< neutron anomalous magnetic moment

#define KAPPA_L -0.613    //!< lambda anomalous magnetic moment

#define KAPPA_S0 0.79     //!< sigma^0 anomalous magnetic moment

#define KAPPA_SP 1.458    //!< sigma^+ anomalous magnetic moment

#define KAPPA_SM -0.160   //!< sigma^- anomalous magnetic moment

/* Transition magnetic moment */
#define KAPPA_S0L 1.61    //!< sigma^0 - lambda transition magnetic moment


////////////////////////////////////////
// A number of user-defined constants //
////////////////////////////////////////

#define PARTICLEMAX 4 //!< max. number of particles in a class

#define CLASSMAX 21 //!< max. number of classes

#define CLASSMAX_NONCGLN 18 //!< max. number of classes in the non-CGLN decomposed calculations

#define CLASSMAX_CGLN CLASSMAX //!< max. number of classes in the CGLN decomposed calculations

#define ISOMAX 13 //!< number of different isospin channels + 1 !!!

#define MAXNRFVAR 5 //!< maximal number of free parameters for 1 particle

#define MAXEXP 20 //!< max number of sets of experimental data.

#define EXPNAME 256 //!< max number of characters in an experiment name

#define MAX_S_TRAJ 14  //!< Number of s-channel trajectories

#define MAX_LOCATION_STRING 4096 //!< Maximum length of string containing location input/output folder

#define DUMP_SIZE 100 //! Size of char dump[DUMP_SIZE] (used for temporary char-arrays)

#define BUFFER_SIZE 1024 //!< Size of char buffer[BUFFER_SIZE] (used for temporary char-arrays)

/*! \brief Underflow in the strangecalc code
 *
 * double smaller than underflow should be treated as being zero 
 */
#define STRANGEUFLOW 1.0e-9


/*!
 * Function to compare two doubles, a and b, within StrangeCalc's underflow
 */
namespace
{
  inline bool compareDoubleUFlow(double a, double b)
  {
    return fabs(1 - a/b) < STRANGEUFLOW;
  }
}


////////////////
// Datatypes //
///////////////

/*!
 * Structure corresponding to 1 particular (intermediate) particle, 
 * with its properties read from the mass.in and coupl.iso files 
 */
typedef struct
{
  char nickname[4]; //!< intermediate particle nickname, eg N12
  double mass;
  double width;
  double spin;
  double E; //!< electric charge divided by |e_electron|
  double G; //!< g_KYN, mult. by factor
  double H; //!< kappa, mult. by factor
  double I;
  double J;
  double X;
  double Y;
  double Z;
  double r_kappa_n_p; //!< ratio of trans. mag. moms k_{nN*}/k_{pN*} (1/2 N*)
  double r_kappa_1_n_p; //!< k^(1)_{nN*}/k^(1)_{pN*} (3/2 N*)
  double r_kappa_2_n_p; //!< k^(2)_{nN*}/k^(2)_{pN*} (3/2 N*)
  short regge_phase; //!< constant (0) or rotating (1) phase, positive (2) or negative (3) signature.
  short regge_phase_nonbase; //!< constant (0) or rotating (1) phase, positive (2) or negative (3) signature.
  short long_coupling; //!< for spin-1/2 resonances: longitudinal coupling or not?
  int gic; //!< use gauge-invariant couplings

  FormFactor* formfactorE; //!< Formfactor for coupling constant E
  FormFactor* formfactorG; //!< Formfactor for coupling constant G
  FormFactor* formfactorH; //!< Formfactor for coupling constant H
  FormFactor* formfactorI; //!< Formfactor for coupling constant I

} Properties;



typedef struct 
{
  short nr;
  double mass;
  short spin_shift;
  short sign;
  short rel_sign;
  double r_slope;
  double i_slope;
  double i_intercept;
} Trajectory;

/*! Structure corresponding to a class of particles (eg S,T,U,A,...) */
typedef struct
{
  int particount; //!< number of particles per class
  int trajectcount;
  Properties partic[PARTICLEMAX]; //!< array of Properties structures, 1 for each particle in the class
  Trajectory traject;

} Class;
typedef Class ClassDiagram;

/*! Structure containing the fit options, as chosen by the user */
typedef struct
{
  short polweight; //!< weight factor for polarization observables? default 1
  short isoWeight[ISOMAX]; //!< weight factor for different isospin channels
  short narrow_grid; //!< narrow or wide grid? 1 or 0
  int seed; //!< integer random number
  double starttemp; //!< start temperature (should be comparable to chisq)
  double tempdecl; //!< declination of temperature (eg 0.90)
  int iter; //!< number of iterations
  char fit_nr[50]; //!< nr to be attached to log files etc

  /*! database_info[iso][exp] is string containing name of 
   * the dataset "exp" corresponding to isospin channel "iso" */
  char database_info[ISOMAX][MAXEXP][EXPNAME];
  short nr_experiments[ISOMAX]; //!< for each isospin: number of datasets

  short photo_diffcs[ISOMAX]; //!< photoprod - for each isospin: calc diffcs?
  short photo_recpol[ISOMAX]; //!< photoprod - for each isospin: calc recoil asym ?
  short photo_phopol[ISOMAX]; //!< photoprod - for each isospin: calc photon pol ?
  short photo_tarpol[ISOMAX]; //!< photoprod - for each isospin: calc target asym ?
  short photo_totcs[ISOMAX];  //!< photoprod - for each isospin: calc totcs ?
  short photo_brpol[ISOMAX];  //!< photoprod - for each isospin: calc beam-recoil asym?
  short elec_diffcs[ISOMAX];  //!< electroprod - for each isospin: calc diffcs?

  short kaoncap_diffcs[ISOMAX]; //!< kaoncap - for each isospin: calc diffcs ?
  short kaoncap_totcs[ISOMAX];  //!< kaoncap - for each isospin: calc totcs ?
  short kaoncap_branrat[ISOMAX]; //!< kaoncap - for each isospin: calc branrat ?

  short restart; //!< restart from log file from previous fit? 1 or 0
  short randomstart; //!< random cc's? 1 or 0
  char input[MAX_LOCATION_STRING]; //!< input file for restarting from previous fit

} Fitt;

/*! Structure containing the EM form factor specifications */
typedef struct
{
  short gauge_gr; //!< Gross-Riska gauge restoration procedure ???
  short gauge_modif_ff;
  short nuc_eq_k; //!< in Regge model: p and K+ ff equal or not?
  short nuc_gk; //!< nucleon ff's with Gari-Kruempelmann parameters.
  short nuc_lo; //!< nucleon ff's with Lomon parameters.
  short nuc_gm; //!< nucleon ff's with Guidal parameters.
  short hyp_nuc;  /* L,S0 elastic ff's equal to n ff's */
  short hyp_bonn; /* "    "       "    dipoles calculated in Bonn model */
  short kplus_david;    /* K+ elastic ff acc. to David */
  short kplus_gm;       /* "  "       "  "    "  Guidal (monopole) */
  short kplus_monopole; /* "  "       " monopole with cutoff to be chosen */
  double kplus_monopole_cutoff; /* cutoff for K+ ff, to be chosen by user */
  short kstar_david;    /* F_{K,K*} acc. to David */
  short kstar_munz;     /* "        "    " Munz */
  short kstar_gm;       /* "        "    " Guidal (monopole) */
  short kstar_monopole; /* "        "    " monopole with cutoff to be chosen */
  double kstar_monopole_cutoff; /* cutoff for F_{KK*}, to be chosen by user */

  short kstar2_gm;    /* Guidal monopole */
  short kstar2_monopole; /* "        "    " monopole with cutoff to be chosen */
  double kstar2_monopole_cutoff; /* cutoff for F_{KK*}, to be chosen by user */

  double k1_monopole_cutoff; /* cutoff for F_{KK1}, to be chosen by user */
  short n_y_res_dipole; /* dipole or Bonn parameteriz. for F2_{YY*,NN*}: 
			   1 or 0 */
  double n_res_dipole_cutoff; /* cutoff for nucleon res. ff's F2_{NN*} */
  double y_res_dipole_cutoff; /* cutoff for hyp. res. ff's F2_{YY*} */
} Emfftypes;

/*! Specification of baryon polarization type in electroproduction case */
typedef struct
{
  short tarpol; /* target pol. */
  short recpol; /* recoil pol. */
  short x_pol; /* spin quantization axis */
  short y_pol; /* " */
  short z_pol; /* " */

} Elec_Bar_Pol;

/*! User options and kinematic variables in electroproduction */
typedef struct
{
  short long_coupl; /* longitudinal (electric) couplings for the
		   * spin-1/2 resonances? 1 or 0 */
  short ds_dt; /* virtual photon cs: dsigma/dt instead of 
		  dsigma_{..}/d(omega_K*) */
  short ds_du; /* dito dsigma/du instead of dsigma_{..}/d(omega_K*) */
  short dqq_dv; /* 0: d sigma /(d e' d omega_e'); 1: d sigma / d Q^2 d v*/

  Elec_Bar_Pol bar_pol; /* baryon pol. specifications */

  short baryon_pol; /* polarized baryon? 1 or 0 */
  short elec_pol; /* polarized electron? 1 or 0 */
  short fullelectro; /* full elctroproduction diff cs is computed */
  
  short beam_ener_input; /* 1 if beam energy is given, 0 if epsilon is 
			  * given instead */
  double epsilon; /* degree of transverse pol. (stijn notes 99), can be chosen
		     by user in case of full electroproduction */
  double e_beam_ener; /* the beam energy */

  short cs_convention; /* response functions front-factor convention */

  double phi_k; /* phi_K_star */
  double qq_max; /* Q^2_max (for fixed Q^2, set to 0) in GeV^2 */
  double qq_min; /* Q^2_min (for fixed Q^2, this is the Q^2 value) in GeV^2 */
  double iw_max; /* W_max (for fixed W, set to 0) in GeV */
  double iw_min; /* W_min (for fixed W, this is the W value) in GeV */
  double mandel_t_max; /* t (<0) with |t| max, ie t_min 
			* (for fixed t, irrelevant) in GeV^2 */
  double mandel_t_min; /* t (<0) with |t| min, ie t_max 
			* (for fixed t, this is the t value) in GeV^2 */
  double costheta; /* cos_theta_K_star */

  short t_var; /* cross-sec as function of a varying t: t_max and t_min as 
		* chosen by user */
  short t_fixed; /* cross-section at fixed t as chosen by user */
  short cos_var; /* cross-sec as function of a varying costheta_K* = -1..1 */
  short cos_fixed; /* cross-section at fixed costheta_K* as chosen by user */

  Emfftypes emff; /* EM formfac specifications */

} Elec;

/*! Double pol. asym. specifications in photoproduction */
typedef struct
{
  short beamtar; /* beam-target? 1 or 0 */
  short beamrec; /* beam-recoil? 1 or 0 */
  short circbeam; /* circ polarized photon beam? 1 or 0 */
  short linbeam; /* lin polarized photon beam? 1 or 0 */
  short x_barvec; /* baryon polarized along x-axis? 1 or 0 */
  short z_barvec; /* baryon polarized along z-axis? 1 or 0 */
  short beamhel; /* beam helicity */

} Double;

/*! Single pol. asym. specifications in photoproduction */
typedef struct
{
  short phopol; /* photon polarization? 1 or 0 */
  short x_phovec; /* epsilon along x */
  short y_phovec; /* epsilon along y */
  short tarpol; /* target asym? 1 or 0 */
  short recpol; /* recoil asym? 1 or 0 */

} Single;

/*! Specification of polarization type (none, single, double) in photoprod. */
typedef struct
{
  short nopol; /* unpol. cross-sec */
  short spol; /* single pol. */
  short dpol; /* double pol. */
  Single sinpol;
  Double doubpol;

} Polar;

/*! User options and kinematic variables in photoproduction */
typedef struct
{
  short ds_dt; /* dsigma/dt instead of dsigma/d(omega) */
  short ds_du; /* dsigma/du instead of dsigma/d(omega) */
  short onepoint; /* for one energy and one angle */
  short func_e_a;  /* for all energies and angles */
  short funcenergy; /* for all energies and one angle */
  double cosang; /* cos_theta_K_star */
  double minenergy; /* wlab_min (MeV) */
  double maxenergy; /* wlab_max (MeV) */
  short funcangle; /* for all angles and one energy */
  double labenergy; /* w_lab */
  double labmomentum; /* p_K_star */
  short funcmomentum; /* for all p_K_star's */
  double minmomentum; /* min p_K_star */
  double maxmomentum; /* max p_K_star */
  short funct; /* as a function of t instead of cos_theta_K_star */
  double t_min; /* abs(min value of t (<0)) = max value of -t */
  short funcu; /* as a function of u instead of cos_theta_K_star */
  double u_min; /* max value of abs(u) */
  short anglestep; /* number of steps in integration for computation
		      of sigma_tot and legendre coefficients */

} Kinetic;

/*! Cross-section type */
typedef struct
{
  short legendre; /* legendre decomp? 1 or 0 */
  short legendre_max; /* max. l value */
  short differential; /* diff cs: 1 or 0 */
  short total; /* tot cs: 1 or 0  */
  Polar pol; /* polarization specifications */
  Kinetic kin; /* kinetic variables and other user options */

} Pho;


typedef struct
{
  short stoppedkaon;
  short differential;
  short total;
  short partialtotcs;
  Polar pol;
  Kinetic kin;
} Kaon;


/*! Isospin specifications */
typedef struct
{
  short isospin; /* isospin channel (1...6) in integer form */
  short iso_channel[ISOMAX]; /* array with isospin channel numbers of processes
			      * to be included in fit or chisq calc (1...6) */
  short nr_iso_channels; /* number of isospin channels: can be >1 in fit or 
			  * chisq calc, is 1 in observable calculation */
  short iso_base; /* 1 or 2 if the base isospin channel is K+,L resp K+,S0 */

} Iso;

/*! Hadronic form factor specifications */
typedef struct
{
  double born_cutoff; /* global cutoff for Born terms */
  double res_cutoff; /* cutoff for all resonances */
  double fraction; 
  short ohta_fhat; /* choice of gauge restoration procedure */
  short haberzettl_fhat; /* " */
  short davidson_fhat;  /* " */
  short gauss_strongff; /* gauss or dipole strong ff's:  0=dipole; 1=gaussian; 2=lorentz ; 3=lorentz*gaussian */

} Formfac;

/*! Regge trajectory specifications */
typedef struct
{
  short asymp; /* asymptotic or exact behaviour? */
  short kaon_rot; /* rotating or constant phase for degenerate K traj.? */
  short kstar_rot; /* "       "  "        "     "   "          K*(892) traj.? */
  short k1_rot; /* "       "  "        "     "   "          K1(1400) traj.? */
  short kst2_rot; /* "       "  "        "     "   "          K1(1400) traj.? */
  short lam_rot; /*   "       "  "        "     "   "          L traj.? */
  short sig_rot; /*   "       "  "        "     "   "          S traj.? */
  short sigstar_rot; /* "       "  "        "     "   "        S* traj.? */
  short rot; /* rotating or constant phase for degenerate trajectories? */
  short iso3_born; /* Add electric part of Born terms in case of
		    * t-ch. reggeization in iso=3? yes/no */
  short s_modif; /* 1: replace s by (s-u)/2  to account for 
		  * low-energy behaviour; else 0 */
  short sat_traj; /* saturating traj's? 1 or 0 */
  short dual_corr; /* duality corrections? 1 or 0 */
  short s_traj_type[MAX_S_TRAJ];
  short res_on_bg; /* resonances superimposed on Regge background? 
		    * 1(yes) or 0(no) */ 
  short t_channel; /* t-channel=1 (forward angles) or u-channel=0 (backward 
		     * angles) reggeization? */
  short t_and_u_channel; /* calculation where t-(u-) channel reggeization is 
			  * used for cos>(<)0 points, 
			  * starting from the same coupl.iso file */
  short spinshift_guidal; /* recipe alpha->alpha-n of Guidal (replacement
			   * only in Gamma() and s^exp) or Stijn? 1 or 0 */
} Regge;


/*! Complete description of the kind of observable to be computed, as chosen by the user */
typedef struct
{
  short chisquare; /* chi-squared  compuation? 1 or 0 */
  short fitting; /* fitting procedure? 1 or 0 */
  Fitt fit;
  short electroprod; /* electroprod? 1 or 0 */
  Elec elec;
  short photoprod; /* photoprod.? 1 or 0 */
  Pho photo;
  short kaoncapture; /* kaon capture? 1 or 0 */
  Kaon kaoncap;
  short hadronformfac; /* hadronic form factors? 1 or 0 */
  Formfac ffac;
  short regge; /* regge model? 1 or 0 */
  Regge reg; 
  Iso iso; 
  short backgr_width; /* 0 or finite width in t and u channel? (choose 0) */
  short pvcoupling;
  short cgln; /* Calculation method: 0 use TCalculateDiagram method, 
		1 use CGLN decomposition method */
  short quant_axis_thesis; /* choice of quant. axis in comp. of asymmetries 
			    * with definite baryon spins */
  short quant_axis_drechsel; /* " */
  
  char inFolder[MAX_LOCATION_STRING]; /* Holds the location of the 
				       * input folder */
  char outFolder[MAX_LOCATION_STRING]; /* Holds the location of the 
					* output folder */
  char dataFolder[MAX_LOCATION_STRING]; /* Holds the location of the
					 * data folder */
  char modelType[MAX_LOCATION_STRING]; /* holds the model name */

} Observable;

/*! Array with numerical values of electroproduction response functions */
typedef struct
{
  double L;
  double T;
  double c_TT;
  double s_TT;
  double c_TL;
  double s_TL;
  double TT_pol;
  double c_TL_pol;
  double s_TL_pol;
} Response_func;


#endif
