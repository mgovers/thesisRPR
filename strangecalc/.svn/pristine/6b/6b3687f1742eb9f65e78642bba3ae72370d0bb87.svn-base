/*!
 * \file TCalcInfo.cpp
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 * \author Jannes Nys <jannes.nys@ugent.be>

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

#include <cstring>
#include <iostream>
#include <cstdio>
using std::printf; using std::scanf;
#include <cstdlib>
using std::exit;
#include "TCalcInfo.h"
#include <io_specific.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TCalcInfo                                                                 //
//                                                                           //
// ROOT object that wraps the Data struct defined in fitting.h. This class   //
// does not have ROOT I/O.                                                   //
//                                                                           //
// When doing a model calculation, one obviously needs to tell strangecalc   //
// what it should calculate. This information resides in the Data* argument  //
// we pass to TStrangeCalc::GetCalcpoint(..). This structure is defined in   //
// 'fitting.h'. Historically, all C-type structs begin with a capital letter.//
// This however brings about a problem when we work in an interactive ROOT   //
// environment, because ROOT interprets such types as classes that are know  //
// to the CINT dictonary, which isn't the case. This calls for a ROOT object //
// that wraps the Data struct, hence TCalcInfo.                              //
//                                                                           //
// The user has 2 options to create TCalcInfo objects:                       //
// * with the named constructor CreateCalcInfo().                            //
// * with the explicit constructors TCalcInfo(..).                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*!
 * \class TCalcInfo
 *
 * See http://rprmodel.ugent.be/api/root/TCalcInfo.html
 *
 */

ClassImp(TCalcInfo)


//_____________________________________________________________________
TCalcInfo::TCalcInfo(EInitVerbosity verbosity)
: TObject(), fData(0), fVerbosity(verbosity)
{
  // Default constructor
  static Data emptydata = {0};  // declare a empty Data-struct once
  emptydata.photo.label   = -1;
  emptydata.elec.label	  = -1;
  emptydata.kaoncap.label = -1;
  fData = new Data(emptydata) ;
}

//_____________________________________________________________________
TCalcInfo::TCalcInfo(EReaction reaction, int isospin, const char *observable,
		     EInitVerbosity verbosity)
  : TObject(), fData(0), fVerbosity(verbosity)
{
  // Explicit constructor
  //
  // This constructor provide the user total control. This can lead to errors
  // however, so pay attention. TCalcInfo doesn't do exhaustive checks whether
  // the settings provided by the user make sense.
  // TCalcInfo needs 4 pieces of information to be complete:
  // - reaction channel: kPhoto -> photoproduction
  //                     kElectro -> electroproduction
  //                     kCapture -> kaon capture
  // - isospin channel:  the different final states have their usual labels
  // - observable:       each observable has a label (see SetObservable(TString))
  // - jacobian:         the observables can be differential to Omega, t or u.
  //                     The user should provide this info with the appropriate setter.
  //                     By default all observables are differential to Omega.
  //
  // When calculating electroproduction observables, 2 additional pieces of information
  // need to be provided:
  // - beam energy or epsilon
  // - cross section convention: see SetElectroConvention(int).

  static Data emptydata = {0};  // declare a empty Data-struct once
  emptydata.photo.label   = -1;
  emptydata.elec.label	  = -1;
  emptydata.kaoncap.label = -1;
  fData = new Data(emptydata) ;

  switch(reaction) {
  case kPhoto: fData->photo_prod = 1;; break;
  case kElectro: fData->electro_prod = 1; break;
  case kCapture: fData->kaoncapture = 1; break;
  default: Error("TCalcInfo","incorrect parameter value.");
  }

  SetIsospin(isospin);
  SetObservable(observable);
}

//_____________________________________________________________________
TCalcInfo::TCalcInfo(const Data &data)
  : TObject(), fData(new Data(data))
{
  // Cast Data -> TCalcInfo
}

//_____________________________________________________________________
TCalcInfo::TCalcInfo(const TCalcInfo& toCopy)
  : TObject(toCopy), fData(0)
{
  // Copy constructor
  fData = new Data(*toCopy.fData);
}

//_____________________________________________________________________
TCalcInfo* TCalcInfo::Clone(const char *newname) const
{
  // Virtual copy constructor
  // The argument is ignored.
  return new TCalcInfo(*this);
}

//_____________________________________________________________________
TCalcInfo::~TCalcInfo()
{
  // Destructor
  delete fData;
}

//_____________________________________________________________________
TCalcInfo& TCalcInfo::operator=(const TCalcInfo& toCopy)
{
  // Assignment
  if(this != &toCopy) // avoid self-assignment
    {
      *this->fData = *toCopy.fData;
    }

  return *this;
}

//_____________________________________________________________________
void TCalcInfo::SetIsospin(short iso)
{
  // The usual isospin channels:
  //      (1)   ===>   g + p   ->  K+ + L0
  //      (2)   ===>   g + p   ->  K+ + S0
  //      (3)   ===>   g + p   ->  K0 + S+
  //      (4)   ===>   g + n   ->  K0 + L0
  //      (5)   ===>   g + n   ->  K0 + S0
  //      (6)   ===>   g + n   ->  K+ + S-
  //      (7)   ===>   g + p   ->  P0 + p
  //      (8)   ===>   g + p   ->  P+ + n
  //      (9)   ===>   g + n   ->  P0 + n
  //      (10)  ===>   g + n   ->  P- + p
  //      (11)  ===>   K- + p  ->  g + L0
  //      (12)  ===>   K- + p  ->  g + S0

  fData->iso = iso;
  fData->elec.iso = iso;
  fData->photo.iso = iso;
  fData->kaoncap.iso = iso;
}

//_____________________________________________________________________
void TCalcInfo::SetPhotoproduction()
{
  fData->photo_prod = 1;
  fData->electro_prod = 0;
  fData->kaoncapture = 0;

  if(fVerbosity == kVerbose)
    std::cout << "WARNING in TCalcInfo::SetPhotoproduction(): "
	      << "The observable may not be set.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetElectroproduction()
{
  fData->photo_prod = 0;
  fData->electro_prod = 1;
  fData->kaoncapture = 0;

  if(fVerbosity == kVerbose)
    std::cout << "WARNING in TCalcInfo::SetElectroproduction(): "
	      << "The observable may not be set.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetKaoncapture()
{
  fData->photo_prod = 0;
  fData->electro_prod = 0;
  fData->kaoncapture = 1;

  if(fVerbosity == kVerbose)
    std::cout << "WARNING in TCalcInfo::SetKaoncapture(): "
	      << "The observable may not be set.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetObservable(TString observable)
{
  // List of observables
  //
  // Photoproduction         dcs             differential cross section
  //                         tcs             total cross section
  //                         rec             recoil(y) polarization
  //                         pho             photon polarization
  //                         tar             target(y) polarization
  //                         c_x             beam(circ)-recoil(x) polarization
  //                         c_xp            beam(circ)-recoil(x') polarization
  //                         c_z             beam(circ)-recoil(z) polarization
  //                         c_zp            beam(circ)-recoil(z') polarization
  //                         o_x             beam(lin)-recoil(x) polarization
  //                         o_xp            beam(lin)-recoil(x') polarization
  //                         o_z             beam(lin)-recoil(z) polarization
  //                         o_zp            beam(lin)-recoil(z') polarization
  //
  //                         b{i}^2          transversity amplitudes (PRC 87, 055205 (2013))
  //                         b{i}_R            {i} = 1, 2, 3, or 4 (e.g. b1_R, b4^2)
  //                         b{i}_I            ^2: norm, _R: Re part, _I: Im part
  //
  //                         a{i}^2          normalized transversity amplitudes
  //                         a{i}_R            {i} = 1, 2, 3, or 4 (e.g. a1_R, a4^2)
  //                         a{i}_I            ^2: norm, _R: Re part, _I: Im part
  //
  //                         H{i}^2          helicity amplitudes (PRC 87, 055205 (2013))
  //                         H{i}_R            {i} = 1, 2, 3, or 4 (e.g. H1_R, H4^2)
  //                         H{i}_I            ^2: norm, _R: Re part, _I: Im part
  //
  //                         h{i}^2          normalized helicity amplitudes
  //                         h{i}_R            {i} = 1, 2, 3, or 4 (e.g. h1_R, h4^2)
  //                         h{i}_I            ^2: norm, _R: Re part, _I: Im part
  //
  //                                         transversity amplitudes implementation
  //                                         --------------------------------------
  //                         diffcs          differential cross section
  //                         totcs           total cross section
  //                         S               beam asymmetry
  //                         T               target asymmetry
  //                         P               recoil asymmetry
  //                         C_x             beam(c)-recoil(x) asymmetry
  //                         C_xp            beam(c)-recoil(x') asymmetry
  //                         C_z             beam(c)-recoil(z) asymmetry
  //                         C_zp            beam(c)-recoil(z') asymmetry
  //                         O_x             beam(l)-recoil(x) asymmetry
  //                         O_xp            beam(l)-recoil(x') asymmetry
  //                         O_z             beam(l)-recoil(z) asymmetry
  //                         O_zp            beam(l)-recoil(z') asymmetry
  //                         E               beam(c)-target(z) asymmetry
  //                         F               beam(c)-target(x) asymmetry
  //                         G               beam(l)-target(z) asymmetry
  //                         H               beam(l)-target(x) asymmetry
  //                         T_x             target(x)-recoil(x) asymmetry
  //                         T_xp            target(x)-recoil(x') asymmetry
  //                         T_z             target(x)-recoil(z) asymmetry
  //                         T_zp            target(x)-recoil(z') asymmetry
  //                         L_x             target(z)-recoil(x) asymmetry
  //                         L_xp            target(z)-recoil(x') asymmetry
  //                         L_z             target(z)-recoil(z) asymmetry
  //                         L_zp            target(z)-recoil(z') asymmetry
  // Electroproduction       diff_phi        virtual photon cross section (d{sigma}/d{t}d{phi})
  //                         diff_l          longitudinal response
  //                         diff_t          transverse response
  //                         diff_t+l        dsigma_t + epsilon dsigma_l
  //                         diff_r_lt       ratio longitudinal / transverse
  //                         diff_tt_unpol   TT response (with costhkcm dependence)
  //                         diff_tl_unpol   TL response (with costhkcm dependence)
  //                         diff_tl_epol    TL' response (with sinthkcm dependence)
  //                         diff_tt_epol    TT' response
  //                         a_lu            azimuthal moment associated with beam SSA
  //                         induced_pol_y   induced polarization (recoil along y)
  //                         induced_pol_yp  induced polarization (recoil along y')
  //                         induced_pol_n   induced polarization (recoil along n)
  //                         induced_pol_yh  induced polarization (recoil along h)
  //                         transf_pol_x    transfered polarization (recoil along x)
  //                         transf_pol_xp   transfered polarization (recoil along x')
  //                         transf_pol_t    transfered polarization (recoil along t)
  //                         transf_pol_xh   tranfered polarization (recoil along xh)
  //                         transf_pol_z    transfered polarization (recoil along z)
  //                         transf_pol_zp   transfered polarization (recoil along z')
  //                         transf_pol_l    transfered polarization (recoil along l)
  //                         transf_pol_zh   transfered polarization (recoil along zh)
  //                         M{L}{Lp}{Ly}_R  real part of the hadronic part of the amplitude
  //                                         see e.g. M in thesis Lesley (eq. 2.25)
  //                                         L is the photon helicity (+, 0 ,-)
  //                                         Lp is the proton helicity (+,-)
  //                                         Ly is the hyperon helicity (+,-)
  //                         M{L}{Lp}{Ly}_I  imaginary part of the hadronic part of the amplitude
  //                         M{L}{Lp}{Ly}^2  squared norm of the hadronic part of the amplitude
  // Radiative kaon capture  dcs             differential cross section
  //                         tcs             total cross section
  //                         ptcs            partial total cross section
  //                         bran            branching ratio for stopped kaons
  //
  //
  // We have defined these reference frames:
  // (x,y,z)         z along photon
  //                 y perpendicular to electron plane
  // (x',y',z')      z' along kaon
  //                 y' perpendicular to hadron plane
  // (l,t,n)         l along hyperon
  //                 n perpendicular to hadron plane
  // (xh,yh,zh)      zh along photon
  //                 yh perpendicular to hadron plane

  if(observable.Length() > 14) {
    std::cerr << "ERROR in TCalcInfo::SetObservable(TString): "
	      << "string to long!\n";
    exit(1);
  }

  if(fData->photo_prod)
    std::strcpy(fData->photo.observable,observable);

  else if(fData->electro_prod)
    std::strcpy(fData->elec.observable,observable);

  else if(fData->kaoncapture)
    std::strcpy(fData->kaoncap.observable,observable);

}

//_____________________________________________________________________
const char *TCalcInfo::GetObservable() const
{
  if(fData->photo_prod)
    return fData->photo.observable;

  else if(fData->electro_prod)
    return fData->elec.observable;

  return fData->kaoncap.observable;
}

//_____________________________________________________________________
void TCalcInfo::SetDsDt()
{
  // BEGIN_LATEX
  // #frac{d#sigma}{dt}
  // END_LATEX
  if(fData->photo_prod) {
    fData->photo.ds_dt = 1;
    fData->photo.ds_du = 0;
  }

  if(fData->electro_prod) {
    fData->elec.ds_dt = 1;
  }

  if(fData->kaoncapture)
    if(fVerbosity == kVerbose)
      std::cout << "WARNING in TCalcInfo::SetDsDt(): "
		<< "This option is not valid for kaon capture.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetDsDu()
{
  // BEGIN_LATEX
  // #frac{d#sigma}{du}
  // END_LATEX
  if(fData->photo_prod) {
    fData->photo.ds_dt = 0;
    fData->photo.ds_du = 1;
  }

  if(fData->electro_prod)
    if(fVerbosity == kVerbose)
      std::cout << "WARNING in TCalcInfo::SetDsDu(): "
		<< "This option is not valid for electroproduction.\n";

  if(fData->kaoncapture)
    if(fVerbosity == kVerbose)
      std::cout << "WARNING in TCalcInfo::SetDsDu(): "
		<< "This option is not valid for kaon capture.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetDsDomega()
{
  // BEGIN_LATEX
  // #frac{d#sigma}{d#Omega}
  // END_LATEX
  if(fData->photo_prod) {
    fData->photo.ds_dt = 0;
    fData->photo.ds_du = 0;
  }

  else if(fData->electro_prod) {
    fData->elec.ds_dt = 0;
  }

  else if(fData->kaoncapture)
    if(fVerbosity == kVerbose)
      std::cout << "WARNING in TCalcInfo::SetDsDomega(): "
		<< "This option is not valid for kaon capture.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetElectroConvention(int convention)
{
  // This specifies the electroproduction cross section convention as
  // defined in strangecalc:
  //
  // 0 -> e ds_L, sqrt(e(e+1)) ds_TL(')
  // 1 -> e ds_L, sqrt(2e(e+1)) ds_TL(')
  // 2 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega^2 / Q^2
  // 3 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega_lab^2 / Q^2

  if(convention>=0 && convention<4)
    fData->elec.cs_convention = convention;

  else if(fVerbosity == kVerbose)
    std::cerr << "WARNING in TCalcInfo::SetElectroConvention(int): "
	      << "Convention no." << convention << " is undefined.\n";
}

//_____________________________________________________________________
void TCalcInfo::SetElectroBeamEnergy(double energy)
{
  // Beam energy in MeV.
  if(energy>0.0) {
    fData->elec.e_beam_ener = energy;
    fData->elec.beam_ener_input = 1;
  }

  else if(fVerbosity == kVerbose)
    std::cerr << "WARNING in TCalcInfo::SetElectroBeamEnergy(double): "
	      << "Beam energy should be positive!\n";
}

//_____________________________________________________________________
void TCalcInfo::SetElectroEpsilon(double epsilon)
{
  if(epsilon>0.0) {
    fData->elec.eps = epsilon;
    fData->elec.beam_ener_input = 0;
  }

  else if(fVerbosity == kVerbose)
    std::cerr << "WARNING in TCalcInfo::SetElectroEpsilon(double): "
	      << "Epsilon should be positive!\n";
}

//_____________________________________________________________________
void TCalcInfo::Help()
{
  std::cout
    << "\n******************************\n"
    << "* Instructions for TCalcInfo *\n"
    << "******************************\n\n"

    << "What:\n-----\n"
    << "ROOT object that wraps the Data struct defined in fitting.h\n\n"

    << "Motivation:\n-----------\n"
    << "When doing a model calculation, one obviously needs to tell strangecalc what it should calculate.\nThis information resides in the Data* argument we pass to TStrangeCalc::GetCalcpoint(..). This structure\nis defined in 'fitting.h'. Historically, all C-type structs begin with a capital letter. This however\nbrings about a problem when we work in an interactive ROOT environment, because ROOT interprets such\ntypes as classes that are know to the CINT dictonary, which isn't the case. This calls for a ROOT object\nthat wraps the Data struct, hence TCalcInfo.\n\n"

    << "Constructor:\n------------\n"
    << "The user has 2 options to create TCalcInfo objects:\n"
    << "* The named constructor TCalcInfo::CreateCalcInfo() prompts the user for info (interacively) and returns\n  a pointer to a TCalcInfo object. The number of observables accessible to the user is limited with this\n  named constructor.\n  CAUTION! The user owns this object and is responsible for deleting it.\n"
    << "* The explicit constructors provide the user total control. This can lead to errors however, so pay\n  attention. TCalcInfo doesn't do exhaustive checks whether the settings provided by the user make sense.\n  TCalcInfo needs 4 pieces of information to be complete:\n"
    << "  - reaction channel: 1 -> photoproduction\n"
    << "                      2 -> electroproduction\n"
    << "                      3 -> kaon capture\n"
    << "  - isospin channel:  the different final states have their usual labels\n"
    << "  - observable:       each observable has a label (see below)\n"
    << "  - jacobian:         the observables can be differential to Omega, t or u. The user should provide this\n                      info with the appropriate setter. By default all observables are differential to Omega.\n"
    << "  When calculating electroproduction observables, 2 additional pieces of information need to be provided:\n"
    << "  - beam energy or epsilon\n"
    << "  - cross section convention: 0 -> e ds_L, sqrt(e(e+1)) ds_TL('), e_L = e k_lab^2 / Q^2\n"
    << "                              1 -> e ds_L, sqrt(2e(e+1)) ds_TL('), e_L = e k_lab^2 / Q^2\n"
    << "                              2 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega^2 / Q^2\n"
    << "                              3 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega_lab^2 / Q^2\n\n"

    << "List of observables:\n--------------------\n"
    << "Photoproduction\t\tdcs\t\tdifferential cross section\n"
    << "\t\t\ttcs\t\ttotal cross section\n"
    << "\t\t\trec\t\trecoil(y) polarization\n"
    << "\t\t\tpho\t\tphoton polarization\n"
    << "\t\t\ttar\t\ttarget(y) polarization\n"
    << "\t\t\tc_x\t\tbeam(circ)-recoil(x) polarization\n"
    << "\t\t\tc_xp\t\tbeam(circ)-recoil(x') polarization\n"
    << "\t\t\tc_z\t\tbeam(circ)-recoil(z) polarization\n"
    << "\t\t\tc_zp\t\tbeam(circ)-recoil(z') polarization\n"
    << "\t\t\to_x\t\tbeam(lin)-recoil(x) polarization\n"
    << "\t\t\to_xp\t\tbeam(lin)-recoil(x') polarization\n"
    << "\t\t\to_z\t\tbeam(lin)-recoil(z) polarization\n"
    << "\t\t\to_zp\t\tbeam(lin)-recoil(z') polarization\n"
    << "\t\t\tb1^2\t\ttransverse amplitudes (defined in PRC75(2007)024002)\n"
    << "\t\t\tb2^2\t\twith - linear photon polarization\n"
    << "\t\t\tb3^2\t\t     - hyperon and nucleon polarized\n"
    << "\t\t\tb4^2\t\t       along y axis.\n"
    << "Electroproduction\tdiff_l\t\tlongitudinal response\n"
    << "\t\t\tdiff_t\t\ttransverse response\n"
    << "\t\t\tdiff_t+l\tdsigma_t + epsilon dsigma_l\n"
    << "\t\t\tdiff_tt_unpol\tTT response (with costhkcm dependence)\n"
    << "\t\t\tdiff_tl_unpol\tTL response (with costhkcm dependence)\n"
    << "\t\t\tdiff_tl_epol\tTL' response (with sinthkcm dependence)\n"
    << "\t\t\tdiff_tt_epol\tTT' response\n"
    << "\t\t\tinduced_pol_y\tinduced polarization (recoil along y)\n"
    << "\t\t\tinduced_pol_yp\tinduced polarization (recoil along y')\n"
    << "\t\t\tinduced_pol_n\tinduced polarization (recoil along n)\n"
    << "\t\t\tinduced_pol_yh\tinduced polarization (recoil along h)\n"
    << "\t\t\ttransf_pol_x\ttransfered polarization (recoil along x)\n"
    << "\t\t\ttransf_pol_xp\ttransfered polarization (recoil along x')\n"
    << "\t\t\ttransf_pol_t\ttransfered polarization (recoil along t)\n"
    << "\t\t\ttransf_pol_xh\ttranfered polarization (recoil along xh)\n"
    << "\t\t\ttransf_pol_z\ttransfered polarization (recoil along z)\n"
    << "\t\t\ttransf_pol_zp\ttransfered polarization (recoil along z')\n"
    << "\t\t\ttransf_pol_l\ttransfered polarization (recoil along l)\n"
    << "\t\t\ttransf_pol_zh\ttransfered polarization (recoil along zh)\n"
    << "\t\t\tM{L}{Lp}{Ly}{_R/_I/^2} hadronic part of the amplitude\n"
    << "\t\t\tsee e.g. M in thesis Lesley (eq. 2.25)\n"
    << "Radiative kaon capture\tdcs\t\tdifferential cross section\n"
    << "\t\t\ttcs\t\ttotal cross section\n"
    << "\t\t\tptcs\t\tpartial total cross section\n"
    << "\t\t\tbran\t\tbranching ratio for stopped kaons\n\n"

    << "Reference frames:\n-----------------\n"
    << "(x,y,z)\t\tz along photon\n"
    << "\t\ty perpendicular to electron plane\n"
    << "(x',y',z')\tz' along kaon\n"
    << "\t\ty' perpendicular to hadron plane\n"
    << "(l,t,n)\t\tl along hyperon\n"
    << "\t\tn perpendicular to hadron plane\n"
    << "(xh,yh,zh)\tzh along photon\n"
    << "\t\tyh perpendicular to hadron plane\n"
    << std::endl;
}
