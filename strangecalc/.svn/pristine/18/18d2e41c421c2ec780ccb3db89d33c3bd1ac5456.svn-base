/*!
 * \file TStrangeModel.cpp
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
using std::cout; using std::cerr; using std::endl;
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstdlib>
using std::atoi;
#include <cmath>
#include <cassert>
#include <TString.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TDirectory.h>
#include "TDataset.h"
#include "TStrangeModel.h"
#include <FourVector.h>
#include <TMatrixElement.h>
#include <FormFactorSpecification.h>
#include <FormFactorParametrization.h>
#include <io_specific.h>
#include <fitting.h>

//////////////////////////////////////////////////////////////////////////////////////
//
// TStrangeModel
//
// A TStrangeModel wraps the main strangecalc object TStrangeCalc which is defined
// in src/TStrangeCalc.h.
// 
// A TStrangeModel object can be written to ROOT files.
//
// A strangecalc model can be initialized with SetStrangeModel.
//
// This class provides a number of members to perform calculations or produce graphs:
// * GetCalcpoint
// * GetCalcpoints
// * MakeGraph
//
// After initializing a model, the current operator can be accessed with GetCurrent.
//
//////////////////////////////////////////////////////////////////////////////////////

/*!
 * \class TStrangeModel
 *
 * ROOT object representing a strangecalc model.
 * See http://rprmodel.ugent.be/api/root/TStrangeModel.html
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 * 
 */

ClassImp(TStrangeModel)

//_____________________________________________________________________
TStrangeModel::TStrangeModel(const char* name, const char* title)
: TNamed::TNamed(name,title), fCalc(0)
{ 
  // Constructor
  // Don't forget to initialize the model with SetStrangeModel(const char*)
}

//_____________________________________________________________________
TStrangeModel::TStrangeModel(TRootIOCtor* rio)
  : TNamed(), fCalc(0)
{
  // ROOT I/O Constructor
}
//_____________________________________________________________________
TStrangeModel::TStrangeModel(const TStrangeModel& toCopy)
  : TNamed(toCopy), fCalc( toCopy.fCalc ? new TStrangeCalc(*toCopy.fCalc) : 0)
{
  // Copy constructor
}

//_____________________________________________________________________
TStrangeModel& TStrangeModel::operator=(const TStrangeModel& toCopy)
{
  // Assignment
  if( this!=&toCopy ) { // avoid self-assignment
    TNamed::operator=(toCopy);
    if(toCopy.fCalc) {
      if(fCalc) *fCalc = *toCopy.fCalc;
      else fCalc = new TStrangeCalc(*toCopy.fCalc);
    } else {
      delete fCalc;
      fCalc = 0;
    }
  }

  return *this;
}

//_____________________________________________________________________
TStrangeModel* TStrangeModel::Clone(const char *newname) const
{
  // Virtual copy constructor
  TStrangeModel *clone = new TStrangeModel(*this);
  if (newname && std::strlen(newname)) clone->SetName(newname);
  return clone;
}

//_____________________________________________________________________
TStrangeModel::~TStrangeModel()
{
  // Destructor
  delete fCalc;
}

//_____________________________________________________________________
void TStrangeModel::SetStrangeModel(const char* folder, EInitVerbosity silenceOutput)
{ 
  // Initialize a strangecalc model with a fit_specification file. 
  //
  // Examples of fit_specification files can be found in the strangecalc-models
  // folder in the svn repository.
  //
  // If a strangecalc model has been initialized previously, it will be deleted
  // before the new model is loaded.
  //
  // By setting the setting the 2nd argument silenceOutput to kSilent (not default),
  // no initialization information will be printed to the screen. Error messages
  // are never silenced.
  delete fCalc;
  
  // Redirect C-style and C++-style standard output when needed
  fpos_t pos;
  int fd;
  std::streambuf *coutBuffer = cout.rdbuf();
  static std::ofstream devNull( "/dev/null" ); // Will never be closed

  if( silenceOutput==kSilent ) {
    
    fflush(stdout);
    cout.flush();
    
    fgetpos(stdout,&pos);
    fd = dup(fileno(stdout));
    freopen("/dev/null","w", stdout);

    cout.rdbuf( devNull.rdbuf() );
  }
  
  // Create TStrangeCalc object
  fCalc = new TStrangeCalc(stdout,folder,TDataset::GetDataFolder());

  // Initialize the model
  InitialiseStrangeModel();

  // Restore C-style and C++-style standard output when needed
  if( silenceOutput==kSilent ) {

    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout,&pos);

    cout.rdbuf( coutBuffer );
  }  
}

//_____________________________________________________________________
void TStrangeModel::InitialiseStrangeModel()
{ 

  // Strangecalc models are initialized for both photo- and electroprouduction
  // calculations. It is also initialized for calculations in the base isospin
  // final state as well as all the related isospin final states.
  //
  // The kaon capture isospin channels are related to the photo- and electro-
  // production channels, but still essentially different. Loading a model of
  // one type will not initialize the TStrangeModel for the other type of 
  // calculation.
  //

  // Get the Observable structure
  Observable* observ = fCalc->GetObserv();

  // Make sure that the EM form factors for the base isospin channel are initialized
  if( !observ->fit.elec_diffcs[observ->iso.iso_base])
  {
    observ->electroprod = 1;
    em_formfactor_specification(observ,fCalc->GetParticles(observ->iso.iso_base));
    observ->electroprod = 0;
  }

  // Make sure that the EM form factors are initialized for all the isospin channels
  // initialised in the fit_specification
  for (int i=0; i < observ->iso.nr_iso_channels; i++)
  {
    // Only proceed when:
    // - it isn't the base channel
    // - it was only initialised for photoproduction
    // - not kaon capture
    if( observ->iso.iso_channel[i] != observ->iso.iso_base
	&& observ->iso.iso_channel[i] <=10 
	&& !observ->fit.elec_diffcs[observ->iso.iso_channel[i]])
    {
      observ->electroprod = 1;
      em_formfactor_specification(observ,fCalc->GetParticles(observ->iso.iso_channel[i]));
      observ->electroprod = 0;
    }
  }

  // fCalc was initialized for one or more isospin channels
  // here we want to make sure it's initialized for all brother
  // and sister channels of the base isospin channel (but not for kaon capture!)
  for(int iso=1; iso<=10; ++iso) {
    if( (getisospinbasechannel(iso) == observ->iso.iso_base) &&
	// make sure it has not yet been initialised already:
	(fCalc->GetParticles(iso)[0].partic[0].mass==0.0)   ) {
      
      observ->iso.isospin = iso;
      copyparticles(fCalc->GetParticles(iso),
		    fCalc->GetParticles(observ->iso.iso_base),observ);
      modify_nicknames(fCalc->GetParticles(iso),observ);
      import_particle_properties(fCalc->GetParticles(iso),observ);
      if(iso > 3) 
	magnetic_transition_ratio(fCalc->GetParticles(iso),observ);
      
      if( observ->hadronformfac )
	strong_formfactor_specification(observ,fCalc->GetParticles(iso));
      
      observ->electroprod = 1;
      em_formfactor_specification(observ,fCalc->GetParticles(iso));
      observ->electroprod = 0;
      
    }
  }
  
}

//_____________________________________________________________________
bool TStrangeModel::DoCalculation(const TKinematics& kinematics, const TCalcInfo* calcInfo, double& result)
{
  // Calculate the observable specified in the 'calcInfo' object for the 
  // kinematics specified in the 'kinematics' object. The result is stored
  // in 'result'.
  //
  // Beware that in case of electroproduction observables the situation is
  // somewhat more complicated, because the epsilon or beam energy specified in the 'calcInfo'
  // object also has an influence on what is or isn't a physical state. This function will 
  // check whether the kinematics are compatible with the electron beam energy.
  //
  // The return value indicates whether the provided kinematics are physical:
  // - return true = physical point
  // - return false = unphysical point (result is random number)

  // Check whether elementary kinematics are physical
  if( !kinematics.IsPhysical() ) return false;

  // For electroproduction, check if kinematics are physical,
  // i.e. if |costheta_e| <= 1
  if( calcInfo->IsElectroproduction() ) {
    double costheta_e = 2.;
    if( calcInfo->fData->elec.beam_ener_input ) {
	double e_beam_ener = calcInfo->fData->elec.e_beam_ener;
	costheta_e = 1.0 - kinematics.GetQsquared() 
	  / ( 2. * e_beam_ener * (e_beam_ener - kinematics.GetWlab()));
      }
      else {
	double epsMin1 = 1./calcInfo->fData->elec.eps;
	double tmp = 2.*kinematics.GetKlab()*kinematics.GetKlab()
	  / kinematics.GetQsquared();
	costheta_e = (tmp-epsMin1+1.) / (tmp+epsMin1-1.);
      }

    if( std::fabs(costheta_e)-1. > 0. ) return false;
  }

  // Do the actual calculation
  result 
    = GetTStrangeCalc(kinematics.GetIsospin())
    .GetCalcpoint(calcInfo->fData,kinematics.GetWcm(),
		  kinematics.GetKcm(),kinematics.GetCosthkcm(),
		  kinematics.GetPk(),kinematics.GetS(),
		  kinematics.GetPhi(),kinematics.GetPhiMin());
  
  return true;
}

//_____________________________________________________________________
double* TStrangeModel::GetCalcpoints(const TKinematics& const_kinematics, const TCalcInfo* calcInfo)
{
  // Perform calculations of the observable specified in the 'calcInfo' object for the 
  // kinematics specified in the 'kinematics' object.
  //
  // Calculations in unphysical points will result in random numbers. The user is informed
  // when calculations in unphysical points are performed.
  //
  // The returned array of doubles is of size TKinematics::GetNumberOfSteps() and is
  // owned by the user.

  // Loop over all kinematic points that need to be calculated
  // --------------------------------------------------------
  TKinematics kinematics = const_kinematics;
  double *calcpoints = new double[kinematics.GetNumberOfSteps()]; // result
  int     nrOfUnphysicalPoints = 0;
  
  for(int i=0; i<kinematics.GetNumberOfSteps(); ++i) {
    kinematics.GoTo(i);
    if( !DoCalculation(kinematics,calcInfo,calcpoints[i]) )
      ++nrOfUnphysicalPoints;
  }
  
  if( nrOfUnphysicalPoints )
    cerr << "WARNING in TStrangeModel::GetCalcpoints(const TKinematics&,const TCalcInfo*): "
	      << "Requested " << nrOfUnphysicalPoints <<" calculation(s) "
	      << "outside the physical plane.\n";

  return calcpoints;
}

//_____________________________________________________________________
double* TStrangeModel::GetCalcpoints(const TDataset& dataset)
{
  // Perform calculations for the  type of observable and the kinematics relevant
  // to the data in the 'dataset'.
  //
  // The returned array of doubles is of size dataset.GetKinematics().GetNumberOfSteps()
  // and is owned by the user.
  //
  // See TDataset's documentation on how to set the kinematics of a TDataset.
  return GetCalcpoints(*dataset.GetKinematics(),dataset.GetSpecification());
}

//_____________________________________________________________________
double TStrangeModel::GetCalcpoint(const TKinematics& kinematics, const TCalcInfo* calcInfo)
{
  // Perform calculations of the observable specified in the 'calcInfo' object for the 
  // current kinematic point specified in the 'kinematics' object.
  //
  // Calculations in an unphysical point will result in a random number. The user is informed
  // when a calculation in an unphysical point is performed.

  double calcpoint = 0;

  if( !DoCalculation(kinematics,calcInfo,calcpoint) )
    cerr << "WARNING in TStrangeModel::GetCalcpoint(const TKinematics&,const TCalcInfo*): "
	 << "Requested calculation outside the physical plane.\n";
  
  return calcpoint;
}

//_____________________________________________________________________
double* TStrangeModel::GetPhysicalCalcpoints(const TKinematics& const_kinematics, const TCalcInfo* calcInfo)
{  
  // Perform calculations of the observable specified in the 'calcInfo' object for the 
  // kinematics specified in the 'kinematics' object.
  //
  // Calculations in unphysical points are skipped. For electroproduction the situation is
  // somewhat more complicated, because the epsilon or beam energy specified in the 'calcInfo'
  // object also has an influence on what is or isn't a physical state. This function will 
  // check whether the kinematics are compatible with the electron beam energy.
  // If for one of the steps, their exists an unphysical situation, a warning will be printed
  // and the returned result will be undefined.
  //
  // The returned array of doubles is of size TKinematics::GetNumberOfPhysicalSteps() and is
  // owned by the user.

  // Loop over all kinematic points that need to be calculated
  // ---------------------------------------------------------
  TKinematics kinematics = const_kinematics;
  double *calcpoints = new double[kinematics.GetNumberOfPhysicalSteps()]; // result
  int     nrOfPhysicalPoints = 0; // = physical elementary kinematics + electroproduction kinematics
  
  for(int i=0,calcpointsIndex=0; i<kinematics.GetNumberOfSteps(); ++i) {
    kinematics.GoTo(i);
    if( kinematics.IsPhysical() )
      if( DoCalculation(kinematics,calcInfo,calcpoints[calcpointsIndex++]) )
	++nrOfPhysicalPoints;
  }

  if( nrOfPhysicalPoints != kinematics.GetNumberOfPhysicalSteps() )
    cerr << "WARNING in TStrangeModel::GetPhysicalCalcpoints(const TKinematics&,const TCalcInfo*): "
	 << "Requested " 
	 << kinematics.GetNumberOfPhysicalSteps() - nrOfPhysicalPoints 
	 << " calculation(s) with unphysical electron kinematics.\n";
 
  return calcpoints;
}

//_____________________________________________________________________
double* TStrangeModel::GetPhysicalCalcpoints(const TDataset& dataset)
{  
  // Perform calculations for the  type of observable and the physical kinematics
  // relevant to the data in the 'dataset'.
  //
  // The returned array of doubles is of size dataset.GetKinematics().GetNumberOfPhysicalSteps()
  // and is owned by the user.
  //
  // See TDataset's documentation on how to set the kinematics of a TDataset.
  // See TStrangeModel::GetPhysicalCalcpoints(TKinematics&,TCalcInfo*) for more info.
  return GetPhysicalCalcpoints(*dataset.GetKinematics(),dataset.GetSpecification());
}

//_____________________________________________________________________
TGraph* TStrangeModel::MakeGraph(const TKinematics& kinematics, const TCalcInfo* calcInfo, 
				 const TString& variable,
				 double xunit, double yunit)
{
  // Create a graph of the observable specified in the 'calcInfo' object in the physical  
  // points specified in the 'kinematics' object as a function of the kinematic variable
  // 'xvariable'.
  //
  // For a list of possibilities for 'xvariable' see the documentation of 
  // TKinematics::GetVarArray(const TString&).
  // 
  // The x- and y-values are scaled by 'xunit' and 'yunit' before being drawn.
  //
  // The user owns the returned TGraph.
  if( variable.IsNull() ) {
    cerr << "ERROR in TStrangeModel::MakeGraph(kinematics,Data*,TString,double=1.0,double=1.0): "
	 << "Specify a kinematic variable for the x-axis.\n";
    return 0;
  }

  if( kinematics.GetNumberOfVariables() != 1 ) {
    cerr << "ERROR in TStrangeModel::MakeGraph(kinematics,Data*,TString,double=1.0,double=1.0): "
	 << "Only one kinematics variable can have a range.\n";
    return 0;
  }

  double *x = kinematics.GetPhysicalVarArray(variable);
  double *y = GetPhysicalCalcpoints(kinematics,calcInfo);

  if(x==0 || y==0) return 0;

  // Scale the independent variable with the 'xunit', and
  // the dependent variable with the 'yunit'.
  for (int i = 0; i < kinematics.GetNumberOfPhysicalSteps(); i++){
    x[i] /= xunit;
    y[i] /= yunit;
  }

  TGraph *graph = new TGraph(kinematics.GetNumberOfPhysicalSteps(),x,y);

  delete[] x;
  delete[] y;

  return graph;
}

//_____________________________________________________________________
TGraph* TStrangeModel::MakeGraph(const TDataset& dataset, const TString& variable,
				 double xunit, double yunit)
{
  // Create a graph for the  type of observable and the kinematics relevant
  // to the data in the 'dataset' as a function of the kinematic variable
  // 'xvariable'.
  //
  // See TStrangeModel::MakeGraph(TKinematics&,TCalcInfo*,const TString&,double,double)
  // for more info.
  //
  // See TDataset's documentation on how to set the kinematics of a TDataset.
  return MakeGraph(*dataset.GetKinematics(),dataset.GetSpecification(),variable,
		   xunit, yunit);
}

//_____________________________________________________________________
FourVector<GammaStructure> TStrangeModel::GetCurrent(int iso, 
						     const TKinematics& tk) const
{
  // This getter returns the amputated current for this reaction.
  // The current is evaluated in the centre-of-mass frame of the photon and the
  // struck nucleon. The z-axis is taken along the photon's momentum and the x-axis
  // is chosen so that the azimuthal angle of the kaon is zero. 
  //
  // iso = the isospin channel
  // tk = current will be evaluated in the current step of tk.
  //
  // Requesting a calculation in an unphysical point will result in a premature
  // exit with an error message.

  if(!tk.IsPhysical()) {
    cerr << "ERROR in TStrangeModel::GetCurrent(..): "
	 << "Requested a calculation outside the physical plane.\n"; 
    exit(1);
  }
  return GetTStrangeCalc(iso).GetCurrent(iso,tk.GetWcm(),
							tk.GetKcm(),
							tk.GetCosthkcm(),
							tk.GetPk(),
							tk.GetNucleonMass(),
							tk.GetMesonMass(),
							tk.GetHyperonMass());
}

//_____________________________________________________________________
void TStrangeModel::Help()
{
  cout
    << "\n\n**********************************\n"
    << "* Instructions for TStrangeModel *\n"
    << "**********************************\n\n"
    << "Constructor:\n";
}

//_____________________________________________________________________
void TStrangeModel::SetParticleProperty(int iso,
					const TString& particlename,
					const TString& property,
					double value)
{
  // Set the value of a particle's property. Known properties are (case sensitive):
  // - mass:          particle's mass
  // - width:         particle's total decay width (should be 0 for all except s-channel resonances)
  // - G:             EM c.c. (kappa) for diagrams S,U,A,B,C,D,E,F,G
  //                  EM c.c. (kappa_1) for diagrams H,I,J,L
  // - H:             Strong c.c. for diagrams S,T,U,A,D,E,F,G
  //                  Strong vector c.c. for diagrams B,C
  //                  EM c.c. (kappa_2) for diagrams H,I,J,L
  // - I:             Strong tensor c.c. for diagrams B,C
  //                  Strong c.c. for diagrams H,I,J,L
  // - J:             (none)
  // - X:             off-shell parameter for diagrams H,I,J,L
  // - Y:             off-shell parameter for diagrams H,I,J,L
  // - Z:             off-shell parameter for diagrams H,I,J,L
  // - r_kappa_n_p:   ratio of EM c.c. (kappa) neutron/proton for diagrams D,E,F,G
  // - r_kappa_1_n_p: ratio of EM c.c. (kappa_1) neutron/proton for diagrams H,I,J,L
  // - r_kappa_2_n_p: ratio of EM c.c. (kappa_2) neutron/proton for diagrams H,I,J,L
  //
  // - cutoff:		the hadronic formfactor's cutoff value.
  //                   Note that particlename in this case is either "res" (resonances)
  //                   or "born" (born diagrams)
  
  GetTStrangeCalc(iso).SetParticleProperty(particlename,property,value);
}

//_____________________________________________________________________
void TStrangeModel::SetModelType(const int iso, const char* modelType, const bool cgln) 
{
  // Set the model and implementation type of a TStrangeModel.
  // - modelType: 	theory, e.g. offshell or consistent
  // - cgln: 		implementation type, cgln decomposition (cgln=true) 
  // 			or the explicit diagram summation (cgln=false).
  GetTStrangeCalc(iso).SetModelType(modelType,cgln);
}

//_____________________________________________________________________
void TStrangeModel::SetModelImplementation(const int iso, const bool cgln) 
{
  // Set the implementation type of the underlying model:
  //  - cgln=true:  CGLN decomposition
  //  - cgln=false: explicit diagram summation
  GetTStrangeCalc(iso).SetModelImplementation(cgln);
}


//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b)
{
  // Stream an object of class TKinematics.
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    TNamed::Streamer(R__b);
    bool readTStrangeCalc;
    R__b >> readTStrangeCalc;
    if( readTStrangeCalc ) {
      fCalc = new TStrangeCalc();
      Streamer( R__b , R__v, fCalc);
    } else {
      fCalc = 0;
    }
    R__b.CheckByteCount(R__s, R__c, TStrangeModel::IsA());

  } else {
    R__c = R__b.WriteVersion(TStrangeModel::IsA(), kTRUE);
    TNamed::Streamer(R__b);
    R__b << (bool)(fCalc);
    if( fCalc ) Streamer( R__b , 0 , fCalc );
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, Version_t R__v, TStrangeCalc* calc)
{
  // Stream an object of type TStrangeCalc (defined in src/TStrangeCalc.h)
  if (R__b.IsReading()) {
    R__b >> calc->fDim;
    R__b >> calc->fChiSquared;
    R__b.ReadFastArray( calc->fVertex, MAXNDIM );
    for(int i=0; i<MAXNDIM; ++i) Streamer( R__b , calc->fLimits[i] );
    for(int i=0; i<MAXNDIM; ++i) R__b.ReadFastArray( calc->fName[i], calc->kNameLength );
    Streamer( R__b, calc->observ,R__v);
    calc->datapoints = new Data*[ISOMAX];
    for (int i=0; i<ISOMAX; ++i ) {
      calc->datapoints[i] = new Data[DATAMAX];
      if( R__v < 2 ) {
	// In version 1, DATAMAX was equal to 5000
	for (int j=0; j<5000; ++j ) Streamer( R__b, calc->datapoints[i][j], R__v);
	static Data emptydata = {0};
	for (int j=5000; j<DATAMAX; ++j ) calc->datapoints[i][j] = emptydata;
      }
      else {
	for (int j=0; j<DATAMAX; ++j ) Streamer( R__b, calc->datapoints[i][j],R__v);
      }
    }
    R__b.ReadFastArray( calc->datacount, ISOMAX );
    calc->particles = new ClassDiagram*[ISOMAX];
    int classMax= CLASSMAX;
    static ClassDiagram emptydiagram= {0};
    for (int i=0; i<ISOMAX; ++i ) {
      calc->particles[i] = new ClassDiagram[CLASSMAX];
      // The spin 5/2 particles were not implemented before v3
      if (R__v<3)
      {
	classMax = 14;
	for(int j=classMax; j<CLASSMAX; ++j ) calc->particles[i][j] = emptydiagram;
      }
      for(int j=0; j<classMax; ++j ) Streamer( R__b, calc->particles[i][j],R__v);
    }
    calc->printparticles = new ClassDiagram[CLASSMAX];
    for(int i=0; i<classMax; ++i ) Streamer( R__b, calc->printparticles[i],R__v);
    calc->varinfo = new Varinfo[CLASSMAX+1];
    for(int i=0; i<(classMax+1); ++i ) Streamer( R__b, calc->varinfo[i] );

    // The spin 5/2 particles were not implemented before v3
    if (R__v<3)
      {
	for(int i=classMax; i<CLASSMAX; ++i ) calc->printparticles[i] = emptydiagram;
	static Varinfo emptyvarinfo = {{{{0}}}};
	for(int i=classMax+1; i<(CLASSMAX+1); ++i ) calc->varinfo[i] = emptyvarinfo;
      }
    R__b >> calc->pho_diff_setup;
    R__b >> calc->pho_rec_setup;
    R__b >> calc->pho_tar_setup;
    R__b >> calc->elec_diff_setup;
    R__b >> calc->isospin;

  } else {

    R__b << calc->fDim;
    R__b << calc->fChiSquared;
    R__b.WriteFastArray( calc->fVertex, MAXNDIM );
    for(int i=0; i<MAXNDIM; ++i) Streamer( R__b , calc->fLimits[i] );
    for(int i=0; i<MAXNDIM; ++i) R__b.WriteFastArray( calc->fName[i], calc->kNameLength );
    Streamer( R__b, calc->observ,kVersion);
    for (int i=0; i<ISOMAX; ++i )
      for (int j=0; j<DATAMAX; ++j ) Streamer( R__b, calc->datapoints[i][j],kVersion );
    R__b.WriteFastArray( calc->datacount, ISOMAX );
    for (int i=0; i<ISOMAX; ++i )
      for(int j=0; j<CLASSMAX; ++j ) Streamer( R__b, calc->particles[i][j] ,kVersion);
    for(int i=0; i<CLASSMAX; ++i ) Streamer( R__b, calc->printparticles[i] ,kVersion);
    for(int i=0; i<(CLASSMAX+1); ++i ) Streamer( R__b, calc->varinfo[i] );
    R__b << calc->pho_diff_setup;
    R__b << calc->pho_rec_setup;
    R__b << calc->pho_tar_setup;
    R__b << calc->elec_diff_setup;
    R__b << calc->isospin;
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, Limits& limits)
{
  // Stream an object of type Limits (defined in strangecalc/fitting.h)
  if (R__b.IsReading()) {
    R__b >> limits.up;
    R__b >> limits.low;
    R__b >> limits.bound;

  } else {
    R__b << limits.up;
    R__b << limits.low;
    R__b << limits.bound;
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, Observable& obs, Version_t R__v)
{
  // Stream an object of type Observable (defined in strangecalc/classes/Structures.h)
  if (R__b.IsReading()) {
    R__b >> obs.chisquare;
    R__b >> obs.fitting;
    R__b >> obs.fit.polweight;
    R__b.ReadFastArray( obs.fit.isoWeight, ISOMAX );
    R__b >> obs.fit.narrow_grid;
    R__b >> obs.fit.seed;
    R__b >> obs.fit.starttemp;
    R__b >> obs.fit.tempdecl;
    R__b >> obs.fit.iter;
    R__b.ReadFastArray( obs.fit.fit_nr, 50 );
    for(int i=0; i<ISOMAX; ++i)
      for(int j=0; j<MAXEXP; ++j)
	R__b.ReadFastArray( obs.fit.database_info[i][j], EXPNAME );
    R__b.ReadFastArray( obs.fit.nr_experiments, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_diffcs, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_recpol, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_phopol, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_tarpol, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_totcs, ISOMAX );
    R__b.ReadFastArray( obs.fit.photo_brpol, ISOMAX );
    R__b.ReadFastArray( obs.fit.elec_diffcs, ISOMAX );
    R__b.ReadFastArray( obs.fit.kaoncap_diffcs, ISOMAX );
    R__b.ReadFastArray( obs.fit.kaoncap_totcs, ISOMAX );
    R__b.ReadFastArray( obs.fit.kaoncap_branrat, ISOMAX );
    R__b >> obs.fit.restart;
    R__b >> obs.fit.randomstart;
    R__b.ReadFastArray( obs.fit.input, MAX_LOCATION_STRING );
    R__b >> obs.electroprod;
    R__b >> obs.elec.long_coupl;
    R__b >> obs.elec.ds_dt;
    R__b >> obs.elec.ds_du;
    R__b >> obs.elec.dqq_dv;
    R__b >> obs.elec.bar_pol.tarpol;
    R__b >> obs.elec.bar_pol.recpol;
    R__b >> obs.elec.bar_pol.x_pol;
    R__b >> obs.elec.bar_pol.y_pol;
    R__b >> obs.elec.bar_pol.z_pol;
    R__b >> obs.elec.baryon_pol;
    R__b >> obs.elec.elec_pol;
    R__b >> obs.elec.fullelectro;
    R__b >> obs.elec.beam_ener_input;
    R__b >> obs.elec.epsilon;
    R__b >> obs.elec.e_beam_ener;
    R__b >> obs.elec.cs_convention;
    R__b >> obs.elec.phi_k;
    R__b >> obs.elec.qq_max;
    R__b >> obs.elec.qq_min;
    R__b >> obs.elec.iw_max;
    R__b >> obs.elec.iw_min;
    R__b >> obs.elec.mandel_t_max;
    R__b >> obs.elec.mandel_t_min;
    R__b >> obs.elec.costheta;
    R__b >> obs.elec.t_var;
    R__b >> obs.elec.t_fixed;
    R__b >> obs.elec.cos_var;
    R__b >> obs.elec.cos_fixed;
    R__b >> obs.elec.emff.gauge_gr;
    R__b >> obs.elec.emff.gauge_modif_ff;
    R__b >> obs.elec.emff.nuc_eq_k;
    R__b >> obs.elec.emff.nuc_gk;
    R__b >> obs.elec.emff.nuc_lo;
    R__b >> obs.elec.emff.nuc_gm;
    R__b >> obs.elec.emff.hyp_nuc;
    R__b >> obs.elec.emff.hyp_bonn;
    R__b >> obs.elec.emff.kplus_david;
    R__b >> obs.elec.emff.kplus_gm;
    R__b >> obs.elec.emff.kplus_monopole;
    R__b >> obs.elec.emff.kplus_monopole_cutoff;
    R__b >> obs.elec.emff.kstar_david;
    R__b >> obs.elec.emff.kstar_munz;
    R__b >> obs.elec.emff.kstar_gm;
    R__b >> obs.elec.emff.kstar_monopole;
    R__b >> obs.elec.emff.kstar_monopole_cutoff;
    R__b >> obs.elec.emff.kstar2_gm;
    R__b >> obs.elec.emff.kstar2_monopole;
    R__b >> obs.elec.emff.kstar2_monopole_cutoff;
    R__b >> obs.elec.emff.k1_monopole_cutoff;
    R__b >> obs.elec.emff.n_y_res_dipole;
    R__b >> obs.elec.emff.n_res_dipole_cutoff;
    R__b >> obs.elec.emff.y_res_dipole_cutoff;
    R__b >> obs.photoprod;
    R__b >> obs.photo.legendre;
    R__b >> obs.photo.legendre_max;
    R__b >> obs.photo.differential;
    R__b >> obs.photo.total;
    R__b >> obs.photo.pol.nopol;
    R__b >> obs.photo.pol.spol;
    R__b >> obs.photo.pol.dpol;
    R__b >> obs.photo.pol.sinpol.phopol;
    R__b >> obs.photo.pol.sinpol.x_phovec;
    R__b >> obs.photo.pol.sinpol.y_phovec;
    R__b >> obs.photo.pol.sinpol.tarpol;
    R__b >> obs.photo.pol.sinpol.recpol;
    R__b >> obs.photo.pol.doubpol.beamtar;
    R__b >> obs.photo.pol.doubpol.beamrec;
    R__b >> obs.photo.pol.doubpol.circbeam;
    R__b >> obs.photo.pol.doubpol.linbeam;
    R__b >> obs.photo.pol.doubpol.x_barvec;
    R__b >> obs.photo.pol.doubpol.z_barvec;
    R__b >> obs.photo.pol.doubpol.beamhel;
    R__b >> obs.photo.kin.ds_dt;
    R__b >> obs.photo.kin.ds_du;
    R__b >> obs.photo.kin.onepoint;
    R__b >> obs.photo.kin.func_e_a;
    R__b >> obs.photo.kin.funcenergy;
    R__b >> obs.photo.kin.cosang;
    R__b >> obs.photo.kin.minenergy;
    R__b >> obs.photo.kin.maxenergy;
    R__b >> obs.photo.kin.funcangle;
    R__b >> obs.photo.kin.labenergy;
    R__b >> obs.photo.kin.labmomentum;
    R__b >> obs.photo.kin.funcmomentum;
    R__b >> obs.photo.kin.minmomentum;
    R__b >> obs.photo.kin.maxmomentum;
    R__b >> obs.photo.kin.funct;
    R__b >> obs.photo.kin.t_min;
    R__b >> obs.photo.kin.funcu;
    R__b >> obs.photo.kin.u_min;
    R__b >> obs.photo.kin.anglestep;
    R__b >> obs.kaoncapture;
    R__b >> obs.kaoncap.stoppedkaon;
    R__b >> obs.kaoncap.differential;
    R__b >> obs.kaoncap.total;
    R__b >> obs.kaoncap.partialtotcs;
    R__b >> obs.kaoncap.pol.nopol;
    R__b >> obs.kaoncap.pol.spol;
    R__b >> obs.kaoncap.pol.dpol;
    R__b >> obs.kaoncap.pol.sinpol.phopol;
    R__b >> obs.kaoncap.pol.sinpol.x_phovec;
    R__b >> obs.kaoncap.pol.sinpol.y_phovec;
    R__b >> obs.kaoncap.pol.sinpol.tarpol;
    R__b >> obs.kaoncap.pol.sinpol.recpol;
    R__b >> obs.kaoncap.pol.doubpol.beamtar;
    R__b >> obs.kaoncap.pol.doubpol.beamrec;
    R__b >> obs.kaoncap.pol.doubpol.circbeam;
    R__b >> obs.kaoncap.pol.doubpol.linbeam;
    R__b >> obs.kaoncap.pol.doubpol.x_barvec;
    R__b >> obs.kaoncap.pol.doubpol.z_barvec;
    R__b >> obs.kaoncap.pol.doubpol.beamhel;
    R__b >> obs.kaoncap.kin.ds_dt;
    R__b >> obs.kaoncap.kin.ds_du;
    R__b >> obs.kaoncap.kin.onepoint;
    R__b >> obs.kaoncap.kin.func_e_a;
    R__b >> obs.kaoncap.kin.funcenergy;
    R__b >> obs.kaoncap.kin.cosang;
    R__b >> obs.kaoncap.kin.minenergy;
    R__b >> obs.kaoncap.kin.maxenergy;
    R__b >> obs.kaoncap.kin.funcangle;
    R__b >> obs.kaoncap.kin.labenergy;
    R__b >> obs.kaoncap.kin.labmomentum;
    R__b >> obs.kaoncap.kin.funcmomentum;
    R__b >> obs.kaoncap.kin.minmomentum;
    R__b >> obs.kaoncap.kin.maxmomentum;
    R__b >> obs.kaoncap.kin.funct;
    R__b >> obs.kaoncap.kin.t_min;
    R__b >> obs.kaoncap.kin.funcu;
    R__b >> obs.kaoncap.kin.u_min;
    R__b >> obs.kaoncap.kin.anglestep;
    R__b >> obs.hadronformfac;
    R__b >> obs.ffac.born_cutoff;
    R__b >> obs.ffac.res_cutoff;
    R__b >> obs.ffac.fraction;
    R__b >> obs.ffac.ohta_fhat;
    R__b >> obs.ffac.haberzettl_fhat;
    R__b >> obs.ffac.davidson_fhat;
    R__b >> obs.ffac.gauss_strongff;
    R__b >> obs.regge;
    R__b >> obs.reg.asymp;
    R__b >> obs.reg.kaon_rot;
    R__b >> obs.reg.kstar_rot;
    R__b >> obs.reg.k1_rot;
    R__b >> obs.reg.kst2_rot;
    R__b >> obs.reg.lam_rot;
    R__b >> obs.reg.sig_rot;
    R__b >> obs.reg.sigstar_rot;
    R__b >> obs.reg.rot;
    R__b >> obs.reg.iso3_born;
    R__b >> obs.reg.s_modif;
    R__b >> obs.reg.sat_traj;
    R__b >> obs.reg.dual_corr;
    R__b.ReadFastArray( obs.reg.s_traj_type, MAX_S_TRAJ );
    R__b >> obs.reg.res_on_bg;
    R__b >> obs.reg.t_channel;
    R__b >> obs.reg.t_and_u_channel;
    R__b >> obs.reg.spinshift_guidal;
    R__b >> obs.iso.isospin;
    R__b.ReadFastArray( obs.iso.iso_channel, ISOMAX );
    R__b >> obs.iso.nr_iso_channels;
    R__b >> obs.iso.iso_base;
    R__b >> obs.backgr_width;
    R__b >> obs.pvcoupling;
    // full backwards compatibility
    if (R__v <3)
    {
      obs.cgln=0;
      short newnorm;
      R__b >> newnorm;
      if (newnorm)
	    strcpy(obs.modelType, "consistent");
      else
	    strcpy(obs.modelType, "offshell");
    }
    else
    {
      R__b >> obs.cgln;
    }
    R__b >> obs.quant_axis_thesis;
    R__b >> obs.quant_axis_drechsel;
    R__b.ReadFastArray( obs.inFolder, MAX_LOCATION_STRING );
    R__b.ReadFastArray( obs.outFolder, MAX_LOCATION_STRING );
    R__b.ReadFastArray( obs.dataFolder, MAX_LOCATION_STRING );
    if(R__v>=3)
      R__b.ReadFastArray( obs.modelType, MAX_LOCATION_STRING );
  } else {
    R__b << obs.chisquare;
    R__b << obs.fitting;
    R__b << obs.fit.polweight;
    R__b.WriteFastArray( obs.fit.isoWeight, ISOMAX );
    R__b << obs.fit.narrow_grid;
    R__b << obs.fit.seed;
    R__b << obs.fit.starttemp;
    R__b << obs.fit.tempdecl;
    R__b << obs.fit.iter;
    R__b.WriteFastArray( obs.fit.fit_nr, 50 );
    for(int i=0; i<ISOMAX; ++i)
      for(int j=0; j<MAXEXP; ++j)
	R__b.WriteFastArray( obs.fit.database_info[i][j], EXPNAME );
    R__b.WriteFastArray( obs.fit.nr_experiments, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_diffcs, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_recpol, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_phopol, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_tarpol, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_totcs, ISOMAX );
    R__b.WriteFastArray( obs.fit.photo_brpol, ISOMAX );
    R__b.WriteFastArray( obs.fit.elec_diffcs, ISOMAX );
    R__b.WriteFastArray( obs.fit.kaoncap_diffcs, ISOMAX );
    R__b.WriteFastArray( obs.fit.kaoncap_totcs, ISOMAX );
    R__b.WriteFastArray( obs.fit.kaoncap_branrat, ISOMAX );
    R__b << obs.fit.restart;
    R__b << obs.fit.randomstart;
    R__b.WriteFastArray( obs.fit.input, MAX_LOCATION_STRING );
    R__b << obs.electroprod;
    R__b << obs.elec.long_coupl;
    R__b << obs.elec.ds_dt;
    R__b << obs.elec.ds_du;
    R__b << obs.elec.dqq_dv;
    R__b << obs.elec.bar_pol.tarpol;
    R__b << obs.elec.bar_pol.recpol;
    R__b << obs.elec.bar_pol.x_pol;
    R__b << obs.elec.bar_pol.y_pol;
    R__b << obs.elec.bar_pol.z_pol;
    R__b << obs.elec.baryon_pol;
    R__b << obs.elec.elec_pol;
    R__b << obs.elec.fullelectro;
    R__b << obs.elec.beam_ener_input;
    R__b << obs.elec.epsilon;
    R__b << obs.elec.e_beam_ener;
    R__b << obs.elec.cs_convention;
    R__b << obs.elec.phi_k;
    R__b << obs.elec.qq_max;
    R__b << obs.elec.qq_min;
    R__b << obs.elec.iw_max;
    R__b << obs.elec.iw_min;
    R__b << obs.elec.mandel_t_max;
    R__b << obs.elec.mandel_t_min;
    R__b << obs.elec.costheta;
    R__b << obs.elec.t_var;
    R__b << obs.elec.t_fixed;
    R__b << obs.elec.cos_var;
    R__b << obs.elec.cos_fixed;
    R__b << obs.elec.emff.gauge_gr;
    R__b << obs.elec.emff.gauge_modif_ff;
    R__b << obs.elec.emff.nuc_eq_k;
    R__b << obs.elec.emff.nuc_gk;
    R__b << obs.elec.emff.nuc_lo;
    R__b << obs.elec.emff.nuc_gm;
    R__b << obs.elec.emff.hyp_nuc;
    R__b << obs.elec.emff.hyp_bonn;
    R__b << obs.elec.emff.kplus_david;
    R__b << obs.elec.emff.kplus_gm;
    R__b << obs.elec.emff.kplus_monopole;
    R__b << obs.elec.emff.kplus_monopole_cutoff;
    R__b << obs.elec.emff.kstar_david;
    R__b << obs.elec.emff.kstar_munz;
    R__b << obs.elec.emff.kstar_gm;
    R__b << obs.elec.emff.kstar_monopole;
    R__b << obs.elec.emff.kstar_monopole_cutoff;
    R__b << obs.elec.emff.kstar2_gm;
    R__b << obs.elec.emff.kstar2_monopole;
    R__b << obs.elec.emff.kstar2_monopole_cutoff;
    R__b << obs.elec.emff.k1_monopole_cutoff;
    R__b << obs.elec.emff.n_y_res_dipole;
    R__b << obs.elec.emff.n_res_dipole_cutoff;
    R__b << obs.elec.emff.y_res_dipole_cutoff;
    R__b << obs.photoprod;
    R__b << obs.photo.legendre;
    R__b << obs.photo.legendre_max;
    R__b << obs.photo.differential;
    R__b << obs.photo.total;
    R__b << obs.photo.pol.nopol;
    R__b << obs.photo.pol.spol;
    R__b << obs.photo.pol.dpol;
    R__b << obs.photo.pol.sinpol.phopol;
    R__b << obs.photo.pol.sinpol.x_phovec;
    R__b << obs.photo.pol.sinpol.y_phovec;
    R__b << obs.photo.pol.sinpol.tarpol;
    R__b << obs.photo.pol.sinpol.recpol;
    R__b << obs.photo.pol.doubpol.beamtar;
    R__b << obs.photo.pol.doubpol.beamrec;
    R__b << obs.photo.pol.doubpol.circbeam;
    R__b << obs.photo.pol.doubpol.linbeam;
    R__b << obs.photo.pol.doubpol.x_barvec;
    R__b << obs.photo.pol.doubpol.z_barvec;
    R__b << obs.photo.pol.doubpol.beamhel;
    R__b << obs.photo.kin.ds_dt;
    R__b << obs.photo.kin.ds_du;
    R__b << obs.photo.kin.onepoint;
    R__b << obs.photo.kin.func_e_a;
    R__b << obs.photo.kin.funcenergy;
    R__b << obs.photo.kin.cosang;
    R__b << obs.photo.kin.minenergy;
    R__b << obs.photo.kin.maxenergy;
    R__b << obs.photo.kin.funcangle;
    R__b << obs.photo.kin.labenergy;
    R__b << obs.photo.kin.labmomentum;
    R__b << obs.photo.kin.funcmomentum;
    R__b << obs.photo.kin.minmomentum;
    R__b << obs.photo.kin.maxmomentum;
    R__b << obs.photo.kin.funct;
    R__b << obs.photo.kin.t_min;
    R__b << obs.photo.kin.funcu;
    R__b << obs.photo.kin.u_min;
    R__b << obs.photo.kin.anglestep;
    R__b << obs.kaoncapture;
    R__b << obs.kaoncap.stoppedkaon;
    R__b << obs.kaoncap.differential;
    R__b << obs.kaoncap.total;
    R__b << obs.kaoncap.partialtotcs;
    R__b << obs.kaoncap.pol.nopol;
    R__b << obs.kaoncap.pol.spol;
    R__b << obs.kaoncap.pol.dpol;
    R__b << obs.kaoncap.pol.sinpol.phopol;
    R__b << obs.kaoncap.pol.sinpol.x_phovec;
    R__b << obs.kaoncap.pol.sinpol.y_phovec;
    R__b << obs.kaoncap.pol.sinpol.tarpol;
    R__b << obs.kaoncap.pol.sinpol.recpol;
    R__b << obs.kaoncap.pol.doubpol.beamtar;
    R__b << obs.kaoncap.pol.doubpol.beamrec;
    R__b << obs.kaoncap.pol.doubpol.circbeam;
    R__b << obs.kaoncap.pol.doubpol.linbeam;
    R__b << obs.kaoncap.pol.doubpol.x_barvec;
    R__b << obs.kaoncap.pol.doubpol.z_barvec;
    R__b << obs.kaoncap.pol.doubpol.beamhel;
    R__b << obs.kaoncap.kin.ds_dt;
    R__b << obs.kaoncap.kin.ds_du;
    R__b << obs.kaoncap.kin.onepoint;
    R__b << obs.kaoncap.kin.func_e_a;
    R__b << obs.kaoncap.kin.funcenergy;
    R__b << obs.kaoncap.kin.cosang;
    R__b << obs.kaoncap.kin.minenergy;
    R__b << obs.kaoncap.kin.maxenergy;
    R__b << obs.kaoncap.kin.funcangle;
    R__b << obs.kaoncap.kin.labenergy;
    R__b << obs.kaoncap.kin.labmomentum;
    R__b << obs.kaoncap.kin.funcmomentum;
    R__b << obs.kaoncap.kin.minmomentum;
    R__b << obs.kaoncap.kin.maxmomentum;
    R__b << obs.kaoncap.kin.funct;
    R__b << obs.kaoncap.kin.t_min;
    R__b << obs.kaoncap.kin.funcu;
    R__b << obs.kaoncap.kin.u_min;
    R__b << obs.kaoncap.kin.anglestep;
    R__b << obs.hadronformfac;
    R__b << obs.ffac.born_cutoff;
    R__b << obs.ffac.res_cutoff;
    R__b << obs.ffac.fraction;
    R__b << obs.ffac.ohta_fhat;
    R__b << obs.ffac.haberzettl_fhat;
    R__b << obs.ffac.davidson_fhat;
    R__b << obs.ffac.gauss_strongff;
    R__b << obs.regge;
    R__b << obs.reg.asymp;
    R__b << obs.reg.kaon_rot;
    R__b << obs.reg.kstar_rot;
    R__b << obs.reg.k1_rot;
    R__b << obs.reg.kst2_rot;
    R__b << obs.reg.lam_rot;
    R__b << obs.reg.sig_rot;
    R__b << obs.reg.sigstar_rot;
    R__b << obs.reg.rot;
    R__b << obs.reg.iso3_born;
    R__b << obs.reg.s_modif;
    R__b << obs.reg.sat_traj;
    R__b << obs.reg.dual_corr;
    R__b.WriteFastArray( obs.reg.s_traj_type, MAX_S_TRAJ );
    R__b << obs.reg.res_on_bg;
    R__b << obs.reg.t_channel;
    R__b << obs.reg.t_and_u_channel;
    R__b << obs.reg.spinshift_guidal;
    R__b << obs.iso.isospin;
    R__b.WriteFastArray( obs.iso.iso_channel, ISOMAX );
    R__b << obs.iso.nr_iso_channels;
    R__b << obs.iso.iso_base;
    R__b << obs.backgr_width;
    R__b << obs.pvcoupling;
    R__b << obs.cgln;
    R__b << obs.quant_axis_thesis;
    R__b << obs.quant_axis_drechsel;
    R__b.WriteFastArray( obs.inFolder, MAX_LOCATION_STRING );
    R__b.WriteFastArray( obs.outFolder, MAX_LOCATION_STRING );
    R__b.WriteFastArray( obs.dataFolder, MAX_LOCATION_STRING );
    R__b.WriteFastArray( obs.modelType, MAX_LOCATION_STRING );
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, Data& data,Version_t R__v)
{
  // Stream an object of type Data (defined in strangecalc/fitting.h)
  if (R__b.IsReading()) {
    R__b >> data.iso;
    R__b >> data.photo_prod;
    R__b >> data.electro_prod;
    R__b >> data.kaoncapture;
    R__b >> data.photo.iso;
    R__b.ReadFastArray( data.photo.observable, 15 );
    R__b >> data.photo.ds_dt;
    R__b >> data.photo.ds_du;
    R__b >> data.photo.emin;
    R__b >> data.photo.emax;
    R__b >> data.photo.cos;
    R__b >> data.photo.ampli;
    R__b >> data.photo.error;
    if (R__v >=3)
	R__b >> data.photo.label;
    R__b >> data.photo.tch;
    R__b >> data.elec.iso;
    R__b.ReadFastArray( data.elec.observable, 15 );
    R__b >> data.elec.ds_dt;
    R__b >> data.elec.qsquared;
    R__b >> data.elec.s;
    R__b >> data.elec.cos_ang;
    R__b >> data.elec.t;
    R__b >> data.elec.cos;
    R__b >> data.elec.ampli;
    R__b >> data.elec.error;
    if (R__v >=3)
      R__b >> data.elec.label;
    R__b >> data.elec.tch;
    R__b >> data.elec.beam_ener_input;
    R__b >> data.elec.e_beam_ener;
    R__b >> data.elec.eps;
    R__b >> data.elec.cs_convention;
    R__b >> data.kaoncap.pk;
    R__b >> data.kaoncap.err_pk;
    R__b >> data.kaoncap.cos;
    R__b >> data.kaoncap.err_cos;
    R__b >> data.kaoncap.cosmin;
    R__b >> data.kaoncap.cosmax;
    R__b >> data.kaoncap.ampli;
    R__b >> data.kaoncap.ratio;
    R__b >> data.kaoncap.error;
    if (R__v >=3)
      R__b >> data.kaoncap.label;
    R__b >> data.kaoncap.iso;
    R__b.ReadFastArray( data.kaoncap.observable, 15 );

  } else {
    R__b << data.iso;
    R__b << data.photo_prod;
    R__b << data.electro_prod;
    R__b << data.kaoncapture;
    R__b << data.photo.iso;
    R__b.WriteFastArray( data.photo.observable, 15 );
    R__b << data.photo.ds_dt;
    R__b << data.photo.ds_du;
    R__b << data.photo.emin;
    R__b << data.photo.emax;
    R__b << data.photo.cos;
    R__b << data.photo.ampli;
    R__b << data.photo.error;
    R__b << data.photo.label;
    R__b << data.photo.tch;
    R__b << data.elec.iso;
    R__b.WriteFastArray( data.elec.observable, 15 );
    R__b << data.elec.ds_dt;
    R__b << data.elec.qsquared;
    R__b << data.elec.s;
    R__b << data.elec.cos_ang;
    R__b << data.elec.t;
    R__b << data.elec.cos;
    R__b << data.elec.ampli;
    R__b << data.elec.error;
    R__b << data.elec.label;
    R__b << data.elec.tch;
    R__b << data.elec.beam_ener_input;
    R__b << data.elec.e_beam_ener;
    R__b << data.elec.eps;
    R__b << data.elec.cs_convention;
    R__b << data.kaoncap.pk;
    R__b << data.kaoncap.err_pk;
    R__b << data.kaoncap.cos;
    R__b << data.kaoncap.err_cos;
    R__b << data.kaoncap.cosmin;
    R__b << data.kaoncap.cosmax;
    R__b << data.kaoncap.ampli;
    R__b << data.kaoncap.ratio;
    R__b << data.kaoncap.error;
    R__b << data.kaoncap.label;
    R__b << data.kaoncap.iso;
    R__b.WriteFastArray( data.kaoncap.observable, 15 );
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, ClassDiagram& diagram, Version_t R__v)
{
  // Stream an object of type Class (defined in strangecalc/classes/Structures.h)
    if (R__b.IsReading()) {
      R__b >> diagram.particount;
      R__b >> diagram.trajectcount;
      for(int i=0; i<PARTICLEMAX; ++i) {
      	R__b.ReadFastArray( diagram.partic[i].nickname, 3 );
      	R__b >> diagram.partic[i].mass;
      	R__b >> diagram.partic[i].width;
      	R__b >> diagram.partic[i].E;
      	R__b >> diagram.partic[i].G;
      	R__b >> diagram.partic[i].H;
      	R__b >> diagram.partic[i].I;
      	R__b >> diagram.partic[i].J;
      	R__b >> diagram.partic[i].X;
      	R__b >> diagram.partic[i].Y;
      	R__b >> diagram.partic[i].Z;
      	R__b >> diagram.partic[i].r_kappa_n_p;
      	R__b >> diagram.partic[i].r_kappa_1_n_p;
      	R__b >> diagram.partic[i].r_kappa_2_n_p;
      	R__b >> diagram.partic[i].regge_phase;
	R__b >> diagram.partic[i].regge_phase_nonbase;
      	R__b >> diagram.partic[i].long_coupling;
      	if(R__v >=3)
	  R__b >> diagram.partic[i].gic;
	else
	  diagram.partic[i].gic=0;
	Streamer( R__b, diagram.partic[i].formfactorE, R__v,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorG, R__v,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorH, R__v,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorI, R__v,diagram.partic[i].gic);
      }
      R__b >> diagram.traject.nr;
      R__b >> diagram.traject.mass;
      R__b >> diagram.traject.spin_shift;
      R__b >> diagram.traject.sign;
      R__b >> diagram.traject.rel_sign;
      R__b >> diagram.traject.r_slope;
      R__b >> diagram.traject.i_slope;
      R__b >> diagram.traject.i_intercept;
      
  } else {
      R__b << diagram.particount;
      R__b << diagram.trajectcount;
      for(int i=0; i<PARTICLEMAX; ++i) {
      	R__b.WriteFastArray( diagram.partic[i].nickname, 3 );
      	R__b << diagram.partic[i].mass;
      	R__b << diagram.partic[i].width;
      	R__b << diagram.partic[i].E;
      	R__b << diagram.partic[i].G;
      	R__b << diagram.partic[i].H;
      	R__b << diagram.partic[i].I;
      	R__b << diagram.partic[i].J;
      	R__b << diagram.partic[i].X;
      	R__b << diagram.partic[i].Y;
      	R__b << diagram.partic[i].Z;
      	R__b << diagram.partic[i].r_kappa_n_p;
      	R__b << diagram.partic[i].r_kappa_1_n_p;
      	R__b << diagram.partic[i].r_kappa_2_n_p;
      	R__b << diagram.partic[i].regge_phase;
	R__b << diagram.partic[i].regge_phase_nonbase;
      	R__b << diagram.partic[i].long_coupling;
	R__b << diagram.partic[i].gic;
	Streamer( R__b, diagram.partic[i].formfactorE, kVersion ,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorG, kVersion ,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorH, kVersion ,diagram.partic[i].gic);
	Streamer( R__b, diagram.partic[i].formfactorI, kVersion ,diagram.partic[i].gic);
      }
      R__b << diagram.traject.nr;
      R__b << diagram.traject.mass;
      R__b << diagram.traject.spin_shift;
      R__b << diagram.traject.sign;
      R__b << diagram.traject.rel_sign;
      R__b << diagram.traject.r_slope;
      R__b << diagram.traject.i_slope;
      R__b << diagram.traject.i_intercept;
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, Varinfo& varinfo)
{
  // Stream an object of type Varinfo (defined in strangecalc/fitting.h)
  if (R__b.IsReading()) {
    for(int i=0; i<PARTICLEMAX; ++i) {
      for(int j=0; j<MAXNRFVAR; ++j) {
	R__b >> varinfo.partic[i][j].low;
	R__b >> varinfo.partic[i][j].up;
	R__b >> varinfo.partic[i][j].var;
	R__b >> varinfo.partic[i][j].bound;
      }
    }

  } else {
    for(int i=0; i<PARTICLEMAX; ++i) {
      for(int j=0; j<MAXNRFVAR; ++j) {
	R__b << varinfo.partic[i][j].low;
	R__b << varinfo.partic[i][j].up;
	R__b << varinfo.partic[i][j].var;
	R__b << varinfo.partic[i][j].bound;
      }
    }
  }
}

//_____________________________________________________________________
void TStrangeModel::Streamer(TBuffer& R__b, FormFactor*& formfactor, Version_t R__v, int gic)
{
  // [BROKEN] Stream an object of class FormFactor (defined in strangecalc/classes/FormFactor.h)

  cerr << "ERROR in TStrangeModel::Streamer(TBufer&,FormFactor*&,Version_t,int): "
       << "streaming is broken." << endl;
  assert(1==0);

  // You will encounter some pretty ugly code here. This is because of the
  // many global variables used in the context of form factors.
  // It wouldn't hurt to rewrite the whole shebang.

  // Before we read/write the formfactor's parameters, we always read/write a 
  // bit that indicates whether the formfactor exists.
  int formfactorExists=0;

  if (R__b.IsReading()) {
    R__b >> formfactorExists;
    if( formfactorExists ) {
      int i,j;
      R__b >> i;
      R__b >> j;
      ffCalculator parametrization;
      if( i < NrOfParametrizations ) {
	if(gic)
	  parametrization = gic_ff_parametrizations[i][j];
	else
	  parametrization = ff_parametrizations[i][j];
      } else {
	if( j==0 ) parametrization = strong_dipole;
	else if( j==1 ) parametrization = strong_gaussian;
      }	
      double parameters[15];
      R__b.ReadFastArray( parameters, 15 );
      formfactor = new FormFactor( parametrization,
				   parameters[0],parameters[1],parameters[2],
				   parameters[3],parameters[4],parameters[5],
				   parameters[6],parameters[7],parameters[8],
				   parameters[9],parameters[10],parameters[11],
				   parameters[12],parameters[13],parameters[14]);
    } else {
      formfactor = 0;
    }
      
    
  } else {
    if( formfactor ) formfactorExists = 1;
    
    R__b << formfactorExists;

    if( formfactorExists ) {
      // Look for the coordinates of the formfactor parametrization in the
      // global ff_parametrizations array defined in 
      // strangecalc/FormFactorParametrization.h
      int i=0,j=0; bool found=false;
      for(i=0; (i<NrOfParametrizations) && !found; ++i) {
        for(j=0; (j<NrOfCC) && !found; ++j) {
	  if(gic)
	  {
	    if(formfactor->GetParametrization() == gic_ff_parametrizations[i][j])
	      found=true;
	  }
	  else
	  {
	    if(formfactor->GetParametrization() == ff_parametrizations[i][j])
	      found=true;
	  }
	}
	if(found) { --j; break; }
      }
      if( !found ) {
	// The strong formfactors are not in the ff_parametrizations array
	// We need to check for them by hand.
	if( formfactor->GetParametrization() == strong_dipole ) {
	  i = NrOfParametrizations;
	  j = 0;
	} else if( formfactor->GetParametrization() == strong_gaussian ) {
	  i = NrOfParametrizations;
	  j = 1;
	} else {
	  cerr << "ERROR in Streamer(TBuffer&,FormFactor*): "
	       << "cannot find the formfactor parametrization.\n";
	  assert(1==0);
	}
      }	
      R__b << i;
      R__b << j;
      double parameters[15];
      for(int i=0; i<15; ++i) parameters[i] = formfactor->GetParameter(i+1);
      R__b.WriteFastArray( parameters, 15 );
    } 
  }
}

//_____________________________________________________________________
const TStrangeCalc& TStrangeModel::GetTStrangeCalc(int iso) const
{
  // Retrieve the underlying TStrangeCalc object for isospin channel 'iso'.
  if( !fCalc ) {
    cerr << "ERROR in TStrangeModel::GetTStrangeCalc(int): "
	 << "First initialize a model with "
	 << "TStrangeModel::SetStrangeModel(const char*).\n";
    exit(1);
  }

  if( getisospinbasechannel(iso) != fCalc->observ.iso.iso_base ) {
    cerr << "ERROR in TStrangeModel::GetTStrangeCalc(int): "
	 << "TStrangeModel hasn't been initialized for isospin channel "
	 << iso << "." << endl;
    exit(1);
  }

  return *fCalc;
}

//_____________________________________________________________________
TStrangeCalc& TStrangeModel::GetTStrangeCalc(int iso)
{
  // Retrieve the underlying TStrangeCalc object for isospin channel 'iso'.
  if( !fCalc ) {
    cerr << "ERROR in TStrangeModel::GetTStrangeCalc(int): "
	 << "First initialize a model with "
	 << "TStrangeModel::SetStrangeModel(const char*).\n";
    exit(1);
  }

  if( getisospinbasechannel(iso) != fCalc->observ.iso.iso_base ) {
    cerr << "ERROR in TStrangeModel::GetTStrangeCalc(int): "
	 << "TStrangeModel hasn't been initialized for isospin channel "
	 << iso << "." << endl;
    exit(1);
  }

  return *fCalc;
}

//_____________________________________________________________________
double TStrangeModel::ChiSquared(const TString& isospins, const TString& datasets, const short polarizationWeight) const
{
  // Calculate the reduced chi-squared for any number of datasets
  //
  // The first argument is a comma-separated list of isospin channels.
  // The different channels can belong to different base channels. Make
  // sure that the TStrangeModel has been initialized for all these
  // channels.
  //
  // If the second argument is empty (default), the user will be prompted
  // interactively to select a list of datasets for each isospin channel.
  // Alternatively, the datasets can be specified in the second argument.
  // For each isospin channel, the method expects a comma-separated list.
  // The lists for the different isospin channels should be separated by
  // semi-colons.
  //
  // As a third optional argument, the user can specify a weight factor 
  // for polarization observables.
  //
  // An example: ChiSquared("1,2","101,102;301,302")
  // will give datasets 101 and 102 for isospin channel 1 and datasets
  // 301 and 302 for isospin channel 2.

  // The user needs to chose datasets. He/she has 2 options:
  // * specify datasets in the second argument
  // * or interactively
  // We recycle databasespecification() defined in io_specific.*.c
  // To be able to do this, we need to declare the C-struct Observable
  // defined in Structures.h
  Observable observ = {0};
  std::strcpy(observ.dataFolder,TDataset::GetDataFolder());
  short dummy;

  // Parse the list of isospin channels
  TObjArray *isoList = isospins.Tokenize(",");
  observ.iso.nr_iso_channels = isoList->GetEntries();
  for(int i=0; i<observ.iso.nr_iso_channels; ++i)
    observ.iso.iso_channel[i] = atoi( dynamic_cast<TObjString*>(isoList->At(i))->GetName() );
  
  // Determine whether we will work interactively
  FILE* input; // input we feed to databasespecification
  FILE* output; // stream where output of databasespecification will be directed
  if(datasets.IsNull()) { 
    // interactif
    input = stdin;
    output = stdout;
  }
  else {                
    // user given
    output = fopen( "/dev/null", "a");
    TObjArray *dataList = datasets.Tokenize(";");
    if( dataList->GetEntries() != observ.iso.nr_iso_channels ) {
      cerr << "WARNING in TStrangeModel::ChiSquared(const TString&,const TString&): "
	   << "list of datasets incompatible with list of isospin channels." << endl;
      return -1.;
    }

    // create a temporary stream
    input = tmpfile();
    for(int i=0; i<observ.iso.nr_iso_channels; ++i)
      fprintf(input,"%s\n",dynamic_cast<TObjString*>(dataList->At(i))->GetName());
    fseek(input,0,SEEK_SET);

    // Clean up
    delete dataList;
  }

  // Read the datasets
  databasespecification(&observ,&dummy,&dummy,&dummy,input,output,output);
  
  // Close the input & output stream when it's a temporary one
  if(!datasets.IsNull()) {
    fclose(input);
    fclose(output);
  }
  
  // We are now ready to read the data. It will be stored in another C-struct
  // Data defined in fitting.h
  int datacount[ISOMAX];
  static Data emptydata = {0};
  Data **datapoints = new Data*[ISOMAX];
  for(int i=0; i<ISOMAX; ++i) {
    datacount[i] = 0;
    datapoints[i] = new Data[DATAMAX];
    for(int j=0; j<DATAMAX; ++j) {
      datapoints[i][j] = emptydata;
      datapoints[i][j].photo.label = -1;
      datapoints[i][j].elec.label = -1;
      datapoints[i][j].kaoncap.label = -1;
    }
  }
  int nrOfDatapoints=0;

  // We recycle getdatastructure() defined in fitting.c, to do the conversion.
  for(int i=0; i<observ.iso.nr_iso_channels; ++i) {
    int iso = observ.iso.iso_channel[i];
    observ.iso.isospin = iso;
    getdatastructure(datapoints[iso],&datacount[iso],&observ,0);
  }

  // Calculate the total chi-squared (for each isospin channel separate)
  double chisquared = 0.;
  double weightModifiedDatacountTotal=0.;
  for(int i=0; i<observ.iso.nr_iso_channels; ++i) {
    int iso = observ.iso.iso_channel[i];

    // Create a temporary TStrangeCalc and init it correctly
    TStrangeCalc calc(GetTStrangeCalc(iso));
    calc.observ.iso.nr_iso_channels = 1;
    calc.observ.iso.iso_channel[0] = iso;
    calc.observ.fit.polweight = polarizationWeight;
    calc.datacount[iso] = datacount[iso];

    for(int j=0; j<datacount[iso]; ++j) // copy all the data
      calc.datapoints[iso][j] = datapoints[iso][j];

    // Calculate the (reduced) chi-squared for this isospin channel
    calc.SetChiSquared();
    
    // Obtain the number of data points (including weight factor)
    double weightModifiedDatacount = calc.GetNDF();
    weightModifiedDatacountTotal += weightModifiedDatacount;

    // Unreduce the chi-squared and add to total
    chisquared += calc.GetChiSquared() * weightModifiedDatacount;
  }

  // Find the reduced total chi-squared
  chisquared /= weightModifiedDatacountTotal;

  // Clean up
  for(int i=ISOMAX; i>0; --i) delete[] datapoints[i-1];
  delete[] datapoints;
  delete isoList;

  return chisquared;
}
