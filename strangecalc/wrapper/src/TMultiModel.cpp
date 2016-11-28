/*!
 * \file TMultiModel.cpp
 * \ingroup wrapper
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \date started on November 18, 2010
 
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

#include "TMultiModel.h"
#include "strangecalc_path.h"
#include <map>
#include <vector>
#include <iostream>
using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

//////////////////////////////////////////////////////////////////////////////
//
// TMultiModel
//
// A set of TStrangeModel's for multiple isospin channels.
//
//////////////////////////////////////////////////////////////////////////////

/*!
 * \class TMultiModel
 *
 * A set of TStrangeModel's for multiple isospin channels.
 * See http://rprmodel.ugent.be/api/root/TMultiModel.html
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 *
 */

ClassImp(TMultiModel)

//_____________________________________________________________________
TMultiModel::TMultiModel(const char* name, const char* title)
  : TStrangeModel(name,title), fModels()
{
}

//_____________________________________________________________________
TMultiModel::TMultiModel(TRootIOCtor* rio)
  : TStrangeModel(rio), fModels()
{
  // ROOT I/O Constructor
}

//_____________________________________________________________________
TMultiModel::TMultiModel(const TMultiModel& rhs)
  : TStrangeModel(rhs), fModels()
{
  // Copy constructor
  map<int,TStrangeCalc*>::const_iterator it;
  for ( it=rhs.fModels.begin() ; it != rhs.fModels.end(); it++ )
    fModels[it->first] = CopyTStrangeCalc(it->second);
}

//_____________________________________________________________________
TMultiModel& TMultiModel::operator=(const TMultiModel& rhs)
{
  // Assignment
  if( this!=&rhs ) {
    TStrangeModel::operator=(rhs);
    
    // Clear all models
    map<int,TStrangeCalc*>::const_iterator it;
    for ( it=fModels.begin() ; it != fModels.end(); it++ )
      delete it->second;
    fModels.clear();

    // Copy the right-hand side models
    for ( it=rhs.fModels.begin() ; it != rhs.fModels.end(); it++ )
      fModels[it->first] = CopyTStrangeCalc(it->second);
  }

  return *this;
}

//_____________________________________________________________________
TMultiModel* TMultiModel::Clone(const char *newname) const
{
  // Virtual copy constructor
  TMultiModel *clone = new TMultiModel(*this);
  if (newname && std::strlen(newname)) clone->SetName(newname);
  return clone;
}

//_____________________________________________________________________
TMultiModel::~TMultiModel()
{
  // Destructor
  map<int,TStrangeCalc*>::const_iterator it;
  for ( it=fModels.begin() ; it != fModels.end(); it++ )
    delete it->second;
}

//_____________________________________________________________________
void TMultiModel::Streamer(TBuffer& R__b)
{
  // Stream an object of class TMultiModel.
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    TStrangeModel::Streamer(R__b);
    int fModelsSize;
    R__b >> fModelsSize;
    for ( int i=0; i<fModelsSize; ++i) {
      int iso;
      TStrangeCalc *calc = new TStrangeCalc();
      R__b >> iso;
      TStrangeModel::Streamer( R__b , 2 , calc );
      StripData(calc);
      fModels[iso] = calc;
    }
    R__b.CheckByteCount(R__s, R__c, TMultiModel::IsA());

  } else {
    R__c = R__b.WriteVersion(TMultiModel::IsA(), kTRUE);
    TStrangeModel::Streamer(R__b);
    R__b << fModels.size();
    map<int,TStrangeCalc*>::const_iterator it;
    for ( it=fModels.begin() ; it != fModels.end(); it++ ) {
      R__b << it->first;
      for (int i=0; i<ISOMAX; i++ ) it->second->datapoints[i] = new Data[DATAMAX];
      TStrangeModel::Streamer( R__b , 0 , it->second );
      StripData(it->second);
    }
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________
TStrangeCalc& TMultiModel::GetTStrangeCalc(int iso)
{
  // Retrieve the underlying TStrangeCalc object for isospin channel 'iso'.
  //
  // Attention!
  // The method TStrangeCalc::GetDatapoints() will always return an empty
  // array.
  map<int,TStrangeCalc*>::iterator model = fModels.find(iso);
  if( model==fModels.end() ) {
    cerr << "ERROR in TMultiModel::GetTStrangeCalc(int): "
	 << "TMultiModel hasn't been initialized for isospin channel "
	 << iso << "." << endl;
    exit(1);
  }

  return *(model->second);
}

//_____________________________________________________________________
const TStrangeCalc& TMultiModel::GetTStrangeCalc(int iso) const
{
  // Retrieve the underlying TStrangeCalc object for isospin channel 'iso'.
  //
  // Attention!
  // The method TStrangeCalc::GetDatapoints() will always return an empty
  // array.
  map<int,TStrangeCalc*>::const_iterator model = fModels.find(iso);
  if( model==fModels.end() ) {
    cerr << "ERROR in TMultiModel::GetTStrangeCalc(int): "
	 << "TMultiModel hasn't been initialized for isospin channel "
	 << iso << "." << endl;
    exit(1);
  }

  return *(model->second);
}

//_____________________________________________________________________
void TMultiModel::SetStrangeModel(const char *fit_spec_file, EInitVerbosity silenceOutput)
{
  // Initialize a strangecalc model with a fit_specification file. 
  //
  // Examples of fit_specification files can be found in the strangecalc-models
  // folder in the svn repository.
  //
  // Strangecalc models will be initialized for both photo- and electroproduction
  // calculations. They will be initialized for calculations in the base isospin
  // final state as well as all the related isospin final states.
  //
  // The kaon capture isospin channels are related to the photo- and electro-
  // production channels, but still essentially different. Loading a model of
  // one type will not initialize the TStrangeModel for the other type of 
  // calculation.
  //
  // If a strangecalc model has been initialized previously, it will be deleted
  // before the new model is loaded.
  //
  // By setting the setting the 3rd argument silenceOutput to kSilent (not default),
  // no initialization information will be printed to the screen. Error messages
  // are never silenced.

  // Initialize a model using TStrangeModel's method.
  // The model will be temporarily stored in TStrangeModel::fCalc.
  TStrangeModel::SetStrangeModel(fit_spec_file,silenceOutput);
  StripData(fCalc);

  // Load the model for its base channel
  SetStrangeModel(fCalc->observ.iso.iso_base,*fCalc);

  // Loop over all non-base channels and at the model when necessary
  for(int iso=3; iso<=(ISOMAX-1); ++iso) {
    if( getisospinbasechannel(iso) == fCalc->observ.iso.iso_base )
      SetStrangeModel(iso,*fCalc);
  }
  
  delete fCalc;
  fCalc=0;
}

//_____________________________________________________________________
void TMultiModel::SetStrangeModel(int iso, const char *fit_spec_file, EInitVerbosity silenceOutput)
{
  // Initialize a strangecalc model with a fit_specification file for a 
  // specific isospin channel.
  //
  // Examples of fit_specification files can be found in the strangecalc-models
  // folder in the svn repository.
  //
  // Strangecalc models will be initialized for both photo- and electroproduction
  // calculations. 
  //
  // If a strangecalc model has been initialized previously, it will be deleted
  // before the new model is loaded.
  //
  // By setting the setting the 3rd argument silenceOutput to kSilent (not default),
  // no initialization information will be printed to the screen. Error messages
  // are never silenced.
  
  // Clear the existing model if necessary
  map<int,TStrangeCalc*>::iterator model = fModels.find(iso);
  if( model!=fModels.end() ) {
    delete model->second;
  }

  // Initialize a model using TStrangeModel's method.
  // The model will be temporarily stored in TStrangeModel::fCalc.
  TStrangeModel::SetStrangeModel(fit_spec_file,silenceOutput);
  
  // Check if the model matches the isospin channel
  if( getisospinbasechannel(iso) != fCalc->observ.iso.iso_base ) {
    cerr << "WARNING in TMultiModel::SetStrangeModel(int,const char*): "
	 << "fit_specification file cannot be used to initialize "
	 << "isospin channel " << iso << ". "
	 << "This isospin channel is uninitialized for now." << endl;
    return;
  }

  // Add the model to the list
  StripData(fCalc);
  fModels[iso] = fCalc;
  fCalc = 0;
}

//_____________________________________________________________________
void TMultiModel::SetStrangeModel(int iso, const TStrangeCalc& calc)
{
  // Initialize a strangecalc model with a TStrangeCalc object for a 
  // specific isospin channel. The model will be copied.
  //
  // If a strangecalc model has been initialized previously, it will be deleted
  // before the new model is loaded.
  
  // Check if the model matches the isospin channel
  if( getisospinbasechannel(iso) != fCalc->observ.iso.iso_base ) {
    cerr << "WARNING in TMultiModel::SetStrangeModel(int,const char*): "
	 << "fit_specification file cannot be used to initialize "
	 << "isospin channel " << iso << "." << endl;
    return;
  }

  // Clear the existing model if necessary
  map<int,TStrangeCalc*>::iterator model = fModels.find(iso);
  if( model!=fModels.end() ) {
    delete model->second;
  }

  // Add the model to the list
  fModels[iso] = CopyTStrangeCalc(&calc);
}

//_____________________________________________________________________
void TMultiModel::SetStrangeModel(const TMultiModel::EPublishedModel& model, EInitVerbosity silenceOutput)
{
  // Initialize a strangecalc model corresponding to a published model.
  //
  // - kRPR2011:
  //   The RPR-2011 model that is published in arXiv:1111.6511. It is only
  //   initialized for KL channels (isospin 1 & 4).
  //
  // - kRegge2011:
  //   The background model corresponding to kRPR2011.
  //
  // - kRPR2007:
  //   The RPR model published in T.Corthals' Phd thesis.
  //   For the KL channels this corresponds to RPR2 with a missing D_13(1900).
  //   For the K+S channels it is called RPR3.
  //
  // - kRegge2007:
  //   The background model corresponding to kRegge2007.
  //
  // - kRPR2007alt:
  //   Same as kRPR2007, except the K0S production channels,
  //   i.e. isospin 3 and 5, are those corresponding to RPR-B presented in
  //   section 3.3 of P.Vancraeyveld's Phd thesis.
  //
  // - kRegge2007alt:
  //   The background model corresponding to kRegge2007alt.
  //
  // If a strangecalc model has been initialized previously, it will be deleted
  // before the new model is loaded.
  //
  // By setting the setting the 2nd argument silenceOutput to kSilent (not default),
  // no initialization information will be printed to the screen. Error messages
  // are never silenced.
  
  // Reset the TMultiModel
  map<int,TStrangeCalc*>::const_iterator it;
  for ( it=fModels.begin() ; it != fModels.end(); it++ )
    delete it->second;
  fModels.clear();

  // Set the location of the strangecalc models
  TString modelDir = MODELS_PATH;

  // Initialize the TMultiModel
  switch( model ) {

  case kRPR2011:
    SetStrangeModel(modelDir+TString("/rpr-2011/iso1+4/init/fit_specification"),silenceOutput);
    break;

  case kRPR2007:
    SetStrangeModel(modelDir+TString("/rpr-2007/iso1+4/init/fit_specification"),silenceOutput);
    SetStrangeModel(modelDir+TString("/rpr-2007/iso2+6/init/fit_specification"),silenceOutput);
    break;

  case kRPR2007alt:
    SetStrangeModel(modelDir+TString("/rpr-2007/iso1+4/init/fit_specification"),silenceOutput);
    SetStrangeModel(modelDir+TString("/rpr-2007/iso2+6/init/fit_specification"),silenceOutput);
    SetStrangeModel(3,modelDir+TString("/rpr-2007/iso3+5/init/fit_specification"),silenceOutput);
    SetStrangeModel(5,modelDir+TString("/rpr-2007/iso3+5/init/fit_specification"),silenceOutput);
    break;
   
  case kRegge2011:
    SetStrangeModel(modelDir+TString("/regge-2011/iso1+4/init/fit_specification"),silenceOutput);
    break;

  case kRegge2007:
    SetStrangeModel(modelDir+TString("/regge-2007/iso1+4/init/fit_specification"),silenceOutput);
    SetStrangeModel(modelDir+TString("/regge-2007/iso2+6/init/fit_specification"),silenceOutput);
    break;

  case kRegge2007alt:
    SetStrangeModel(modelDir+TString("/regge-2007/iso1+4/init/fit_specification"),silenceOutput);
    SetStrangeModel(modelDir+TString("/regge-2007/iso2+6/init/fit_specification"),silenceOutput);
    SetStrangeModel(3,modelDir+TString("/regge-2007/iso3+5/init/fit_specification"),silenceOutput);
    SetStrangeModel(5,modelDir+TString("/regge-2007/iso3+5/init/fit_specification"),silenceOutput);
    break;

  default:
    cerr << "WARNING in TMultiModel::SetStrangeModel(const EPublishedModel&): "
	 << "unknown model." << endl;
  }
}

//_____________________________________________________________________
void TMultiModel::StripData(TStrangeCalc* calc)
{
  for (int i=0; i<ISOMAX; i++ ) {
    delete[] calc->datapoints[i];
    calc->datapoints[i] = 0;
  }
}

//_____________________________________________________________________
TStrangeCalc* TMultiModel::CopyTStrangeCalc(const TStrangeCalc* rhs)
{

  // Create a non-const copy
  TStrangeCalc *toCopy = const_cast<TStrangeCalc*>(rhs);

  // Create datapoints for the right-hand side
  for (int i=0; i<ISOMAX; i++ ) toCopy->datapoints[i] = new Data[DATAMAX];
  
  // Create the copy
  TStrangeCalc *copy = new TStrangeCalc(*toCopy);

  // Strip the data on both objects
  StripData(toCopy);
  StripData(copy);

  return copy;
}
//_____________________________________________________________________
double TMultiModel::ChiSquared(const TString& isospins, const TString& datasets, const short polarizationWeight) const
{
  // Calculate the reduced chi-squared for any number of datasets.
  //
  // See TStrangeModel::ChiSquared for more info.

  // We need to override the original implementation because it assumes the underlying
  // TStrangeCalc objects have their data structs correctly initialized. This is not
  // true for TMultiModel objects. We strip their data to save memory. Here we will
  // temporarely allocate memory for data structs in the isospin channels we need.
  
  // Parse the list of isospin channels, and allocate data structs for each channel
  vector<TStrangeCalc*> calcs;
  TObjArray *isoList = isospins.Tokenize(",");
  for(int i=0; i<isoList->GetEntries(); ++i) {
    const int isospin = atoi( dynamic_cast<TObjString*>(isoList->At(i))->GetName() );
    calcs.push_back( const_cast<TStrangeCalc*>( &GetTStrangeCalc(isospin) ) );
    for (int i=0; i<ISOMAX; i++ ) calcs.back()->datapoints[i] = new Data[DATAMAX];
  }

  // Do the actual chi-squared calculation
  const double chisquared = TStrangeModel::ChiSquared(isospins,datasets,polarizationWeight);

  // Clean up the data memory
  for(int i=0; i<calcs.size(); ++i)
    StripData(calcs[i]);
  delete isoList;

  return chisquared;
}
