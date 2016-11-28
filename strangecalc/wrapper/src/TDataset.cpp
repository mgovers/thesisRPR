/*!
 * \file TDataset.cpp
 * \ingroup wrapper
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 
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
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <set>
#include <vector>
#include <TEventList.h>
#include <TDirectory.h>
#include "TDataset.h"
#include "strangecalc_path.h"
#include <io_specific.h>

//////////////////////////////////////////////////////////////////////////////////////
//
// TDataset
//
// A ROOT object that holds a strangecalc dataset.
//
// A TDataset can contain only a single dataset. A dataset is imported by calling
// ImportDataset(int,const TString&). Make sure that the location of the strangecalc
// datasets has been set correctly with TDataset::SetDataFolder(const TString&).
//
// A TDataset is a TTree. Every data point is stored as an entry in this tree 
// with 16 branches:
// * observable (char[15])
// * wlab       (double)
// * wcm        (double)
// * klab       (double)
// * kcm        (double)
// * pk         (double)
// * pklab      (double)
// * coshtkcm   (double)
// * qsquared   (double)
// * w          (double)
// * s          (double)
// * t          (double)
// * u          (double)
// * isPhysical (int)
// * amplitude  (double)
// * error      (double)
//
// A user can create subsets of the dataset, coined selections. Conceptually, selections
// represent a set of data points of the same observable and/or for one or more fixed
// kinematical variables. A particular selection can be activated with SetSelection(int).
// The existing selections can be inspected with ViewSelections(). The data points of the
// current selection can be displayed with Scan. The total dataset is given selection nr.0.
//
// Every selection has two attributes: 
// * a kinematics: the type of kinematics associated with a current selection can be set 
//                 by the user with SetKinematics(TKinematics*). Giving this TKinematics
//                 object a range is important for MakeGraph to function correctly.
// * a description: a selection is given a default description, but it can be changed by
//                  the user with SetDescription for the current selection.
//
// There are 2 ways of creating selections. Note that selections are always created with
// respect to the current selection.
// * AddSelection(const char*,TKinematics*): creates a selection where the first argument
//                                           is a string that specifies a specfic subset
//                                           of the current selection.
// * MakeSelection(const TString&): the argument should be a colon-separated list of branch
//                                  names. The function will create selections for every
//                                  set of data points were the value of these branches
//                                  is the same.
//
// The user can produce a graphic representation of the currect selection with MakeGraph.
//
//////////////////////////////////////////////////////////////////////////////////////

/*!
 * \class TDataset
 *
 * \brief ROOT object that can hold a strangecalc dataset.
 *
 * Data is internally stored as a TTree.
 * See http://rprmodel.ugent.be/api/root/TDataset.html
 *
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

ClassImp(TDataset)

//_____________________________________________________________________
TString TDataset::fgDataFolder(DATA_PATH);

//_____________________________________________________________________
TDataset::TDataset(const char* name, const char* title)
  : TTree(name,title), fSpecification(0), fIsospin(0),
    fElectro_prod(0),fPhoto_prod(0),fKaoncapture(0),fDs_dt(0),fDs_du(0),
    fBeam_ener_input(0),fE_beam_ener(0.0),fEps(0.0),fCs_convention(0), 
    fCosmin(0), fCosmax(0), fSelection(-1),
    fSelecDescription(0),fSelecEventLists(0),fSelecKinematics(0)
{
  // Constructor
  SetScanField(0);
}

//_____________________________________________________________________
TDataset::TDataset(const TDataset& toCopy)
  : fSpecification(toCopy.fSpecification), fIsospin(toCopy.fIsospin),
    fElectro_prod(toCopy.fElectro_prod), fPhoto_prod(toCopy.fPhoto_prod),
    fKaoncapture(toCopy.fKaoncapture), fDs_dt(toCopy.fDs_dt), fDs_du(toCopy.fDs_du),
    fBeam_ener_input(toCopy.fBeam_ener_input), fE_beam_ener(toCopy.fE_beam_ener),
    fEps(toCopy.fEps), fCs_convention(toCopy.fCs_convention),
    fCosmin(toCopy.fCosmin), fCosmax(toCopy.fCosmax), fSelection(toCopy.fSelection),
    fSelecDescription(0),fSelecEventLists(0),fSelecKinematics(0)
{
  // Private copy constructor. Please use Clone()
  std::cerr << "You should not have called "
	    << "TDataset::TDataset(const TDataset&)! "
	    << "Use TDataset::Clone(). Goodbye.\n";
  exit(1);
}

//_____________________________________________________________________
TDataset& TDataset::operator=(const TDataset& toCopy)
{
  // Private assignment operator. Please use Clone()
  std::cerr << "You should not have called "
	    << "TDataset::operator=(const TDataset&)! "
	    << "Use TDataset::Clone(). Goodbye.\n";
  exit(1);
}

//_____________________________________________________________________
TDataset::~TDataset()
{
  // Destructor
  delete fSpecification;
  delete fSelecDescription;
  delete fSelecEventLists;
  delete fSelecKinematics;
}

//_____________________________________________________________________
void TDataset::Streamer(TBuffer &R__b)
{
  // Stream an object of class TDataset.

  // This function was generated by rootcint with #pragma link c++ class TDataset+
  // We implement it ourselfs, because we need to call TDataset::SetSelection(int)
  // at the end of the 'read' process.
  // See ROOT user guide chpt.11 -> Streamers for more info

  // Since january 2009 there exist multiple versions of this class.
  // We need to maintain schema evolution manually.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TTree::Streamer(R__b);
      R__b >> fIsospin;
      R__b >> fElectro_prod;
      R__b >> fPhoto_prod;
      R__b >> fDs_dt;
      R__b >> fDs_du;
      R__b >> fBeam_ener_input;
      R__b >> fE_beam_ener;
      R__b >> fEps;
      R__b >> fCs_convention;
      R__b >> fSelection;
      R__b >> fSelecDescription;
      R__b >> fSelecEventLists;
      R__b >> fSelecKinematics;
    
      // Version specific code
      if( R__v > 1) {               // VERSION 2
	R__b >> fKaoncapture;
	R__b >> fCosmin;
	R__b >> fCosmax;
      }
      else if ( R__v < 2 ) {        // VERSION 1
	fKaoncapture = 0;
	fCosmin = fCosmax = 0.0;
	
	// Version 1 does not have pklab branch
	// Create the branch and set branch addresses for independent variables
	double wcm=0., costhkcm=0., qsquared=0., pklab=0.;
	TBranch *pklabBranch = Branch("pklab",&pklab,"pklab/D");
	SetBranchAddress("wcm",&wcm);
	SetBranchAddress("costhkcm",&costhkcm);
	SetBranchAddress("qsquared",&qsquared);
	
	// desactivate other branches
	SetBranchStatus("*",0);
	SetBranchStatus("wcm",1);
	SetBranchStatus("pklab",1);
	SetBranchStatus("costhkcm",1);
	SetBranchStatus("qsquared",1);

	// Create TKinematics object for calculations
	TKinematics tk("tk","kaon capture kinematics",fIsospin,
		       "wcm:costhkcm:qsquared",0.0,0.0,0.0);

	// loop over all entries and set pklab
	for(int i=0; i<GetEntries(); ++i) {
	  GetEntry(i);
	  tk.SetVar(1,wcm);
	  tk.SetVar(2,costhkcm);
	  tk.SetVar(3,qsquared);
	  pklab = tk.GetPklab();
	  pklabBranch->Fill();
	}

	// reactivate all branches and release branch addresses
	SetBranchStatus("*",1);
	ResetBranchAddresses();
      }
      
      R__b.CheckByteCount(R__s, R__c, TDataset::IsA());
      SetSelection(fSelection); // Added by PVC

   } else {
      R__c = R__b.WriteVersion(TDataset::IsA(), kTRUE);
      TTree::Streamer(R__b);
      R__b << fIsospin;
      R__b << fElectro_prod;
      R__b << fPhoto_prod;
      R__b << fDs_dt;
      R__b << fDs_du;
      R__b << fBeam_ener_input;
      R__b << fE_beam_ener;
      R__b << fEps;
      R__b << fCs_convention;
      R__b << fSelection;
      R__b << fSelecDescription;
      R__b << fSelecEventLists;
      R__b << fSelecKinematics;
      R__b << fKaoncapture;
      R__b << fCosmin;
      R__b << fCosmax;

      R__b.SetByteCount(R__c, kTRUE);
   }
}

//_____________________________________________________________________
void TDataset::SetDataFolder(const TString& dataFolder)
{
  // Specify the location of the strangecalc-data folder. It is advisable
  // to use absolute paths. Although relative paths should work in principle.
  fgDataFolder = dataFolder;
  
  if( !dataFolder.EndsWith("/") )
    fgDataFolder += "/";
}

//_____________________________________________________________________
void TDataset::SetSelection(int selection)
{
  // Activate a selection. See ViewSelections() for possible selections.
  if( fSelection < 0 ) {
    std::cerr << "WARNING in TDataset::SetSelection(int): "
	      << "First import a dataset!\n";
    return;
  }

  else if( selection >= fSelecEventLists->GetEntries() ) {
    std::cerr << "WARNING in TDataset::SetSelection(int):\n"
	      << "Selection no." << selection << " does not exist. "
	      << "See TDataset::ViewSelections() for options.\n";
    fSelection = 0;
  }

  else {
    fSelection = selection;
  }
   
  TEventList *eventList = (TEventList*)fSelecEventLists->At(fSelection);
  SetEventList(eventList);
}

//_____________________________________________________________________
void TDataset::SetDescription(int selection, const TString& description)
{
  // Private member to set the description of a selection.
  // TDataset::SetDescription(TString) is the public version.

  if( fSelection<0 ) {
    std::cerr << "WARNING in TDataset::SetDescription(int,TString): "
	      << "First import a dataset!\n";
    return;
  }
  
  if( selection >= 0 && selection<fSelecEventLists->GetEntries() ) {
    delete fSelecDescription->RemoveAt(selection);
    fSelecDescription->AddAt(new TObjString(description),selection);
  }
  else
    std::cerr << "WARNING in TDataset::SetDescription(int,TString): "
	      << "Selection no." << selection << " does not exist!\n";
}

//_____________________________________________________________________
void TDataset::SetDescription(const TString& description)
{
  // Set a description for the current active selection.
  SetDescription(fSelection,description);
}

//_____________________________________________________________________
TString TDataset::GetDescription(int selection) const
{
  // Private member to retrieve the description of a selection.
  // TDataset::GetDescription() is the public version.
  
  TString description = "";
  description += ( (TObjString*)fSelecDescription->At(0) )->GetString();

  if( selection > 0 ) {
    description += " ("; 
    description += ((TObjString*)fSelecDescription->At(selection))->GetString();
    description += ")";
  }

  return description.Data();
}

//_____________________________________________________________________
TString TDataset::GetDescription() const
{
  // Retrieve the description of the current active selection.
  if( fSelection <0 ) {
    std::cerr << "WARNING in TDataset::GetDescription(): "
	      << "First import a dataset!\n";
    return "";
  }
  
  return GetDescription(fSelection);
}

//_____________________________________________________________________
void TDataset::ViewSelections() const
{
  // Prints the list of selections of the dataset to the screen. The current
  // active selection is marked with an asterisk.
  if( fSelection <0 )
    std::cerr << "WARNING in TDataset::ViewSelections(): "
	      << "First import a dataset!\n";
   
  else {
    std::cout << " " << std::setw(3) << "No." << std::setw(9) << "Entries"
	      << "\t" << "Description\n"
	      << " " << std::setw(3) << "---" << std::setw(9) << "-------"
	      << "\t" << "-----------\n";

    for(int i=0; i<fSelecEventLists->GetEntries(); ++i) 
      std::cout << ( fSelection==i? "*" : " ")
		<< std::setw(3) << i << std::setw(9) 
		<< ((TEventList*)fSelecEventLists->At(i))->GetN()
		<< "\t" 
		<< GetDescription(i)
		<< "\n";
  }
}

//_____________________________________________________________________
void TDataset::SetKinematics(const TKinematics* kinematics)
{
  // Link a TKinematics object to the current active selection. The argument
  // will be cloned, meaning the user keeps ownership over the argument.
  if( fSelection >= 0 ) {
    if( fIsospin == kinematics->GetIsospin() ) {
      delete fSelecKinematics->RemoveAt(fSelection);
      fSelecKinematics->AddAtAndExpand(kinematics->Clone(),fSelection);
    }
    else
      std::cerr << "WARNING in TDataset::SetKinematics(TKinematics*): "
		<< "Isospin channels do not match!\n";
    
  }
  else 
    std::cerr << "WARNING in TDataset::SetKinematics(TKinematics*): "
	      << "First import a dataset!\n";
}

//_____________________________________________________________________
TKinematics* TDataset::GetKinematics() const
{
  // Retrieve the TKinematics object associated with the current
  // active selection. Do NOT delete it.
  if( fSelection <0 ) {
    std::cerr << "WARNING in TDataset::GetKinematics(): "
	      << "First import a dataset!\n";
    return 0;
  }

  TKinematics *kinematics = ((TKinematics*) fSelecKinematics->At(fSelection));

  if( !kinematics)
    std::cerr << "WARNING in TDataset::GetKinematics(): "
	      << "No kinematics specified for current selection.\n";

  return kinematics;
}

//_____________________________________________________________________
const TCalcInfo* TDataset::GetSpecification(Long64_t entry) const
{
  // Retrieve the correct TCalcInfo object corresponding to data point 'entry'
  // of the current active selection.
  //
  // Do NOT delete the object.

  TDataset *nonconst_this = const_cast<TDataset*>(this);

  // Fill the Data-struct when it hasn't been initialized
  if( !nonconst_this->fSpecification ) {
    nonconst_this->fSpecification = new TCalcInfo();
    nonconst_this->fSpecification->fData->iso = fIsospin;
    nonconst_this->fSpecification->fData->electro_prod = fElectro_prod;
    nonconst_this->fSpecification->fData->photo_prod = fPhoto_prod;
    nonconst_this->fSpecification->fData->kaoncapture = fKaoncapture;
    if( fElectro_prod ) {
      nonconst_this->fSpecification->fData->elec.iso = fIsospin;
      nonconst_this->fSpecification->fData->elec.ds_dt = fDs_dt;
      nonconst_this->fSpecification->fData->elec.beam_ener_input = fBeam_ener_input;
      nonconst_this->fSpecification->fData->elec.e_beam_ener = fE_beam_ener;
      nonconst_this->fSpecification->fData->elec.eps = fEps;
      nonconst_this->fSpecification->fData->elec.cs_convention = fCs_convention;
    }
    else if ( fPhoto_prod ) {
      nonconst_this->fSpecification->fData->photo.iso = fIsospin;
      nonconst_this->fSpecification->fData->photo.ds_dt = fDs_dt;
      nonconst_this->fSpecification->fData->photo.ds_du = fDs_du;
    }
    else {
      nonconst_this->fSpecification->fData->kaoncap.iso = fIsospin;
      nonconst_this->fSpecification->fData->kaoncap.cosmin = fCosmin;
      nonconst_this->fSpecification->fData->kaoncap.cosmax = fCosmax;
    }
  }
  
  // Update the observable
  // We do this by giving an address to branch 'observable' and reading
  // it for 'entry'.
  // We must be carefull however that the user doesn't notice we have been
  // doing something with the branch addresses.
  TObjArray *branches = nonconst_this->GetListOfBranches();
  
  // Deactivate all branches
  bool activeBranches[branches->GetEntries()];
  for(int i=0; i<branches->GetEntries(); ++i) {
    TBranch *branch = (TBranch*)(branches->At(i));
    activeBranches[i] = GetBranchStatus(branch->GetName());
    nonconst_this->SetBranchStatus(branch->GetName(),0);
  }

  // Set the 'observable' branch
  nonconst_this->SetBranchStatus("observable",1);
  char observable[15];
  TBranch *addressBranch = nonconst_this->GetBranch("observable");
  char *observableAddress = (char*) addressBranch->GetAddress();
  if( !observableAddress ) {
    std::strcpy(observable,"(empty)(empty)");
    observableAddress = observable;
    nonconst_this->SetBranchAddress("observable",observableAddress);
  }

  // Get the observable
  nonconst_this->GoTo(entry);

  // Store the observable in the TCalcInfo object
  if( fElectro_prod )
  {
    std::strcpy(nonconst_this->fSpecification->fData->elec.observable,"(empty)(empty)");
    std::strcpy(nonconst_this->fSpecification->fData->elec.observable, observableAddress);
  }
  else if( fPhoto_prod )
  {
    std::strcpy(nonconst_this->fSpecification->fData->photo.observable,"(empty)(empty)" );
    std::strcpy(nonconst_this->fSpecification->fData->photo.observable, observableAddress);
  }
  else
  {
    std::strcpy(nonconst_this->fSpecification->fData->photo.observable,"(empty)(empty)" );
    std::strcpy(nonconst_this->fSpecification->fData->kaoncap.observable, observableAddress);
  }

  // Re-activate all branches as before
  if( observableAddress == observable ) nonconst_this->ResetBranchAddress(addressBranch);
  for(int i=0; i<branches->GetEntries(); ++i) {
    TBranch *branch = (TBranch*)(branches->At(i));
    nonconst_this->SetBranchStatus(branch->GetName(),activeBranches[i]);
  }

  return fSpecification;
}

//_____________________________________________________________________
void TDataset::ImportDataset(int isospin, const TString& datasets)
{
  // Fill the tree with a strangecalc dataset. The user should specifiy the
  // isospin channel of the data.
  // If the second argument is empty (default), the user will be prompted
  // interactively to select a dataset. The number of the dataset to be
  // imported can be given as the second argument.
  if( fSelection >= 0 ) {
    std::cerr << "WARNING in TDataset::ImportDataset(int):\n"
	      << "A TDataset object can contain only one(!) dataset!\n";
    return;
  }
  
  // The user can choose a dataset. He/she has 2 options:
  // * specify a dataset in the second argument
  // * or interactively
  // We recycle databasespecification() defined in io_specific.*.c
  // To be able to do this, we need to declare the C-struct Observable
  // defined in Structures.h
  Observable observ = {0};
  std::strcpy(observ.dataFolder,fgDataFolder);
  observ.iso.nr_iso_channels = 1;
  observ.iso.iso_channel[0] = isospin;
  observ.iso.isospin = isospin;
  short dummy;

  FILE* input; // input we feed to databasespecification
  if(datasets.IsNull()) // interactif
    input = stdin;
  else {                // user given
    input = tmpfile();
    fprintf(input,"%s",(const char*)datasets);
    fseek(input,0,SEEK_SET);
  }

  databasespecification(&observ,&dummy,&dummy,&dummy,input,stdout,stdout);
  
  if(!datasets.IsNull()) fclose(input);

  if( observ.fit.nr_experiments[isospin] != 1 ) {
    std::cerr << "ERROR in TDataset::ImportDataset(int):\n"
	      << "A TDataset object can contain only one(!) dataset!\n";
    exit(1);
  }


  // After choosing the dataset, we store it in another C-struct
  // Data defined in fitting.*.h.
  // we recycle getdatastructure() defined in fitting.*.c, to
  // do the conversion.
  Data datapoints[DATAMAX];
  int nrOfDatapoints=0;

  getdatastructure(datapoints,&nrOfDatapoints,&observ,0);

  // Store general specifications of the dataset
  // -------------------------------------------
  fIsospin = isospin;
  fElectro_prod = datapoints[0].electro_prod;
  fPhoto_prod = datapoints[0].photo_prod;
  fKaoncapture = datapoints[0].kaoncapture;
  if( fElectro_prod ) { // electroproduction specific
    fDs_dt = datapoints[0].elec.ds_dt;
    fBeam_ener_input = datapoints[0].elec.beam_ener_input;
    fE_beam_ener = datapoints[0].elec.e_beam_ener;
    fEps = datapoints[0].elec.eps;
    fCs_convention = datapoints[0].elec.cs_convention;
  }
  else if( fPhoto_prod ) { // photoproduction specific
    fDs_dt = datapoints[0].photo.ds_dt;
    fDs_du = datapoints[0].photo.ds_du;
  }
  else {
    fCosmin = datapoints[0].kaoncap.cosmin;
    fCosmax = datapoints[0].kaoncap.cosmax;
  }

  // Create branches to hold the dataset
  // -----------------------------------
  char *observable = new char[15]; std::strcpy(observable,"(empty)(empty)");
  Branch("observable",observable,"observable/C");

  double wlab=0.0, wcm=0.0, klab=0.0, kcm=0.0, pk=0.0, pklab=0.0, 
    costhkcm=0.0, qsquared=0.0, w=0, s=0.0, t=0.0, u=0.0;
  int isPhysical=0;
  Branch("wlab",&wlab,"wlab/D");
  Branch("wcm",&wcm,"wcm/D");
  Branch("klab",&klab,"klab/D");
  Branch("kcm",&kcm,"kcm/D");
  Branch("pk",&pk,"pk/D");
  Branch("pklab",&pklab,"pklab/D");
  Branch("costhkcm",&costhkcm,"costhkcm/D");
  Branch("qsquared",&qsquared,"qsquared/D");
  Branch("w",&w,"w/D");
  Branch("s",&s,"s/D");
  Branch("t",&t,"t/D");
  Branch("u",&u,"u/D");
  Branch("isPhysical",&isPhysical,"isPhysical/I");

  double amplitude = 0.0, error = 0.0;
  Branch("amplitude",&amplitude,"amplitude/D");
  Branch("error",&error,"error/D");

  // Fill TTree
  // ----------
  // Declare TKinematics object for conversions
  TKinematics *tk = 0;

  // construct the correct tk
  if( fPhoto_prod )
    tk = new TKinematics("photoprod","photoprod kinematics",isospin,
			 "wlab:costhkcm:qsquared",0.0,0.0,0.0);
  else if( fElectro_prod ) {
    if( datapoints[0].elec.cos_ang )
      tk = new TKinematics("electroProdC",
			   "electroprod kinematics with cos",isospin,
			   "s:costhkcm:qsquared",0.0,0.0,0.0);
    else
      tk = new TKinematics("electroProdT",
			   "electroprod kinematics with t",isospin,
			   "s:t:qsquared",0.0,0.0,0.0);
  }
  else
    tk = new TKinematics("kaoncap","kaon capture kinematics",isospin,
			 "pklab:costhkcm:qsquared",0.0,0.0,0.0);

  // Walk throught the dataset
  for(int i=0; i<nrOfDatapoints; ++i) {
    if( fPhoto_prod ) {
      std::strcpy(observable,"(empty)(empty)");
      std::strcpy(observable,datapoints[i].photo.observable);
    
      wlab = (datapoints[i].photo.emin+datapoints[i].photo.emax)/2.0;
      costhkcm = datapoints[i].photo.cos;
      tk->SetVar(1,wlab);
      tk->SetVar(2,costhkcm);
      wcm = tk->GetWcm();
      klab = tk->GetKlab();
      kcm = tk->GetKcm();
      pk = tk->GetPk();
      pklab = tk->GetPklab();
      w = tk->GetW();
      s = tk->GetS();
      t = tk->GetT();
      u = tk->GetU();
      isPhysical = ( tk->IsPhysical() ? 1 : 0 );
      
      amplitude = datapoints[i].photo.ampli;
      error = datapoints[i].photo.error;
    }
    else if( fElectro_prod ) {
      std::strcpy(observable,datapoints[i].elec.observable);

      s = datapoints[i].elec.s;
      qsquared = datapoints[i].elec.qsquared;
      if( datapoints[i].elec.cos_ang ) {
	costhkcm = datapoints[i].elec.cos;
	tk->SetVar(2,costhkcm);
      }
      else {
	t = datapoints[i].elec.t;
	tk->SetVar(2,t);
      }
      tk->SetVar(1,s);
      tk->SetVar(3,qsquared);
      wlab = tk->GetWlab();
      wcm = tk->GetWcm();
      klab = tk->GetKlab();
      kcm = tk->GetKcm();
      pk = tk->GetPk();
      pklab = tk->GetPklab();
      costhkcm = tk->GetCosthkcm();
      w = tk->GetW();
      t = tk->GetT();
      u = tk->GetU();
      isPhysical = ( tk->IsPhysical() ? 1 : 0 );

      amplitude = datapoints[i].elec.ampli;
      error = datapoints[i].elec.error;
    }
    else {
      std::strcpy(observable,datapoints[i].kaoncap.observable);

      pklab = datapoints[i].kaoncap.pk;
      costhkcm = datapoints[i].kaoncap.cos;
      tk->SetVar(1,pklab);
      tk->SetVar(2,costhkcm);
      wlab = tk->GetWlab();
      wcm = tk->GetWcm();
      klab = tk->GetKlab();
      kcm = tk->GetKcm();
      pk = tk->GetPk();
      w = tk->GetW();
      s = tk->GetS();
      t = tk->GetT();
      u = tk->GetU();
      isPhysical = ( tk->IsPhysical() ? 1 : 0 );

      amplitude = datapoints[i].kaoncap.ampli;
      error = datapoints[i].kaoncap.error;
    }
    
    Fill();
  
  } // end of dataset

  delete tk; // release TKinematics		 
  ResetBranchAddresses(); // necessary because variables will go out-of-scope

  // Show the user the contents of the TDataset
  // ------------------------------------------
  double nrOfEntries =Scan("observable:wlab:wcm:klab:kcm:pk:pklab:costhkcm:w:s:t:u");

  if( nrOfEntries != nrOfDatapoints )
    std::cerr << "WARNING in TDataset::ImportDataset(int): "
	      << "Not all datapoints were imported!\n";
  else
    std::cout << "INFO in TDataset::ImportDataset(int): "
	      << "Imported " << nrOfEntries << " datapoints.\n";


  // Initialize the selection arrays
  // -------------------------------
  fSelection = 0; // first entry is full selection
  
  // Create selection lists
  fSelecDescription = new TObjArray(1);
  fSelecEventLists = new TObjArray(1);
  fSelecKinematics = new TObjArray(1);
  
  // Make them owner of their entries
  fSelecDescription->SetOwner();
  fSelecEventLists->SetOwner();
  fSelecKinematics->SetOwner();

  AddSelection("",0);
  SetDescription(observ.fit.database_info[isospin][0]);

  // Clean up
  delete[] observable;
  
}

//_____________________________________________________________________
int TDataset::AddSelection(const char* selectionString, TKinematics* kinematics)
{
  // Create a new selection with the data points of the current active selection
  // that satisfy the selectionString. 
  //
  // selectionString is a boolean expression with a combination of the columns.
  // In a selectionString all the C++ operators are authorized. If the result is true,
  // the data point is added to the new selection.
  // Example: "pklab<klab && acos(costhkcm)>1.5"
  //
  // With the second argument, the user can associate a TKinematics object to
  // the new selection. See SetKinematics for more info.
  //
  // Returns the number of the new selection.

  if( fSelection <0 ) {
    std::cerr << "WARNING in TDataset::AddSelection(char*,TKinematics*): "
	      << "First import a dataset!\n";
    return -1;
  }

  // Put selection in the next slot
  int entry = fSelecEventLists->GetEntries();
  
  // Create a TEventList
  Draw(">> eventList",selectionString);
  TEventList *selecEventList = 
    (TEventList*)gDirectory->GetList()->FindObject("eventList");
  selecEventList->SetDirectory(0); // I want to own the TEventList
  fSelecEventLists->AddAtAndExpand(selecEventList,entry);
  
  // Add a description
  TString selectionStringFull = "";
  if( fSelection != 0 ) {
    selectionStringFull += 
      ((TObjString*) fSelecDescription->At(fSelection))->GetString();
    selectionStringFull += " && ";
  }
  selectionStringFull += selectionString;
  fSelecDescription->AddAtAndExpand(new TObjString(selectionStringFull),entry);
  
  // Add kinematics
  if(kinematics)
    fSelecKinematics->AddAtAndExpand(kinematics->Clone(),entry);
  else
    fSelecKinematics->AddAtAndExpand(kinematics,entry);
  
  // return position of selection
  return ( fSelecEventLists->GetEntries() - 1 );  
}

//_____________________________________________________________________
void TDataset::RemoveSelection(int selection)
{
  // Delete the specified selection. Selection nr.0 can obviously not be deleted.
   if( fSelection < 0 )
    std::cerr << "WARNING in TDataset::RemoveSelection(int): "
	      << "First import a dataset!\n";

   else if( selection == 0 )
     std::cerr << "WARNING in TDataset::RemoveSelection(int): "
	       << "Selection no.0 can NOT be removed.\n";

   else if( selection >= fSelecEventLists->GetEntries() )
     std::cerr << "WARNING in TDataset::RemoveSelection(int):\n"
	       << "Selection no." << selection << " does not exist. "
	       << "See TDataset::ViewSelections() for options.\n";

   else {
     // Remove the entries
     delete fSelecDescription->RemoveAt(selection);
     delete fSelecEventLists->RemoveAt(selection);
     delete fSelecKinematics->RemoveAt(selection);
     
     // Reconfigure the selection lists
     int nrOfEntries = fSelecEventLists->GetEntries();
     for(int i=selection; i<nrOfEntries; ++i) {
       fSelecDescription->AddAt(fSelecDescription->At(i+1),i);
       fSelecEventLists->AddAt(fSelecEventLists->At(i+1),i);
       fSelecKinematics->AddAt(fSelecKinematics->At(i+1),i);
     }
     fSelecDescription->RemoveAt(nrOfEntries);
     fSelecEventLists->RemoveAt(nrOfEntries);
     fSelecKinematics->RemoveAt(nrOfEntries);

     // Reset the position
     if( fSelection == selection) SetSelection(0);
     else if( fSelection > selection ) SetSelection(fSelection-1);
   }
}

//_____________________________________________________________________
int TDataset::RemoveEmptySelections()
{
  // Delete selections that don't have entries.
  // Returns the number of deleted selections.
  int currentSelection = fSelection;
  int removed = 0;

  for( int i=1; i<fSelecEventLists->GetEntries(); ++i) {
    SetSelection(i);
    if( GetEntries("") == 0 ) {
      RemoveSelection(i--);
      removed++;
    }
  }

  SetSelection(currentSelection);

  return removed;
}

//_____________________________________________________________________
int TDataset::MakeSelections(const TString& listOfVariables)
{
  // list is a colon-separated list of branch names. This function will 
  // create a selection for every combination of values of the specified
  // branches and fill them with data points.
  // 
  // Returns the number of newly created selections.

  // Tokenize the list of variables
  TObjArray *selectors = listOfVariables.Tokenize(":");
  int nrOfSelectors = selectors->GetEntries();
  
  if( nrOfSelectors == 0 ) {
    std::cerr << "ERROR in TDataset::MakeSelections(TString):\n"
	      << "Select a least one variable!\n";
    exit(1);
  }
  
  // Make a list of values of the selection variables
  // and the corresponding selection expressions
  // selectorStrings are the actual selection strings
  // prettySelectorStrings are the selection descriptions that
  // are printed on screen.
  std::set<double> values[nrOfSelectors];
  std::set<TString> svalues[nrOfSelectors];
  std::vector<TString> selectorStrings[nrOfSelectors];
  std::vector<TString> prettySelectorStrings[nrOfSelectors];
    
  // Determine selections
  // --------------------
  // We will loop over all selectors. Everytime we encounter a
  // new value for the selector, we build a selection string and
  // a selection description. To keep track of the number of 
  // different selector values, we use a set<double>.
  // When the selector is the "observable", it requires a special
  // treatment, since it's a char-array, whereas all other selectors
  // are doubles.
  TString selector;              // name of selection variable
  TString selectionString;       // selection string
  TString prettySelectionString; // selection description
  TBranch *selectorBranch;       // points to branch of selection variable
  double value=0.0;              // temporary value of selection variable
  char observable[15];     // temporary char array
  
  for(int i=0; i<nrOfSelectors; ++i) { // evaluate all selectors
    
    selector = ((TObjString*)selectors->At(i))->GetString();
    selectorBranch = FindBranch(selector);  // Grap the correct branch
    
    if( !selectorBranch ) { // does the branch exist?
      std::cerr << "ERROR in TDataset::MakeSelections(TString):\n"
		<< "Branch \"" << selector << "\" does not exist!\n";
      exit(1);
    }
    
    // Give an address to the branch
    if( selector == "observable" ) // requires special treatment
      selectorBranch->SetAddress(observable);
    else 
      selectorBranch->SetAddress(&value);
    
    // Loop over all entries
    for(int j=0; j<selectorBranch->GetEntries(); ++j) {
      std::strcpy(observable,"(empty)(empty)");
      selectorBranch->GetEntry(j);
      bool second;
      // Evaluate the 'value' of the selector
      if( selector == "observable" ) // requires special treatment
	second = (svalues[i].insert(observable)).second;
      else
	second= (values[i].insert(value)).second;

      // When it's a new selector 'value', build the
      // selection string and selection description.
      // The selection description differs in the way doubles
      // are stored.
      if(second) 
	{
	  selectionString = selector;
	  selectionString += "==";
	  if( selector == "observable" ) { // requires special treatment
	    selectionString += "\""; 
	    selectionString += observable;
	    selectionString += "\"";
	    prettySelectionString = selectionString;
	  }
	  else {
	    prettySelectionString = selectionString;
	    selectionString += TString::Format("%.16e",value);
	    prettySelectionString += TString::Format("%g",value);
	  }
	  selectorStrings[i].push_back(selectionString);
	  prettySelectorStrings[i].push_back(prettySelectionString);
	}
    } // end loop over all branch entries
    
    selectorBranch->ResetAddress();

  } // end loop over selectors
  

  // Construct full selection expressions and add selections
  // -------------------------------------------------------
  // This piece of code looks like hell, but it works
  // and that's all that matters! ;)
  selectionString = "";
  int index[nrOfSelectors+1];
  int maxIndex[nrOfSelectors+1];
  maxIndex[0] = 1;
  for(int i=0; i<nrOfSelectors; ++i)
    maxIndex[i+1] = maxIndex[i] * selectorStrings[i].size();
  
  for(int i=0; i<maxIndex[nrOfSelectors]; ++i) {
    for(int j=nrOfSelectors; j>0; --j) {
      index[j-1] = i;
      for(int k=j; k<nrOfSelectors; ++k)
	index[j-1] -= index[k] * maxIndex[k]; 
      index[j-1] /= maxIndex[j-1];
    }
    
    selectionString = selectorStrings[0][index[0]];
    if( fSelection == 0 ) 
      prettySelectionString = "";
    else {
      prettySelectionString = 
	((TObjString*) fSelecDescription->At(fSelection))->GetString();
      prettySelectionString += " && ";
    }
    prettySelectionString += prettySelectorStrings[0][index[0]];
    for(int j=1; j<nrOfSelectors; ++j) {
      selectionString += " && ";
      selectionString += selectorStrings[j][index[j]];
      prettySelectionString += " && ";
      prettySelectionString += prettySelectorStrings[j][index[j]];
    }
    
    int position = AddSelection(selectionString,0);
    SetDescription(position,prettySelectionString);    
  }

  int removedSelections = RemoveEmptySelections();

  // Print some info
  std::cout << "INFO in TDataset::MakeSelections(TString): "
	    << "Created " << maxIndex[nrOfSelectors] - removedSelections
	    << " selections.\n";

  return maxIndex[nrOfSelectors] - removedSelections;
}

//_____________________________________________________________________
Int_t TDataset::GoTo(Long64_t entry, Int_t getall)
{
  // Read all branches of entry nr. 'entry' of the current active selection
  // and return the total number of bytes read.
  // See TTree::GetEntry for more info.
  return TTree::GetEntry(GetEntryNumber(entry),getall);
}

//_____________________________________________________________________
TGraphErrors* TDataset::MakeGraph(const TString& variable,
				  double xunit, double yunit)
{
  // Create a TGraphErrors of the data points in the current active
  // selection. The data points are given as a function of "xvariable".
  //
  // xvariable can be any branch name. "-t" is also allowed.
  //
  // the x- and y-axis will be scaled with 'xunit' and 'yunit'.
  //
  // The user owns the returned graph and is responsible for deleting it.
  if( variable.IsNull() ) {
    std::cerr << "ERROR in TDataset::MakeGraph(TString,double=1.0,double=1.0): "
	      << "Specify a kinematic variable for the x-axis.\n";
    return 0;
  }

  TString varexp;
  if( std::strcmp(variable,"-t") )
    varexp = variable;
  else
    varexp = "t";
  varexp += ":amplitude:error";

  TTree::Draw(varexp,"","goff");

  double *x = TTree::GetV1();
  double *y = TTree::GetV2();
  double *z = TTree::GetV3();

  if(!std::strcmp(variable,"-t"))
    for(int i=0; i<TTree::GetSelectedRows(); ++i)
      x[i] *= -1.;

  // Scale the independent variable with the 'xunit', and
  // the amplitude+error with the 'yunit'.
  for(int i=0; i<TTree::GetSelectedRows(); ++i){
    x[i] /= xunit;
    y[i] /= yunit;
    z[i] /= yunit;
  }

  return new TGraphErrors(TTree::GetSelectedRows(),x,y,0,z);  
}

//_____________________________________________________________________
Long64_t TDataset::Scan(const char* varexp, const char* selection,
			Option_t* option, Long64_t nentries,
			Long64_t firstentry)
{
  // Loop over tree entries and print entries passing selection.
  // 
  // The default varexp is "observable:qsquared:wlab:w:pklab:costhkcm:amplitude:error"
  // See TTree::Scan for more info.
  TString varexp2;
  if( !strcmp(varexp,"") )
    varexp2 = "observable:qsquared:wlab:w:pklab:costhkcm:amplitude:error";
  else
    varexp2 = varexp;

  return TTree::Scan(varexp2,selection,option,nentries,firstentry);
}

//_____________________________________________________________________
void TDataset::Help()
{
  std::cout << "You seek help, young Padawan?\n"
            << "1) Open emacs\n"
            << "2) type Alt+x\n"
            << "3) type doctor\n"
            << "4) press return\n";
}
