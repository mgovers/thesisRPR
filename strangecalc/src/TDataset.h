/*!
 * \file TDataset.h
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

/*
 * TDataset.h
 *
 * ROOT object that can hold a strangecalc dataset.
 * Data is internally stored as a TTree.
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

#ifndef TDATASET_H
#define TDATASET_H

#include <TTree.h>
#include <TString.h>
#include <TGraphErrors.h>
#include "TKinematics.h"
#include "TCalcInfo.h"

class TDataset : public TTree
{
 public:
  // Constructors
  TDataset(const char* name="",const char* title="");
  
  // Destructor
  virtual ~TDataset();

  // Methods
  void               ImportDataset(int isospin,const TString& datasets ="");
  int                AddSelection(const char *selectionString,TKinematics* =0); // function will -copy- TKinematics
  int                MakeSelections(const TString& list); // list separated by :
  void               RemoveSelection(int);
  void               ViewSelections() const;
  Int_t              GoTo(Long64_t entry, Int_t getall =0);
  TGraphErrors*      MakeGraph(const TString& xvariable, double xunit= 1.0, double yunit= 1.0); // user owns TGraphErrors
  virtual Long64_t   Scan(const char* ="",const char* ="", Option_t* ="",
			  Long64_t=1000000000, Long64_t=0);

  // Getters

  static const char *GetDataFolder()          { return fgDataFolder; }
  int                GetSelection() const     { return fSelection; }
  TString            GetDescription() const;
  TKinematics       *GetKinematics() const;   // do not delete!
  const TCalcInfo   *GetSpecification(Long64_t entry=0) const;// do not delete!
  const TCalcInfo   *GetCalcInfo(Long64_t entry=0) const { return GetSpecification(entry); } // do not delete!
  Long64_t           GetNrOfEntries()         { return GetEntries(""); } 
  
  // Setters
  static void        SetDataFolder(const TString&);
  void               SetDescription(const TString&);
  void               SetKinematics(const TKinematics*); // function will -copy- TKinematics
  void               SetSelection(int);

  // Static Methods
  static void        Help();
 

 private:
  // copy constructor and assignment are private
  // the user should use TDataset::Clone()
  TDataset(const TDataset&);
  TDataset& operator=(const TDataset&);

  // private methods
  TString            GetDescription(int) const;
  void               SetDescription(int,const TString&);
  int                RemoveEmptySelections();


  // Data members
  // ------------
  static TString fgDataFolder;   // Location of datasets

  // Dataset specifications
  TCalcInfo  *fSpecification;    //! full specification
  int         fIsospin;          // isospin channel of dataset
  short       fElectro_prod;     // electroproduction data?
  short       fPhoto_prod;       // photoproduction data?
  short       fKaoncapture;      // kaon capture data?
  short       fDs_dt;            // observables function of t (only electro)
  short       fDs_du;            // observables function of u (only electro)
  short       fBeam_ener_input;  // is beam energy given (only electro)
  double      fE_beam_ener;      // beam energy (only electro)
  double      fEps;              // epsilon prefactor (only electro)
  short       fCs_convention;    // cross section convention (only electro)
  double      fCosmin;           // lower cosine limit (ptcs kaoncap only)
  double      fCosmax;           // upper cosine limit (ptcs kaoncap only)

  // Selections
  int         fSelection;        // current selection (0 for full dataset)
  TObjArray  *fSelecDescription; //-> List of selection descriptions
  TObjArray  *fSelecEventLists;  //-> List of selection event lists
  TObjArray  *fSelecKinematics;  //-> List of selection kinematics (owns entries)


  ClassDef(TDataset,2)           // Dataset for strangeViewer
};

#endif
