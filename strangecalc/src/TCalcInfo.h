/*!
 * \file TCalcInfo.h
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

#ifndef TCALCINFO_H
#define TCALCINFO_H

#include <TObject.h>
#include <TString.h>
#include <fitting.h>

class TCalcInfo : public TObject
{
  // TDataset and TStrangeModel are friends
  friend class TDataset;
  friend class TStrangeModel;

 public:
  enum EReaction { kPhoto = 1,
		   kElectro = 2,
		   kCapture = 3 };

  enum EInitVerbosity {
    // Specifies whether warnings will be printed to the screen.
    kVerbose = 1,
    kSilent = 2
  };
  
 public:
  // Constructors
  TCalcInfo(EInitVerbosity=kVerbose);
  TCalcInfo(EReaction reaction,int isospin,const char *observable, EInitVerbosity=kVerbose);
  // cast TCalcInfo* to Data*
  operator Data*() const       { return fData; }
  TCalcInfo(const Data&);
  TCalcInfo(const TCalcInfo&);
  virtual TCalcInfo* Clone(const char *newname ="") const;
  virtual ~TCalcInfo();
  TCalcInfo&        operator=(const TCalcInfo&);

  static void       Help();

  void SetIsospin(short iso);
  void SetPhotoproduction();
  void SetElectroproduction();
  void SetKaoncapture();
  void SetObservable(TString observable);
  void SetDsDt();
  void SetDsDu();
  void SetDsDomega();
  void SetElectroConvention(int convention);
  void SetElectroBeamEnergy(double energy); // in [MeV]
  void SetElectroEpsilon(double epsilon);

  short       GetIsospin() const           { return fData->iso; }
  bool        IsPhotoproduction() const    { return fData->photo_prod; }
  bool        IsElectroproduction() const  { return fData->electro_prod; }
  bool        IsKaoncapture() const        { return fData->kaoncapture; }
  int         GetElectroConvention() const { return fData->elec.cs_convention; }
  double      GetElectroBeamEnergy() const { return fData->elec.e_beam_ener; }
  double      GetElectroEpsilon() const    { return fData->elec.eps; }
  const char *GetObservable() const;

  	      
  
 private:
  Data  *fData;              //! The internal Data struct
  EInitVerbosity fVerbosity;  // verbosity mode

  ClassDef(TCalcInfo,0)      // Data struct wrapper
};

#endif
