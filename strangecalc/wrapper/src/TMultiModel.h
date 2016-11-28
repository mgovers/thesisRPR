/*!
 * \file TMultiModel.h
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

#ifndef TMULTIMODEL_H
#define TMULTIMODEL_H

#include "TStrangeModel.h"
#include <TStrangeCalc.h>
#include <map>

//////////////////////////////////////////////////////////////////////////////
//
// TMultiModel
//
// A set of TStrangeModel's for multiple isospin channels.
//
//////////////////////////////////////////////////////////////////////////////

class TMultiModel : public TStrangeModel
{
 public:
  enum EPublishedModel { 
    // Constants to initialize published RPR models.
    // See documentation of TMultiModel::SetStrangeModel(const EPublishedModel&)
    // for more information.
    kRPR2011 = 1,
    kRPR2007 = 2,
    kRPR2007alt = 5,
    kRegge2011 = 3,
    kRegge2007 = 4,
    kRegge2007alt = 6
  };

 public:
  TMultiModel(const char *name="", const char *title="");
  TMultiModel(TRootIOCtor* rio); // ROOT I/O Constructor
  TMultiModel(const TMultiModel&);
  TMultiModel& operator=(const TMultiModel&);
  virtual TMultiModel* Clone(const char* newname="") const;
  virtual ~TMultiModel();

  virtual void SetStrangeModel(const char *fit_spec_file, EInitVerbosity verbosity=kVerbose);
  virtual void SetStrangeModel(int iso,const char *fit_spec_file, EInitVerbosity verbosity=kVerbose);
  virtual void SetStrangeModel(int iso,const TStrangeCalc&); // model will be copied
  virtual void SetStrangeModel(const EPublishedModel&, EInitVerbosity verbosity=kVerbose);
  virtual double  ChiSquared(const TString& isospins, const TString& datasets ="", const short polarizationWeight=1) const;

  virtual const TStrangeCalc& GetTStrangeCalc(int iso) const;
  virtual TStrangeCalc& GetTStrangeCalc(int iso);

 private:
  static TStrangeCalc* CopyTStrangeCalc(const TStrangeCalc*);
  static void StripData(TStrangeCalc*);

 private:
  std::map<int,TStrangeCalc*> fModels; // map of models

  ClassDef(TMultiModel,1); // A set of TStrangeModel's

}; // class TMultiModel

#endif // TMULTIMODEL_H
