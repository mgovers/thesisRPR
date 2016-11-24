/*!
 * \file TStrangeModel.h
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

#ifndef TSTRANGEMODEL_H
#define TSTRANGEMODEL_H

#include <TNamed.h>
#include <TGraph.h>
#include <FourVector.h>
#include <GammaStructure.h>
#include <TStrangeCalc.h>
#include "TKinematics.h"
#include "TDataset.h"
#include "TCalcInfo.h"
#include <Structures.h>
#include <fitting.h>
class TString;

class TStrangeModel : public TNamed
{
 public:
  enum EInitVerbosity {
    // Specifies whether initialization information will be printed
    // to the screen.
    kVerbose = 1,
    kSilent = 2
  };

 public:
  // Constructors
  TStrangeModel(const char *name="", const char *title="");
  TStrangeModel(TRootIOCtor* rio); // ROOT I/O Constructor
  TStrangeModel(const TStrangeModel&);
  TStrangeModel& operator=(const TStrangeModel&);
  virtual TStrangeModel* Clone(const char* newname = "") const;

  // Destructor
  virtual ~TStrangeModel();

  // Static Methods
  static void   Help();

  // Methods
  TGraph* MakeGraph(const TKinematics&, const TCalcInfo*, const TString& xvariable, double xunit= 1.0, double yunit= 1.0); // user owns TGraph*
  TGraph* MakeGraph(const TDataset&, const TString& xvariable, double xunit= 1.0, double yunit= 1.0); // user owns TGraph*
  virtual double  ChiSquared(const TString& isospins, const TString& datasets ="", const short polarizationWeight=1) const;

  // Getters
  double* GetCalcpoints(const TKinematics&, const TCalcInfo*); // user owns double*
  double* GetCalcpoints(const TDataset&); // user owns double*
  double GetCalcpoint(const TKinematics&, const TCalcInfo*);
  double* GetPhysicalCalcpoints(const TKinematics&, const TCalcInfo*); // user owns double*
  double* GetPhysicalCalcpoints(const TDataset&); // user owns double*
  virtual const TStrangeCalc& GetTStrangeCalc(int iso) const;
  virtual TStrangeCalc& GetTStrangeCalc(int iso);
  virtual FourVector<GammaStructure> GetCurrent(int iso,const TKinematics& tk) const;

  // Setters
  virtual void SetStrangeModel (const char *fit_spec_file, EInitVerbosity verbosity=kVerbose);
  void SetModelType(const int iso,const char* modelType, const bool cgln);
  void SetModelImplementation(const int iso,const bool cgln);
  void SetParticleProperty(int iso,const TString& particlename,const TString& property,double value);

 protected:
  void Streamer(TBuffer&, Version_t, TStrangeCalc*);
  void Streamer(TBuffer&, Limits&);
  void Streamer(TBuffer& R__b, Observable& obs, Version_t R__v=2);
  void Streamer(TBuffer& R__b, Data& data, Version_t R__v); 
  void Streamer(TBuffer& R__b, ClassDiagram& diagram, Version_t R__v); 
  void Streamer(TBuffer&, Varinfo&);
  void Streamer(TBuffer& R__b, FormFactor*& formfactor, Version_t R__v, int gic);
  void InitialiseStrangeModel();
  bool DoCalculation(const TKinematics&, const TCalcInfo*, double&);
  
  // Data members
  // ------------
  TStrangeCalc *fCalc;          // interface to strangecalc
  const static int kVersion=3;

  ClassDef(TStrangeModel,kVersion);    // strangecalc model


};


#endif
