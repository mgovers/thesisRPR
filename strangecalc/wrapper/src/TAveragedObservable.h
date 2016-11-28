/*!
 * \file TAveragedObservable.h
 * \ingroup wrapper
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \date Work started on March, 1 2012
 
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

#ifndef TAVERAGEDOBSERVABLE_H
#define TAVERAGEDOBSERVABLE_H

#include <TObject.h>
#include <Math/IFunction.h>
#include <Math/IntegrationTypes.h>
#include "TKinematics.h"
#include "TCalcInfo.h"
class TStrangeModel;
class TGraph;
class TGraphErrors;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TAveragedObservable                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TAveragedObservable : public TObject
{
 public:
  class Integrand; // forward declaration

 public:
  TAveragedObservable(const TCalcInfo&, const char* formatString, double var1, double var2, double var3);
  TAveragedObservable(TRootIOCtor*); // ROOT I/O constructor
  TAveragedObservable(const TAveragedObservable&);
  virtual ~TAveragedObservable();
  TAveragedObservable& operator=(const TAveragedObservable&);
  virtual TAveragedObservable* Clone(const char *newname="") const;

  const TKinematics& GetKinematics() const { return *fCenteredKinematics; }
  void SetVar(int varNum, double value);
  void SetVarRange(int varNum, double lower, double upper, int steps);
  void SetBinningRange(int varNum, double lower, double upper);
  int GetNumberOfSteps() const;
  int GetNumberOfPhysicalSteps() const;

  const TCalcInfo& GetCalcInfo() const { return *fCalcInfo; }
  void SetCalcInfo(const TCalcInfo&);

  TGraph* MakeBinCenteredGraph(TStrangeModel&, const TString& xvariable, double xunit=1.0, double yunit=1.0) const; // user owns TGraph*
  TGraph* MakeBinAveragedGraph(TStrangeModel&, const TString& xvariable, double xunit=1.0, double yunit=1.0) const; // user owns TGraph*
  TGraphErrors* MakeBinAveragedGraphWithErrors(TStrangeModel&, const TString& xvariable, double xunit=1.0, double yunit=1.0) const; // user owns TGraphErrors*

 protected:
  static bool IsValidVariableNumber(int varNum);
  unsigned int NumberOfBinnedVariables() const;
  TGraph* MakeGraph(double,TStrangeModel&, const TString& xvariable, double xunit=1.0, double yunit=1.0) const; // user owns TGraphErrors*
  double Integrate(const Integrand&, double* lower, double* upper) const;

 private:
  enum VarState {
    kFixed = 1,
    kVariable = 2,
    kBinned = 3
  }; // State of independant kinematics variables
  
  VarState fVarState[3]; //[3] Keeps track of the state of each kinematic variable
  double fLowerBounds[3]; //[3] Lower bound of kinematic bins
  double fUpperBounds[3]; //[3] Upper bound of kinematic bins
  TKinematics *fCenteredKinematics; // Bin-centered kinematics
  TCalcInfo *fCalcInfo; // Specification of observable

  const static double kAbsoluteTolerance; // absolute precision of numerical integrals
  const static double kRelativeTolerance; // relative precision of numerical integrals
  const static unsigned int kMultiDimMaximumPoints; // Maximum number of function evaluations for multi-dimensional integral
  const static unsigned int kWorkingArraySize; // number of sub-division used for calculating numerical integrals
  const static ROOT::Math::IntegrationOneDim::Type kOneDimIntegrationType; // integration type for the 1-dimensional integrals
  

  ClassDef(TAveragedObservable,1); // Observables averaged over binned kinematic variables

  ///////////////////////////////////////////////////////////////////////////////
  //                                                                           //
  // TAveragedObservable::Integrand                                            //
  //                                                                           //
  ///////////////////////////////////////////////////////////////////////////////

 public:
  class Integrand : public ROOT::Math::IBaseFunctionMultiDim
    {
    // The parent class is a friend because it needs access to the constructor
    friend class TAveragedObservable;

    private:
    Integrand(const TAveragedObservable*, TStrangeModel&, double power);

    public:
      virtual Integrand* Clone() const { return 0; }
      virtual unsigned int NDim() const { return fParent->NumberOfBinnedVariables(); }

    private:
      virtual double DoEval(const double* x) const;

    private: // hide copy constructor and assignment
      Integrand(const Integrand&); // hidden
      Integrand& operator=(const Integrand&); // hidden
      
    protected:
      const TAveragedObservable* fParent; // observable to be integrated
      TStrangeModel& fModel; // the RPR model used to calculate observable
      const double fPower; // power of the observable to be integrated

      ClassDef(Integrand,0); // Helper class to integrate observables 

    }; // class TAveragedObservable::Integrand

}; // class TAveragedObservable

#endif // TAVERAGEDOBSERVABLE_H
