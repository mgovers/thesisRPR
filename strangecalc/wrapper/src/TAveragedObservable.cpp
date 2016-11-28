/*!
 * \file TAveragedObservable.cpp
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

#include "TStrangeModel.h"
#include "TAveragedObservable.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <RVersion.h>
#include <Math/AdaptiveIntegratorMultiDim.h>
#include <Math/GSLIntegrator.h>
#include <cmath>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////
//
// OneDimWrapperMultiDim
//
// Wraps an object of type IBaseFunctionMultiDim and implements the interface
// of a IBaseFunctionOneDim ABC. This is usefull to pass integrands of type
// TAveragedObservable::Integrand to a 1-dimensional fitter in
// TAveragedObservable::Integrate.
//
///////////////////////////////////////////////////////////////////////////////

class OneDimWrapperMultiDim : public ROOT::Math::IBaseFunctionOneDim
{
public:
  OneDimWrapperMultiDim(const ROOT::Math::IBaseFunctionMultiDim& multiDim) : fParent(multiDim) { assert(multiDim.NDim()==1); }
  virtual OneDimWrapperMultiDim* Clone() const { return 0; }

private:
  virtual double DoEval(double x) const { return fParent(&x); }

private:
  const ROOT::Math::IBaseFunctionMultiDim& fParent; // reference to the multi-dimensional function

}; // class OneDimWrapperMultiDim


///////////////////////////////////////////////////////////////////////////////
//
// TAveragedObservable
//
// This class addresses the issue where calculations are needed for observables
// that are binned over a set of kinematic variables.
//
// The type of observable and the kinematics of the observable are initially 
// specified in the constructor. In a next step, the user can:
// - change the value of an independent kinematic variable with SetVar.
// - set a range for an independent kinematic variable with SetVarRange.
// - set a binning range for an independent kinematic variable with SetBinningRange.
//
// Results can be obtained with 3 different methods:
// - MakeBinCenteredGraph:
//   Calculates observables in the bin centers of the binned variables for all 
//   kinematical points that have been given a range.
// - MakeBinAveragedGraph:
//   Calculates observables averaged over the binned variables for all kinematical
//   points that have been given a range.
// - MakeBinAveragedGraphWithErrors:
//   Same as MakeBinAveragedGraph, but in addition each point has an error bar that
//   corresponds to the standard deviation of the observable in its bin.
//
// The parameters for the numerical integrals are set at compile time with the constants:
// - kAbsoluteTolerance: absolute precision of numerical integrals
// - kRelativeTolerance: relative precision of numerical integrals
// - kMultiDimMaximumPoints: maximum number of function evaluations for multi-dimensional integral
// - kWorkingArraySize: number of sub-division used for calculating numerical integrals
// - kOneDimIntegrationType: integration type for the 1-dimensional integrals
//
///////////////////////////////////////////////////////////////////////////////

/*! 
 * \class TAveragedObservable
 * \brief See http://rprmodel.ugent.be/api/root/TAveragedObservable.html
 */

ClassImp(TAveragedObservable)

//_____________________________________________________________________
// Absolute and relative precision of numerical integrals.
// Because of the bug in AdaptiveIntegratorMultiDim the absolute precision
// should always be zero.
const double TAveragedObservable::kAbsoluteTolerance = 0.;
const double TAveragedObservable::kRelativeTolerance = .5e-3;

//_____________________________________________________________________
// Parameters for the integration
// - Maximum number of function evaluations for multi-dimensional integral
// - Size of working array, i.e. number of sub-division used for calculating the integral.
const unsigned int TAveragedObservable::kMultiDimMaximumPoints = 1000000;
const unsigned int TAveragedObservable::kWorkingArraySize = 100;

//_____________________________________________________________________
// The integration type for the 1-dimensional integrals
const ROOT::Math::IntegrationOneDim::Type TAveragedObservable::kOneDimIntegrationType = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;

//_____________________________________________________________________
TAveragedObservable::TAveragedObservable(const TCalcInfo& calcInfo, const char* formatString, double var1, double var2, double var3)
  : TObject()
  , fVarState()
  , fLowerBounds()
  , fUpperBounds()
  , fCenteredKinematics( new TKinematics("","",calcInfo.GetIsospin(),formatString,var1,var2,var3) )
  , fCalcInfo( calcInfo.Clone() )
{
  // Constructor
  //
  // The type of observable is specified by the 'calcInfo' argument. The kinematics are
  // fixed by providing a 'formatString' and values for the 3 independent kinematic variables
  // as argument (similar to the contructors of TKinematics). See the documentation of
  // TKinematics::SetFormat for info about what constitutes a valid format string.
  //
  // The isospin channel of the observable is taken from the 'calcInfo' argument.
  for(int i=0; i<3; ++i) {
    fVarState[i] = kFixed;
    fLowerBounds[i] = 0.;
    fUpperBounds[i] = 0.;
  }
}

//_____________________________________________________________________
TAveragedObservable::TAveragedObservable(TRootIOCtor *rio)
  : TObject()
  , fVarState()
  , fLowerBounds()
  , fUpperBounds()
  , fCenteredKinematics(0)
  , fCalcInfo(0)
{
  // ROOT I/O Constructor
  for(int i=0; i<3; ++i) {
    fVarState[i] = kFixed;
    fLowerBounds[i] = 0.;
    fUpperBounds[i] = 0.;
  }
}

//_____________________________________________________________________
TAveragedObservable::TAveragedObservable(const TAveragedObservable& rhs)
  : TObject(rhs)
  , fVarState()
  , fLowerBounds()
  , fUpperBounds()
  , fCenteredKinematics( rhs.fCenteredKinematics->Clone() )
  , fCalcInfo( rhs.fCalcInfo->Clone() )
{
  // Copy constructor
  for(int i=0; i<3; ++i) {
    fVarState[i] = rhs.fVarState[i];
    fLowerBounds[i] = rhs.fLowerBounds[i];
    fUpperBounds[i] = rhs.fUpperBounds[i];
  }
}

//_____________________________________________________________________
TAveragedObservable::~TAveragedObservable()
{
  // Destructor
  delete fCenteredKinematics;
  delete fCalcInfo;
}

//_____________________________________________________________________
TAveragedObservable& TAveragedObservable::operator=(const TAveragedObservable& rhs)
{
  // Assignment
  if( this!=&rhs ) { // avoid self-assignment
    TObject::operator=(rhs);
    for(int i=0; i<3; ++i) {
      fVarState[i] = rhs.fVarState[i];
      fLowerBounds[i] = rhs.fLowerBounds[i];
      fUpperBounds[i] = rhs.fUpperBounds[i];
    }
    *fCenteredKinematics = *rhs.fCenteredKinematics;
    *fCalcInfo = *rhs.fCalcInfo;
  }
  
  return *this;
}

//_____________________________________________________________________
TAveragedObservable* TAveragedObservable::Clone(const char *newname) const
{
  // Virtual copy constructor
  return new TAveragedObservable(*this);
}

//_____________________________________________________________________
void TAveragedObservable::SetCalcInfo(const TCalcInfo& calcInfo)
{
  // Specify a new type of observable
  *fCalcInfo = calcInfo;
  fCenteredKinematics->SetIsospin( fCalcInfo->GetIsospin() );
}

//_____________________________________________________________________
void TAveragedObservable::SetVar(int varNum, double value)
{
  // Fixes kinematic variable number 'varNum' to a single value
  if( IsValidVariableNumber(varNum) ) {
    fVarState[varNum-1] = kFixed;

    fCenteredKinematics->FixVariable(varNum);
    fCenteredKinematics->SetVar(varNum,value);
  }
}

//_____________________________________________________________________
void TAveragedObservable::SetVarRange(int varNum, double lower, double upper, int steps)
{
  // Gives kinematic variable number 'varNum' a range in the interval [lower,upper] taking
  // 'steps' steps.
  if( IsValidVariableNumber(varNum) ) {
    fVarState[varNum-1] = kVariable;

    fCenteredKinematics->SetVarRange(varNum,lower,upper,steps);
  }
}

//_____________________________________________________________________
void TAveragedObservable::SetBinningRange(int varNum, double lower, double upper)
{
  // Specifies a binning interval [lower,upper] for kinematic variable number 'varNum'.
  if( IsValidVariableNumber(varNum) ) {
    fVarState[varNum-1] = kBinned;

    fLowerBounds[varNum-1] = lower;
    fUpperBounds[varNum-1] = upper;
    
    fCenteredKinematics->FixVariable(varNum);
    fCenteredKinematics->SetVar(varNum,(upper+lower)/2.);
  }
}

//_____________________________________________________________________
unsigned int TAveragedObservable::NumberOfBinnedVariables() const
{
  // Returns the number of independent kinematic variables that are given a binning range ([0,3])
  unsigned int count = 0;

  for(int varNum=0; varNum<3; ++varNum) {
    if( fVarState[varNum] == TAveragedObservable::kBinned )
      ++count;
  }

  return count;
}

//_____________________________________________________________________
bool TAveragedObservable::IsValidVariableNumber(int varNum)
{
  // Checks whether 'varNum' is a valid variable number, i.e. [1,3]. 
  // If not a warning is printed.
  if( varNum>0 && varNum<=3 ) return true;

  // else
  std::cerr << "WARNING in TAveragedObservable::IsValidVariableNumber(int): "
	    << "invalid variable number."
	    << std::endl;

  return false;
}

//_____________________________________________________________________
int TAveragedObservable::GetNumberOfSteps() const
{
  // Number of points in the grid of kinematic variables
  return fCenteredKinematics->GetNumberOfSteps();
}

//_____________________________________________________________________
int TAveragedObservable::GetNumberOfPhysicalSteps() const
{
  // Number of physical points in the grid of kinematic variables. Physical
  // points are evaluated for binned variables at their bin center.
  return fCenteredKinematics->GetNumberOfPhysicalSteps();
}

//_____________________________________________________________________
TGraph* TAveragedObservable::MakeBinCenteredGraph(TStrangeModel& model, const TString& xvariable, double xunit, double yunit) const
{
  // Create a graph of the observable using 'model' in the physical kinematic points specified
  // by setting ranges with SetVarRange as a function of the kinematic variable 'xvariable'.
  // For binned kinematic variables, calculations are performed at the bin center.
  //
  // For a list of possibilities for 'xvariable' see the documentation of 
  // TKinematics::GetVarArray(const TString&).
  // 
  // The x- and y-values are scaled by 'xunit' and 'yunit' before being drawn.
  //
  // The user owns the returned TGraph.
  return model.MakeGraph(*fCenteredKinematics,fCalcInfo,xvariable,xunit,yunit);
}

//_____________________________________________________________________
double TAveragedObservable::Integrate(const Integrand& integrand, double* lower, double* upper) const
{
  // [Protected]
  // Integrates a given integrand over the bounds given in the argument.
  // This method deals with the hassle over integrating 1- or multi-dimensional integrands.

  // For 1-dimensional integrands
  if( NumberOfBinnedVariables()==1 ) {
    ROOT::Math::GSLIntegrator onedimAlgorithm(kOneDimIntegrationType,
					      kAbsoluteTolerance,
					      kRelativeTolerance,
					      kWorkingArraySize);
    
    return onedimAlgorithm.Integral(OneDimWrapperMultiDim(integrand),lower[0],upper[0]);
  }

  // else for multi-dimensional integrands
#if ROOT_SVN_REVISION >= 36764 // Deal with ROOT API change
  ROOT::Math::AdaptiveIntegratorMultiDim algorithm(integrand,
						   kAbsoluteTolerance,
						   kRelativeTolerance,
						   kMultiDimMaximumPoints,
						   kWorkingArraySize);
#else
  ROOT::Math::AdaptiveIntegratorMultiDim algorithm(integrand,
						   kAbsoluteTolerance,
						   kRelativeTolerance,
						   kMultiDimMaximumPoints);
#endif
						   
  return algorithm.Integral(lower,upper);
}

//_____________________________________________________________________
TGraph* TAveragedObservable::MakeGraph(double power, TStrangeModel& model, const TString& xvariable, double xunit, double yunit) const
{
  // [Protected]
  // Creates a graph of the observable raised to the power 'power' using 'model' in the physical
  // kinematic points specified by setting ranges with SetVarRange as a function of the
  // kinematic variable 'xvariable'. For binned kinematic variables, calculations are integrated
  // over the bins.
  //
  // For a list of possibilities for 'xvariable' see the documentation of 
  // TKinematics::GetVarArray(const TString&).
  // 
  // The x- and y-values are scaled by 'xunit' and 'yunit' before being drawn.
  //
  // The user owns the returned TGraph.
  
  // Make a non-const version of the kinematics
  TKinematics* kinematics = const_cast<TKinematics*>( fCenteredKinematics );
  const int step_position = kinematics->GetStep();

  // Initialize the integrand
  Integrand integrand(this,model,power);
  
  // Prepare the integration bounds and volume
  double volume = 1.;
  double *lowerBounds = new double[NumberOfBinnedVariables()];
  double *upperBounds = new double[NumberOfBinnedVariables()];
  for(int varNum=0, index=0; varNum<3; ++varNum) {
    if( fVarState[varNum] == TAveragedObservable::kBinned ) {
      lowerBounds[index] = fLowerBounds[varNum];
      upperBounds[index] = fUpperBounds[varNum];
      volume *= ( upperBounds[index] - lowerBounds[index] );
      ++index;
    }
  }

  // Arrays for x- and y-values
  double *x = fCenteredKinematics->GetVarArray( xvariable );
  double *y = new double[GetNumberOfPhysicalSteps()];

  // Loop over all non-binned variables and integrate
  for(int step=0, index=0; step<GetNumberOfSteps(); ++step) {
    kinematics->GoTo(step);
    if( kinematics->IsPhysical() )
      y[index++] = Integrate(integrand,lowerBounds,upperBounds) / volume;
  }

  // Stuff everything into a graph and clean up
  TGraph *g = new TGraph( GetNumberOfPhysicalSteps() , x , y );
  delete[] x;
  delete[] y;
  delete[] lowerBounds;
  delete[] upperBounds;

  // Return the kinematics to their old state
  kinematics->GoTo(step_position);

  return g;
}

//_____________________________________________________________________
TGraph* TAveragedObservable::MakeBinAveragedGraph(TStrangeModel& model, const TString& xvariable, double xunit, double yunit) const
{
  // Create a graph of the observable using 'model' in the physical kinematic points specified
  // by setting ranges with SetVarRange as a function of the kinematic variable 'xvariable'.
  // For binned kinematic variables, calculations are averaged over the bins.
  //
  // For a list of possibilities for 'xvariable' see the documentation of 
  // TKinematics::GetVarArray(const TString&).
  // 
  // The x- and y-values are scaled by 'xunit' and 'yunit' before being drawn.
  //
  // The user owns the returned TGraph.
  return MakeGraph(1.,model,xvariable,xunit,yunit);
}

//_____________________________________________________________________
TGraphErrors* TAveragedObservable::MakeBinAveragedGraphWithErrors(TStrangeModel& model, const TString& xvariable, double xunit, double yunit) const
{
  // Create a graph of the observable using 'model' in the physical kinematic points specified
  // by setting ranges with SetVarRange as a function of the kinematic variable 'xvariable'.
  // For binned kinematic variables, calculations are averaged over the bins and are a given
  // an error corresponding to the standard deviation in that bin.
  //
  // For a list of possibilities for 'xvariable' see the documentation of 
  // TKinematics::GetVarArray(const TString&).
  // 
  // The x- and y-values are scaled by 'xunit' and 'yunit' before being drawn.
  //
  // The user owns the returned TGraph.

  // Integrated observable
  TGraph *obs = MakeGraph(1.,model,xvariable,xunit,yunit);

  // Integrated observable^2
  TGraph *obs2 = MakeGraph(2.,model,xvariable,xunit,yunit);

  // Combine both
  TGraphErrors *g = new TGraphErrors(obs->GetN(),obs->GetX(),obs->GetY());
  for(int i=0; i<g->GetN(); ++i)
    g->SetPointError(i, 0., std::sqrt( obs2->GetY()[i] - obs->GetY()[i]*obs->GetY()[i] ));

  // Clean up
  delete obs;
  delete obs2;

  return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// TAveragedObservable::Integrand
//
// Helper class for TAveragedObservable by giving it the interface of a
// multi-dimensional function that can be given to a integrator.
//
///////////////////////////////////////////////////////////////////////////////

/*! 
 * \class TAveragedObservable::Integrand
 * \brief See http://rprmodel.ugent.be/api/root/TAveragedObservable__Integrand.html
 */

ClassImp(TAveragedObservable::Integrand)

//_____________________________________________________________________
TAveragedObservable::Integrand::Integrand(const TAveragedObservable* parent, TStrangeModel& model, double power)
  : ROOT::Math::IBaseFunctionMultiDim()
  , fParent( parent )
  , fModel( model )
  , fPower( power )
{
  // Constructor
  //
  // Arguments:
  // - parent to give access to the kinematics and observable
  // - an RPR model
  // - power of the observable
}

//_____________________________________________________________________
double TAveragedObservable::Integrand::DoEval(const double* x) const
{
  // Calculate the observable with 'x' fixing the kinematics of the
  // binned variables at the current value of the running variables
  // in fParent->fCenteredKinematics.
  
  // Prepare the kinematics
  TKinematics kinematics( *fParent->fCenteredKinematics );
  kinematics.FixVariables();
  for(int varNum=0, index=0; varNum<3; ++varNum) {
    if( fParent->fVarState[varNum] == TAveragedObservable::kBinned )
      kinematics.SetVar(varNum+1,x[index++]);
  }

  // Calculate
  return pow( fModel.GetCalcpoint(kinematics,fParent->fCalcInfo) , fPower );
}
