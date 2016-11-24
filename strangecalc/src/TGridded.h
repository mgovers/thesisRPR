/*!
 * \file TGridded.h
 * \ingroup wrapper
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \date started on December 30, 2010
 
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

#ifndef TGRIDDED_H
#define TGRIDDED_H

#include "TStrangeModel.h"
#include "TMultiModel.h"
#include <FourVector.h>
#include <map>

class GammaStructure;
class TKinematics;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TGridded<Model>                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

template <typename Model>
class TGridded : public Model
{
 public:
  TGridded(const char* name="", const char* title="");
  TGridded(TRootIOCtor* rio); // ROOT I/O Constructor
  TGridded(const TGridded&);
  TGridded& operator=(const TGridded&);
  virtual TGridded* Clone(const char* newname="") const;
  virtual ~TGridded();

  virtual FourVector<GammaStructure> GetCurrent(int iso,const TKinematics& tk) const;
  void SetMaximumMemory(unsigned int maxMem);
  void ClearCache();

 private:
  static UInt_t IndexCosthkcm(double costhkcm);
  static UInt_t IndexW(double w);
  static UInt_t Index(int iso, double w, double costhkcm);

 private:
  mutable std::map<UInt_t,GammaStructure*> fCurrents; // pre-calculated currents
  unsigned int fMaxSize; // Maximum size of fCurrents

  static const double kWStep; // Step size for invariant mass caching
  static const double kWThreshold; // Smallest invariant mass (Pi0 P threshold)
  static const int kCosGridSize; // Grid size for scattering angle caching
  static const double kCosThreshold; // scattering angle

  double fQsquared; // Q^2 of last calculated currents

  ClassDef(TGridded,1); // strangecalc model whose current is calculated in a grid

}; // class TGridded<Model>

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TGridded<TStrangeModel>                                               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

typedef TGridded<TStrangeModel> TGriddedStrangeModel; // TStrangeModel whose current is calculated in a grid

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TGridded<TMultiModel>                                                 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

typedef TGridded<TMultiModel> TGriddedMultiModel; // TMultiModel whose current is calculated in a grid

#endif // TGRIDDED_H
