/*!
 * \file TGridded.cpp
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

#include "TGridded.h"
#include <Structures.h>
#include <GammaStructure.h>
#include <TKinematics.h>
#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
using std::map;
using std::fabs;
using std::complex;
using std::size_t;

///////////////////////////////////////////////////////////////////////////
//
// TGridded<Model>
//
// This template class is a derived class of classes of type TStrangeModel
// where the GetCurrent method is overridden.
//
// This class keeps a cache of calculated currents in memory using a grid in
// the invariant mass 'W' and the CM scattering angle 'costhkcm'. The size
// of the grid for 'costhkcm' is fixed by kCosGridSize (at compile time). 
// For the 'W' grid we only specify the smallest step in 'W' for which we
// expect only negligable changes in the current. This is set with kWStep (at
// compile time).
//
// The current is obviously also function of Q^2. We don't cache in this 
// variable however. Each time a calculation is requested for a new Q^2,
// the cache is emptied.
//
// The calculated currents consume quite some memory. The user can specify
// the maximum memory each instance of this class can use with SetMaximumMemory.
// 
// The user can force a reset at any time with ClearCache.
//
///////////////////////////////////////////////////////////////////////////

/*!
 * \class TGridded
 *
 * See http://rprmodel.ugent.be/api/root/TGridded_TStrangeModel_.html or http://rprmodel.ugent.be/api/root/TGridded_TMultiModel_.html
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \date started on December 30, 2010
 *
 */

templateClassImp(TGridded)

//______________________________________________________________________________
// Step size for invariant mass caching
template <typename Model>
const double TGridded<Model>::kWStep = 1.; 

//______________________________________________________________________________
// Smallest invariant mass (Pi0 P threshold)
template <typename Model>
const double TGridded<Model>::kWThreshold = 1073.24;

//______________________________________________________________________________
// Grid size for scattering angle caching
template <typename Model>
const int TGridded<Model>::kCosGridSize = 201; 

//______________________________________________________________________________
// Smallest scattering angle
template <typename Model>
const double TGridded<Model>::kCosThreshold = -1.;


//______________________________________________________________________________
template <typename Model>
TGridded<Model>::TGridded(const char* name, const char* title)
: Model(name,title), fCurrents(), fMaxSize(640000), fQsquared(0.)
{
  // Constructor
  //
  // By default, the maximum memory usage is ~820MB. This can be changed
  // using SetMaximumMemory.
}

//______________________________________________________________________________
template <typename Model>
TGridded<Model>::TGridded(TRootIOCtor* rio)
  : Model(rio), fCurrents(), fMaxSize(0), fQsquared(0.)
{
  // ROOT I/O Constructor
}

//______________________________________________________________________________
template <typename Model>
TGridded<Model>::TGridded(const TGridded& rhs)
  : Model(rhs), fCurrents(), fMaxSize(rhs.fMaxSize), fQsquared(rhs.fQsquared)
{
  // Copy Constructor
  //
  // The entire cache will be copied. Make sure you have enough memory.
  map<UInt_t,GammaStructure*>::const_iterator it;
  for(it=rhs.fCurrents.begin(); it!=rhs.fCurrents.end(); ++it) {
    fCurrents[ it->first ] = new GammaStructure[4];
    for(int i=0; i<4; ++i) 
      fCurrents[ it->first ][i] = it->second[i];
  }
}

//______________________________________________________________________________
template <typename Model>
TGridded<Model>& TGridded<Model>::operator=(const TGridded& rhs)
{
  // Assignment
  //
  // The entire cache will be copied. Make sure you have enough memory.
  if( this!=&rhs ) {
    Model::operator=(rhs);
    fMaxSize = rhs.fMaxSize;
    fQsquared = rhs.fQsquared;

    // empty the map
    ClearCache();

    // and refill it
    map<UInt_t,GammaStructure*>::const_iterator cit;
    for(cit=rhs.fCurrents.begin(); cit!=rhs.fCurrents.end(); ++cit) {
      fCurrents[ cit->first ] = new GammaStructure[4];
      for(int i=0; i<4; ++i) 
	fCurrents[ cit->first ][i] = cit->second[i];
    }
  }

  return *this;
}

//______________________________________________________________________________
template <typename Model>
TGridded<Model>* TGridded<Model>::Clone(const char* newname) const
{
  // Virtual copy constructor
  //
  // The entire cache will be copied. Make sure you have enough memory.
  TGridded *clone = new TGridded(*this);
  if (newname && std::strlen(newname)) clone->SetName(newname);
  return clone;
}

//______________________________________________________________________________
template <typename Model>
TGridded<Model>::~TGridded()
{
  // Destructor
  ClearCache();
}

//______________________________________________________________________________
template <typename Model>
void TGridded<Model>::ClearCache()
{
  // Remove all previously calculated currents
  map<UInt_t,GammaStructure*>::iterator it;
  for(it=fCurrents.begin(); it!=fCurrents.end(); ++it)
    delete[] it->second;
  fCurrents.clear();
}

//______________________________________________________________________________
template <typename Model>
void TGridded<Model>::SetMaximumMemory(unsigned int maxMem)
{
  // Specify the maximum amount of memory this object can use (in MB).
  // Keep in mind this is an upper bound. The cacke can be emptied before
  // the maximum memory is reached.
  const size_t sizeGammaStructure 
    = 16 * ( sizeof(complex<double>) + sizeof(int) );

  fMaxSize = maxMem*1000000 / (4*sizeGammaStructure);
}

//______________________________________________________________________________
template <typename Model>
UInt_t TGridded<Model>::IndexCosthkcm(double cos)
{
  // Given the CM scattering angle of the kaon returns a unique key

  // Calculate the atomic step for our cosine-grid
  static const double cosStep = 2./(kCosGridSize-1.);

  return static_cast<int>( floor( ((cos-kCosThreshold)/cosStep) +.5 ) );
}

//______________________________________________________________________________
template <typename Model>
UInt_t TGridded<Model>::IndexW(double w)
{
  // Given the invariant mass returns a unique key
  return static_cast<int>( floor( ((w-kWThreshold)/kWStep) +.5 ) );
}

//______________________________________________________________________________
template <typename Model>
UInt_t TGridded<Model>::Index(int iso, double w, double cos)
{
  // Given the isospin channel and the independent kinematic variables
  // returns a unique key to store the current in the cache
  return iso + (ISOMAX-1)*(IndexCosthkcm(cos) + kCosGridSize*IndexW(w));
}

//______________________________________________________________________________
template <typename Model>
FourVector<GammaStructure> TGridded<Model>::GetCurrent(int iso,
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

  const double w = tk.GetW();
  const double costhkcm = tk.GetCosthkcm();
  const UInt_t index = Index(iso,w,costhkcm);

  // Empty the cache when 
  // - it's taking to much memory
  // - we are calculating in a different Q^2 point
  if( (fCurrents.size() > fMaxSize) || 
      (fabs(tk.GetQsquared()-fQsquared)>STRANGEUFLOW) ) {
    TGridded* nonconst_copy = const_cast<TGridded*>(this);
    nonconst_copy->ClearCache();
    const_cast<double&>(fQsquared) = tk.GetQsquared();
  }

  // Look for our index
  map<UInt_t,GammaStructure*>::const_iterator lookup = fCurrents.find(index);

  // When it is in our cache, return
  if( lookup!=fCurrents.end() ) {
    return FourVector<GammaStructure>( lookup->second[0],
				       lookup->second[1],
				       lookup->second[2],
				       lookup->second[3] );
  }

  // If not, calculate the current..
  static const double cosStep = 2./(kCosGridSize-1.);
  TKinematics kinematics("","",iso,"w:costhkcm:qsquared",
			 IndexW(w)*kWStep+kWThreshold,
			 IndexCosthkcm(costhkcm)*cosStep+kCosThreshold,
			 fQsquared);
  FourVector<GammaStructure> current = Model::GetCurrent(iso,kinematics);

  // ..and store it in the cache
  GammaStructure *newCurrent = new GammaStructure[4];
  for(int i=0; i<4; ++i)
    newCurrent[i] = current[i];
  fCurrents[ index ] = newCurrent;

  // Finally return the current
  return current;
}

//______________________________________________________________________________
// Define the template specialiazations we want
template class TGridded<TStrangeModel>;
template class TGridded<TMultiModel>;
