/*!
 * \file TCalculateDiagram.cpp
 * \ingroup wrapper
 * \date started on February,22 2010
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

#include "TCalculateDiagram.h"

//________________________________________________________________________
/*! \brief Constructor
 * \param func A function which calculates a particle-exchange diagram contribution.
 *             These are declared in Lagrangian.h.
 * \param memoizeDiagram whether or not to cache func.
 */
TCalculateDiagram::TCalculateDiagram(diagramBuilder func,
				     bool memoizeDiagram)
: fMemoizeDiagram(memoizeDiagram), fCachemap(TCachemap<FourVector<GammaStructure> >()),
  fCalculateDiagram_func(func)
{
  // empty body
}

//________________________________________________________________________
/*! \brief Evaluate the diagram contribution for a set of coupling constants and
 *        kinematics.
 * \param particle Properties of the exchanged particle.
 * \param k4vect Photon four-momentum in CM frame.
 * \param p4vect Proton four-momentum in CM frame.
 * \param pK4vect Kaon four-momentum in CM frame.
 * \param pY4vect Hyperon four-momentum in CM frame.
 */
FourVector<GammaStructure> 
TCalculateDiagram::operator()(const Properties& particle,
			      const FourVector<double>& k4vect,
			      const FourVector<double>& p4vect,
			      const FourVector<double>& pK4vect,
			      const FourVector<double>& pY4vect,
			      const bool kaoncapture) const
{
    if(!fMemoizeDiagram)
      //no memoizing: simply evaluate the function.
      return (*fCalculateDiagram_func)(particle,k4vect,p4vect,pK4vect,pY4vect,kaoncapture);
    
    else {
      // Variable formfactor parameters are res or born cutoff in strongFF.
      // Diagrams S,T,U, A,D, E,F,G: formfactorH->getp1()
      // Diagrams B,C: Regge terms
      // Diagrams H,I,J,L: formfactorI->getp1()
      
      double strongff_parameter = 4294967295;
      if(particle.formfactorH != NULL)
	strongff_parameter  = particle.formfactorH->Getp2();
      if(particle.formfactorI != NULL)
	strongff_parameter  = particle.formfactorI->Getp2();
      
      // Other parameters needed for this diagam are: 
      // G, H, I, X, Y, Z, mass & width of particle. (=8 doubles)
      // 4 kinematical vectors (=16 doubles),
      
      float argumentArray[25]= {(float)particle.G,(float)particle.H,(float)particle.I,(float)particle.X,(float)particle.Y,(float)particle.Z,
				(float)particle.mass, (float)particle.width, (float)strongff_parameter,
				(float)k4vect[0],(float)p4vect[0],(float)pK4vect[0],(float)pY4vect[0],(float)k4vect[1],(float)p4vect[1],(float)pK4vect[1],(float)pY4vect[1],
				(float)k4vect[2],(float)p4vect[2],(float)pK4vect[2],(float)pY4vect[2], (float)k4vect[3],(float)p4vect[3],(float)pK4vect[3],(float)pY4vect[3]};
      
      // Create unique key for cacheMap lookups
      char cacheKey[sizeof(float)*25*2 + 1];
      floatstochar(25, argumentArray, cacheKey);
      
      // Create a temp object to hold the result
      FourVector<GammaStructure> result;

      // Look up cachekey in cacheMap and return its associated value if found.
      if( fCachemap.find(cacheKey,result) ) return result;
      
      // if cacheKey was not found in the map, we evaluate the function
      result = (*fCalculateDiagram_func)(particle,k4vect,p4vect,pK4vect,pY4vect,kaoncapture);
      
      // and store its result in cacheMap under cacheKey.
      fCachemap.insert(cacheKey,result);
      return result;

    }
}
