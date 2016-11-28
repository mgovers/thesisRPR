/*!
 * \file TCalculateDiagram.h
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

#ifndef __TCALCULATEDIAGRAM_H__
#define __TCALCULATEDIAGRAM_H__

#include <Structures.h>
#include <TCachemap.h>
#include <FourVector.h>
#include <numtoa.h>
#include <GammaStructure.h>

//________________________________________________________________________
/*! \class TCalculateDiagram
 * \brief Wrapper for calculation of particle-exchange diagrams.
 *
 * In Lagrangian.h, we declare function to evaluate the contribution of a
 * type of particle-exchange diagram to the matrix element.
 * To allow speed-up, we can cache this contribution, hence this wrapper.
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 */
class TCalculateDiagram
{
 private:
  /*! \typedef diagramBuilder
   * Definition of a particle-exchange diagram-contribution calculation function.
   */
  typedef FourVector<GammaStructure> (*diagramBuilder)(const Properties&,
						       const FourVector<double>&,
						       const FourVector<double>&,
						       const FourVector<double>&,
						       const FourVector<double>&,
						       const bool);

 public:
  TCalculateDiagram(diagramBuilder, bool);
  FourVector<GammaStructure> operator()(const Properties&, const FourVector<double>&, 
					const FourVector<double>&, const FourVector<double>&,
					const FourVector<double>&, const bool=false) const;
  
 private:
  bool fMemoizeDiagram; //!< Will we be caching this diagram?
  mutable TCachemap <FourVector < GammaStructure > >  fCachemap; //!< the cache container
  diagramBuilder fCalculateDiagram_func; //!< function to cache
  
}; // class TCalculateDiagram

#endif // __TCALCULATEDIAGRAM_H__
