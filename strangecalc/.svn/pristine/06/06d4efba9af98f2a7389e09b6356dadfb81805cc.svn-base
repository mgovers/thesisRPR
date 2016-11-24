/*!
 * \file Lagrangian.h
 * \ingroup wrapper
 *
 * This header declares all types of Lagrangians and propagators
 * necessary to construct the amputated current.
 *
 * We declare a number of functions (calculateDiagram*()) that put
 * together the necessary lagrangians and propagators to determine 
 * the contribution of a specific diagram to the amputated current.
 *
 * To improve speed these calculateDiagram*() functions are wrapped
 * in a TCalculateDiagram class that can cache the results if
 * memoizing is turned on (on a per-diagram basis, hard-coded).
 *
 * In this header we declare an array of TCalculateDiagram objects,
 * one for each type of diagram. We recomment users to call the 
 * operator() member function of one of these objects to calculate
 * a diagrams contribution to the matrix element.
 *
 * All references are to Stijn's notes
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 *
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

#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include "Structures.h"
#include "FourVector.h"
#include "FourVectorHandle.h"
#include "GammaStructure.h"
#include "Tensor.h"
#include <complex>
#include "TCalculateDiagram.h"


/* ************** *
 * EM Lagrangians *
 * ************** */

FourVector<GammaStructure> EMvertex_12p_12p(const Properties&,
					    const FourVector<double>&);

FourVector<GammaStructure> EMvertex_12p_12m(const Properties&,
					    const FourVector<double>&);

FourVector< std::complex<double> > EMvertex_0m_0m(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&);

TensorRank2 EMvertex_1m_0m(const Properties&,
			   const FourVector<double>&,
			   const FourVector<double>&);

TensorRank2 EMvertex_1p_0m(const Properties&,
			   const FourVector<double>&,
			   const FourVector<double>&);

TensorRank2 EMvertex_12p_32p(const Properties&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     bool =false);

TensorRank2 EMvertex_12p_32m(const Properties&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     bool =false);

TensorRank3 EMvertex_12p_52p(const Properties&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     const FourVector<double>&);

TensorRank3 EMvertex_12p_52m(const Properties&,
			     const FourVector<double>&,
			     const FourVector<double>&,
			     const FourVector<double>&);

/* ****************** *
 * Strong Lagrangians *
 * ****************** */

GammaStructure StrongVertex_12p_0m_12p(const Properties&, double);

GammaStructure StrongVertex_12p_0m_12m(const Properties&, double);

FourVector<GammaStructure> StrongVertex_12p_1m_12p(const Properties&,
						   const FourVector<double>&,
						   double);

FourVector<GammaStructure> StrongVertex_12p_1p_12p(const Properties&,
						   const FourVector<double>&,
						   double);

FourVector<GammaStructure> StrongVertex_32p_0m_12p(const Properties&,
						   const FourVector<double>&,
						   const FourVector<double>&,
						   double);

FourVector<GammaStructure> StrongVertex_32m_0m_12p(const Properties&,
						   const FourVector<double>&,
						   const FourVector<double>&,
						   double);

TensorRank2 StrongVertex_52p_0m_12p(const Properties&,
				    const FourVector<double>&,
				    const FourVector<double>&,
				    double);

TensorRank2 StrongVertex_52m_0m_12p(const Properties&,
				    const FourVector<double>&,
				    const FourVector<double>&,
				    double);

/* *********** *
 * PROPAGATORS *
 * *********** */

std::complex<double> propagatorSpin0(const Properties&,
				const FourVector<double>&,
				const bool=false);

GammaStructure propagatorSpin12(const Properties&,
				const FourVector<double>&,
				const bool=false);

TensorRank2 propagatorSpin1(const Properties&,
			    const FourVector<double>&,
				const bool=false);

TensorRank2 propagatorSpin32(const Properties&,
			     const FourVector<double>&,
				const bool=false);

TensorRank4 propagatorSpin52(const Properties&,
			     const FourVector<double>&,
			     const bool=false);

std::complex<double> propagatorRegge(const Properties&, double, double,
				     double, double, double, const Observable*);


/* ******************************* *
 * Calculate Diagram contributions *
 * ******************************* */


FourVector<GammaStructure> calculateDiagramS(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramT(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramU(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);


FourVector<GammaStructure> calculateDiagramA(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramB(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramC(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramD(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramE(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramF(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramG(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramH(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramI(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramJ(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramL(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramM(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramN(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramO(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

FourVector<GammaStructure> calculateDiagramQ(const Properties&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const FourVector<double>&,
					     const bool=false);

/*!
 * \brief Array containing diagram building functions
 *
 * Array of TCalculateDiagrams calculateDiagram[],
 * which wraps the calculateDiagram*() functions declared above.
 * Each element of the array corresponds to the correct
 * function to determine the contribution of a specific
 * diagram to the amputated current.
 * The corresponce is as follows:
 *
 *\verbatim
  (0)    S ----- p         -  s Born term, proton exchange
  (1)    T ----- K         -  t Born term, kaon exchange
  (2)    U ----- Y         -  u Born term, hyperon exchange
  (3)    A ----- Y'        -  hyperon exchange
  (4)    B ----- K*(-)     -  vector meson exchange
  (5)    C ----- K1(+)     -  vector meson exchange
  (6)    D ----- N*(1/2+)  -  nucleon resonance exchange
  (7)    E ----- N*(1/2-)  -  nucleon resonance exchange
  (8)    F ----- Y*(1/2+)  -  hyperon resonance exchange
  (9)    G ----- Y*(1/2-)  -  hyperon resonance exchange
  (10)   H ----- N*(3/2+)  -  nucleon resonance exchange
  (11)   I ----- N*(3/2-)  -  nucleon resonance exchange
  (12)   J ----- Y*(3/2+)  -  hyperon resonance exchange
  (13)   L ----- Y*(3/2-)  -  hyperon resonance exchange
  (14)   M ----- N*(5/2+)  -  nucleon resonance exchange
  (15)   N ----- N*(5/2-)  -  nucleon resonance exchange
  (16)   O ----- Y*(5/2+)  -  hyperon resonance exchange
  (17)   Q ----- Y*(5/2-)  -  hyperon resonance exchange
 \endverbatim
 */
const TCalculateDiagram calculateDiagram[CLASSMAX_NONCGLN] 
= { TCalculateDiagram(&calculateDiagramS,true),
    TCalculateDiagram(&calculateDiagramT,true),
    TCalculateDiagram(&calculateDiagramU,true),
    TCalculateDiagram(&calculateDiagramA,true),
    TCalculateDiagram(&calculateDiagramB,true),
    TCalculateDiagram(&calculateDiagramC,true),
    TCalculateDiagram(&calculateDiagramD,true),
    TCalculateDiagram(&calculateDiagramE,true),
    TCalculateDiagram(&calculateDiagramF,true),
    TCalculateDiagram(&calculateDiagramG,true),
    TCalculateDiagram(&calculateDiagramH,true),
    TCalculateDiagram(&calculateDiagramI,true),
    TCalculateDiagram(&calculateDiagramJ,true),
    TCalculateDiagram(&calculateDiagramL,true),
    TCalculateDiagram(&calculateDiagramM,true),
    TCalculateDiagram(&calculateDiagramN,true),
    TCalculateDiagram(&calculateDiagramO,true),
    TCalculateDiagram(&calculateDiagramQ,true) };

#endif
