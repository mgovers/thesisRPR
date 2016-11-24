/*!
 * \file calcMatrixElement.h
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

#include <Structures.h>
#include "TMatrixElement.h"

#ifndef CALCMATRIXELEMENT_H
#define CALCMATRIXELEMENT_H


/*! \typedef Response_func_calc
 * This structure tells calc_response_func()
 * which response functions need to be determined.
 */
typedef struct
{
  short L;
  short T;
  short c_TT;
  short s_TT;
  short c_TL;
  short s_TL;
  short TT_pol;
  short c_TL_pol;
  short s_TL_pol;
} Response_func_calc;


double determine_M2_photo(double, double,double, double, 
			  Class[], Observable*, int label);

double determine_M2_photo_withME( TMatrixElement*, Observable*);

int calc_response_func(Response_func*,Response_func_calc*,
		       double,double,double,
		       double,double,double,Class*,Observable*,int label);

TMatrixElement* construct_matrixelement(double,double,double,
					double,Class*,Observable*,int label);

void get_transversity_amplitude( TMatrixElement*, double*, int);

void get_helicity_amplitude( TMatrixElement*, double*, int);

int release_matrixelement( TMatrixElement*);

int updateSpinDependencies_matrixelement( TMatrixElement*, int* label);

void set_datasize(int size);
       
#endif
