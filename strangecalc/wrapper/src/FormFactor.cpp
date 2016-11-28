/*!
 * \file FormFactor.cpp
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

#include <iostream>
#include <cstdlib>
#include "FormFactor.h"

using namespace std;


/* Constructor
 */
FormFactor::FormFactor ( double ( *function ) ( double, double, double, double,
                         double, double, double, double,
                         double, double, double, double,
                         double, double, double, double ),
                         double i1,double i2,double i3,
                         double i4,double i5,double i6,double i7,
                         double i8,double i9,double i10,double i11,
                         double i12,double i13,double i14,double i15 )
      : p1 ( i1 ), p2 ( i2 ),  p3 ( i3 ),   p4 ( i4 ),   p5 ( i5 ),   p6 ( i6 ),   p7 ( i7 ), p8 ( i8 ),
      p9 ( i9 ), p10 ( i10 ), p11 ( i11 ), p12 ( i12 ), p13 ( i13 ), p14 ( i14 ), p15 ( i15 )

{
   ffParametrization = function;
}
FormFactor::FormFactor( const FormFactor& toCopy ):
  p1 ( toCopy.p1 ), p2 ( toCopy.p2 ),  p3 ( toCopy.p3 ),   p4 ( toCopy.p4 ),   p5 ( toCopy.p5 ),   p6 ( toCopy.p6 ),   p7 ( toCopy.p7 ), p8 ( toCopy.p8 ), p9 ( toCopy.p9 ), p10 ( toCopy.p10 ), p11 ( toCopy.p11 ), p12 ( toCopy.p12 ), p13 ( toCopy.p13 ), p14 ( toCopy.p14 ), p15 ( toCopy.p15 ), ffParametrization (toCopy.ffParametrization)
// copy constructor 
{
}

FormFactor& FormFactor::operator= ( const FormFactor& toCopy ) // explicit assignment operator (cf. TStrangeCalc)
{
   if ( this != &toCopy )
   {
      p1 = toCopy.p1;
      p2 = toCopy.p2;
      p3 = toCopy.p3;
      p4 = toCopy.p4;
      p5 = toCopy.p5;
      p6 = toCopy.p6;
      p7 = toCopy.p7;
      p8 = toCopy.p8;
      p9 = toCopy.p9;
      p10 = toCopy.p10;
      p11 = toCopy.p11;
      p12 = toCopy.p12;
      p13 = toCopy.p13;
      p14 = toCopy.p14;
      p15 = toCopy.p15;

      ffParametrization = toCopy.ffParametrization;
   }
   return *this;
}

/*! Given Q^2 and the coupling constant cc in the real photon point
 * this function determines the coupling strength
 */
double FormFactor::value ( double cc, double Q2 , double par3, double par4, double par5, double par6 )
{
   // calculate the value of the formfactor
  double ff = ( *ffParametrization ) ( Q2,p1,p2,par3,par4,par5,par6,p7,p8,
                                        p9,p10,p11,p12,p13,p14,p15 );

   /* Depending on the value of the coupling constant in the real photon
    * point, we return the formfactor or its product with the coupling
    * constant */

   // when the coupling constant is non-zero, the answer is simple
   if ( cc!=0.0 )
      return cc*ff;
   /* when cc==0.0, we want to return ff
    * but first we want to make sure that ff in the limit Q2->0
    * is ==0 */
   else if ( Q2!=0.0 )
   {
      if ( value ( cc,0.0 ) ==0.0 )
         return ff;
      else
      {
         cout << "Inconsistency with form factors!" << endl
         << "The limit Q2->0 does not match "
         << "real-photon point coupling constants"
         << endl << endl;
         exit ( 1 );
      }
   }
   else
      return cc;
}


/*! Straight forward function to set the parameters 1->15 to a different value
 */
void FormFactor::setParameters ( double i1,double i2,double i3,
                                 double i4,double i5,double i6,double i7,
                                 double i8,double i9,double i10,double i11,
                                 double i12,double i13,double i14,double i15 )
{
   p1 = i1;
   p2 = i2;
   p3 = i3;
   p4 = i4;
   p5 = i5;
   p6 = i6;
   p7 = i7;
   p8 = i8;
   p9 = i9;
   p10 = i10;
   p11 = i11;
   p12 = i12;
   p13 = i13;
   p14 = i14;
   p15 = i15;
}

/*! Access the individual parameters
 */
double FormFactor::GetParameter(int p)
{
  switch(p) {
  case 1: return p1; break;
  case 2: return p2; break;
  case 3: return p3; break;
  case 4: return p4; break;
  case 5: return p5; break;
  case 6: return p6; break;
  case 7: return p7; break;
  case 8: return p8; break;
  case 9: return p9; break;
  case 10: return p10; break;
  case 11: return p11; break;
  case 12: return p12; break;
  case 13: return p13; break;
  case 14: return p14; break;
  case 15: return p15; break;
  default: 
    cerr << "ERROR in FormFactor::GetParameter(int): "
	 << "Parameter index out-of-bounds." << endl;
    exit(1);
  }

  // dummy return
  return 0.;
}

