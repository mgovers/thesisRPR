/*!
 * \file FormFactor.h
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

#ifndef FORMFACTOR_H
#define FORMFACTOR_H


/*! Type definition of a function pointer to a formfactor calculating function
 */
typedef double (*ffCalculator)(double,double,double,double,
			       double,double,double,double,
			       double,double,double,double,
			       double,double,double,double);

/*! \class FormFactor
 *   \brief The FormFactor class represents a general formfactor.
 * It's data members are:
 * - a pointer to the function containting the parametrization of
 *   the formfactor
 * - The formfactor parameters (max.15)
 */
class FormFactor
{
 public:
  FormFactor(double (*)(double, double, double, double,
			double, double, double, double,
			double, double, double, double,
			double, double, double, double),
	     double=0.0,double=0.0,double=0.0,
	     double=0.0,double=0.0,double=0.0,double=0.0,
	     double=0.0,double=0.0,double=0.0,double=0.0,
	     double=0.0,double=0.0,double=0.0,double=0.0); // constructor
  FormFactor& operator= (const FormFactor& toCopy);    // explicit assignment operator (cf. TStrangeCalc)
  FormFactor( const FormFactor& toCopy );
  double value(double, double, double = 0.0, double = 0.0, double = 0.0, double = 0.0); // Calculate formfactor for specific Q^2

  void setParameters(double=0.0,double=0.0,double=0.0,
		     double=0.0,double=0.0,double=0.0,double=0.0,
		     double=0.0,double=0.0,double=0.0,double=0.0,
		     double=0.0,double=0.0,double=0.0,double=0.0); // set new parameters
  double GetParameter(int);
  double Getp2() {return GetParameter(2);} // return second parameter (first parameter = mass!)
  ffCalculator GetParametrization() const { return ffParametrization; }

 private:
  // p1 -> p15 are the formfactor parameters
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  double p7;
  double p8;
  double p9;
  double p10;
  double p11;
  double p12;
  double p13;
  double p14;
  double p15;

  // Pointer to the function containing the formfactor parametrization. 
  ffCalculator ffParametrization; //DO NOT DELETE!!
};

#endif
