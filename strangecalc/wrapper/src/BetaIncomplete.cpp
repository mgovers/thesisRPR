/*!
 * \file BetaIncomplete.cpp
 * \ingroup wrapper
 *
 * Implementation of the Generalized Incomplete Beta Function.
 * 
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

#include "BetaIncomplete.h"
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

using namespace std;

#define pi 4*atan(1)
#define x_limit 0.5
#define n_max 500


/*! 
 * Removes numerical noise from x which results from adding or subtracting doubles. 
 * This method is only effective when x should equal zero but instead equals a multiple 
 * or an inverse multiple of epsilon, with epsilon being the difference between 1 and
 * the smallest value greater than 1 for double objects.
 */
double remove_noise(double x)
{
  double epsilon = numeric_limits<double>::epsilon(), n = x/epsilon;

  return ((n == (int)n && n < n_max) || (1/n == (int)(1/n) && 1/n < n_max) ? 0 : x);
}


/*! Checks whether x is an integer or not.
 */
bool is_integer(double x)
{
  return remove_noise(x - round(x)) == 0;
}


/*! Computes the value of exp(i*theta).
 */
complex<double> phase(double theta)
{
  complex<double> result;

  if(is_integer(theta))
    result = pow(-1,round(theta));
  else
    result = exp(complex<double>(0,pi*theta));

  return result;
}


/*! Computes the value of u^v for real values of u and v.
 */
complex<double> cpow(double u, double v)
{
  return pow(abs(u),v)*(u < 0 ? phase(v) : 1);
}


/*! Computes the difference in order of magnitude between u and v.
 */
int magnitude_difference(double u, double v)
{
  double d;
  int n, m;

  d = log(abs(u))/log(10);
  n = (is_integer(d) ? round(d) : (int)d) - (abs(u) < 1 ? 1 : 0);

  d = log(abs(v))/log(10);
  m = (is_integer(d) ? round(d) : (int)d) - (abs(v) < 1 ? 1 : 0);

  return n - m;
}


/*! Computes the value of (x)_n/n!.
 */
double Pochhammer_factorial(double x, int n)
{
  double result = 1.;

  for(int i = 0; i < n; i++)
    result *= (x + i)/(i + 1);

  return result;
}


/*! Error exit for the combination of the values x, a, and b.
 */
void error_exit(string error_message, double x, double a, double b)
{
  cerr << "Error in beta(" << (x != 1 ? "x," : "") << "a,b):" << endl
       << error_message << endl;
  if(x != 1) 
    cerr << "x = " << x << endl;
  cerr << "a = " << a << endl
       << "b = " << b << endl;

  exit(1);
}


/*! Verifies the analyticity of B_x(a,b) for the combination of the values x, a, and b.
 */
void verify_analyticity(double x, double a, double b)
{
  if(is_integer(a) && a <= 0 && !(is_integer(b) && b > 0 && a + b <= 0))
     error_exit("This function is divergent for nonpositive, integer a,\nexcept for the case a + b <= 0 for positive, integer b.",x,a,b);

  if(x == 0 && a < 0)
    error_exit("This function is divergent for x = 0 and a < 0.",x,a,b);

  if(x == 1 && is_integer(b) && b <= 0 && !(is_integer(a) && a > 0 && a + b <= 0))
    error_exit("This function is divergent for nonpositive, integer b,\nexcept for the case a + b <= 0 for positiver, integer a.",x,a,b);
}


/*! Calculates the value of B(a,b).
 */
double beta_(double a, double b)
{
  if(a == -b)
    return (is_integer(a) ? pow(-1,a)/abs(a) : 0);
  if(is_integer(a) && is_integer(b))
    if(a < 0 && b > 0 && a + b < 0)
      return 1/Pochhammer_factorial(a,b)/b;
    else if(a > 0 && b < 0 && a + b < 0)
      return 1/Pochhammer_factorial(b,a)/a;

  return gsl_sf_beta(a,b);
}


/*!
 * Calculates the value of B_x(a,b) for -1 < x < 1 through a series expansion.
 * If a is a non-positive integer, the term for n = -a is left out from the summation.
 * The summation is stopped when double precision is acquired or when the number
 * of terms in the summation exceeds n_max.
 */
complex<double> beta_series(double x, double a, double b, bool *converged_)
{
  int n = 0;
  complex<double> result = 0, result_previous = 1;
  bool converged = false;

  while(!converged && n < n_max)
    {
      result += (a == -n ? 0 : Pochhammer_factorial(1 - b,n)/(a + n)*cpow(x,n + a));
      
      converged = (result == result_previous && a != -n && n > 0);
      
      result_previous = result;
      n++;
    }

  *converged_ = converged;
  
  return result;
}


/*!
 * Implementation of the functions \mathcal{B}_1(x,a,n) and \mathcal{B}_2(x,a,n).
 * See documentation for more information.
 */
complex<double> limit(double x, double a, int n, int *mag_diff, short option)
{
  complex<double> result;
  int m = round(a);
  bool a_int = is_integer(a);
  double psi, psi_diff;

  if(option == 0)
    if(!(a_int && m - n <= 0))
      {
	psi = gsl_sf_psi(a - n); psi_diff = psi - gsl_sf_psi(n + 1);

	result = phase(n + 1)/(a*beta_(a - n,n + 1))*
	  (psi_diff + log(complex<double>(1 - x,0)));

	*mag_diff = magnitude_difference(psi,psi_diff);

      }
    else
      result = phase(m)*beta_(m,n - m + 1);

  else if(option == 1)
    if(!(a_int && m - n <= 0))
      {
	psi = gsl_sf_psi(a); psi_diff = gsl_sf_psi(n + 1) - psi;

	result = phase(-a)*phase(n)/(a*beta_(a - n,n + 1))*
	  (psi_diff + log(x) + complex<double>(0,pi));

	*mag_diff = magnitude_difference(psi,psi_diff);
      }
    else
      result = 0;
  
  return result;
}


/*!
 * Calculates B_x(a,b) for all real values of x, a, and b from a certain
 * identity of the incomplete beta function, which is specified by 'option'.
 */
complex<double> compute(double x, double a, double b, int *mag_diff, bool *converged, short option)
{
  complex<double> result, term;
  int n = round(b), m = round(a + b);
  bool b_int = (is_integer(b) && n <= 0), ab_int = (is_integer(a + b) && m > 0);

  *mag_diff = 0;

  if(option == 0)
    result = beta_series(x,a,b,converged);
  if(option == 1)
    {
      term = beta_series(1 - x,b_int ? n : b,a,converged);
      result = (b_int ? limit(x,a,-n,mag_diff,0) : beta_(a,b)) - term;
    }
  else if(option == 2)
    {
      term = phase(b)*beta_series(1/x,ab_int ? -m + 1 : 1 - a - b,b,converged);
      result = (ab_int ? limit(x,a,m - 1,mag_diff,1) : phase(-a)*beta_(1 - a - b,a)) + term;

      if(b == n && n > 0) result = complex<double>( result.real(), 0. );
    }
  else if(option == 3)
    {
      term = -phase(b)*beta_series(1 - 1/x,b_int ? n : b,1 - a - b,converged);
      result = (b_int ? conj(phase(a)*limit(x/(x - 1),a,-n,mag_diff,1)) : beta_(a,b)) + term;
    }
  else if(option == 4)
    result = phase(a)*beta_series(x/(x - 1),a,1 - a - b,converged);
  else if(option == 5)
    {
      term = phase(a)*beta_series(1/(1 - x),ab_int ? -m + 1 : 1 - a - b,a,converged);
      result = phase(a)*(ab_int ? limit(x/(x - 1),a,m - 1,mag_diff,0) :
			 beta_(1 - a - b,a)) - term;
    }

  if(!(*mag_diff > 1) && option != 0 && option != 4)
    *mag_diff = max(magnitude_difference(real(term),real(result)),
		    magnitude_difference(imag(term),imag(result)));
  
  return complex<double>(remove_noise(real(result)),remove_noise(imag(result)));
}


/*! Calculates B_x(a,b) for all real values of x, a, and b.
 */
complex<double> beta(double x, double a, double b)
{
  verify_analyticity(x,a,b);

  complex<double> result;
  bool converged;
  int mag_diff;

  if(x == 0 || (x == -b && a == x - 1))
    result = 0;
  else if(x == 1)
    result = beta_(a,b);
  else if(is_integer(a) && is_integer(b))
    result = compute(x,a,b,&mag_diff,&converged,(b > 0 ? 0 : 1));
  else 
    {
      result = compute
	(x,a,b,&mag_diff,&converged,((x > 0 && x <= x_limit) ? 0 :
				     ((x > x_limit && x < 1) ? 1 :
				      (x > 1 ? (1/x <= x_limit ? 2 : 3) :
				       (x/(x - 1) <= x_limit ? 4 : 5)))));
      
      if(mag_diff > 1 || (x > 1 && x < 2 && is_integer(b) && b <= 0 &&
			  (mag_diff <= 1 && mag_diff >= -1)))
      	{
      	  complex<double> result_ = compute(x,a,b,&mag_diff,&converged,
      					    (x > x_limit && x < 1 ? 0 :
      					     (x > 1 ? (1/x <= x_limit ? 3 : 2) : 4)));

	  result = (converged ? result_ : result);
      	}
    }
  
  return result;
}


/*! Returns the value of B(a,b). This function is for extern use only.
 */
double beta(double a, double b)
{
  return real(beta(1,a,b));
}
