/*!
 * \file FourVectorHandle.cpp
 * \ingroup wrapper
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

#include <complex>
#include <iostream>
#include "FourVectorHandle.h"

using namespace std;

//-------------------------------------------------------------------------

// Fourvector constructors
Fourvector::Fourvector( const FourVector< double >& input) 
  : isGamma(false)    
{
  complexx = FourVector< complex<double> >(input[0],
					   input[1],
					   input[2],
					   input[3]);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::Fourvector( const FourVector< complex<double> >& input) 
  : complexx(input), isGamma(false)
{
  // empty body
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::Fourvector( const FourVector< GammaStructure >& input)
  : gamma(input), isGamma(true)     
{
  // empty body
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::Fourvector( const Fourvector& toCopy )
  : complexx( toCopy.complexx ),
    gamma( toCopy.gamma ),
    isGamma( toCopy.isGamma )
{
  // empty body
}


//-------------------------------------------------------------------------

// Overload assignment operator
Fourvector& Fourvector::operator=(const Fourvector& toCopy)
{
  if(this != &toCopy) { // avoid self-assignment
    isGamma = toCopy.isGamma;
    complexx = toCopy.complexx;
    gamma = toCopy.gamma;
  }
  
  return *this;
}

//-------------------------------------------------------------------------

// Casting
// -------

Fourvector::operator const FourVector< complex<double> >&() const
{
  // When the Fourvector is of type FourVector<GammaStructure>
  // the cast is meaningless
  if( isGamma )
    {
      cerr << "Error in 'Cast from Fourvector to "
	   << "FourVector<complex>':" << endl
	   << "Attempt to cast Fourvector with GammaStructure components "
	   << "to FourVector<complex<double>>!" << endl;
      exit(1);
    }

  return complexx;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::operator const FourVector< GammaStructure >&() const
{
  // When the Fourvector is of type FourVector<complex<double>>
  // the cast is meaningless
  if( !isGamma )
    {
      cerr << "Error in 'Cast from Fourvector to "
	   << "FourVector<GammaStructure>':" << endl
	   << "Attempt to cast Fourvector with complex<double> components "
	   << "to FourVector<GammaStructure>!" << endl;
      exit(1);
    }

  return gamma;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::operator FourVector< complex<double> >&()
{
  // When the Fourvector is of type FourVector<GammaStructure>
  // the cast is meaningless
  if( isGamma )
    {
      cerr << "Error in 'Cast from Fourvector to "
	   << "FourVector<complex>':" << endl
	   << "Attempt to cast Fourvector with GammaStructure components "
	   << "to FourVector<complex<double>>!" << endl;
      exit(1);
    }

  return complexx;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector::operator FourVector< GammaStructure >&()
{
  // When the Fourvector is of type FourVector<complex<double>>
  // the cast is meaningless
  if( !isGamma )
    {
      cerr << "Error in 'Cast from Fourvector to "
	   << "FourVector<GammaStructure>':" << endl
	   << "Attempt to cast Fourvector with complex<double> components "
	   << "to FourVector<GammaStructure>!" << endl;
      exit(1);
    }

  return gamma;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Overloaded << operators used as cast operators
FourVector< complex<double> >&
operator<<(FourVector< complex<double> >& left, const Fourvector& right)
{
  // When the Fourvector is of type FourVector<GammaStructure>
  // the cast is meaningless
  if( right.isGamma )
    {
      cerr << "Error in 'FourVector<complex<double>>& "
	   << "operator*(const FourVector<complex<double>>&, "
	   << "const Fourvector&)':" << endl
	   << "Attempt to cast Fourvector with GammaStructure components "
	   << "to FourVector<complex<double>>!" << endl;
      exit(1);
    }

  // Do assignment
  left[0] = right.complexx[0];
  left[1] = right.complexx[1];
  left[2] = right.complexx[2];
  left[3] = right.complexx[3];
 
  return left;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Overloaded << operators used as cast operators
FourVector< GammaStructure >&
operator<<(FourVector< GammaStructure >& left, const Fourvector& right)
{
  // When the Fourvector is of type FourVector<complex<double>>
  // the cast is meaningless
  if( !right.isGamma )
    {
      cerr << "Error in 'FourVector<GammaStructure>& "
	   << "operator*(const FourVector<GammaStructure>&, "
	   << "const Fourvector&)':" << endl
	   << "Attempt to cast Fourvector with complex<double> components "
	   << "to FourVector<GammaStructure>!" << endl;
      exit(1);
    }

  // Do assignment
  left[0] = right.gamma[0];
  left[1] = right.gamma[1];
  left[2] = right.gamma[2];
  left[3] = right.gamma[3];
 
  return left;
}

//-------------------------------------------------------------------------

/// Hermitian conjugate
Fourvector& Fourvector::HermitianConjugate()
{
  if(isGamma)
    for(int i=0; i<4; i++)
      gamma[i] = gamma[i].Hermitian();

  else
    for(int i=0; i<4; i++)
      complexx[i] = real(complexx[i]) - complex<double>(0.0,1.0)*imag(complexx[i]);

  return *this;
}

//-------------------------------------------------------------------------

/// Scalar multiplication
Fourvector& Fourvector::operator*=(double factor)
{
  if(isGamma)
    for(int i=0; i<4; ++i)
      gamma[i] *= factor;

  else
    for(int i=0; i<4; ++i)
      complexx[i] *= factor;
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector& Fourvector::operator*=(const complex<double>& factor)
{
  if(isGamma)
    for(int i=0; i<4; ++i)
      gamma[i] *= factor;

  else
    for(int i=0; i<4; ++i)
      complexx[i] *= factor;
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector& Fourvector::operator*=(const GammaStructure& factor)
{
  if(isGamma)
    for(int i=0; i<4; ++i)
      gamma[i] *= factor;

  else
    {
      isGamma = true;
      for(int i=0; i<4; ++i)
	gamma[i] = complexx[i] * factor;
    }  

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector Fourvector::operator*(double factor) const
{
  Fourvector product = *this;
  product *= factor;

  return product;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector Fourvector::operator*(const complex<double>& factor) const
{
  Fourvector product = *this;
  product *= factor;

  return product;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector Fourvector::operator*(const GammaStructure& factor) const
{
  Fourvector product = *this;
  product *= factor;

  return product;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector operator*(double factor, const Fourvector& fourvector)
{
  Fourvector product = fourvector;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector operator*(const complex<double>& factor,const Fourvector& fourvector)
{
  Fourvector product = fourvector;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// Scalar multiplication
Fourvector operator*(const GammaStructure& factor,const Fourvector& fourvector)
{
  Fourvector product = fourvector;

  if(product.isGamma)
    for(int i=0; i<4; ++i)
      product.gamma[i] = factor * product.gamma[i];

  else
    {
      product.isGamma = true;
      for(int i=0; i<4; ++i)
	product.gamma[i] = product.complexx[i] * factor;
    }  

  return product;
}

//-------------------------------------------------------------------------

/// Vector Multiplication
GammaStructure Fourvector::operator*(const Fourvector& right) const
{
  if( !this->isGamma && !right.isGamma )
    {
      cerr << "Error in GammaStructure Fourvector::operator*(Fourvector): "
	   << "Multiplying two FourVector<complex>'s can never yield "
	   << "a GammaStructure!" << endl;
      exit(1);
    }

  //( _left->isGamma ? _left->gamma[i] : _left->complexx[i] )
  
  GammaStructure product = 
    ( this->isGamma ? this->gamma[0] : this->complexx[0] ) *
    ( right.isGamma ? right.gamma[0] : right.complexx[0] );

  for(int i=1; i<4; ++i)
    product -= 
      ( this->isGamma ? this->gamma[i] : this->complexx[i] ) *
      ( right.isGamma ? right.gamma[i] : right.complexx[i] );

  return product;
}

//-------------------------------------------------------------------------

// Output
// ------

ostream& operator<<(ostream& stream, const Fourvector& fourvector)
{
  if(fourvector.isGamma)
    stream << fourvector.gamma;
  else
    stream << fourvector.complexx;
  
  return stream;
}
