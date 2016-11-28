/*!
 * \file TensorRank4.cpp
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

#include "TensorRank2.h"
#include "TensorRank3.h"
#include "TensorRank4.h"
using std::complex;
using std::cout; using std::cerr; using std::endl;
using std::ostream;

typedef std::vector<int>::size_type size_type;

/*!
 * \class TensorRank4
 *
 * Represents direct products of FourVector's:
 * rank4 = a^{alpha} b^{beta} c^{gamma} d^{delta}
 *
 * Common operations like +,-,* are provided.
 *
 * When the result of a multiplication is a FourVector, the overloaded
 * operator*() will return a Fourvector (= FourVector<> handle).
 * It is up to the user to cast it to the appropriate type when necessary.
 *
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>

 */

// Constructors
// ------------

TensorRank4::TensorRank4()
{
  array_.reserve(2);
  array_.push_back( TensorRank2(FourVector<double>(),
				FourVector<double>()) );
  array_.push_back( TensorRank2(FourVector<double>(),
				FourVector<double>()) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4::TensorRank4(const TensorRank2& left,
			 const TensorRank2& right)
{
  array_.reserve(2);
  array_.push_back( left );
  array_.push_back( right ); 
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4::TensorRank4(const Fourvector& v1, const Fourvector& v2,
			 const Fourvector& v3, const Fourvector& v4)
{
  array_.reserve(2);
  array_.push_back( TensorRank2(v1,v2) );
  array_.push_back( TensorRank2(v3,v4) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4::TensorRank4(const TensorRank2& left,
			 const Fourvector& v1, const Fourvector& v2)
{
  array_.reserve(2);
  array_.push_back( left );
  array_.push_back( TensorRank2(v1,v2) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4::TensorRank4(const Fourvector& v1, const Fourvector& v2,
			 const TensorRank2& right)
{
  array_.reserve(2);
  array_.push_back( TensorRank2(v1,v2) );
  array_.push_back( right );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4::TensorRank4(const TensorRank4& toCopy)
  : array_(toCopy.array_)
{
  // empty body
}

//----------------------------------------------------------------------

// Destructor
// ----------

TensorRank4::~TensorRank4()
{
  // empty body
}

//----------------------------------------------------------------------

// Assignment
// ----------

TensorRank4& TensorRank4::operator=(const TensorRank4& toCopy)
{
  array_ = toCopy.array_;

  return *this;
}

//----------------------------------------------------------------------

// Addition
// --------
// Addition is done by concatenating both array_'s

TensorRank4& TensorRank4::operator+=(const TensorRank4& right)
{
  array_.reserve( array_.size() + right.array_.size() );

  for(size_type i=0; i<right.array_.size(); ++i)
    array_.push_back( right.array_[i] );
 
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4& TensorRank4::operator-=(const TensorRank4& right)
{
  array_.reserve( array_.size() + right.array_.size() );

  for(size_type i=0; i<right.array_.size(); ++i)
    array_.push_back( ( i%2==0 ? -1.0 : 1.0 ) * right.array_[i] );
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator+(const TensorRank4& right) const
{
  TensorRank4 sum = *this;
  sum += right;

  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator-(const TensorRank4& right) const
{
  TensorRank4 sum = *this;
  sum -= right;

  return sum;
}

//----------------------------------------------------------------------

// Scalar Multiplication
// ---------------------

TensorRank4& TensorRank4::operator*=(double factor)
{
  for(size_type i=0; i<array_.size(); i+=2)
    array_[i] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4& TensorRank4::operator*=(const complex<double>& factor)
{
  for(size_type i=0; i<array_.size(); i+=2)
    array_[i] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4& TensorRank4::operator*=(const GammaStructure& factor)
{
  for(size_type i=1; i<array_.size(); i+=2)
    array_[i] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator*(double factor) const
{
  TensorRank4 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator*(const complex<double>& factor) const
{
  TensorRank4 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator*(const GammaStructure& factor) const
{
  TensorRank4 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 operator*(double factor, const TensorRank4& tensor)
{
  TensorRank4 product = tensor;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 operator*(const complex<double>& factor, const TensorRank4& tensor)
{
  TensorRank4 product = tensor;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 operator*(const GammaStructure& factor, const TensorRank4& tensor)
{
  TensorRank4 product = tensor;
  for(size_type i=0; i<product.array_.size(); i+=2)
    product.array_[i] = factor * product.array_[i];

  return product;
}

//----------------------------------------------------------------------

// Vector Multiplication
// ---------------------

TensorRank3 TensorRank4::operator*(const Fourvector& vector) const
{
  TensorRank3 product(array_[0],array_[1]*vector);

  for(size_type i=2; i<array_.size(); i+=2)
    product += TensorRank3(array_[i],array_[i+1]*vector);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank4::operator*(const FourVector<double>& vector) const
{
  TensorRank3 product(array_[0],array_[1]*vector);

  for(size_type i=2; i<array_.size(); i+=2)
    product += TensorRank3(array_[i],array_[i+1]*vector);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank4::operator*(const FourVector<complex<double> >& vector) const
{
  TensorRank3 product(array_[0],array_[1]*vector);

  for(size_type i=2; i<array_.size(); i+=2)
    product += TensorRank3(array_[i],array_[i+1]*vector);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank4::operator*(const FourVector<GammaStructure>& vector) const
{
  TensorRank3 product(array_[0],array_[1]*vector);

  for(size_type i=2; i<array_.size(); i+=2)
    product += TensorRank3(array_[i],array_[i+1]*vector);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const Fourvector& vector,
		      const TensorRank4& tensor)
{
  TensorRank3 product(vector*tensor.array_[0],tensor.array_[1]);

  for(size_type i=2; i<tensor.array_.size(); i+=2)
    product += TensorRank3(vector*tensor.array_[i],tensor.array_[i+1]);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const FourVector<double>& vector,
		      const TensorRank4& tensor)
{
  TensorRank3 product(vector*tensor.array_[0],tensor.array_[1]);

  for(size_type i=2; i<tensor.array_.size(); i+=2)
    product += TensorRank3(vector*tensor.array_[i],tensor.array_[i+1]);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const FourVector< complex<double> >& vector,
		      const TensorRank4& tensor)
{
  TensorRank3 product(vector*tensor.array_[0],tensor.array_[1]);

  for(size_type i=2; i<tensor.array_.size(); i+=2)
    product += TensorRank3(vector*tensor.array_[i],tensor.array_[i+1]);

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const FourVector<GammaStructure>& vector,
		      const TensorRank4& tensor)
{
  TensorRank3 product(vector*tensor.array_[0],tensor.array_[1]);

  for(size_type i=2; i<tensor.array_.size(); i+=2)
    product += TensorRank3(vector*tensor.array_[i],tensor.array_[i+1]);

  return product;
}

//----------------------------------------------------------------------

// TensorRank2 Multiplication
// --------------------------

TensorRank4& TensorRank4::operator*=(const TensorRank2& rank2)
{
  for(size_type i=1; i<array_.size(); i+=2)
    array_[i] *= rank2;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 TensorRank4::operator*(const TensorRank2& rank2) const
{
  TensorRank4 product = *this;
  product *= rank2;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank4 operator*(const TensorRank2& rank2, const TensorRank4& rank4)
{
  TensorRank4 product = rank4;
  for(size_type i=0; i<product.array_.size(); i+=2)
    product.array_[i] = rank2 * product.array_[i];

  return product;
}

//----------------------------------------------------------------------

// TensorRank2 Contraction
// -----------------------

TensorRank2 TensorRank4::operator%(const TensorRank2& rank2) const
{
  TensorRank2 contraction(FourVector<double>(0.0,0.0,0.0,0.0),FourVector<double>(0.0,0.0,0.0,0.0));
  TensorRank2 term;
  double delta;

  rank2.toMatrix();

  for(size_type i=0; i<array_.size(); i+=2) {
    for(int mu=0; mu<4; ++mu) {
      for(int nu=0; nu<4; ++nu) {
	if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	  delta = -1.0;
	else
	  delta = 1.0;

	array_[i+1].toMatrix();

	term = delta * array_[i] * (array_[i+1]._isGamma ? array_[i+1]._gmatrix[4*mu+nu] : 
	       array_[i+1]._cmatrix[4*mu+nu]) * (rank2._isGamma ? rank2._gmatrix[4*mu+nu] : 
	       rank2._cmatrix[4*mu+nu]);

	if(term._isGamma && !contraction._isGamma)
	  contraction *= GammaStructure(1.0);
	contraction += term;
      }
    }
  }
  
  return contraction;
}

//----------------------------------------------------------------------

// Output
// ------

ostream& operator<<(ostream& stream, const TensorRank4& tensor)
{
  stream << "{{ " << endl;
  for(size_type i=0; i<tensor.array_.size(); ++i) {
    stream << tensor.array_.at(i) << " X ";
    stream << tensor.array_.at(++i)
	   << endl << "+";
  }
  stream << "\r" << "}}";

  return stream;
}
