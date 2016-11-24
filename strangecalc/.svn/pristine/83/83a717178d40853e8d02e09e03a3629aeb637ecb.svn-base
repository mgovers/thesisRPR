/*!
 * \file TensorRank2.cpp
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
#include <vector>
using std::complex;
using std::cout; using std::cerr; using std::endl;
using std::ostream;

typedef std::vector<int>::size_type size_type;

/*!
 * \class TensorRank2
 *
 * Represents direct products of FourVector's:
 * rank2 = a^{alpha} b^{beta}
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*! The default constructor constructs the metric */
TensorRank2::TensorRank2()
  : _isFourVector(false), _isGamma(false),   _left(0),    _right(0),    _gmatrix(0)
{
  // allocate memory from heap
  _cmatrix = new complex<double>[4*4];

  /* We declare 3 static complex doubles
   * to improve performance */
  static complex<double> zero(0.0);
  static complex<double> one(1.0);
  static complex<double> minusOne(-1.0);

  _cmatrix[0*4+0] = one;
  _cmatrix[0*4+1] = zero;
  _cmatrix[0*4+2] = zero;
  _cmatrix[0*4+3] = zero;
  _cmatrix[1*4+0] = zero;
  _cmatrix[1*4+1] = minusOne;
  _cmatrix[1*4+2] = zero;
  _cmatrix[1*4+3] = zero;
  _cmatrix[2*4+0] = zero;
  _cmatrix[2*4+1] = zero;
  _cmatrix[2*4+2] = minusOne;
  _cmatrix[2*4+3] = zero;
  _cmatrix[3*4+0] = zero;
  _cmatrix[3*4+1] = zero;
  _cmatrix[3*4+2] = zero;
  _cmatrix[3*4+3] = minusOne;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const Fourvector& left,
			 const Fourvector& right)
  : _isFourVector(true),
    _cmatrix(0),
    _gmatrix(0)
{
  if(left.isGamma || right.isGamma)
    _isGamma = true;
  else
    _isGamma = false;

  _left = new Fourvector(left);
  _right = new Fourvector(right);  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const Fourvector& left,
			 const FourVector<double>& right)
  : _isFourVector(true),
    _cmatrix(0),
    _gmatrix(0)
{
  if(left.isGamma)
    _isGamma = true;
  else
    _isGamma = false;

  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<double>& left,
			 const Fourvector& right)
  : _isFourVector(true),
    _cmatrix(0),
    _gmatrix(0)
{
  if(right.isGamma)
    _isGamma = true;
  else
    _isGamma = false;

  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const Fourvector& left,
			 const FourVector< complex<double> >& right)
  : _isFourVector(true),
    _cmatrix(0),
    _gmatrix(0)
{
  if(left.isGamma)
    _isGamma = true;

  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector< complex<double> >& left,
			 const Fourvector& right)
  : _isFourVector(true),
    _cmatrix(0),
    _gmatrix(0)
{
  if(right.isGamma)
    _isGamma = true;
  else
    _isGamma = false;

  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const Fourvector& left,
			 const FourVector< GammaStructure >& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector< GammaStructure >& left,
			 const Fourvector& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector< complex<double> >& left,
			 const FourVector< complex<double> >& right)
  : _isFourVector(true),
    _isGamma(false),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<double>& left,
			 const FourVector<double>& right)
  : _isFourVector(true),
    _isGamma(false),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<GammaStructure>& left,
			 const FourVector<GammaStructure>& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<double>& left,
			 const FourVector< complex<double> >& right)
  : _isFourVector(true),
    _isGamma(false),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector< complex<double> >& left,
			 const FourVector<double>& right)
  : _isFourVector(true),
    _isGamma(false),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector< complex<double> >& left,
			 const FourVector<GammaStructure>& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<GammaStructure>& left,
			 const FourVector< complex<double> >& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<double>& left,
			 const FourVector<GammaStructure>& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const FourVector<GammaStructure>& left,
			 const FourVector<double>& right)
  : _isFourVector(true),
    _isGamma(true),
    _cmatrix(0),
    _gmatrix(0)
{
  _left = new Fourvector(left);
  _right = new Fourvector(right);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const complex<double>& c1,
			 const complex<double>& c2, 
			 const complex<double>& c3, 
			 const complex<double>& c4,
			 const complex<double>& c5, 
			 const complex<double>& c6, 
			 const complex<double>& c7, 
			 const complex<double>& c8, 
			 const complex<double>& c9, 
			 const complex<double>& c10, 
			 const complex<double>& c11, 
			 const complex<double>& c12, 
			 const complex<double>& c13, 
			 const complex<double>& c14, 
			 const complex<double>& c15, 
			 const complex<double>& c16 )
  : _isFourVector(false),
    _isGamma(false),
    _left(0),
    _right(0),
    _gmatrix(0)
{
  // allocate memory from heap
  _cmatrix = new complex<double>[4*4];
 
  // initialize
  _cmatrix[0*4+0] = c1;
  _cmatrix[0*4+1] = c2;
  _cmatrix[0*4+2] = c3;
  _cmatrix[0*4+3] = c4;
  _cmatrix[1*4+0] = c5;
  _cmatrix[1*4+1] = c6;
  _cmatrix[1*4+2] = c7;
  _cmatrix[1*4+3] = c8;
  _cmatrix[2*4+0] = c9;
  _cmatrix[2*4+1] = c10;
  _cmatrix[2*4+2] = c11;
  _cmatrix[2*4+3] = c12;
  _cmatrix[3*4+0] = c13;
  _cmatrix[3*4+1] = c14;
  _cmatrix[3*4+2] = c15;
  _cmatrix[3*4+3] = c16;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const GammaStructure& c1,
			 const GammaStructure& c2, 
			 const GammaStructure& c3, 
			 const GammaStructure& c4,
			 const GammaStructure& c5, 
			 const GammaStructure& c6, 
			 const GammaStructure& c7, 
			 const GammaStructure& c8, 
			 const GammaStructure& c9, 
			 const GammaStructure& c10, 
			 const GammaStructure& c11, 
			 const GammaStructure& c12, 
			 const GammaStructure& c13, 
			 const GammaStructure& c14, 
			 const GammaStructure& c15, 
			 const GammaStructure& c16 )
  : _isFourVector(false),
    _isGamma(true),
    _left(0),
    _right(0),
    _cmatrix(0)
{
  // allocate memory from heap
  _gmatrix = new GammaStructure[4*4];
 
  // initialize
  _gmatrix[0*4+0] = c1;
  _gmatrix[0*4+1] = c2;
  _gmatrix[0*4+2] = c3;
  _gmatrix[0*4+3] = c4;
  _gmatrix[1*4+0] = c5;
  _gmatrix[1*4+1] = c6;
  _gmatrix[1*4+2] = c7;
  _gmatrix[1*4+3] = c8;
  _gmatrix[2*4+0] = c9;
  _gmatrix[2*4+1] = c10;
  _gmatrix[2*4+2] = c11;
  _gmatrix[2*4+3] = c12;
  _gmatrix[3*4+0] = c13;
  _gmatrix[3*4+1] = c14;
  _gmatrix[3*4+2] = c15;
  _gmatrix[3*4+3] = c16;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2::TensorRank2(const TensorRank2& toCopy)
  : _isFourVector(toCopy._isFourVector),
    _isGamma(toCopy._isGamma)
{
  if(_isFourVector) 
    {
      _left = new Fourvector(*toCopy._left);
      _right = new Fourvector(*toCopy._right);
      
      _cmatrix = 0;
      _gmatrix = 0;
    }
  else if(_isGamma) 
    {
      _left = 0;
      _right = 0;
      _cmatrix = 0;

      _gmatrix = new GammaStructure[4*4];
      
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _gmatrix[i*4+j] = toCopy._gmatrix[i*4+j];
    }
  else
    {
      _left = 0;
      _right = 0;
      _gmatrix = 0;

      _cmatrix = new complex<double>[4*4];
      
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _cmatrix[i*4+j] = toCopy._cmatrix[i*4+j];
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// private constructor
TensorRank2::TensorRank2(bool isFourVector, bool isGamma,
			 Fourvector* left, Fourvector* right,
			 complex<double>* cmatrix, GammaStructure* gmatrix)
  : _isFourVector( isFourVector ),
    _isGamma( isGamma ),
    _left( left ),
    _right( right ),
    _cmatrix( cmatrix ),
    _gmatrix( gmatrix )
{
  // empty body
}

//----------------------------------------------------------------------

// Destructor
// ----------

TensorRank2::~TensorRank2()
{
  delete _left;
  delete _right;
  delete[] _cmatrix;
  delete[] _gmatrix;
}

//----------------------------------------------------------------------

// Overloaded assignment operator
// ------------------------------

TensorRank2& TensorRank2::operator=(const TensorRank2& toCopy)
{
  if(this != &toCopy) // avoid self-assignment
    {
      if(toCopy._isFourVector) // Right-side stored as FourVector's
	{
	  if(this->_isFourVector) { // left-side too
	    *_left = *toCopy._left;
	    *_right = *toCopy._right;	    
	  }
	  
	  else {                    // left side stored as matrix
	    delete _left;
	    delete _right;
	    _left = new Fourvector(*toCopy._left);
	    _right = new Fourvector(*toCopy._right);

	    if(this->_isGamma) {
	      delete[] _gmatrix;
	      _gmatrix = 0;
	    }
	    else {
	      delete[] _cmatrix;
	      _cmatrix = 0;
	    }
	  }
	}
      else if(toCopy._isGamma) // Right-side stored as Gamma matrix
	{
	  if(this->_isFourVector) { // left-side stored as FourVector's
	    delete _left;
	    delete _right;
	    _left = 0;
	    _right = 0;
	     
	    _gmatrix = new GammaStructure[4*4];
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_gmatrix[i*4+j] = toCopy._gmatrix[i*4+j];
	  }

	  else if(this->_isGamma) { // left-side stored as gamma matrix
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_gmatrix[i*4+j] = toCopy._gmatrix[i*4+j];
	  }
	  else {                    // left-side stored as complex matrix
	    delete[] _cmatrix;
	    _cmatrix = 0;

	    _gmatrix = new GammaStructure[4*4];
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_gmatrix[i*4+j] = toCopy._gmatrix[i*4+j];
	  }
	}
      else                      // Right-side stored as complex matrix
	{
	  if(this->_isFourVector) { // left-side stored as FourVector's
	    delete _left;
	    delete _right;
	    _left = 0;
	    _right = 0;
	     
	    _cmatrix = new complex<double>[4*4];
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_cmatrix[i*4+j] = toCopy._cmatrix[i*4+j];
	  }

	  else if(this->_isGamma) { // left-side stored as gamma matrix
	    delete[] _gmatrix;
	    _gmatrix = 0;

	    _cmatrix = new complex<double>[4*4];
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_cmatrix[i*4+j] = toCopy._cmatrix[i*4+j];
	  }
	  else {                    // left-side stored as complex matrix
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
		_cmatrix[i*4+j] = toCopy._cmatrix[i*4+j];
	  }
	}

      // Finally we copy the bool's
      _isFourVector = toCopy._isFourVector;
      _isGamma = toCopy._isGamma;
    }
  
  return *this;
}

//----------------------------------------------------------------------

/*! \brief Internal conversion
 *
 * A TensorRank2 object is stored internally as two Fourvector's
 * or as a 4x4 matrix.
 * This member function makes sure the object is stored as a
 * matrix.
 */
void TensorRank2::toMatrix() const
{
  TensorRank2 *self = const_cast<TensorRank2*>(this);

  if(self->_isFourVector) // only do something when stored as Fourvector's
    {
      if(self->_isGamma)
	{
	  self->_gmatrix = new GammaStructure[4*4];

	  if( self->_left->isGamma ) {
	    if( self->_right->isGamma )
	      for(int i=0; i<4; i++)
	  	for(int j=0; j<4; j++)
	  	  self->_gmatrix[i*4+j] = self->_left->gamma[i] * self->_right->gamma[j];
	    else
	      for(int i=0; i<4; i++)
	  	for(int j=0; j<4; j++)
	  	  self->_gmatrix[i*4+j] = self->_left->gamma[i] * self->_right->complexx[j];
	    
	  } else { // !_left->isGamma
	    for(int i=0; i<4; i++)
	      for(int j=0; j<4; j++)
	  	self->_gmatrix[i*4+j] = self->_left->complexx[i] * self->_right->gamma[j];
	  } 
	}
      else
	{
	  self->_cmatrix = new complex<double>[4*4];
	  for(int i=0; i<4; i++)
	    for(int j=0; j<4; j++)
	      self->_cmatrix[i*4+j] = self->_left->complexx[i] * self->_right->complexx[j];
	}
      
      self->_isFourVector = false;
    }
}

//----------------------------------------------------------------------

// Addition
// --------

TensorRank2& TensorRank2::operator+=(const TensorRank2& that)
{
  this->toMatrix();
  if(!that._isFourVector)
    {  
      // when both are of type GammaStructure
      if(this->_isGamma && that._isGamma)
	for(int i=0; i<4; ++i)
	  for(int j=0; j<4; ++j)
	    this->_gmatrix[i*4+j] += that._gmatrix[i*4+j];
      
      // when both are complex
      else if(!this->_isGamma && !that._isGamma)
	for(int i=0; i<4; ++i)
	  for(int j=0; j<4; ++j)
	    this->_cmatrix[i*4+j] += that._cmatrix[i*4+j];
	  
      else {
	cerr << "Error in TensorRank2& TensorRank2::operator+=(TensorRank2): "
	     << "Rank-2 tensors are not of same internal type!" << endl;
	exit(1);
      }
    }
  else
    {
      if(this->_isGamma && that._isGamma){
	for(int i=0; i<4; ++i){
	  for(int j=0; j<4; ++j){
	    if(that._left->isGamma){
	      if(that._right->isGamma)
		this->_gmatrix[i*4+j].add(that._left->gamma[i],that._right->gamma[j]);
	      else
		this->_gmatrix[i*4+j].add(that._right->complexx[j], that._left->gamma[i]);
	    }
	    else
	      this->_gmatrix[i*4+j].add(that._left->complexx[i], that._right->gamma[j]);
	  }  
	}
      }
      else if(!this->_isGamma && !that._isGamma)
	for(int i=0; i<4; ++i)
 	  for(int j=0; j<4; ++j)
	    this->_cmatrix[i*4+j] += that._left->complexx[i] * that._right->complexx[j];
      else{
	cerr << "Error in TensorRank2& TensorRank2::operator+=(TensorRank2): "
 	     << "Rank-2 tensors are not of same internal type!" << endl;
 	exit(1);
      }
    }
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2& TensorRank2::operator-=(const TensorRank2& that)
{
  this->toMatrix();
  if(!that._isFourVector)
    {  
      // when both are of type GammaStructure
      if(this->_isGamma && that._isGamma)
	for(int i=0; i<4; ++i)
	  for(int j=0; j<4; ++j)
	    this->_gmatrix[i*4+j] -= that._gmatrix[i*4+j];
      
      // when both are complex
      else if(!this->_isGamma && !that._isGamma)
	for(int i=0; i<4; ++i)
	  for(int j=0; j<4; ++j)
	    this->_cmatrix[i*4+j] -= that._cmatrix[i*4+j];
	  
      else {
	cerr << "Error in TensorRank2& TensorRank2::operator-=(TensorRank2): "
	     << "Rank-2 tensors are not of same internal type!" << endl;
	exit(1);
      }
    }
  else
    {
      if(this->_isGamma && that._isGamma){
	for(int i=0; i<4; ++i){
	  for(int j=0; j<4; ++j){
	    if(that._left->isGamma){
	      if(that._right->isGamma)
		this->_gmatrix[i*4+j].add(-1.0 * that._left->gamma[i],that._right->gamma[j]);
	      else
		this->_gmatrix[i*4+j].add(-1.0 * that._right->complexx[j], that._left->gamma[i]);
	    }
	    else
	      this->_gmatrix[i*4+j].add(-1.0 * that._left->complexx[i], that._right->gamma[j]);
	  }  
	}
      }
      else if(!this->_isGamma && !that._isGamma)
	for(int i=0; i<4; ++i)
 	  for(int j=0; j<4; ++j)
	    this->_cmatrix[i*4+j] -= that._left->complexx[i] * that._right->complexx[j];
      else{
	cerr << "Error in TensorRank2& TensorRank2::operator-=(TensorRank2): "
 	     << "Rank-2 tensors are not of same internal type!" << endl;
 	exit(1);
      }
    }
  return *this;
}

// TensorRank2& TensorRank2::operator-=(const TensorRank2& that)
// {
//   // Change the way both terms are stored internally
//   this->toMatrix();
//   TensorRank2 right = that;
//   right.toMatrix();
// 
//   // when both are of type GammaStructure
//   if(this->_isGamma && right._isGamma)
//     for(int i=0; i<4; ++i)
//       for(int j=0; j<4; ++j)
// 	this->_gmatrix[i*4+j] -= right._gmatrix[i*4+j];
// 
//   // when both are complex
//   else if(!this->_isGamma && !right._isGamma)
//     for(int i=0; i<4; ++i)
//       for(int j=0; j<4; ++j)
// 	this->_cmatrix[i*4+j] -= right._cmatrix[i*4+j];
//   
//   else {
//     cerr << "Error in TensorRank2& TensorRank2::operator-=(TensorRank2): "
// 	 << "Rank-2 tensors are not of same internal type!" << endl;
//     exit(1);
//   }
// 
//   return *this;
// }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank2::operator+(const TensorRank2& right) const
{
  TensorRank2 sum = *this;
  sum += right;

  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank2::operator-(const TensorRank2& right) const
{
  TensorRank2 sum = *this;
  sum -= right;

  return sum;
}

//----------------------------------------------------------------------

// Scalar Multiplication
// ---------------------

TensorRank2& TensorRank2::operator*=(double scalar)
{
  if(_isFourVector)
    {
      *_left *= scalar;
    }
  else if(_isGamma)
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _gmatrix[i*4+j] *= scalar;
    }
  else
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _cmatrix[i*4+j] *= scalar;
    }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2& TensorRank2::operator*=(const complex<double>& scalar)
{
  if(_isFourVector)
    {
      *_left *= scalar;
    }
  else if(_isGamma)
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _gmatrix[i*4+j] *= scalar;
    }
  else
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _cmatrix[i*4+j] *= scalar;
    }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2& TensorRank2::operator*=(const GammaStructure& scalar)
{
  if(_isFourVector)
    {
      _isGamma = true;
      *_right *= scalar;
    }
  else if(_isGamma)
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _gmatrix[i*4+j] *= scalar;
    }
  else
    {
      _isGamma = true;

      _gmatrix = new GammaStructure[4*4];
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  _gmatrix[i*4+j] = _cmatrix[i*4+j] * scalar;

      delete[] _cmatrix;
      _cmatrix = 0;
    }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank2::operator*(double factor) const
{
  TensorRank2 tensor = *this;
  tensor *= factor;

  return tensor;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank2::operator*(const complex<double>& factor) const
{
  TensorRank2 tensor = *this;
  tensor *= factor;

  return tensor;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank2::operator*(const GammaStructure& factor) const
{
  TensorRank2 tensor = *this;
  tensor *= factor;

  return tensor;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(double factor, const TensorRank2& tensor)
{
  return ( tensor * factor );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const complex<double>& factor, const TensorRank2& tensor)
{
  return ( tensor * factor );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const GammaStructure& factor, const TensorRank2& tensor)
{
  TensorRank2 product = tensor;

  if(product._isFourVector)
    {
      product._isGamma = true;
      *(product._left) = factor * *(product._left);
    }
  else if(product._isGamma)
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  product._gmatrix[i*4+j] = factor * product._gmatrix[i*4+j];
    }
  else
    {
      product._isGamma = true;

      product._gmatrix = new GammaStructure[4*4];
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  product._gmatrix[i*4+j] = factor * product._cmatrix[i*4+j];

      delete[] product._cmatrix;
      product._cmatrix = 0;
    }

  return product;
}

//------------------------------------------------------------------------

// Vector multiplication
// ---------------------

// Very very ugly implementation, but it works...
Fourvector TensorRank2::operator*(const Fourvector& vector) const
{
  if(_isFourVector) // result = _left * ( _right*vector )
    {
      Fourvector result = *_left;

      if(_right->isGamma || vector.isGamma)
	result *= ( *_right * vector );
      
      else  
	result *= (_right->complexx * vector.complexx );
                  
      return result;
    }

  else if(_isGamma || vector.isGamma) { // result = _gmatrix * vector
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += 
	  ( _isGamma ? _gmatrix[i*4+j] : _cmatrix[i*4+j] ) *
	  ( vector.isGamma ? 
	    ((vector.gamma).lowerIndex())[j] : 
	    ((vector.complexx).lowerIndex())[j] );
    
    return result;
  }
  else { // result = _cmatrix * vector
    FourVector< complex<double> > result;

    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += ( _cmatrix[i*4+j] * ((vector.complexx).lowerIndex())[j] );

    return result;
  }
  
  // dummy return
  return vector;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector TensorRank2::operator*(const FourVector<double>& vector) const
{
  FourVector< complex<double> > cvector(vector[0],vector[1],
					vector[2],vector[3]);

  return ( *this * cvector );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector TensorRank2::operator*(const FourVector< complex<double> >& vector) const
{
  if(_isFourVector) // result = _left * ( _right*vector )
    {
      Fourvector result = *_left;

      if(_right->isGamma)
	result *= ( _right->gamma * vector );
      
      else  
	result *= (_right->complexx * vector );
                  
      return result;
    }

  else if(_isGamma) { // result = _gmatrix * vector

    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += _gmatrix[i*4+j] * (vector.lowerIndex())[j];

    return result;
  }
  else { // result = result = _cmatrix * vector
    FourVector< complex<double> > result;

    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += _cmatrix[i*4+j] * (vector.lowerIndex())[j];

    return result;
  }
  
  // dummy return
  return vector;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector TensorRank2::operator*(const FourVector<GammaStructure>& vector) const
{
  if(_isFourVector) // result = _left * ( _right*vector )
    {
      Fourvector result = *_left;

      if(_right->isGamma)
	result *= ( _right->gamma * vector );

      else
	result *= ( _right->complexx * vector );
           
      return result;
    }

  else if(_isGamma) { // result = _gmatrix * vector
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += _gmatrix[i*4+j] * (vector.lowerIndex())[j];
    
    return result;
  }
  else { // result = _cmatrix * vector
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += _cmatrix[i*4+j] * (vector.lowerIndex())[j];
    
    return result;
  }
  
  // dummy return
  return vector;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Very very ugly implementation, but it works...
Fourvector operator*(const Fourvector& vector,
		     const TensorRank2& tensor)
{
  if(tensor._isFourVector) // result =(vector* _left) * _right
    {
      Fourvector result = *(tensor._right);

      if(tensor._left->isGamma || vector.isGamma)
	result = ( vector * *(tensor._left) ) * result;
      
      else  
	result = ( vector.complexx * tensor._left->complexx ) * result;
                  
      return result;
    }

  else if(tensor._isGamma || vector.isGamma) { // result = vector * _gmatrix
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += 
	  ( vector.isGamma ? 
	    ((vector.gamma).lowerIndex())[j] : 
	    ((vector.complexx).lowerIndex())[j] ) *
	  ( tensor._isGamma ? tensor._gmatrix[j*4+i] : tensor._cmatrix[j*4+i] );
    
    return result;
  }
  else { // result = _cmatrix * vector
    FourVector< complex<double> > result;

    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += 
	  ((vector.complexx).lowerIndex())[j] * tensor._cmatrix[j*4+i];

    return result;
  }
  
  // dummy return
  return vector;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector operator*(const FourVector<double>& vector,
		     const TensorRank2& tensor)
{
  FourVector< complex<double> > cvector(vector[0],vector[1],
					vector[2],vector[3]);

  return ( cvector * tensor );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector operator*(const FourVector< complex<double> >& vector,
		     const TensorRank2& tensor)
{
  if(tensor._isFourVector) // result =(vector* _left) * _right
    {
      Fourvector result = *(tensor._right);

      if(tensor._left->isGamma)
	result = ( vector * tensor._left->gamma ) * result;

      else
	result = ( vector * tensor._left->complexx ) * result;
           
      return result;
    }

  else if(tensor._isGamma) { // result = vector * _gmatrix
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += (vector.lowerIndex())[j] * tensor._gmatrix[j*4+i];
    
    return result;
  }
  else { // result = result = _cmatrix * vector
    FourVector< complex<double> > result;

    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += (vector.lowerIndex())[j] * tensor._cmatrix[j*4+i];

    return result;
  }
  
  // dummy return
  return vector;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourvector operator*(const FourVector<GammaStructure>& vector,
		     const TensorRank2& tensor)
{
  if(tensor._isFourVector) // result =(vector* _left) * _right
    {
      Fourvector result = *(tensor._right);

      if(tensor._left->isGamma)
	result = ( vector * tensor._left->gamma ) * result;

      else
	result = ( vector * tensor._left->complexx ) * result;
           
      return result;
    }

  else if(tensor._isGamma) { // result = vector * _gmatrix
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += (vector.lowerIndex())[j] * tensor._gmatrix[j*4+i];
    
    return result;
  }
  else { // result = _cmatrix * vector
    FourVector<GammaStructure> result;
    
    for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j)
	result[i] += (vector.lowerIndex())[j] * tensor._cmatrix[j*4+i];
    
    return result;
  }
  
  // dummy return
  return vector;
}

//------------------------------------------------------------------------

// Tensor Multiplication
// ---------------------

// another ugly implementation...
TensorRank2 TensorRank2::operator*(const TensorRank2& right) const
{
  // Both tensors are stored internally as matrices
  // result = this->matrix * metric * right.matrix
  if( !this->_isFourVector && ! right._isFourVector ) {
    
    if( this->_isGamma || right._isGamma ) // at least one GammaStructure
      {
	GammaStructure* gmatrix = new GammaStructure[4*4];
	for(int i=0; i<4; ++i) {
	  for(int j=0; j<4; ++j) {
	    gmatrix[i*4+j] = GammaStructure(0.0);
	    for(int k=0; k<4; ++k)
	      gmatrix[i*4+j] += 
		( this->_isGamma ? this->_gmatrix[i*4+k] : this->_cmatrix[i*4+k] ) * 
		( k>0 ? -1.0 : 1.0 ) * // (=metric)
		( right._isGamma ? right._gmatrix[k*4+j] : right._cmatrix[k*4+j] );
	  }
	}
	return TensorRank2(false,true,0,0,0,gmatrix);
      }
    else  // only complex numbers
      { 
	complex<double>* cmatrix = new complex<double>[4*4];
	for(int i=0; i<4; ++i) {
	  for(int j=0; j<4; ++j) {
	    cmatrix[i*4+j] = 0.0;
	    for(int k=0; k<4; ++k)
	      cmatrix[i*4+j] += 
		this->_cmatrix[i*4+k] * 
		( k>0 ? -1.0 : 1.0 ) * right._cmatrix[k*4+j];
	  }
	}
	return TensorRank2(false,false,0,0,cmatrix,0);
      }
  }

  // At least one tensor is stored as Fourvector's
  else {
    
    TensorRank2 result(true,
		       ( (_isGamma || right._isGamma) ? true : false ),
		       0,0,0,0); // this variable will hold the final result.
                                 // for the time being it's initialized
                                 // with all pointers =NULL
    
    if(this->_isFourVector) {
      if(right._isFourVector) {
	// result = this->_left * (factor) * right._right
	// factor = this->_right * right._left
	result._left = new Fourvector( *(this->_left) );
	if(this->_right->isGamma || right._left->isGamma) {
	  GammaStructure factor = *(this->_right) * *(right._left);
	  result._right = new Fourvector( factor * *(right._right) );
	}
	else {
	  complex<double> factor = 
	    this->_right->complexx * right._left->complexx;
	  result._right = new Fourvector( factor * *(right._right) );
	}
      }
      else {
	// result = this->_left * ( this->_right*right)
	result._left = new Fourvector( *(this->_left) );
	result._right = new Fourvector( *(this->_right) * right );
      }
    }
    else {
      // result = (this*right._left) * right._right
      result._left = new Fourvector( *this * *(right._left) );
      result._right = new Fourvector( *(right._right) );
    }

    return result;
  }
  
  // dummy return
  return TensorRank2(false,false,0,0,0,0);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2& TensorRank2::operator*=(const TensorRank2& right)
{
  *this = *this * right;
  
  return *this;
}

//----------------------------------------------------------------------

// TensorRank3 Contraction
// -----------------------

Fourvector TensorRank2::operator%(const TensorRank3& rank3) const
{
  double delta;
  bool init = false;

  toMatrix();

  if(rank3.tensor_[0]._isGamma || rank3.vector_[0].isGamma || _isGamma)
  {
    FourVector<GammaStructure> contraction;
    FourVector<GammaStructure> term;

    for( size_type i=0; i<rank3.tensorLeft_.size(); ++i) {
      for(int mu=0; mu<4; ++mu) {
	for(int nu=0; nu<4; ++nu) {
	  if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	    delta = -1.0;
	  else
	    delta = 1.0;

	  rank3.tensor_[i].toMatrix();

	  if(rank3.tensorLeft_[i])
	  {
	    term = delta * (_isGamma ? _gmatrix[4*mu+nu] : _cmatrix[4*mu+nu]) * 
	           (rank3.tensor_[i]._isGamma ? rank3.tensor_[i]._gmatrix[4*mu+nu] : 
		   rank3.tensor_[i]._cmatrix[4*mu+nu]) * rank3.vector_[i];
	  }
	  else
          {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * (_isGamma ? _gmatrix[4*mu+nu] : _cmatrix[4*mu+nu]) * 
                          (rank3.vector_[i].isGamma ? rank3.vector_[i].gamma[mu] :
		          rank3.vector_[i].complexx[mu]) * (rank3.tensor_[i]._isGamma ? 
			  rank3.tensor_[i]._gmatrix[4*nu+rho] : rank3.tensor_[i]._cmatrix[4*nu+rho]);
	    }
	  }

	  if(init)
          {
	    contraction += term;
	  }
	  else
          {
	    contraction = term;
	    init = true;
	  }
	}
      }
    }
  
    return contraction; 
  }
  else
  {
    FourVector< complex<double> > contraction;
    FourVector< complex<double> > term;

    for( size_type i=0; i<rank3.tensorLeft_.size(); ++i) {
      for(int mu=0; mu<4; ++mu) {
	for(int nu=0; nu<4; ++nu) {
	  if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	    delta = -1.0;
	  else
	    delta = 1.0;

	  rank3.tensor_[i].toMatrix();

	  if(rank3.tensorLeft_[i])
	  {
	    term = delta * _cmatrix[4*mu+nu] * rank3.tensor_[i]._cmatrix[4*mu+nu] * 
		   rank3.vector_[i];
	  }
	  else
          {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * _cmatrix[4*mu+nu] * rank3.vector_[i].complexx[mu] *
			  rank3.tensor_[i]._cmatrix[4*nu+rho];
	    }
	  }

	  if(init)
          {
	    contraction += term;
	  }
	  else
          {
	    contraction = term;
	    init = true;
	  }
	}
      }
    }
  
    return contraction; 
  }

  return FourVector<double>(0.0,0.0,0.0,0.0);
}

//----------------------------------------------------------------------

// TensorRank4 Contraction
// -----------------------

TensorRank2 TensorRank2::operator%(const TensorRank4& rank4) const
{

  TensorRank2 contraction(FourVector<double>(0.0,0.0,0.0,0.0),FourVector<double>(0.0,0.0,0.0,0.0));
  TensorRank2 term;
  double delta;

  toMatrix();

  for(size_type i=0; i<rank4.array_.size(); i+=2) {
    for(int mu=0; mu<4; ++mu) {
      for(int nu=0; nu<4; ++nu) {
	if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	  delta = -1.0;
	else
	  delta = 1.0;

	rank4.array_[i].toMatrix();

	term = delta * (_isGamma ? _gmatrix[4*mu+nu] : _cmatrix[4*mu+nu]) * 
	       (rank4.array_[i]._isGamma ? rank4.array_[i]._gmatrix[4*mu+nu] : 
	       rank4.array_[i]._cmatrix[4*mu+nu]) * rank4.array_[i+1];

	if(term._isGamma && !contraction._isGamma)
	  contraction *= GammaStructure(1.0);
	contraction += term;
      }
    }
  }
  
  return contraction;
}

//------------------------------------------------------------------------

// Output
// ------

void TensorRank2::print()
{
  cout << *this << endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ostream& operator<<(ostream& output, const TensorRank2& tensor)
{
  if(tensor._isFourVector)
    {
      output << "{{ " << *(tensor._left) << " x "
	     << *(tensor._right) << " }}";
    }
  else if(tensor._isGamma)
    {
      output << "{";
      for(int i=0; i<4; i++)
	{
	  output << "{";
	  for(int j=0; j<4; j++)
	    {
	      output << tensor._gmatrix[i*4+j]
		     << ( j<3 ? "," : "" );
	    }
	  output << "}";
	}
      output << "}";
    }
  else
    {
      output << "{";
      for(int i=0; i<4; i++)
	{
	  output << "{";
	  for(int j=0; j<4; j++)
	    {
	      output << tensor._cmatrix[i*4+j]
		     << ( j<3 ? "," : "" );
	    }
	  output << "}";
	}
      output << "}";
    }

  return output;
}
