/*!
 * \file TensorRank3.cpp
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
#include <GammaStructure.h>
#include <FourVector.h>
#include "FourVectorHandle.h"
#include "TensorRank2.h"
#include "TensorRank3.h"
#include "TensorRank4.h"

using namespace std;

typedef vector<int>::size_type size_type;

/*!
 * \class TensorRank3
 *
 * Represents direct products of FourVector's:
 * rank3 = a^{alpha} b^{beta} c^{gamma}
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

TensorRank3::TensorRank3()
{
  tensor_.push_back( TensorRank2(0.0,0.0,0.0,0.0,
				 0.0,0.0,0.0,0.0,
				 0.0,0.0,0.0,0.0,
				 0.0,0.0,0.0,0.0) );
  vector_.push_back( FourVector<double>() );
  tensorLeft_.push_back( true );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3::TensorRank3(const Fourvector& vector,
			 const TensorRank2& tensor)
{
  tensor_.push_back( tensor );
  vector_.push_back( vector );
  tensorLeft_.push_back( false );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3::TensorRank3(const TensorRank2& tensor,
			 const Fourvector& vector)
{
  tensor_.push_back( tensor );
  vector_.push_back( vector );
  tensorLeft_.push_back( true );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3::TensorRank3(const Fourvector& vector1,
			 const Fourvector& vector2,
			 const Fourvector& vector3)
{
  tensor_.push_back( TensorRank2(vector1,vector2) );
  vector_.push_back( vector3 );
  tensorLeft_.push_back( true );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3::TensorRank3(const TensorRank3& toCopy)
  : tensor_( toCopy.tensor_),
    vector_( toCopy.vector_),
    tensorLeft_( toCopy.tensorLeft_ )
{
  // empty body
}

//----------------------------------------------------------------------

// Destructor
// ----------

TensorRank3::~TensorRank3()
{
  // empty body
}

//----------------------------------------------------------------------

// Assignment
// ----------

TensorRank3& TensorRank3::operator=(const TensorRank3& toCopy)
{
  tensor_ = toCopy.tensor_;
  vector_ = toCopy.vector_;
  tensorLeft_ = toCopy.tensorLeft_;

  return *this;
}

//----------------------------------------------------------------------

// Addition
// --------

TensorRank3& TensorRank3::operator+=(const TensorRank3& right)
{
  tensor_.reserve( tensor_.size() + right.tensor_.size() );
  vector_.reserve( vector_.size() + right.vector_.size() );
  tensorLeft_.reserve( tensorLeft_.size() + right.tensorLeft_.size() );

  for(size_type i=0; i<right.tensor_.size(); ++i) {
    tensor_.push_back( right.tensor_[i] );
    vector_.push_back( right.vector_[i] );
    tensorLeft_.push_back( right.tensorLeft_[i] );
  }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3& TensorRank3::operator-=(const TensorRank3& right)
{
  tensor_.reserve( tensor_.size() + right.tensor_.size() );
  vector_.reserve( vector_.size() + right.vector_.size() );
  tensorLeft_.reserve( tensorLeft_.size() + right.tensorLeft_.size() );

  for(size_type i=0; i<right.tensor_.size(); ++i) {
    tensor_.push_back( right.tensor_[i] );
    vector_.push_back( -1.0 * right.vector_[i] );
    tensorLeft_.push_back( right.tensorLeft_[i] );
  }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator+(const TensorRank3& right) const
{
  TensorRank3 sum = *this;
  sum += right;

  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator-(const TensorRank3& right) const
{
  TensorRank3 sum = *this;
  sum -= right;

  return sum;
}

//----------------------------------------------------------------------

// Scalar Multiplication
// ---------------------

TensorRank3& TensorRank3::operator*=(double factor)
{
  for(size_type i=0; i<vector_.size(); ++i)
    vector_[i] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3& TensorRank3::operator*=(const complex<double>& factor)
{
  for(size_type i=0; i<vector_.size(); ++i)
    vector_[i] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3& TensorRank3::operator*=(const GammaStructure& factor)
{
  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      vector_[i] *= factor;
    else
      tensor_[i] *= factor;
  }

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator*(double factor) const
{
  TensorRank3 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator*(const complex<double>& factor) const
{
  TensorRank3 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator*(const GammaStructure& factor) const
{
  TensorRank3 product = *this;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(double factor, const TensorRank3& tensor)
{
  TensorRank3 product = tensor;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const complex<double>& factor, const TensorRank3& tensor)
{
  TensorRank3 product = tensor;
  product *= factor;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const GammaStructure& factor, const TensorRank3& tensor)
{
  TensorRank3 product = tensor;
  for(size_type i=0; i<product.tensorLeft_.size(); ++i) {
    if(product.tensorLeft_[i])
      product.tensor_[i] = factor * product.tensor_[i];
    else
      product.vector_[i] = factor * product.vector_[i];
  }

  return product;
}

//----------------------------------------------------------------------

// Vector Multiplication
// ---------------------

TensorRank2 TensorRank3::operator*(const Fourvector& right) const
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      {
	if(vector_[i].isGamma || right.isGamma)
	  term = tensor_[i] * (vector_[i] * right);
	else
	  term = tensor_[i] * (vector_[i].complexx * right.complexx);
      }
    else
      {
	term = TensorRank2(vector_[i],tensor_[i]*right);
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank3::operator*(const FourVector<double>& right) const
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      {
	if(vector_[i].isGamma)
	  term = tensor_[i] * (vector_[i].gamma * right);
	else
	  term = tensor_[i] * (vector_[i].complexx * 
			  FourVector<complex<double> >(right[0],right[1],
						       right[2],right[3]));
      }
    else
      {
	term = TensorRank2(vector_[i],tensor_[i]*right);
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank3::operator*(const FourVector< complex<double> >& right) const
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      {
	if(vector_[i].isGamma)
	  term = tensor_[i] * (vector_[i].gamma * right);
	else
	  term = tensor_[i] * (vector_[i].complexx * right);
      }
    else
      {
	term = TensorRank2(vector_[i],tensor_[i]*right);
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 TensorRank3::operator*(const FourVector< GammaStructure >& right) const
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      {
	if(vector_[i].isGamma)
	  term = tensor_[i] * (vector_[i].gamma * right);
	else
	  term = tensor_[i] * (vector_[i].complexx * right);
      }
    else
      {
	term = TensorRank2(vector_[i],tensor_[i]*right);
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const Fourvector& left,
		      const TensorRank3& right)
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<right.tensorLeft_.size(); ++i) {
    if(right.tensorLeft_[i])
      {
	term = TensorRank2(left*right.tensor_[i],right.vector_[i]);
      }
    else
      {
	if(left.isGamma || right.vector_[i].isGamma)
	  term = (left * right.vector_[i]) * right.tensor_[i];
	else
	  term = (left.complexx * right.vector_[i].complexx) * right.tensor_[i];
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const FourVector<double>& left,
		      const TensorRank3& right)
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<right.tensorLeft_.size(); ++i) {
    if(right.tensorLeft_[i])
      {
	term = TensorRank2(left*right.tensor_[i],right.vector_[i]);
      }
    else
      {
	if(right.vector_[i].isGamma)
	  term = (left * right.vector_[i].gamma) * right.tensor_[i];
	else
	  term = (FourVector<complex<double> >(left[0],left[1],left[2],left[3])*
		  right.vector_[i].complexx) * right.tensor_[i];
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const FourVector< complex<double> >& left,
		      const TensorRank3& right)
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<right.tensorLeft_.size(); ++i) {
    if(right.tensorLeft_[i])
      {
	term = TensorRank2(left*right.tensor_[i],right.vector_[i]);
      }
    else
      {
	if(right.vector_[i].isGamma)
	  term = (left * right.vector_[i].gamma) * right.tensor_[i];
	else
	  term = (left * right.vector_[i].complexx) * right.tensor_[i];
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank2 operator*(const FourVector< GammaStructure >& left,
		      const TensorRank3& right)
{
  TensorRank2 product(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  TensorRank2 term;

  for(size_type i=0; i<right.tensorLeft_.size(); ++i) {
    if(right.tensorLeft_[i])
      {
	term = TensorRank2(left*right.tensor_[i],right.vector_[i]);
      }
    else
      {
	if(right.vector_[i].isGamma)
	  term = (left * right.vector_[i].gamma) * right.tensor_[i];
	else
	  term = (left * right.vector_[i].complexx) * right.tensor_[i];
      }
    
    if(term._isGamma && !product._isGamma) 
      product *= GammaStructure(1.0);
    product += term;
  }
  
  return product;
}

//----------------------------------------------------------------------

// TensorRank2 Multiplication
// --------------------------

TensorRank3& TensorRank3::operator*=(const TensorRank2& rank2)
{
  for(size_type i=0; i<tensorLeft_.size(); ++i) {
    if(tensorLeft_[i])
      vector_[i] = vector_[i] * rank2;
    else
      tensor_[i] *= rank2;      
  }
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 TensorRank3::operator*(const TensorRank2& rank2) const
{
  TensorRank3 product = *this;
  product *= rank2;

  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorRank3 operator*(const TensorRank2& rank2, const TensorRank3& rank3)
{
  TensorRank3 product = rank3;
  
  for(size_type i=0; i<product.tensorLeft_.size(); ++i) {
    if(product.tensorLeft_[i])
      product.tensor_[i] = rank2 * product.tensor_[i];
    else
      product.vector_[i] = rank2 * product.vector_[i];
  }

  return product;
}

//----------------------------------------------------------------------

// TensorRank2 Contraction
// -----------------------

Fourvector TensorRank3::operator%(const TensorRank2& right) const
{

  double delta;
  bool init = false;

  right.toMatrix();
  
  if(tensor_[0]._isGamma || vector_[0].isGamma || right._isGamma)
  {
    FourVector<GammaStructure> contraction;
    FourVector<GammaStructure> term;

    for(size_type i=0; i<tensorLeft_.size(); ++i) {
      for(int mu=0; mu<4; ++mu) {
	for(int nu=0; nu<4; ++nu) {
	  if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	    delta = -1.0;
	  else
	    delta = 1.0;
	  
	  tensor_[i].toMatrix();

	  if(tensorLeft_[i])
	  {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * (tensor_[i]._isGamma ? tensor_[i]._gmatrix[4*rho+mu] : 
			  tensor_[i]._cmatrix[4*rho+mu]) * (vector_[i].isGamma ? 
			  vector_[i].gamma[nu] : vector_[i].complexx[nu]) * 
		(right._isGamma ? right._gmatrix[4*mu+nu] : right._cmatrix[4*mu+nu]);
	    }
	  }
	  else
	  {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * (vector_[i].isGamma ? vector_[i].gamma[rho] : 
			  vector_[i].complexx[rho]) * (tensor_[i]._isGamma ? 
			  tensor_[i]._gmatrix[4*mu+nu] : tensor_[i]._cmatrix[4*mu+nu]) * 
		(right._isGamma ? right._gmatrix[4*mu+nu] : right._cmatrix[4*mu+nu]);
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

    for(size_type i=0; i<tensorLeft_.size(); ++i) {
      for(int mu=0; mu<4; ++mu) {
	for(int nu=0; nu<4; ++nu) {
	  if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	    delta = -1.0;
	  else
	    delta = 1.0;
	  
	  tensor_[i].toMatrix();

	  if(tensorLeft_[i])
	  {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * tensor_[i]._cmatrix[4*rho+mu] * vector_[i].complexx[nu] * 
                          right._cmatrix[4*mu+nu];
	    }
	  }
	  else
	  {
	    for(int rho=0; rho<4; ++rho) {
	      term[rho] = delta * vector_[i].complexx[rho] * tensor_[i]._cmatrix[4*mu+nu] * 
                          right._cmatrix[4*mu+nu];
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

TensorRank3 TensorRank3::operator%(const TensorRank4& rank4) const
{
  TensorRank3 contraction;
  TensorRank3 term;	
  bool init = false;
  double delta;

  if(tensor_[0]._isGamma || vector_[0].isGamma || rank4.array_[0]._isGamma)
  {
    FourVector<GammaStructure> term1;

    for(size_type j=0; j<tensorLeft_.size(); ++j) {
      for(size_type i=0; i<rank4.array_.size(); i+=2) {
    	for(int mu=0; mu<4; ++mu) {
	  for(int nu=0; nu<4; ++nu) {
	    if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	      delta = -1.0;
	    else
	      delta = 1.0;

	    tensor_[j].toMatrix();
	    rank4.array_[i].toMatrix();
	    
	    if(tensorLeft_[j])
	    {
	      for(int rho=0; rho<4; ++rho) {
		term1[rho] = delta * (tensor_[j]._isGamma ? tensor_[j]._gmatrix[4*rho+mu] :
			     tensor_[j]._cmatrix[4*rho+mu]) * (vector_[j].isGamma ? 
			     vector_[j].gamma[nu] : vector_[j].complexx[nu]) * 
		             (rank4.array_[i]._isGamma ? rank4.array_[i]._gmatrix[4*mu+nu] :
                             rank4.array_[i]._cmatrix[4*mu+nu]);
	      }
	    }
	    else
	    {
	      term1 = delta * vector_[j] * (tensor_[j]._isGamma ? tensor_[j]._gmatrix[4*mu+nu] :
                      tensor_[j]._cmatrix[4*mu+nu]) * (rank4.array_[i]._isGamma ? 
		      rank4.array_[i]._gmatrix[4*mu+nu] : rank4.array_[i]._cmatrix[4*mu+nu]);
	    }
	
	    term = TensorRank3(term1, rank4.array_[i+1]);
					
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
    }

    return contraction;	
  }
  else
  {
    FourVector< complex<double> > term1;

    for(size_type j=0; j<tensorLeft_.size(); ++j) {
      for(size_type i=0; i<rank4.array_.size(); i+=2) {
    	for(int mu=0; mu<4; ++mu) {
	  for(int nu=0; nu<4; ++nu) {
	    if((mu==0 && nu!=0)||(mu!=0 && nu==0))
	      delta = -1.0;
	    else
	      delta = 1.0;
	    
	    tensor_[j].toMatrix();
	    rank4.array_[i].toMatrix();

	    if(tensorLeft_[j])
	    {
	      for(int rho=0; rho<4; ++rho) {
		term1[rho] = delta * tensor_[j]._cmatrix[4*rho+mu] * vector_[j].complexx[nu] * 
		             rank4.array_[i]._cmatrix[4*mu+nu];
	      }
	    }
	    else
	    {
	      term1 = delta * vector_[j] * tensor_[j]._cmatrix[4*mu+nu] *
		      rank4.array_[i]._cmatrix[4*mu+nu];
	    }
	
	    term = TensorRank3(term1, rank4.array_[i+1]);
					
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
    }

    return contraction;	
  }

  return TensorRank3(FourVector<double>(0.0,0.0,0.0,0.0),
		     TensorRank2(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
				 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0));
}

//----------------------------------------------------------------------

// Output
// ------

ostream& operator<<(ostream& stream, const TensorRank3& tensor)
{
  stream << "{" << endl;
  
  for(size_type i=0; i<tensor.tensorLeft_.size(); ++i) {
    if(tensor.tensorLeft_[i])
      stream << tensor.tensor_[i] << " X " << tensor.vector_[i] 
	     << endl << "+";
    else
      stream << tensor.vector_[i] << " X " << tensor.tensor_[i] 
	     << endl << "+";
  }
  
  stream << "\r" << "}" << endl;

  return stream;
}
