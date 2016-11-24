/*!
 * \file TensorRank2.h
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

#ifndef TENSORRANK2_H
#define TENSORRANK2_H

#include <iostream>
#include <complex>
#include "GammaStructure.h"
#include "FourVector.h"
#include "FourVectorHandle.h"

class TensorRank2
{
  friend class TensorRank3;
  friend class TensorRank4;

 public:  
  // Constructors
  // ------------
  TensorRank2();
  TensorRank2(const Fourvector&, const Fourvector&);
  TensorRank2(const Fourvector&, const FourVector<double>&);
  TensorRank2(const FourVector<double>&, const Fourvector&);
  TensorRank2(const Fourvector&, 
	      const FourVector< std::complex<double> >&);
  TensorRank2(const FourVector< std::complex<double> >&, 
	      const Fourvector&);
  TensorRank2(const Fourvector&, 
	      const FourVector< GammaStructure >&);
  TensorRank2(const FourVector< GammaStructure >&,
	      const Fourvector&);
  TensorRank2(const FourVector< std::complex<double> >&,
	      const FourVector< std::complex<double> >&);
  TensorRank2(const FourVector<double>&,
	      const FourVector<double>&);
  TensorRank2(const FourVector<GammaStructure>&,
	      const FourVector<GammaStructure>&);
  TensorRank2(const FourVector<double>&,
	      const FourVector< std::complex<double> >&);
  TensorRank2(const FourVector< std::complex<double> >&,
	      const FourVector<double>&);
  TensorRank2(const FourVector< std::complex<double> >&,
	      const FourVector<GammaStructure>&);
  TensorRank2(const FourVector<GammaStructure>&,
	      const FourVector< std::complex<double> >&);
  TensorRank2(const FourVector<double>&,
	      const FourVector<GammaStructure>&);
  TensorRank2(const FourVector<GammaStructure>&,
	      const FourVector<double>&);
  TensorRank2(const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&,
	      const std::complex<double>&,const std::complex<double>&);
  TensorRank2(const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&,
	      const GammaStructure&,const GammaStructure&);
  TensorRank2(const TensorRank2&); // copy constructor

  // Destructor
  // ----------
  ~TensorRank2();

  // Overloaded assignment operator
  // ------------------------------
  TensorRank2& operator=(const TensorRank2&);

  // Internal conversion
  // -------------------
  void toMatrix() const;

  // Addition
  // --------
  TensorRank2& operator+=(const TensorRank2&);
  TensorRank2& operator-=(const TensorRank2&);
 
  TensorRank2 operator+(const TensorRank2&) const;
  TensorRank2 operator-(const TensorRank2&) const;

  // Scalar Multiplication
  // ---------------------
  TensorRank2& operator*=(double);
  TensorRank2& operator*=(const std::complex<double>&);
  TensorRank2& operator*=(const GammaStructure&);
  
  TensorRank2 operator*(double) const;
  TensorRank2 operator*(const std::complex<double>&) const;
  TensorRank2 operator*(const GammaStructure&) const;

  friend TensorRank2 operator*(double, const TensorRank2&);
  friend TensorRank2 operator*(const std::complex<double>&, const TensorRank2&);
  friend TensorRank2 operator*(const GammaStructure&, const TensorRank2&);

  // Vector Multiplication
  // ---------------------
  Fourvector operator*(const Fourvector&) const;
  Fourvector operator*(const FourVector<double>&) const;
  Fourvector operator*(const FourVector< std::complex<double> >&) const;
  Fourvector operator*(const FourVector< GammaStructure >&) const;

  friend Fourvector operator*(const Fourvector&,
			      const TensorRank2&);
  friend Fourvector operator*(const FourVector<double>&,
			      const TensorRank2&);
  friend Fourvector operator*(const FourVector< std::complex<double> >&,
			      const TensorRank2&);
  friend Fourvector operator*(const FourVector< GammaStructure >&,
			      const TensorRank2&);

  friend TensorRank2 operator*(const Fourvector&,
			       const TensorRank3&);
  friend TensorRank2 operator*(const FourVector<double>&,
			       const class TensorRank3&);
  friend TensorRank2 operator*(const FourVector< std::complex<double> >&,
			       const class TensorRank3&);
  friend TensorRank2 operator*(const FourVector< GammaStructure >&,
			       const class TensorRank3&);

  // Tensor Multiplication
  // ---------------------
  TensorRank2 operator*(const TensorRank2&) const;
  TensorRank2& operator*=(const TensorRank2&);

  // TensorRank3 Contraction
  // -----------------------
  Fourvector operator%(const TensorRank3&) const;
  
  // TensorRank4 Contraction
  // -----------------------
  TensorRank2 operator%(const TensorRank4&) const;

  // Output
  // ------
  void print();
  friend std::ostream& operator<<(std::ostream&, const TensorRank2&);

 private:
  // data members
  // ------------
  bool _isFourVector;         //!< specifies how tensor is stored internaly
  bool _isGamma;
  
  Fourvector* _left;          //!< the FourVector<>'s that built up
  Fourvector* _right;         //!< the tensor
  
  std::complex<double> * _cmatrix; //!< 4x4 matrix containing the direct product (for complex case)
  GammaStructure * _gmatrix; //!< 4x4 matrix containing the direct product (for gamma matrix case)

  // Private constructor
  // -------------------
  TensorRank2(bool,bool,Fourvector*,Fourvector*,
	      std::complex<double>*,GammaStructure*);

}; // end class TensorRank2

#endif // TENSORRANK2_H
