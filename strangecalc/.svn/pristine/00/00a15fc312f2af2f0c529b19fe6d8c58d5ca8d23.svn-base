/*!
 * \file TensorRank3.h
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

#include <iostream>
#include <vector>
#include <complex>
#include "FourVector.h"
#include "FourVectorHandle.h"
class TensorRank2;

#ifndef TENSORRANK3_H
#define TENSORRANK3_H

class TensorRank3
{
  friend class TensorRank2;

 public:
  // Constructors
  // ------------
  TensorRank3();
  TensorRank3(const Fourvector&, const TensorRank2&);
  TensorRank3(const TensorRank2&, const Fourvector&);
  TensorRank3(const Fourvector&, const Fourvector&, const Fourvector&);
  TensorRank3(const TensorRank3&); // copy constructor

  // Destructor
  // ----------
  ~TensorRank3();

  // Assignment
  // ----------
  TensorRank3& operator=(const TensorRank3&);

  // Addition
  // --------
  TensorRank3& operator+=(const TensorRank3&);
  TensorRank3& operator-=(const TensorRank3&);

  TensorRank3 operator+(const TensorRank3&) const;
  TensorRank3 operator-(const TensorRank3&) const;

  // Scalar Multiplication
  // ---------------------
  TensorRank3& operator*=(double);
  TensorRank3& operator*=(const std::complex<double>&);
  TensorRank3& operator*=(const GammaStructure&);

  TensorRank3 operator*(double) const;
  TensorRank3 operator*(const std::complex<double>&) const;
  TensorRank3 operator*(const GammaStructure&) const;

  friend TensorRank3 operator*(double,const TensorRank3&);
  friend TensorRank3 operator*(const std::complex<double>&, const TensorRank3&);
  friend TensorRank3 operator*(const GammaStructure&, const TensorRank3&);
    
  // Vector Multiplication
  // ---------------------
  TensorRank2 operator*(const Fourvector&) const;
  TensorRank2 operator*(const FourVector<double>&) const;
  TensorRank2 operator*(const FourVector< std::complex<double> >&) const;
  TensorRank2 operator*(const FourVector< GammaStructure >&) const;
  
  friend TensorRank2 operator*(const Fourvector&,
			       const TensorRank3&);
  friend TensorRank2 operator*(const FourVector<double>&,
			       const TensorRank3&);
  friend TensorRank2 operator*(const FourVector< std::complex<double> >&,
			       const TensorRank3&);
  friend TensorRank2 operator*(const FourVector< GammaStructure >&,
			       const TensorRank3&);

  // TensorRank2 Multiplication
  // --------------------------
  TensorRank3& operator*=(const TensorRank2&);
  TensorRank3 operator*(const TensorRank2&) const;
  friend TensorRank3 operator*(const TensorRank2&, const TensorRank3&);

  // TensorRank2 Contraction
  // -----------------------
  Fourvector operator%(const TensorRank2&) const;

  // TensorRank4 Contraction
  // -----------------------
  TensorRank3 operator%(const TensorRank4&) const;
  
  // Output
  // ------
  friend std::ostream& operator<<(std::ostream&, const TensorRank3&);

 private:
  // Data members
  // ------------
  std::vector< TensorRank2 > tensor_;
  std::vector< Fourvector > vector_;
  
  std::vector< bool > tensorLeft_;

}; // end class TensorRank3

#endif // TENSORRANK3_H
