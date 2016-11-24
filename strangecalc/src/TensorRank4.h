/*!
 * \file TensorRank4.h
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
class TensorRank2;

#ifndef TENSORRANK4_H
#define TENSORRANK4_H

class TensorRank4
{
  friend class TensorRank2;
  friend class TensorRank3;

 public:
  // Constructors
  // ------------
  TensorRank4();
  TensorRank4(const TensorRank2&, const TensorRank2&);
  TensorRank4(const Fourvector&, const Fourvector&,
	      const Fourvector&, const Fourvector&);
  TensorRank4(const TensorRank2&, const Fourvector&, const Fourvector&);
  TensorRank4(const Fourvector&, const Fourvector&, const TensorRank2&);
  TensorRank4(const TensorRank4&); // copy constructor

  // Destructor
  // ----------
  ~TensorRank4();

  // Assignment
  // ----------
  TensorRank4& operator=(const TensorRank4&);

  // Addition
  // --------
  TensorRank4& operator+=(const TensorRank4&);
  TensorRank4& operator-=(const TensorRank4&);
 
  TensorRank4 operator+(const TensorRank4&) const;
  TensorRank4 operator-(const TensorRank4&) const;

  // Scalar Multiplication
  // ---------------------
  TensorRank4& operator*=(double);
  TensorRank4& operator*=(const std::complex<double>&);
  TensorRank4& operator*=(const GammaStructure&);

  TensorRank4 operator*(double) const;
  TensorRank4 operator*(const std::complex<double>&) const;
  TensorRank4 operator*(const GammaStructure&) const;

  friend TensorRank4 operator*(double,const TensorRank4&);
  friend TensorRank4 operator*(const std::complex<double>&, const TensorRank4&);
  friend TensorRank4 operator*(const GammaStructure&, const TensorRank4&);

  // Vector Multiplication
  // ---------------------
  TensorRank3 operator*(const Fourvector&) const;
  TensorRank3 operator*(const FourVector<double>&) const;
  TensorRank3 operator*(const FourVector< std::complex<double> >&) const;
  TensorRank3 operator*(const FourVector< GammaStructure >&) const;

  friend TensorRank3
    operator*(const Fourvector&, const TensorRank4&);
  friend TensorRank3
    operator*(const FourVector<double>&, const TensorRank4&);
  friend TensorRank3
    operator*(const FourVector< std::complex<double> >&, const TensorRank4&);
  friend TensorRank3
    operator*(const FourVector< GammaStructure >&, const TensorRank4&);

  // TensorRank2 Multiplication
  // --------------------------
  TensorRank4& operator*=(const TensorRank2&);
  TensorRank4 operator*(const TensorRank2&) const;
  friend TensorRank4 operator*(const TensorRank2&, const TensorRank4&);

  // TensorRank2 Contraction
  // -----------------------
  TensorRank2 operator%(const TensorRank2&) const;

  // Output
  // ------
  friend std::ostream& operator<<(std::ostream&, const TensorRank4&);

  

 private:
  // Data members
  // ------------
  std::vector< TensorRank2 > array_;

}; // end class TensorRank4

#endif // TENSORRANK4_H
