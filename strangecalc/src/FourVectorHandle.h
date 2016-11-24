/*!
 * \file FourVectorHandle.h
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

#ifndef FOURVECTORHANDLE_H
#define FOURVECTORHANDLE_H

#include <iostream>
#include <complex>
#include <FourVector.h>
#include <GammaStructure.h>

/*!
 * \class Fourvector
 *
 * \brief Handle for objects of type FourVector< complex<double> > or FourVector< GammaStructure >.
 *
 * Both types can be stored in the object. Using the << operator
 * they can be cast back to their original type.
 * Cast operators for explicit casting are also provided.

 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * 
 */
class Fourvector
{
  // Objects of type TensorRank2/3/4 and there friends
  // are allowed to see the intimate parts of a Fourvector
  friend class TensorRank2;
  friend class TensorRank3;
  friend class TensorRank4;
  friend Fourvector operator*(const Fourvector&,
			      const class TensorRank2&);
  friend Fourvector operator*(const FourVector<double>&,
			      const class TensorRank2&);
  friend Fourvector operator*(const FourVector< std::complex<double> >&,
			      const class TensorRank2&);
  friend Fourvector operator*(const FourVector< GammaStructure >&,
			      const class TensorRank2&);
  friend TensorRank2 operator*(const Fourvector&,
			       const class TensorRank3&);
  friend TensorRank2 operator*(const FourVector<double>&,
			       const class TensorRank3&);
  friend TensorRank2 operator*(const FourVector< std::complex<double> >&,
			       const class TensorRank3&);
  friend TensorRank2 operator*(const FourVector< GammaStructure >&,
			       const class TensorRank3&);
  friend TensorRank3 operator*(const Fourvector&,
			       const class TensorRank4&);
  friend TensorRank3 operator*(const FourVector<double>&,
			       const class TensorRank4&);
  friend TensorRank3 operator*(const FourVector< std::complex<double> >&,
			       const class TensorRank4&);
  friend TensorRank3 operator*(const FourVector< GammaStructure >&,
			       const class TensorRank4&);

 public:
  // Constructors
  // ------------
  Fourvector( const FourVector< double >& );
  Fourvector( const FourVector< std::complex<double> >& );
  Fourvector( const FourVector< GammaStructure >& );
  Fourvector( const Fourvector& );

  // Assignment
  // ----------
  Fourvector& operator=(const Fourvector&);

  // Casting
  // -------
  operator const FourVector< std::complex<double> >&() const;
  operator const FourVector< GammaStructure >&() const;
  operator FourVector< std::complex<double> >&();
  operator FourVector< GammaStructure >&();

  friend
    FourVector< std::complex<double> >&
    operator<<(FourVector< std::complex<double> >&, const Fourvector&);
  friend
    FourVector< GammaStructure >&
    operator<<(FourVector< GammaStructure >&, const Fourvector&);

  // Hermitian Conjugate
  // -------------------
  Fourvector& HermitianConjugate();

  // Scalar Multiplication
  // ---------------------
  Fourvector& operator*=(double);
  Fourvector& operator*=(const std::complex<double>&);
  Fourvector& operator*=(const GammaStructure&);

  Fourvector operator*(double) const;
  Fourvector operator*(const std::complex<double>&) const;
  Fourvector operator*(const GammaStructure&) const;

  friend
    Fourvector operator*(double, const Fourvector&);
  friend
    Fourvector operator*(const std::complex<double>&, const Fourvector&);
  friend
    Fourvector operator*(const GammaStructure&, const Fourvector&);

  // Vector Multiplication
  // ---------------------
  // Always returns a GammaStructure. Thus * does not always act
  // like a user would expect. !! ATTENTION !!
  GammaStructure operator*(const Fourvector&) const;

  // Output
  // ------
  friend std::ostream& operator<<(std::ostream&, const Fourvector&);

 private:
  // Data members
  // ------------
  FourVector< std::complex<double> > complexx; ///< the actual FourVector<>
  FourVector< GammaStructure  > gamma;

  bool isGamma; ///< boolean to specify how the Fourvector is stored internaly
};

#endif
