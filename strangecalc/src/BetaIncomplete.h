/*!
 * \file BetaIncomplete.h
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

#include <complex>

#ifndef BETAINCOMPLETE_H
#define BETAINCOMPLETE_H

std::complex<double> phase(double);

double beta(double, double);

std::complex<double> beta(double, double, double);

#endif
