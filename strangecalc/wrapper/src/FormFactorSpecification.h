/*!
 * \file FormFactorSpecification.h
 * \ingroup wrapper
 *
 * 2 functions to read in formfactor specifications and store
 * them in the particles[] array
 *
 * 1 function to give a FormFactor different parameters
 *
 * 1 function to delete the dynamically allocated memory
 * for the formfactors
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * 
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

#include <Structures.h>


#ifndef FORMFACTORSPECIFICATION_H
#define FORMFACTORSPECIFICATION_H

int em_formfactor_specification(Observable*,Class*);

int strong_formfactor_specification(Observable*,Class*);

int update_strong_formfactor(Observable*,Class*);

int release_formfactors(Class*);

#endif
