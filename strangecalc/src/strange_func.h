/*!
 * \file strange_func.h
 * \ingroup wrapper
 *
 * \author Tamara Corthals
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Stijn Janssen
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
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

#ifndef STRANGE_FUNC_H
#define STRANGE_FUNC_H
 
// Import declaration of datastructures
#include <Structures.h>
#include <string>

int getisospinbasechannel(int);

int observablespecification(Observable*, FILE*, FILE*, FILE*);

int fittingprocess(Observable*); /* fitting */

int particlespecify(Class[], Observable*, FILE*, FILE*);

int import_baseclass_particles(Class[],Observable*,
			       FILE*,FILE*);

int change_isospin_channel(Class[], Observable*);

int change_nickname(Class*, const char*, const char*);

int modify_nicknames(Class[], Observable*);

int magnetic_transition_ratio(Class[], Observable*);

int import_particle_properties(Class[], Observable*);

int change_coupling_constants(Class[], Observable*, int);

int fprintparticspec(FILE*, Class[], Observable);

int findclassindex(char);

char findclasslabel(int);

int check_born_cc(Class[]);

int allocate_born_masses(double*,double*,double*,int);

double calc_epsilon(double, double, double, double);

double kincoeff_photo(double, double, double);

double kincoeff_kaoncap(double, double, double, double);

double spinaveraging(const Observable*);

int lorentztrans_photo(double, double*, double);

int lorentztrans_electro(double, double*, double, double*, double);

int lorentztrans_kaon(double, double*, double, double); /* kaon capture */

double construct_pk(double, double, double, double, double);

double construct_w(double, double, double, double);

double massrelation(double, double, double, double, double, double);

double domega_to_dt_conversion(double, double);

double domega_to_du_conversion(double, double);

void error_exit(const std::string&);

#endif








