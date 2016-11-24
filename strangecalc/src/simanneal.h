/*!
 * \file simanneal.h
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

#ifndef SIMANNEAL_H
#define SIMANNEAL_H

 
/*! Termination tolerance for the simulated annealing procedure */
#define FTOL 0.005

/* Random generator parameters */
#define IA 16807 //!< Random generator parameter
#define IM 2147483647 //!< Random generator parameter
#define AM (1.0/IM) //!< Random generator parameter
#define IQ 127773 //!< Random generator parameter
#define IR 2836 //!< Random generator parameter
#define NTAB 32 //!< Random generator parameter
#define NDIV (1+(IM-1)/NTAB) //!< Random generator parameter
#define EPS 1.2e-7 //!< Random generator parameter
#define RNMX (1.0-EPS) //!< Random generator parameter

/* #define MAXLOOP 50 */
#define MAXLOOP 200

int makestartvertex(Class[], Varinfo[], Observable, double[][MAXNDIM], 
		    Limits[], int*);
void fixtempsearch(double[][MAXNDIM], double[], int, double[],  double*,
		   double, int*, double,Class**, 
		   Varinfo[], Observable, Data**, int[], Limits[]);
double scalesimplex(double[][MAXNDIM], double[], double[], int, double[], 
		    double*, int, double*, double,Class**, 
		    Varinfo[], Observable, Data**, int[], 
		    Limits[]);
double ran1(long*);
int fprintlog(char*, Observable, Class**, Varinfo[], int);
int printchiset(double[], int);
int fprintchiset(FILE*, double[], int);
int fprintlogstatus(char*, double, double, Class[], Observable);
int give_param_info(FILE*, int*, double, Celinfo);
int fprintfinallog(char*, char*, double, double, Class[], Observable);
int fprinttmpinfo(char*, double[][MAXNDIM], int);

#endif


