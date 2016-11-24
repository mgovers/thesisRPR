/*!
 * \file fitting.h
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
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

#ifndef FITTING_H
#define FITTING_H

#include <Structures.h>
#include <string>
#include <vector>
#include <fstream>

#include <iostream>

#define DATAMAX 16384
#define MAXNDIM 40    //!< Max number of free parameters

#define COSGRID 1.0   /* Two values to  obtain an averaging over the cosin */
#define COSDELTA 0.5  /* bins in the recoil polarization. */


struct Photo
{
  short iso;
  char observable[15]; /* a string that specifies the observable */
  short ds_dt;         /* conversion form d_omega to d_t is needed. */
  short ds_du;         /* conversion form d_omega to d_u is needed. */
  short t_ang;         /* 1 if cos as angular variable, 0 if t. */
  short is_w;          /* 1 if W as angular variable, 0 if E_lab. */
  short is_cos_bin;    /* 1 in case cosine bins are specified, 0 if not */
  double emin, emax;
  double cos, cos_max, t;
  double ampli;
  double error;
  int label;		/* unique label of this datapoint */
  short tch;           /* ampli as function of t (tch=1) or cos(theta) (tch=0) */
};


typedef struct
{
  short iso;
  char observable[15]; /* a string that specifies the observable. */
  //short eps_l;         /* d_sigma_L with 'eps' or 'eps_l' in front */ 
  short ds_dt;         /* conversion form d_omega to d_t is needed. */
  double qsquared;
  double s;
  short cos_ang;       /* 1 if cos as angular variable, 0 if t. */
  short is_w;          /* 1 if w as energy variable, 0 if s. */
  double t;
  double cos;
  //
  short cos_binaverage; // bin average (theta_k bins) with uniform acceptance (1) or
                        // linear acceptance (2), 0 if bin average is not required
  double theta_lower;   // lower theta bin limit
  double theta_upper;   // upper theta bin limit
  //
  double ampli;
  double error;
  int label;		/* unique label of this datapoint */
  short tch;           /* ampli as function of t (tch=1) or cos(theta) (tch=0) */
  short beam_ener_input; // 1 if beam energy is given
                         // 0 if epsilon is given
  double e_beam_ener;    // Beam energy (MeV)
  double eps;            // epsilon (if not given --> 0)
  double phi;
  double phiMin;
  short cs_convention;   // cross section convention (->output_format)
  
} Electro;


typedef struct
{
  double pk;
  double err_pk;
  double cos;
  double err_cos;
  double cosmin;    /* lower limit in case of partial total c.s. */
  double cosmax;    /* upper limit in case of partial total c.s. */
  double ampli;
  double ratio;
  double error;
  int label;	/* unique label of this datapoint */
  short iso;
  char observable[15]; /* a string that specifies the observable. */

} Kaoncap;


typedef struct
{
  short iso;
  short photo_prod;
  short electro_prod;
  short kaoncapture;
  Photo photo;
  Electro elec;
  Kaoncap kaoncap;
} Data;


typedef struct{
  double low;
  double up;
  short var;
  short bound;
} Celinfo; 


typedef struct{
  Celinfo partic[PARTICLEMAX][MAXNRFVAR];
} Varinfo;


typedef struct{
  double up;
  double low;
  short bound;
} Limits;


int getdatastructure(Data[], int*, Observable*, Class[]);
int get_variable_info(Varinfo[], Class[], Observable*);
int get_nr_par(int classindex);
int store_vertex_value(double[], Limits[], int*, double, Celinfo);
double get_gridnum(int narrow_grid, char* observable);
int setlabels(Data** datapoints, int* datacount, Observable*);
double chifunc(double*,Class**, Varinfo[], Observable*, 
	       Data**, int[]);
int insertvertex(Class[], Class[], Varinfo[], Observable*, double*);
int copyparticles(Class[], Class[], Observable*);
int simannealing(Class**, Observable*, Data**, 
		 int[]);
double photo_chi(Photo, Class*, Class*, Observable*,int*);
double electro_chi(Electro, Class*, Class*, Observable*,int*);
void is_physical_electro(double, double, double, double, double, double, short);
double calculate_photo_observable(Photo*,Observable*,Class*,Class*,
				  double,double,double,double);
double calculate_electro_observable(Electro*,Observable*,Class*,Class*,
				    double,double,double,double,double,double,double);
double kaoncap_chi(Kaoncap, Class*, Class*, Observable*,int*);
double calculate_kaoncap_observable(Kaoncap*,Observable*,Class*,Class*,
				    double,double,double,double);

int import_exp_rec_saphir_03_1_2(Data[], int*, char*, int,Observable*);
int import_exp_brookhaven_stopped_kaoncap_11_12(Data[], int*, char*, int,Observable*);
int import_exp_diff_crystalball_kaoncap_11(Data[], int*, char*,Observable*);
int import_exp_tot_crystalball_kaoncap_11(Data[], int*, char*,Observable*);
int import_exp_diff_cb2009_kaoncap_11(Data[], int*, char*, Observable*);
int import_exp_tot_cb2009_kaoncap_11(Data[], int*, char*, Observable*);
int import_exp_diff_cb2009_kaoncap_12(Data[], int*, char*, Observable*);
int import_exp_tot_cb2009_kaoncap_12(Data[], int*, char*, Observable*);
int line_count(FILE*, long);

#endif

