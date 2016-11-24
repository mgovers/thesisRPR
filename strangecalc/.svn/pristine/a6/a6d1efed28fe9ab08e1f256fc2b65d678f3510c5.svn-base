/*!
 * \file fitting.cpp
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 * \author Jannes Nys <jannes.nys@ugent.be>

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

#include "version.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TMatrixElement.h>
#include <TDataset.h>
#include "strange_func.h"
#include "fitting.h"
#include "calcMatrixElement.h"
#include "FormFactorSpecification.h"
#include <gsl/gsl_integration.h>
#include <DataHandler.h>
#include <stdexcept>
#include <complex>

#include <iomanip>

using std::cerr;
using std::string;
using std::ifstream;
using std::stringstream;
using std::vector;

/*!
 * In this function, the available data points are read in and stored
 * in the "diffcrdata[]" structure. The coupling constants are determined
 * by a fit to this data points.
 */
int fittingprocess(Observable* observ)
{
  int datacount[ISOMAX];
  int i, iso;
  double* vertex;
  double chisquare;
  Varinfo* varinfo;
  Data **datapoints;
  datapoints = (Data**)calloc(ISOMAX,sizeof(Data*));
  Class **particles;
  particles = (Class**)calloc(ISOMAX,sizeof(Class*));

  for(i=0; i<ISOMAX; ++i)
  {
    datapoints[i] = (Data*)calloc(DATAMAX,sizeof(Data));
    particles[i] = (Class*)calloc(CLASSMAX,sizeof(Class));
  }


  observ->photoprod = 0;
  observ->electroprod = 0;
  observ->kaoncapture = 0;

  // initialize datacount[] to zero
  for(i=0; i<ISOMAX; i++)
    datacount[i] = 0;

  /* First we construct the
   * particle specifications in the "particles[iso][]" array
   * and import the data points for the base isospin channel.
   */
  iso = observ->iso.iso_base;

  observ->iso.isospin = iso;

  /*
   * The particle properties are read in from the
   * ./input/numinput/coupl.iso.* file.
   */

  FILE* ifp = stdin;
  FILE* ofp = stdout;

  particlespecify(particles[iso], observ, ifp, ofp);


  /*
   * We need to change isospin channel
   * change_isospin_channel() can not be used however,
   * because the coupling constants may not be changed yet!
   * (this would disturb makestartvertex()).
   */

  modify_nicknames(particles[iso],observ);
  import_particle_properties(particles[iso],observ);
  if(iso > 3)
    magnetic_transition_ratio(particles[iso],observ);

  /*
   * Store Form factor info in particles[iso][]
   */
  if( observ->hadronformfac )
    strong_formfactor_specification(observ,particles[iso]);

  if ( observ->fit.elec_diffcs[iso] )
    {
      observ->electroprod = 1;
      em_formfactor_specification(observ,particles[iso]);
      observ->electroprod = 0;
    }


  /* Store the data in an easy accessible structure */

  getdatastructure(datapoints[iso], &datacount[iso], observ,
		   particles[iso]);


  for(i=0; i<observ->iso.nr_iso_channels; i++)
    {
      /*
       * For each of the channels under investigation, construct the
       * particle specifications in the "particles[iso][]" array, import
       * the data points isospin channel per isospin channel.
       */

      iso = observ->iso.iso_channel[i];

      if( iso!=observ->iso.iso_base ) {
	/*
	 * We already initialized the base channel!
	 */

	observ->iso.isospin = iso;

	/*
	 * The particle properties are copied from
	 * the base channel particles
	 */

	copyparticles(particles[iso],
		      particles[observ->iso.iso_base],
		      observ);

	/*
	 * We need to change isospin channel
	 * change_isospin_channel() can not be used however,
	 * because the coupling constants may not be changed yet!
	 * (this would disturb makestartvertex()).
	 */

	modify_nicknames(particles[iso],observ);
	import_particle_properties(particles[iso],observ);
	if(iso > 3)
	  magnetic_transition_ratio(particles[iso],observ);

	/*
	 * Store Form factor info in particles[iso][]
	 */
	if( observ->hadronformfac )
	  strong_formfactor_specification(observ,particles[iso]);

	if ( observ->fit.elec_diffcs[iso] )
	  {
	    observ->electroprod = 1;
	    em_formfactor_specification(observ,particles[iso]);
	    observ->electroprod = 0;
	  }


	/* Store the data in an easy accessible structure */

	getdatastructure(datapoints[iso], &datacount[iso], observ,
			 particles[iso]);
      }
    } // end loop isospin channels

    //set the labels for caching
    setlabels(datapoints, datacount, observ);

  if(observ->fitting)
    {
      /* Start the simulated annealing procedure */

      simannealing(particles, observ, datapoints, datacount);

    }
  else if(observ->chisquare)
    {
      /*
       * Calculate the chi squared
       *
       * Since the vertex[] and varinfo[] has no meaning here,
       * we put them equal to NULL
       */

      vertex = NULL;
      varinfo = NULL;

      chisquare = chifunc(vertex,particles, varinfo, observ,
			  datapoints, datacount);

      printf("chi squared for this set of coupling constants is = %f\n",
	     chisquare);


      /*
       * If there is more than one channel used in the chi squared
       * calculation, print the base set, if not, print the used set.
       */

      if(observ->iso.nr_iso_channels > 1 ||
	 observ->iso.iso_base == observ->iso.isospin)
	{
	  fprintparticspec(stdout,particles[observ->iso.iso_base], *observ);
	}
      else if(observ->iso.iso_base != observ->iso.isospin)
	{

	  /*
	   * ATTENTION: the nicknames are OK for this iso channel but the
	   * (printed) coefficients are still the one from the "base" set!!
	   */

	  fprintparticspec(stdout,particles[observ->iso.isospin], *observ);
	}

    }

  /* The memory that was allocated for the form factors
   * needs to be released
   */
  for(i=0; i<observ->iso.nr_iso_channels; i++)
    {
      iso = observ->iso.iso_channel[i];
      release_formfactors(particles[iso]);
    }

  for(i=ISOMAX; i>0; --i)
  {
    free(datapoints[i-1]);
    free(particles[i-1]);
  }
  free(datapoints);
  free(particles);

  return 0;
}

/*!
 * In this function, the data files in the directory  "./data/iso.#/"
 * are read in and stored in the datapoints[] structure according to
 * the files specified in observ.fit.database_info[][].
 *
 * Remark that for every different isospin channel this function
 * is invoked again. So, all the datapoints are for a specific isospin
 * channel.
 */
int getdatastructure(Data datapoints[], int* datacount, Observable* observ,
		     Class particles[])
{
  int i, iso;
  double mp, mk, my;
  char filename[EXPNAME];

  // Set isospin channel
  iso = observ->iso.isospin;

  // Allocate Born masses
  allocate_born_masses(&mp,&mk,&my,iso);

  *datacount = 0;


  for(i=0; i < observ->fit.nr_experiments[iso]; i++)
    {
      strcpy(filename, observ->fit.database_info[iso][i]);

      if(iso == 11)
	{
	  if(!strcmp(filename, "brookhaven_89"))
	    import_exp_brookhaven_stopped_kaoncap_11_12(datapoints, datacount,
							filename, iso,observ);
	  else if(!strcmp(filename, "diff.crystalball_2001"))
	    import_exp_diff_crystalball_kaoncap_11(datapoints, datacount,
						   filename,observ);
	  else if(!strcmp(filename, "tot.crystalball_2001"))
	    import_exp_tot_crystalball_kaoncap_11(datapoints, datacount,
						  filename,observ);
	  else if(!strcmp(filename, "diff.cbc_ags_2009"))
	    import_exp_diff_cb2009_kaoncap_11(datapoints, datacount,
					      filename,observ);
	  else if(!strcmp(filename, "tot.cbc_ags_2009"))
	    import_exp_tot_cb2009_kaoncap_11(datapoints, datacount,
					     filename,observ);
	}
      else if(iso == 12)
	{
	  if(!strcmp(filename, "brookhaven_89"))
	    import_exp_brookhaven_stopped_kaoncap_11_12(datapoints, datacount,
							filename, iso,observ);
	  else if(!strcmp(filename, "diff.cbc_ags_2009"))
	    import_exp_diff_cb2009_kaoncap_12(datapoints, datacount,
					      filename,observ);
	  else if(!strcmp(filename, "tot.cbc_ags_2009"))
	    import_exp_tot_cb2009_kaoncap_12(datapoints, datacount,
					     filename,observ);
	  else if(!strcmp(filename, "tot.cbc_ags_stanislaus_2009"))
	    import_exp_tot_cb2009_kaoncap_12(datapoints, datacount,
					     filename,observ);
	}
      else
        DataHandler::importData(datapoints, datacount, filename);
    }


  return 0;
}


/*!
 * Determine which variables are allowed to vary and which have
 * to be static. If a variable is allowed to vary, it can be a free
 * parameter or a bound one. For the latter, determine the limits.
 */
int get_variable_info(Varinfo varinfo[], Class particles[],
		      Observable* observ)
{
  int i, j, k, match, nr_par;
  char dump[DUMP_SIZE], classlabel[7], tmpnickname[3], ans[4],
    varfile[MAX_LOCATION_STRING], iso[2];
  FILE* ifp;


  // Open varinfo.* file
  strcpy(varfile, observ->inFolder);
  strcat(varfile, "numinput/varinfo.");
  iso[0] = observ->iso.iso_base + 48;
  iso[1] = '\0';
  strcat(varfile, iso);

  ifp = fopen(varfile, "r");
  if(ifp == NULL)
    error_exit("Error in opening varinfo-file!!\n");


  // Dump the header of the varinfo-file.
  char line[BUFFER_SIZE];
  fgets(line, BUFFER_SIZE, ifp);
  fgets(line, BUFFER_SIZE, ifp);

  // Read the class label
  fscanf(ifp, "%s", classlabel);

  /* Scroll the varinfo file until all info is stored */
  while(strcmp(classlabel, "cutoff"))
    {
      i = findclassindex(classlabel[0]);

      match = 0;

      // read nickname
      fscanf(ifp, "%s", tmpnickname);

      // When particles of this class are present in coupl.iso.*, proceed
      if(particles[i].particount != 0)
	{
	  // Look for the corresponding particle(tmpnickname) in particles[]
	  for(j = 0; j < particles[i].particount && !match; j++)
	    if(!strcmp(particles[i].partic[j].nickname,tmpnickname))
	      {
		match = 1;
		nr_par = get_nr_par(i);


		for(k = 0; k < nr_par; k++)
		  {
		    fscanf(ifp, "%s", ans);
		    // fixed parameter
		    if(!strcmp(ans, "x"))
		      varinfo[i].partic[j][k].var = 0;
		    // unbound parameter
		    else if(!strcmp(ans, "f,n"))
		      {
			varinfo[i].partic[j][k].var = 1;
			varinfo[i].partic[j][k].bound = 0;
			varinfo[i].partic[j][k].low = 0;
			varinfo[i].partic[j][k].up = 0;
		      }
		    // free parameter within limits
		    else if(!strcmp(ans, "f,l"))
		      {
			varinfo[i].partic[j][k].var = 1;
			varinfo[i].partic[j][k].bound = 1;

			fscanf(ifp, "%lf%lf",
			       &varinfo[i].partic[j][k].low,
			       &varinfo[i].partic[j][k].up);
		      }
		    else
		      error_exit("Error in varinfo file.\n");
		  }
	      }

	  // Could not find particle with tmpnickname in particles[i]
	  if(!match)
	    {
	      fprintf(stderr, "Problem with class %c:\n", findclasslabel(i));
	      fprintf(stderr, "match = %d, particount[%d] = %d\nnickname in varinfo: %s\n", match,
		      i, particles[i].particount, tmpnickname);
	      fprintf(stderr, "nickname in coupl.iso: %s\n", particles[i].partic[0].nickname);
	      error_exit("No matching \"coupl.iso.#\" and \"varinfo\" files!!!\n");
	    }

	} // end "particles[i].particount != 0"


      // No particles of this class present in coupl.iso.*
      else
        error_exit("No matching \"coupl.iso.#\" and \"varinfo\" files!\n");

      fscanf(ifp, "%s", classlabel);
    }

  /*
   * After the particle information, the formfactor information is stored
   * in the (CLASSMAX)th cel of the varinfo array.
   */

  i = CLASSMAX;

  for(k=0; k < 2; k++)
    {

      // dump "cutoff" (only second time)
      if(k != 0)
	fscanf(ifp, "%s", dump);

      // dump "born:" or "res:"
      fscanf(ifp, "%s%s", dump, ans);


      // fixed parameter
      if(!strcmp(ans, "x"))
	varinfo[i].partic[0][k].var = 0;
      // unbound parameter
      else if(!strcmp(ans, "f,n"))
	{
	  varinfo[i].partic[0][k].var = 1;
	  varinfo[i].partic[0][k].bound = 0;
	}
      // free parameter within limits
      else if(!strcmp(ans, "f,l"))
	{
	  varinfo[i].partic[0][k].var = 1;
	  varinfo[i].partic[0][k].bound = 1;

	  fscanf(ifp, "%lf%lf",
		 &varinfo[i].partic[0][k].low,
		 &varinfo[i].partic[0][k].up);

	}
      else
	error_exit("Error in varinfo file.\n");
    }

  fclose(ifp);

  return 0;
}


/*!
 * In this function, the vertex[] array and the limits[] structure
 * are generated. This function is called in the makestartvertex()
 * function in simannealing process, as e.g.
 * store_vertex_value(tempvertex, limits, ndim, particles[i].partic[j].G,
 * varinfo[i].partic[j][0]);
 * If the parameter to which "value" corresponds is a free one,
 * "value" is stored in the component vertex[ndim]; the limits
 * array is filled in this function as well.
 */
int store_vertex_value(double vertex[], Limits limits[], int* ndim,
		       double value, Celinfo partic_info)
{
  if(partic_info.var)
    {
      vertex[*ndim] = value;

      // if variable varies within boundaries
      if(partic_info.bound)
	{
	  limits[*ndim].bound = 1;
	  limits[*ndim].low = partic_info.low;
	  limits[*ndim].up = partic_info.up;

	  // check if initial value fits in the boundaries
	  if(value < limits[*ndim].low || value > limits[*ndim].up)
	    {
	      std::stringstream ss;
	      ss << "Value " << *ndim << " out of boundaries:" << value << " "
		 << limits[*ndim].low << " " << limits[*ndim].up << "\n";

	      throw std::logic_error(ss.str());
	    }
	}

      // unbound variable
      else
      {
	limits[*ndim].bound = 0;
	limits[*ndim].low = 0;
	limits[*ndim].up = 0;
      }

      // Increase number of free parameters
      (*ndim)++;
    }

  return 0;
}



/*!
 * This function calculates the chi squared.
 */
double chifunc(double* vertex,Class** particles,
	       Varinfo varinfo[], Observable* observ,
	       Data** datapoints, int datacount[])
{
  int i, j, iso, modifdatacount=0;
  double chisquared=0.0;

  Class* workparticles = NULL;
  Class* workparticles_2 = NULL;

  modifdatacount = 0;
  chisquared = 0;

  /*
   * If chi squared has to be calculated over different isospin
   * channels, this is done first on the base set (represented
   * by the vertex[] array, this is isospin channel 1 or 2) and
   * afterwards over the other channels.
   */

  for(j = 0; j < observ->iso.nr_iso_channels; j++)
    {

      /* Determine the isospin channel under investigation */

      iso = observ->iso.iso_channel[j];
      observ->iso.isospin = iso;

      workparticles = (Class*) calloc(CLASSMAX, sizeof(Class));
      if(workparticles == NULL)
	error_exit("Error during allocation of workparticles\n");

      workparticles_2 = (Class*) calloc(CLASSMAX, sizeof(Class));
      if(workparticles_2 == NULL)
	error_exit("Error during allocation of workparticles\n");


      /*
       * The vertex[] array is coupled in the workparticles[] structure.
       * In the case of a simple chi square calculation, vertex == NULL,
       * particles[] is copied in workparticles[].
       */

      if(vertex != NULL)
	insertvertex(workparticles, particles[iso], varinfo, observ, vertex);
      else
	copyparticles(workparticles, particles[iso],observ);

      /* for now, workparticles = workparticles_2 */
      copyparticles(workparticles_2, workparticles,observ);


      /*
       * Transform coupling constants to correct isospin channel
       */
      if (observ->regge && observ->reg.t_and_u_channel)
	{
	  observ->reg.t_channel = 1;
	  change_coupling_constants(workparticles, observ,observ->iso.isospin);

	  observ->reg.t_channel = 0;
	  change_coupling_constants(workparticles_2, observ,observ->iso.isospin);
	}
      else
	change_coupling_constants(workparticles, observ,observ->iso.isospin);

      //fprintparticspec(stdout,workparticles, *observ);


      modifdatacount += datacount[iso]; /* this variable will include the
					 * weight of the polarization
					 * observables */

      for(i=0; i<datacount[iso]; i++)
	{
	  /*
	   * Here, the selection for a photoproduction or an electroproduction
	   * data point is made.
	   */

	  if(datapoints[iso][i].photo_prod)
	    {
	      chisquared += photo_chi(datapoints[iso][i].photo, workparticles, workparticles_2,
				      observ,&modifdatacount);
	    }
	  else if(datapoints[iso][i].electro_prod)
	    {
	      chisquared += electro_chi(datapoints[iso][i].elec, workparticles,
					workparticles_2,
					observ,&modifdatacount);
	    }
	  else if(datapoints[iso][i].kaoncapture)
	    {
	      /*
	      printf("Do not trust chi-squared calculations with rad.cap. data\n");
	      */
	      chisquared += kaoncap_chi(datapoints[iso][i].kaoncap,workparticles,
					workparticles_2, observ,&modifdatacount);
	    }

	} // end loop over datapoints



      /* Clean up the "workparticles", if necessary. */

      free(workparticles);
      free(workparticles_2);

    } // end loop over isospin channels



  /* Return chi square per degree of freedom (= number of data points) */
  return chisquared / (modifdatacount*1.0);



}

/*!
 * The vertex with the modified coupling constants is together with
 * the unchanged parameters of the particles[] list copied into
 * the workparticles[] structure.
 */
int insertvertex(Class workparticles[], Class particles[], Varinfo varinfo[],
		 Observable* observ, double* vertex)
{
  int i, j, count = 0, born_cc_used = 0, new_cutoff = 0;


  /* ***********************
   * Copy coupling constants
   * from - particles[]
   *      - vertex[]
   * *********************** */

  // loop over diagrams
  for(i=0; i< CLASSMAX; i++)
    {
      // loop over particles of same CLASS
      for(j=0; j<particles[i].particount; j++)
	{

	  // Born channels
	  if(i == 0 || i == 1 || i == 2)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		{
		  // Count does not increase, because if g_NKY is fitted
		  // is needs to be the same for all 3 born channels
		  workparticles[i].partic[j].H = vertex[count];
		  born_cc_used = 1;
		}
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].G = particles[i].partic[j].G;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }

	  // A-channel
	  else if(i == 3)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].G = particles[i].partic[j].G;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }

	  // B- and C-channel
	  else if(i == 4 || i == 5)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].I = vertex[count++];
	      else
		workparticles[i].partic[j].I = particles[i].partic[j].I;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].G = particles[i].partic[j].G;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }

	  // D,E,F,G-channel
	  else if(i >= 6 && i <=9)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].G = particles[i].partic[j].G;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }

	  // H,I,J,L-channel
	  else if(i >= 10 && i<= 13)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].G = vertex[count++];
	      else
		workparticles[i].partic[j].G = particles[i].partic[j].G;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      if(varinfo[i].partic[j][2].var)
		workparticles[i].partic[j].X = vertex[count++];
	      else
		workparticles[i].partic[j].X = particles[i].partic[j].X;

	      if(varinfo[i].partic[j][3].var)
		workparticles[i].partic[j].Y = vertex[count++];
	      else
		workparticles[i].partic[j].Y = particles[i].partic[j].Y;

	      if(varinfo[i].partic[j][4].var)
		workparticles[i].partic[j].Z = vertex[count++];
	      else
		workparticles[i].partic[j].Z = particles[i].partic[j].Z;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].I = particles[i].partic[j].I;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }

	  // M,N,O,Q-channel
	  else if(i >= 14 && i<= 17)
	    {
	      // Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].G = vertex[count++];
	      else
		workparticles[i].partic[j].G = particles[i].partic[j].G;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      if(varinfo[i].partic[j][2].var)
		workparticles[i].partic[j].X = vertex[count++];
	      else
		workparticles[i].partic[j].X = particles[i].partic[j].X;

	      if(varinfo[i].partic[j][3].var)
		workparticles[i].partic[j].Y = vertex[count++];
	      else
		workparticles[i].partic[j].Y = particles[i].partic[j].Y;

	      if(varinfo[i].partic[j][4].var)
		workparticles[i].partic[j].Z = vertex[count++];
	      else
		workparticles[i].partic[j].Z = particles[i].partic[j].Z;

	      // Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].I = particles[i].partic[j].I;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }
	  // R-channel
	  else if(i == 18)
	    {
	      //Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].I = vertex[count++];
	      else
		workparticles[i].partic[j].I = particles[i].partic[j].I;

	      if(varinfo[i].partic[j][2].var)
		workparticles[i].partic[j].X = vertex[count++];
	      else
		workparticles[i].partic[j].X = particles[i].partic[j].X;

	      if(varinfo[i].partic[j][3].var)
		workparticles[i].partic[j].Y = vertex[count++];
	      else
		workparticles[i].partic[j].Y = particles[i].partic[j].Y;

	      if(varinfo[i].partic[j][4].var)
		workparticles[i].partic[j].Z = vertex[count++];
	      else
		workparticles[i].partic[j].Z = particles[i].partic[j].Z;

	      //Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].G = particles[i].partic[j].G;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }
	  // V-channel
	  else if(i == 19)
	    {
	      //Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].G = vertex[count++];
	      else
		workparticles[i].partic[j].G = particles[i].partic[j].G;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      //Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].I = particles[i].partic[j].I;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }
	  // W-channel
	  else if(i == 20)
	    {
	      //Free parameters
	      if(varinfo[i].partic[j][0].var)
		workparticles[i].partic[j].G = vertex[count++];
	      else
		workparticles[i].partic[j].G = particles[i].partic[j].G;

	      if(varinfo[i].partic[j][1].var)
		workparticles[i].partic[j].H = vertex[count++];
	      else
		workparticles[i].partic[j].H = particles[i].partic[j].H;

	      //Fixed parameters
	      workparticles[i].partic[j].E = particles[i].partic[j].E;
	      workparticles[i].partic[j].I = particles[i].partic[j].I;
	      workparticles[i].partic[j].gic = particles[i].partic[j].gic;
	    }


	} // end particle loop


      /*
       * Make sure that the free parameters for the Born terms, if they occur,
       * are  nicely stored in the first 3 channels, if they occur.
       */
       if(i >= 2 && born_cc_used)
	{
	  count++;
	  born_cc_used = 0;
	}

    } // end loop over diagrams


  /* ***********************************
   * Copy form factor cutoff information
   * *********************************** */

  /* Form factor info is stored in the (CLASSMAX)th cel */
  if(observ->hadronformfac)
    {
      i = CLASSMAX;
      j = 0;

      // Update strong form factor cutoffs
      if(varinfo[i].partic[j][0].var)
	{
	  observ->ffac.born_cutoff = vertex[count++];
	  new_cutoff = 1;
	}

      if(varinfo[i].partic[j][1].var)
	{
	  observ->ffac.res_cutoff = vertex[count++];
	  new_cutoff = 1;
	}

      /* Update FormFactor objects in particles[],
       * if necessary */
      if(new_cutoff)
	{
	  update_strong_formfactor(observ,particles);
	  new_cutoff = 0;
	}

    }


  /* *************************
   * Copy all the obvious info
   * ************************* */

  // loop over all diagrams
  for(i=0; i< CLASSMAX; i++)
    {
      // nr. of particles
      workparticles[i].particount = particles[i].particount;

      // loop over all particles
      for(j=0; j<particles[i].particount; j++)
	{

	  // name of particle
	  strcpy(workparticles[i].partic[j].nickname,
		 particles[i].partic[j].nickname);

	  // properties
	  workparticles[i].partic[j].mass = particles[i].partic[j].mass;
	  workparticles[i].partic[j].width = particles[i].partic[j].width;
	  workparticles[i].partic[j].spin = particles[i].partic[j].spin;

	  // EM transition coefficients
	  workparticles[i].partic[j].r_kappa_n_p =
	    particles[i].partic[j].r_kappa_n_p;
	  workparticles[i].partic[j].r_kappa_1_n_p =
	    particles[i].partic[j].r_kappa_1_n_p;
	  workparticles[i].partic[j].r_kappa_2_n_p =
	    particles[i].partic[j].r_kappa_2_n_p;

	  // gauge-invariant couplings
	  workparticles[i].partic[j].gic =
	    particles[i].partic[j].gic;

	  // regge phase
	  workparticles[i].partic[j].regge_phase =
	    particles[i].partic[j].regge_phase;
	  workparticles[i].partic[j].regge_phase_nonbase =
	    particles[i].partic[j].regge_phase_nonbase;

	  // longitudinal coupling for resonances?
	  workparticles[i].partic[j].long_coupling =
	    particles[i].partic[j].long_coupling;

	  // formfactor info
	  workparticles[i].partic[j].formfactorE =
	    particles[i].partic[j].formfactorE;
	  workparticles[i].partic[j].formfactorG =
	    particles[i].partic[j].formfactorG;
	  workparticles[i].partic[j].formfactorH =
	    particles[i].partic[j].formfactorH;
	  workparticles[i].partic[j].formfactorI =
	    particles[i].partic[j].formfactorI;
	}

      // Make sure that the masses/widths/charges of the three born terms are copied!
      if(i <= 2 &&  particles[i].particount == 0)
	{
	  strcpy(workparticles[i].partic[0].nickname,
		 particles[i].partic[0].nickname);
	  workparticles[i].partic[0].mass = particles[i].partic[0].mass;
	  workparticles[i].partic[0].width = particles[i].partic[0].width;
	  workparticles[i].partic[0].E = particles[i].partic[0].E;
	}


      /* If duality corrections to the Regge propagator are included,
       * make sure that the s-channel trajectory properties are copied.*/
      /* Obsolete:
	if(observ->reg.dual_corr)
	{
	  workparticles[i].traject.mass = particles[i].traject.mass;
	  workparticles[i].traject.sign = particles[i].traject.sign;
	  workparticles[i].traject.rel_sign = particles[i].traject.rel_sign;
	  workparticles[i].traject.r_slope = particles[i].traject.r_slope;
	  workparticles[i].traject.i_slope = particles[i].traject.i_slope;
	  workparticles[i].traject.i_intercept =
	    particles[i].traject.i_intercept;
	  if(i == 0)
	    workparticles[i].trajectcount = particles[i].trajectcount;
	}
       */

    } // end diagram loop

  return 0;
}




int copyparticles(Class workparticles[], Class particles[], Observable* observ)
{


  int i,j;


  for(i=0; i< CLASSMAX; i++)
    {
      workparticles[i].particount = particles[i].particount;

      for(j=0; j< PARTICLEMAX; j++)
	{
	  // name of particle
	  strcpy(workparticles[i].partic[j].nickname,
		 particles[i].partic[j].nickname);

	  // properties
	  workparticles[i].partic[j].mass = particles[i].partic[j].mass;
	  workparticles[i].partic[j].width = particles[i].partic[j].width;
	  workparticles[i].partic[j].spin = particles[i].partic[j].spin;

	  // coupling constants
	  workparticles[i].partic[j].E = particles[i].partic[j].E;
	  workparticles[i].partic[j].G = particles[i].partic[j].G;
	  workparticles[i].partic[j].H = particles[i].partic[j].H;
	  workparticles[i].partic[j].I = particles[i].partic[j].I;
	  workparticles[i].partic[j].J = particles[i].partic[j].J;
	  workparticles[i].partic[j].X = particles[i].partic[j].X;
	  workparticles[i].partic[j].Y = particles[i].partic[j].Y;
	  workparticles[i].partic[j].Z = particles[i].partic[j].Z;

	  // EM transition coefficients
	  workparticles[i].partic[j].r_kappa_n_p =
	    particles[i].partic[j].r_kappa_n_p;
	  workparticles[i].partic[j].r_kappa_1_n_p =
	    particles[i].partic[j].r_kappa_1_n_p;
	  workparticles[i].partic[j].r_kappa_2_n_p =
	    particles[i].partic[j].r_kappa_2_n_p;

	  // gauge-invariant couplings
	  workparticles[i].partic[j].gic =
	    particles[i].partic[j].gic;

	  // regge phase
	  workparticles[i].partic[j].regge_phase =
	    particles[i].partic[j].regge_phase;
	  workparticles[i].partic[j].regge_phase_nonbase =
	    particles[i].partic[j].regge_phase_nonbase;

	  // longitudinal coupling for resonances?
	  workparticles[i].partic[j].long_coupling =
	    particles[i].partic[j].long_coupling;

	  // formfactor info
	  workparticles[i].partic[j].formfactorE =
	    particles[i].partic[j].formfactorE;
	  workparticles[i].partic[j].formfactorG =
	    particles[i].partic[j].formfactorG;
	  workparticles[i].partic[j].formfactorH =
	    particles[i].partic[j].formfactorH;
	  workparticles[i].partic[j].formfactorI =
	    particles[i].partic[j].formfactorI;
	}

      /*
       * Make sure that the masses of the three born terms are
       * copied!
       */

      if(i <= 2 &&  particles[i].particount == 0)
	workparticles[i].partic[0].mass = particles[i].partic[0].mass;



      /*
       * If duality corrections to the Regge propagator are included,
       * make sure that the s-channel trajectory properties are copied.
       */


      if(observ->reg.dual_corr)
	{
	  workparticles[i].traject.mass = particles[i].traject.mass;
	  workparticles[i].traject.sign = particles[i].traject.sign;
	  workparticles[i].traject.rel_sign = particles[i].traject.rel_sign;
	  workparticles[i].traject.r_slope = particles[i].traject.r_slope;
	  workparticles[i].traject.i_slope = particles[i].traject.i_slope;
	  workparticles[i].traject.i_intercept =
	    particles[i].traject.i_intercept;

	  if(i == 0)
	    workparticles[i].trajectcount = particles[i].trajectcount;
	}

    }

  return 0;
}

/*!
 * This function returns the chi^2 if "datapoint" is a photoproduction point.
 */
double photo_chi(Photo datapoint, Class* workparticles, Class* workparticles_2,
		 Observable* observ, int* modifdatacount)
{
  int weight;
  double energy, w, pk;
  double ampli=0.0, tmpampli=0.0, sumampli=0.0, chisquared=0.0;
  double gridnum, realgridnum=0.0, e_delta;
  double mp, mk, my;
  double costhkcm;

  /* Allocate the Born masses */
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  gridnum = get_gridnum(observ->fit.narrow_grid, datapoint.observable);

  /* If there is only one energy, emax was put = emin in the
   * getdatastructure() function. In which case we need
   * to omit the for-loop construction. */
  if(datapoint.emax >  datapoint.emin)  // Determine the step size
    e_delta = (datapoint.emax - datapoint.emin) / gridnum; //#calculations = gridnum+1
  else {  // no averaging will be done!
    e_delta = 10.;
  }

  // Loop over energy inside bin
  for(energy = datapoint.emin;
      energy < datapoint.emax + 2.0; energy += e_delta)
    {
      // Conversion from lab to c.m. energy (necessary in calculations)
      if(datapoint.is_w)
        w = (energy*energy - mp*mp)/(2.*energy);
      else
        lorentztrans_photo(energy, &w, mp);

      // Calculation of the (three vector) momentum of the kaon.
      pk = construct_pk(w, w, mp, mk, my);

      if(datapoint.t_ang)
        costhkcm = (datapoint.t - mk*mk + 2*sqrt(pk*pk + mk*mk)*w)/(2.*pk*w);
      else
        if(datapoint.is_cos_bin)
          costhkcm = (datapoint.cos + datapoint.cos_max)/2.;
        else
          costhkcm = datapoint.cos;

      // Test whether pk and datapoint.cos are physical
      if(pk > 0);
      else if(fabs(costhkcm) <= 1);
      else if(compareDoubleUFlow(fabs(costhkcm), 1))
        costhkcm = costhkcm > 0 ? 1 : -1;
      else
        throw std::logic_error("A datapoint of type '" + string(datapoint.observable) +
                               "' lies outside the physical region!\n\n");

      // Calculate the appropriate observable
      tmpampli = calculate_photo_observable(&datapoint,observ,workparticles,
                                            workparticles_2,w,w,costhkcm,pk);

      // Add observable to bin-total
      sumampli += tmpampli;

      // Count gridpoints (realgridnum can be different from gridnum!)
      realgridnum += 1;

      if (datapoint.label>=0)  datapoint.label++; // increment it either way, because we have already taken into consideration gridnum points

    } // end loop in energy-bin

  // Average observable
  ampli = sumampli / realgridnum;


  // General check for unphysical behavior. !
  if(ampli < 0.0 && (!strcmp(datapoint.observable, "diffcs") ||
		     !strcmp(datapoint.observable, "totcs")))
    {
      fprintparticspec(stderr, workparticles, *observ);
      std::stringstream ss;

      ss << "Error: neg. modulus !!\n(energy = " << datapoint.emin
         << "; costheta = " << datapoint.cos << ")\n";

      throw std::logic_error(ss.str());
    }

  weight = 1;

  /* weight procedure according to
     H.Thom, Phys.Rev. 151, 1322 (1966). */

  // Polarization observables
  if( !strcmp ( datapoint.observable, "P" ) ||
  !strcmp ( datapoint.observable, "S" ) ||
  !strcmp ( datapoint.observable, "T" ) ||
  !strcmp ( datapoint.observable, "C_x" ) ||
      !strcmp ( datapoint.observable, "C_xp") ||
      !strcmp ( datapoint.observable, "C_z" ) ||
      !strcmp ( datapoint.observable, "C_zp") ||
      !strcmp ( datapoint.observable, "O_x" ) ||
      !strcmp ( datapoint.observable, "O_xp") ||
      !strcmp ( datapoint.observable, "O_z" ) ||
      !strcmp ( datapoint.observable, "O_zp") )
    {
      if(observ->fit.polweight != 0)
	weight *= observ->fit.polweight;
      else
        error_exit("error in polarization weight!!\n");
    }

  // Multiple isospin channels
  if(observ->iso.nr_iso_channels > 1)
    weight *= observ->fit.isoWeight[datapoint.iso];

  // Calculate the chi^2:
  chisquared = weight* pow(((ampli - datapoint.ampli) / datapoint.error), 2);

  /* According to the weight procedure, the total number of degrees of freedom
   * has to increase. */
  (*modifdatacount) += weight - 1;

  observ->photoprod = 0;

  return chisquared;
}

/*!
 * Calculates one photoproduction observable in one kinematic point
 * based on the the info in the Data-struct and the Observable-struct.
 *
 * We can calculate different observables:
 * \verbatim
 datapoint->observable | Actual observable
 --------------------- | ------------------------------------------
 dcs                   | differential cross section
 tcs                   | total cross section
 rec                   | recoil polarization (P) (along y')
 pho                   | photon polarization (Sigma)
 tar                   | target polarization (T) (along y)
 c_x                   | beam(circ)-recoil polarization (along x)
 c_xp                  | beam(circ)-recoil polarization (along x')
 c_z                   | beam(circ)-recoil polarization (along z)
 c_zp                  | beam(circ)-recoil polarization (along z')
 o_x                   | beam(lin)-recoil polarization (along x)
 o_xp                  | beam(lin)-recoil polarization (along x')
 o_z                   | beam(lin)-recoil polarization (along z)
 o_zp                  | beam(lin)-recoil polarization (along z')
 b1^2                  | transversity amplitudes
 b2^2                  |        with - linear photon polarization
 b3^2                  |             - hyperon and nucleon polarized along y axis
 b4^2                  |
 b1_R                  |
 b2_R                  |          ^2: norm
 b3_R                  |          _R: real part
 b4_R                  |          _I: imaginary part
 b1_I                  |
 b2_I                  |
 b3_I                  |
 b4_I                  |
 H1^2                  | helicity amplitudes
 H2^2                  |
 H3^2                  |          ^2: norm
 H4^2                  |          _R: real part
 H1_R                  |          _I: imaginary part
 H2_R                  |
 H3_R                  |
 H4_R                  |
 H1_I                  |
 H2_I                  |
 H3_I                  |
 H4_I                  |
 \endverbatim
 *
 * Transversity amplitudes implementation:
 * \verbatim
 datapoint->observable | Actual observable
 --------------------- | ------------------
 diffcs                | differential cross section
 totcs                 | total cross section
 S                     | beam asymmetry
 T                     | target asymmetry
 P                     | recoil asymmetry
 C_x                   | beam(c)-recoil(x) asymmetry
 C_xp                  | beam(c)-recoil(x') asymmetry
 C_z                   | beam(c)-recoil(z) asymmetry
 C_zp                  | beam(c)-recoil(z') asymmetry
 O_x                   | beam(l)-recoil(x) asymmetry
 O_xp                  | beam(l)-recoil(x') asymmetry
 O_z                   | beam(l)-recoil(z) asymmetry
 O_zp                  | beam(l)-recoil(z') asymmetry
 E                     | beam(c)-target(z) asymmetry
 F                     | beam(c)-target(x) asymmetry
 G                     | beam(l)-target(z) asymmetry
 H                     | beam(l)-target(x) asymmetry
 T_x                   | target(x)-recoil(x) asymmetry
 T_xp                  | target(x)-recoil(x') asymmetry
 T_z                   | target(x)-recoil(z) asymmetry
 T_zp                  | target(x)-recoil(z') asymmetry
 L_x                   | target(z)-recoil(x) asymmetry
 L_xp                  | target(z)-recoil(x') asymmetry
 L_z                   | target(z)-recoil(z) asymmetry
 L_zp                  | target(z)-recoil(z') asymmetry
 \endverbatim
 *
 * Note: The z' axis is chosen along the kaon momentum (drechsel convention).
 */
double calculate_photo_observable(Photo *datapoint, Observable *observ, Class *particles,
				  Class *particles_2, double w, double k, double costhk, double pk)
{
  double amplitude = 0.0;
  double asym = 0.0;
  double asym_x=0.0;
  double asym_z= 0.0;
  double unpol= 0.0;
  double mp =0.0;
  double mk =0.0;
  double my =0.0;
  int iang=0;
  double kincoeff, spinaver;
  double a_delta =0.0;
  TMatrixElement* matrixelement;

  /* Allocate the Born masses */
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  // Set the observable specification for photoproduction !!
  observ->photoprod = 1;
  observ->electroprod = 0;
  observ->kaoncapture = 0;

  // Coeff for conversion from MeV^-2 to mubarn
  double convertcoeff = 197.32 * 197.32 * 1e4;


  /* Calculate observable
   * -------------------- */

  // Fit to differential cross section data
  // --------------------------------------
  if(!strcmp(datapoint->observable, "dcs"))
    {
      observ->photo.pol.nopol = 1;          // no polarization
      spinaver = spinaveraging(observ);     // spinaveraging
      kincoeff = kincoeff_photo(w, pk, mp); // kinematic coefficient

      // Determine amplitude
      if (!observ->reg.t_and_u_channel)
	amplitude = determine_M2_photo(w, w,costhk, pk,particles, observ, datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  amplitude = determine_M2_photo(w, w, costhk, pk, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  amplitude = determine_M2_photo(w, w, costhk, pk, particles_2, observ, datapoint->label);
	}
      }

      amplitude *= kincoeff * convertcoeff * spinaver;

      // Conversion from dsigma / dOmega  --> dsigma / dt(u)
      if(datapoint->ds_dt || datapoint->ds_du) {
	// dOmega / dt(u) = 2*pi / (2*pk*k)
	// Factor 1e-6 since d sigma / dt(u) in (mub * GeV^-2).
	amplitude *= PI / ( pk * w * 1e-6 );
      }

      observ->photo.pol.nopol = 0;
    }

  // Fit to total cross section
  // --------------------------
  else if(!strcmp(datapoint->observable, "tcs"))
    {
      observ->photo.pol.nopol = 1;          // no polarization
      spinaver = spinaveraging(observ);     // spinaveraging
      kincoeff = kincoeff_photo(w, pk, mp); // kinematic coefficient

      amplitude = 0; // reset amplitude

      // Set the stepsize for the integration
      observ->photo.kin.anglestep = 3;

      // Loop over angle
      for(iang=0; iang <= 180; iang += observ->photo.kin.anglestep)// #calculations = fl(180/anglestep)+1 = 61
	{
	  if(iang == 0 || iang >= 180) // fully fwd and bckwd bin have only half the width
	    a_delta = observ->photo.kin.anglestep / 2.0;
	  else
	    a_delta = observ->photo.kin.anglestep;


	  if (!observ->reg.t_and_u_channel)
	    amplitude += determine_M2_photo(w, w, cos(iang*PI/180.0), pk,
					    particles, observ, datapoint->label)
	      * (2 * PI) * sin(iang*PI/180.0) * a_delta * PI/180.0;

	  else {
	    if (cos(iang*PI/180.0) > 0.0 || datapoint->tch) {
	      observ->reg.t_channel = 1;
	      amplitude += determine_M2_photo(w, w, cos(iang*PI/180.0), pk,
					      particles, observ, datapoint->label)
		* (2 * PI) * sin(iang*PI/180.0) * a_delta * PI/180.0;
	    }
	    else {
	      observ->reg.t_channel = 0;
	      amplitude += determine_M2_photo(w, w, cos(iang*PI/180.0), pk,
					      particles_2, observ, datapoint->label)
		* (2 * PI) * sin(iang*PI/180.0) * a_delta * PI/180.0;
	    }
	  }
	  if (datapoint->label>=0) datapoint->label++;
	} // end loop over angle

      amplitude *= kincoeff * convertcoeff * spinaver;

      observ->photo.pol.nopol = 0;
    }

  // Fit to recoil polarization data
  // -------------------------------
  else if(!strcmp(datapoint->observable,"rec"))
    {
	observ->photo.pol.spol = 1;          // single polarization
	observ->photo.pol.sinpol.recpol = 1; // recoil polarization

	// Construct matrixelement for recoil polarization
	if (!observ->reg.t_and_u_channel)
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	else {
	  if (costhk > 0.0 || datapoint->tch) {
	    observ->reg.t_channel = 1;
	    matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	  }
	  else {
	    observ->reg.t_channel = 0;
	    matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	  }
	}

	// Calculate <g_5>
	spinaver =  spinaveraging(observ); // spinaveraging
	asym = determine_M2_photo_withME(matrixelement,observ);
	asym *= spinaver;   //  *kincoeff*convertcoeff;

	// Reset polarizations
	observ->photo.pol.spol = 0;
	observ->photo.pol.sinpol.recpol = 0;
	observ->photo.pol.nopol = 1;         // no polarization

	// Reconstruct matrixelement for unpolarized case
	updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

	// Calculate unpolarized cross section
	spinaver =  spinaveraging(observ);
	unpol = determine_M2_photo_withME(matrixelement,observ);
	unpol *= spinaver; //   *kincoeff*convertcoeff;

	// Determine asymmetry
	amplitude += 2*asym / unpol;

	// Reset everything
	observ->photo.pol.nopol = 0;
	release_matrixelement(matrixelement);

    }

  // Fit to photon polarization data
  // -------------------------------
  else if(!strcmp(datapoint->observable,"pho"))
    {
      observ->photo.pol.spol = 1;          // single polarization
      observ->photo.pol.sinpol.phopol = 1; // photon polarization

      // Construct matrixelement for photon polarization
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5>
      spinaver = spinaveraging(observ);
      asym = determine_M2_photo_withME(matrixelement, observ);
      asym *= spinaver; // *kincoeff*convertcoeff

      // Reset polarizations
      observ->photo.pol.spol = 0;
      observ->photo.pol.sinpol.phopol = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));
      // Calculate unpolarized cross section
      spinaver = spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement, observ);
      unpol *= spinaver; // *kincoeff*convertcoeff

      // Determine asymmetry
      amplitude = asym / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to target polarization data
  // -------------------------------
  else if(!strcmp(datapoint->observable, "tar"))
    {
      observ->photo.pol.spol = 1;          // single polarization
      observ->photo.pol.sinpol.tarpol = 1; // target polarization

      // Construct matrixelement for target polarization
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5>
      spinaver =  spinaveraging(observ);
      asym = determine_M2_photo_withME(matrixelement,observ);
      asym *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.spol = 0;
      observ->photo.pol.sinpol.tarpol = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry
      amplitude = asym / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along x)
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "c_x"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.circbeam = 1; // circular photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Set polarization to hyperon along z'
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.doubpol.z_barvec = 1;

      // Reconstruct matrixelement for hyperon polarized along z'
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.circbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry C_x = C_x' costh + C_z' sinth
      amplitude =
	2 * (asym_x * costhk + asym_z * sqrt(1.0 - costhk*costhk)) / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along z)
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "c_z"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.circbeam = 1; // circular photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Set polarization to hyperon along z'
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.doubpol.z_barvec = 1;

      // Reconstruct matrixelement for hyperon polarized along z'
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.circbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry C_z = -C_x' sinth + C_z' costh
      amplitude =
	2 * (-1.0 * asym_x * sqrt(1.0 - costhk*costhk) + asym_z * costhk) / unpol;

      asym_x *= kincoeff_photo(w, pk, mp)*convertcoeff;
      asym_z *= kincoeff_photo(w, pk, mp)*convertcoeff;
      unpol *= kincoeff_photo(w, pk, mp)*convertcoeff;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along x')
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "c_xp"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.circbeam = 1; // circular photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.circbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry
      amplitude = 2 * asym_x  / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along z')
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "c_zp"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.circbeam = 1; // circular photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.z_barvec = 1; // hyperon pol. along z'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along z'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.circbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.z_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry
      amplitude = 2 * asym_z  / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along x)
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "o_x"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.linbeam = 1;  // linear photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Set polarization to hyperon along z'
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.doubpol.z_barvec = 1;

      // Reconstruct matrixelement for hyperon polarized along z'
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.linbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry O_x = O_x' costh + O_z' sinth
      amplitude =
	2 * (asym_x * costhk + asym_z * sqrt(1.0 - costhk*costhk)) / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along z)
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "o_z"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.linbeam = 1;  // linear photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Set polarization to hyperon along z'
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.doubpol.z_barvec = 1;

      // Reconstruct matrixelement for hyperon polarized along z'
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.linbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry O_z = -O_x' sinth + O_z' costh
      amplitude =
	2 * (-1.0 * asym_x * sqrt(1.0 - costhk*costhk) + asym_z * costhk) / unpol;

      asym_x *= kincoeff_photo(w, pk, mp)*convertcoeff;
      asym_z *= kincoeff_photo(w, pk, mp)*convertcoeff;
      unpol *= kincoeff_photo(w, pk, mp)*convertcoeff;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along x')
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "o_xp"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.linbeam = 1;  // linear photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.x_barvec = 1; // hyperon pol. along x'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along x'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^x+>
      spinaver =  spinaveraging(observ);
      asym_x = determine_M2_photo_withME(matrixelement,observ);
      asym_x *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.linbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.x_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry
      amplitude = 2 * asym_x  / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Fit to beam-recoil polarization data (along z')
  // -----------------------------------------------
  else if(!strcmp(datapoint->observable, "o_zp"))
    {
      observ->photo.pol.dpol = 1;             // double polarization
      observ->photo.pol.doubpol.beamrec = 1;  // beam-recoil pol.
      observ->photo.pol.doubpol.linbeam = 1;  // linear photons
      observ->photo.pol.doubpol.beamhel = 1;  // positive beam helicity
      observ->photo.pol.doubpol.z_barvec = 1; // hyperon pol. along z'
      observ->quant_axis_drechsel = 1;        // z' along kaon momentum
      observ->quant_axis_thesis = 0;

      // Construct matrixelement for beam-recoil polarization along z'
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      // Calculate <g_5^z+>
      spinaver =  spinaveraging(observ);
      asym_z = determine_M2_photo_withME(matrixelement,observ);
      asym_z *= spinaver;   // *kincoeff*convertcoeff;

      // Reset polarizations
      observ->photo.pol.dpol = 0;
      observ->photo.pol.doubpol.beamrec = 0;
      observ->photo.pol.doubpol.linbeam = 0;
      observ->photo.pol.doubpol.beamhel = 0;
      observ->photo.pol.doubpol.z_barvec = 0;
      observ->photo.pol.nopol = 1;         // no polarization

      // Reconstruct matrixelement for unpolarized case
      updateSpinDependencies_matrixelement(matrixelement, &(datapoint->label));

      // Calculate unpolarized cross section
      spinaver =  spinaveraging(observ);
      unpol = determine_M2_photo_withME(matrixelement,observ);
      unpol *= spinaver; // *kincoeff*convertcoeff;

      // Determine asymmetry
      amplitude = 2 * asym_z  / unpol;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Transversity amplitudes b_i and a_i (normalized)
  // -----------------------------------------------------
  else if((datapoint->observable[0] == 'b' || datapoint->observable[0] == 'a')
	  && string(datapoint->observable).length() == 4)
    {
      int index = ((int)datapoint->observable[1]) - 48;

      if(index >= 1 && index <= 4)
	{
	  int option;

	  if(datapoint->observable[2] == '^' && datapoint->observable[3] == '2')
	    option = 0;
	  else if(datapoint->observable[2] == '_' && datapoint->observable[3] == 'R')
	    option = 1;
	  else if(datapoint->observable[2] == '_' && datapoint->observable[3] == 'I')
	    option = 2;
	  else
	    goto error;

	  observ->photo.pol.nopol = 1; // no polarization

	  // Determine matrixelement
	  if (!observ->reg.t_and_u_channel)
	    matrixelement =
	      construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	  else {
	    if (costhk > 0.0 || datapoint->tch) {
	      observ->reg.t_channel = 1;
	      matrixelement =
		construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	    }
	    else {
	      observ->reg.t_channel = 0;
	      matrixelement =
		construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	    }
	  }

	  double amp[2];

	  if(datapoint->observable[0] == 'b')
	    get_transversity_amplitude(matrixelement, amp, index);
	  else
	    {

	      double a[4][2], norm = 0.;

	      for(int i = 0; i < 4; i++)
		{
		  get_transversity_amplitude(matrixelement, a[i], i + 1);
		  norm += (pow(a[i][0], 2) + pow(a[i][1], 2));
		}

	      norm = sqrt(norm);

	      amp[0] = a[index - 1][0]/norm;
	      amp[1] = a[index - 1][1]/norm;
	    }

	  // Determine amplitude
	  if(option == 0)
	    amplitude = amp[0]*amp[0] + amp[1]*amp[1];
	  else if(option == 1)
	    amplitude = amp[0];
	  else if(option == 2)
	    amplitude = amp[1];

	  // Reset everything
	  observ->photo.pol.nopol = 0;
	  release_matrixelement(matrixelement);
	}
      else
	goto error;
    }

  // Helicity amplitudes H_i
  // -----------------------------------------------------
  else if((datapoint->observable[0] == 'H' || datapoint->observable[0] == 'h')
	  && string(datapoint->observable).length() == 4)
    {
      int index = ((int)datapoint->observable[1]) - 48;

      if(index >= 1 && index <= 4)
	{
	  int option;

	  if(datapoint->observable[2] == '^' && datapoint->observable[3] == '2')
	    option = 0;
	  else if(datapoint->observable[2] == '_' && datapoint->observable[3] == 'R')
	    option = 1;
	  else if(datapoint->observable[2] == '_' && datapoint->observable[3] == 'I')
	    option = 2;
	  else
	    goto error;

	  observ->photo.pol.dpol = 1;
	  observ->photo.pol.doubpol.z_barvec = 1;
	  observ->quant_axis_thesis = 1;
	  observ->photo.pol.doubpol.circbeam = 1;

	  // Determine matrixelement
	  if (!observ->reg.t_and_u_channel)
	    matrixelement =
	      construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	  else {
	    if (costhk > 0.0 || datapoint->tch) {
	      observ->reg.t_channel = 1;
	      matrixelement =
		construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	    }
	    else {
	      observ->reg.t_channel = 0;
	      matrixelement =
		construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	    }
	  }

	  double amp[2];

	  if(datapoint->observable[0] == 'H')
	    get_helicity_amplitude(matrixelement, amp, index);
	  else
	    {

	      double h[4][2], norm = 0.;

	      for(int i = 0; i < 4; i++)
		{
		  get_helicity_amplitude(matrixelement, h[i], i + 1);
		  norm += (pow(h[i][0], 2) + pow(h[i][1], 2));
		}

	      norm = sqrt(norm);

	      amp[0] = h[index - 1][0]/norm;
	      amp[1] = h[index - 1][1]/norm;
	    }

	  // Determine amplitude
	  if(option == 0)
	    amplitude = amp[0]*amp[0] + amp[1]*amp[1];
	  else if(option == 1)
	    amplitude = amp[0];
	  else if(option == 2)
	    amplitude = amp[1];

	  // Reset everything
	  observ->photo.pol.nopol = 0;
	  observ->photo.pol.dpol = 0;
	  observ->photo.pol.doubpol.z_barvec = 0;
	  observ->quant_axis_thesis = 0;
	  observ->photo.pol.doubpol.circbeam = 0;
	  release_matrixelement(matrixelement);
	}
      else
	goto error;
    }

  // Differential cross section
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "diffcs"))
    {
      observ->photo.pol.nopol = 1; // no polarization
      kincoeff = kincoeff_photo(w, pk, mp);

      int i;
      double amp[4][2];

      amplitude = 0.;

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  amplitude += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude *= kincoeff * convertcoeff / 4.;

      // Conversion from dsigma / dOmega  --> dsigma / dt(u)
      if(datapoint->ds_dt || datapoint->ds_du) {
	// dOmega / dt(u) = 2*pi / (2*pk*k)
	// Factor 1e-6 since d sigma / dt(u) in (mub * GeV^-2).
	amplitude *= PI / ( pk * w * 1e-6 );
      }

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Total cross section
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "totcs"))
    {
      observ->photo.pol.nopol = 1;          // no polarization
      spinaver = spinaveraging(observ);     // spinaveraging
      kincoeff = kincoeff_photo(w, pk, mp); // kinematic coefficient

      amplitude = 0; // reset amplitude

      // Set the stepsize for the integration
      observ->photo.kin.anglestep = 3;

      int i;
      double amp[4][2];
      amplitude = 0.;

      // Loop over angle
      for(iang=0; iang <= 180; iang += observ->photo.kin.anglestep)// #calculations = fl(180/anglestep)+1 = 61
	{
	  if(iang == 0 || iang >= 180) // fully fwd and bckwd bin have only half the width
	    a_delta = observ->photo.kin.anglestep / 2.0;
	  else
	    a_delta = observ->photo.kin.anglestep;

	  if (!observ->reg.t_and_u_channel)
	    matrixelement = construct_matrixelement(w,w,cos(iang*PI/180.0),pk,particles,observ,datapoint->label);

	  else {
	    if (cos(iang*PI/180.0) > 0.0 || datapoint->tch) {
	      observ->reg.t_channel = 1;
	      matrixelement = construct_matrixelement(w,w,cos(iang*PI/180.0),pk,particles,observ,datapoint->label);
	    }
	    else {
	      observ->reg.t_channel = 0;
	      matrixelement = construct_matrixelement(w,w,cos(iang*PI/180.0),pk,particles_2,observ,datapoint->label);
	    }
	  }

	  for(i = 0; i < 4; i++)
	    {
	      get_transversity_amplitude(matrixelement, amp[i], i + 1);
	      amplitude += (pow(amp[i][0],2) + pow(amp[i][1],2))
		*2*PI*sin(iang*PI/180.0)*a_delta*PI/180.0;
	    }


	  if (datapoint->label>=0) datapoint->label++;
	} // end loop over angle

      amplitude *= kincoeff * convertcoeff * spinaver;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Single asymmetries (S,T,P)
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "S") ||
	  !strcmp(datapoint->observable, "T") ||
	  !strcmp(datapoint->observable, "P"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      int
	obs = ((int)datapoint->observable[0]) - 80,
	index = 5./12.*obs*obs - 23./12.*obs + 2.,
	posIndex[3][2] = {{0,1},{0,3},{0,2}};

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  amplitude += (i == posIndex[index][0] || i == posIndex[index][1] ? 1. : -1.)
	    *(pow(amp[i][0],2) + pow(amp[i][1],2));
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-recoil asymmetry C_x(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "C_x") ||
	  !strcmp(datapoint->observable, "C_xp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // C_x
		   2.*(amp[0][0]*amp[3][1] - amp[0][1]*amp[3][0] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0])
		   :
		   // C_x'
		   costhk*
		   2.*(amp[0][0]*amp[3][1] - amp[0][1]*amp[3][0] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0]) -
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] -
		       amp[1][0]*amp[2][0] - amp[1][1]*amp[2][1]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-recoil asymmetry C_z(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "C_z") ||
	  !strcmp(datapoint->observable, "C_zp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // C_z
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] -
		      amp[1][0]*amp[2][0] - amp[1][1]*amp[2][1])
		   :
		   // C_z'
		   costhk*
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] -
		       amp[1][0]*amp[2][0] - amp[1][1]*amp[2][1]) +
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[3][1] - amp[0][1]*amp[3][0] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-recoil asymmetry O_x(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "O_x") ||
	  !strcmp(datapoint->observable, "O_xp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // O_x
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] +
		       amp[1][0]*amp[2][0] + amp[1][1]*amp[2][1])
		   :
		   // O_x'
		   costhk*
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] +
		       amp[1][0]*amp[2][0] + amp[1][1]*amp[2][1]) -
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][1]*amp[3][0] - amp[0][0]*amp[3][1] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-recoil asymmetry O_z(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "O_z") ||
	  !strcmp(datapoint->observable, "O_zp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // O_z
		   2.*(amp[0][1]*amp[3][0] - amp[0][0]*amp[3][1] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0])
		   :
		   // O_z'
		   costhk*
		   2.*(amp[0][1]*amp[3][0] - amp[0][0]*amp[3][1] +
		       amp[1][0]*amp[2][1] - amp[1][1]*amp[2][0]) +
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[3][0] + amp[0][1]*amp[3][1] +
		       amp[1][0]*amp[2][0] + amp[1][1]*amp[2][1]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-target asymmetry E
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "E"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = 2.*(amp[0][0]*amp[2][0] + amp[0][1]*amp[2][1] -
		       amp[1][0]*amp[3][0] - amp[1][1]*amp[3][1]);

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-target asymmetry F
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "F"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = 2.*(amp[0][0]*amp[2][1] - amp[0][1]*amp[2][0] +
		      amp[1][0]*amp[3][1] - amp[1][1]*amp[3][0]);

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-target asymmetry G
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "G"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = -2.*(amp[0][1]*amp[2][0] - amp[0][0]*amp[2][1] +
		       amp[1][0]*amp[3][1] - amp[1][1]*amp[3][0]);

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Beam-target asymmetry H
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "H"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = 2.*(amp[0][0]*amp[2][0] + amp[0][1]*amp[2][1] +
		      amp[1][0]*amp[3][0] + amp[1][1]*amp[3][1]);

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Target-recoil asymmetry T_x(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "T_x") ||
	  !strcmp(datapoint->observable, "T_xp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // T_x
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] +
		       amp[2][0]*amp[3][0] + amp[2][1]*amp[3][1])
		   :
		   // T_x'
		   costhk*
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] +
		       amp[2][0]*amp[3][0] + amp[2][1]*amp[3][1]) -
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][1]*amp[1][0] - amp[0][0]*amp[1][1] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1]));


      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Target-recoil asymmetry T_z(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "T_z") ||
	  !strcmp(datapoint->observable, "T_zp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // T_z
		   2.*(amp[0][1]*amp[1][0] - amp[0][0]*amp[1][1] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1])
		   :
		   // T_z'
		   costhk*
		   2.*(amp[0][1]*amp[1][0] - amp[0][0]*amp[1][1] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1]) +
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] +
		       amp[2][0]*amp[3][0] + amp[2][1]*amp[3][1]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Target-recoil asymmetry L_x(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "L_x") ||
	  !strcmp(datapoint->observable, "L_xp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // L_x
		   2.*(amp[0][0]*amp[1][1] - amp[0][1]*amp[1][0] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1])
		   :
		   // L_x'
		   costhk*
		   2.*(amp[0][0]*amp[1][1] - amp[0][1]*amp[1][0] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1]) -
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] -
		       amp[2][0]*amp[3][0] - amp[2][1]*amp[3][1]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Target-recoil asymmetry L_z(')
  // -----------------------------------------------------
  else if(!strcmp(datapoint->observable, "L_z") ||
	  !strcmp(datapoint->observable, "L_zp"))
    {
      observ->photo.pol.nopol = 1; // no polarization

      // Determine matrixelement
      if (!observ->reg.t_and_u_channel)
	matrixelement =
	  construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
      else {
	if (costhk > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles,observ,datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  matrixelement =
	    construct_matrixelement(w,w,costhk,pk,particles_2,observ,datapoint->label);
	}
      }

      int i;
      double amp[4][2], norm = 0.;

      amplitude = 0.;

      for(i = 0; i < 4; i++)
	{
	  get_transversity_amplitude(matrixelement, amp[i], i + 1);
	  norm += (pow(amp[i][0],2) + pow(amp[i][1],2));
	}

      amplitude = (datapoint->observable[3] != 'p' ?
		   // L_z
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] -
		       amp[2][0]*amp[3][0] - amp[2][1]*amp[3][1])
		   :
		   // L_z'
		   costhk*
		   2.*(amp[0][0]*amp[1][0] + amp[0][1]*amp[1][1] -
		       amp[2][0]*amp[3][0] - amp[2][1]*amp[3][1]) +
		   sqrt(1. - costhk*costhk)*
		   2.*(amp[0][0]*amp[1][1] - amp[0][1]*amp[1][0] +
		       amp[2][1]*amp[3][0] - amp[2][0]*amp[3][1]));

      amplitude /= norm;

      // Reset everything
      observ->photo.pol.nopol = 0;
      release_matrixelement(matrixelement);
    }

  // Error, datapoint->observable has not been implemented yet!!
  else {
  error:
    fprintf(stderr,"Observable '%s' is not defined in chifunc()!\n",
	    datapoint->observable);
  }

  return amplitude;
}


struct Electro_observable
{
  Electro* datapoint;
  Observable* observ;
  Class* workparticles;
  Class* workparticles_2;
  double w;
  double k;
  double pk;
  double s;
};


double calc_elec_obs_cos(double theta_k, void *params)
{
  struct Electro_observable *elec_obs = (struct Electro_observable *) params;
  double
    theta_f = elec_obs->datapoint->theta_upper,
    theta_i = elec_obs->datapoint->theta_lower,
    theta_av = acos(elec_obs->datapoint->cos),
    acceptance;

  if(elec_obs->datapoint->cos_binaverage == 1)
    acceptance = 1./(theta_f - theta_i);
  else if(elec_obs->datapoint->cos_binaverage == 2)
    acceptance = (6.*(2.*theta_av - theta_f - theta_i)*theta_k
		  + 4.*(theta_f*theta_f + theta_i*theta_i + theta_f*theta_i) -
		  6.*theta_av*(theta_f + theta_i))/pow(theta_f - theta_i,3);

  return calculate_electro_observable(elec_obs->datapoint, elec_obs->observ,
				      elec_obs->workparticles, elec_obs->workparticles_2,
				      elec_obs->w, elec_obs->k, cos(theta_k), elec_obs->pk,
				      elec_obs->s,elec_obs->datapoint->phi,
				      elec_obs->datapoint->phiMin) * acceptance;
}

/*!
 * This function returns the chi^2 if "datapoint" is an electroproduction
 * point.
 */
double electro_chi(Electro datapoint, Class* workparticles, Class* workparticles_2, Observable* observ, int* modifdatacount)
{
  int weight;
  double mp, mk, my;
  double w, wlab, k, klab, pk, costheta_k;
  double s, t, qsquared;
  double calcpoint, chisquared;

  // Allocate the Born masses
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  // Define kinematic variables of the datapoint

  s = datapoint.s; // Energy
  qsquared = datapoint.qsquared; // Q^2

  if(datapoint.is_w)
    s *= s;

  if(datapoint.cos_ang) // Angle
    costheta_k = datapoint.cos;
  else
    t = datapoint.t;

  // The virtual photon energy and momentum are calculated in lab frame
  wlab = (s + qsquared - mp*mp) / (2 * mp);
  klab = sqrt(qsquared + wlab*wlab);

  // Transform virtual photon energy and momentum from lab to c.m. frame.
  lorentztrans_electro(wlab, &w, klab, &k, mp);

  // "pk" is constructed out of the mass relation in the c.m. frame
  pk = construct_pk(w, k, mp, mk, my);

  // Calculate "costheta_k" or "t"
  if (datapoint.cos_ang)
    t = mk*mk - qsquared - 2 * (sqrt(mk*mk + pk*pk) * w - pk * k * costheta_k);
  else
    costheta_k = (t + qsquared + 2 * sqrt(mk*mk + pk*pk) * w - mk*mk) /
      ( 2 * pk * k);

  // Check if "pk" and "costheta_k" are physical
  is_physical_electro(s, t, qsquared, costheta_k, pk, datapoint.ampli, datapoint.cos_ang);

  // Calculate electro observable
  if(datapoint.cos_binaverage > 0)
    {
      double result, error;
      gsl_function F;
      gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
      struct Electro_observable elec_obs = {&datapoint, observ, workparticles,
					    workparticles_2, w, k, pk, s};
      F.function = &calc_elec_obs_cos;
      F.params = &elec_obs;

      gsl_integration_qag(&F, datapoint.theta_lower, datapoint.theta_upper, 0,
		      1e-3, 1000, 1, ws, &result, &error);
      gsl_integration_workspace_free(ws);

      calcpoint = result;
    }
  else
    calcpoint = calculate_electro_observable(&datapoint,observ,workparticles,
					     workparticles_2,w,k,costheta_k,
					     pk,s,datapoint.phi,datapoint.phiMin);

  /* weight procedure according to
     H.Thom, Phys.Rev. 151, 1322 (1966). */
  weight = 1;

  // Polarization observables
  if( !strcmp ( datapoint.observable, "P" ) ||
      !strcmp ( datapoint.observable, "S" ) ||
      !strcmp ( datapoint.observable, "T" ))
    {
      if(observ->fit.polweight != 0)
	weight *= observ->fit.polweight;
      else
	error_exit("error in polarization weight!!\n");
    }

  // Multiple isospin channels
  if(observ->iso.nr_iso_channels > 1)
    weight *= observ->fit.isoWeight[datapoint.iso];

  // Calculate the chi^2:
  chisquared = weight*pow(((calcpoint - datapoint.ampli) /
                           datapoint.error), 2);

  // According to the weight procedure, the total number of degrees
  // of freedom has to increase.
  (*modifdatacount) += weight - 1;

  // reset observ
  observ->electroprod = 0;

  return chisquared;
}

void is_physical_electro(double s, double t, double qsquared, double costheta_k, double pk, double ampli, short cos_ang)
{
  // Check if pk is physical
  if(pk < 0)
    {
      if (cos_ang)
  	printf("s=%lf  cos=%lf  qq=%lf   amp=%lf\n", s, costheta_k, qsquared, ampli);
      else
	printf("s=%lf  t=%lf  qq=%lf   amp=%lf\n", s, t, qsquared, ampli);
      error_exit("Data point not in the physical plane! (pk < 0)");
    }

  // Check if costheta_k is physical
  if(costheta_k > 1. || costheta_k < -1.)
    {
      if (cos_ang)
	printf("s=%lf  cos=%lf  qq=%lf   amp=%lf\n", s, costheta_k, qsquared, ampli);
      else
	printf("s=%lf  t=%lf  qq=%lf   amp=%lf\n", s, t, qsquared, ampli);
      error_exit("Data point not in the physical plane! (|cos| > 1");
    }
}

/*!
 * Calculates one electroproduction observable in one kinematic point
 * based on the the info in the Data-struct and the Observable-struct.
 *
 * We can calculate different observables:
 * \verbatim
 datapoint->observable = diff_l
                         diff_t
                         diff_r_lt
                         diff_t+l
                         diff_tt_unpol
                         diff_tl_unpol
                         diff_tl_epol
                         diff_tt_epol
                         diff_phi
                         diff_int
                         a_lu
                         induced_pol_y
                         induced_pol_yp
                         induced_pol_n
                         induced_pol_yh
                         transf_pol_x
                         transf_pol_xp
                         transf_pol_t
                         transf_pol_xh
                         transf_pol_z
                         transf_pol_zp
                         transf_pol_l
                         transf_pol_zh
 \endverbatim
 *
 * In case of a polarized hadron, one can choose different
 * reference frames. Here we follow the notation of D.Carman's
 * notes:
 * \verbatim
 (x,y,z):    z along photon
             y perpendicular to electron plane
 (xp,yp,zp): zp along kaon
             yp perpendicular to hadron plane
 (l,t,n):    l along hyperon
             n perpendicular to hadron plane
 (xh,yh,zh): zh along photon
             yh perpendicular to hadron plane
 \endverbatim
 *
 */
double calculate_electro_observable(Electro *datapoint, Observable *observ, Class *particles,
				    Class *particles_2, double w, double k, double costheta_k,
				    double pk, double s, double phi, double phiMin)
{
  double calcpoint, sigma0;
  double spinaver, e_el_conv_fac, eps;
  double mp, mk, my;
  short e_L, TL_2;
  Response_func resp = {0};
  Response_func_calc respcalc;

  // Allocate the Born masses
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  // Set the observable specification for electroproduction !!
  observ->electroprod = 1;
  observ->photoprod = 0;
  observ->kaoncapture = 0;

  // Determine the spinaveraging factor
  spinaver = spinaveraging(observ);

  // Some kinematics
  double qsquared = k*k - w*w;
  double wlab = (s + qsquared - mp*mp) / (2 * mp);
  double klab = sqrt(qsquared + wlab*wlab);

  /* determine front-factor convention
   * ---------------------------------
   * 0 -> e ds_L, sqrt(e(e+1)) ds_TL('), e_L = e k_lab^2 / Q^2
   * 1 -> e ds_L, sqrt(2e(e+1)) ds_TL('), e_L = e k_lab^2 / Q^2
   * 2 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega^2 / Q^2
   * 3 -> e_L ds_L, sqrt(2e_L(e+1)) ds_TL('), e_L = e omega_lab^2 / Q^2
   */
  switch(datapoint->cs_convention) {
  case 0:
    e_L = 0;
    e_el_conv_fac = klab * klab / qsquared;
    TL_2 = 0;
    break;
  case 1:
    e_L = 0;
    e_el_conv_fac = klab * klab / qsquared;
    TL_2 = 1;
    break;
  case 2:
    e_L = 1;
    e_el_conv_fac = wlab * wlab / qsquared;
    TL_2 = 1;
    break;
  case 3:
    e_L = 1;
    e_el_conv_fac = w * w / qsquared;
    TL_2 = 1;
    break;
  default:
    error_exit("cs_convention for some datapoint not implemented!");
    break;
  }


  /* Calculate observable
   * --------------------
   * For each observable it suffices to determine
   * a selection of response functions.
   * This info is passed to calc_response_func()
   * through 'respcalc'.
   *
   * The response functions are stored in 'resp'. Inside
   * this function the response functions are combined
   * to form 'calcpoint', which is then used to determine
   * the chisquared of this datapoint
   */
  respcalc.L = 0;
  respcalc.T = 0;
  respcalc.c_TT = 0;
  respcalc.c_TL = 0;
  respcalc.s_TT = 0;
  respcalc.s_TL = 0;
  respcalc.TT_pol = 0;
  respcalc.s_TL_pol = 0;
  respcalc.c_TL_pol = 0;


  if((datapoint->observable[0] == 'M')
                    && string(datapoint->observable).length() == 6)
    {
    /*-----------------------------------
    * Construction of M    (dimensionless)
    *-----------------------------------*/
    int L, Ly, Lp;

    // get the photon helicity
    if(datapoint->observable[1] == '-') L = 0;
    else if(datapoint->observable[1] == '0') L = 1;
    else if(datapoint->observable[1] == '+') L = 2;
    else goto error;

    // get the proton helicity
    if(datapoint->observable[2] == '-') Lp = 0;
    else if(datapoint->observable[2] == '+') Lp = 1;
    else goto error;

    // get the hyperon helicity
    if(datapoint->observable[3] == '-') Ly = 0;
    else if(datapoint->observable[3] == '+') Ly = 1;
    else goto error;

    // get the part we need
    int option;
    if(datapoint->observable[4] == '^' && datapoint->observable[5] == '2')
        option = 0;
    else if(datapoint->observable[4] == '_' && datapoint->observable[5] == 'R')
        option = 1;
    else if(datapoint->observable[4] == '_' && datapoint->observable[5] == 'I')
        option = 2;
    else
        goto error;

    observ->elec.elec_pol = 1;
//    std::cerr << "Setting polarization to 1\n";

    /* Construct the matrixelement
     * ***************************/
    TMatrixElement* me;
    if (!observ->reg.t_and_u_channel)
          me = TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles, observ, datapoint->label);
    else{
        if(costheta_k > 0.0 || datapoint->tch) {
            observ->reg.t_channel = 1;
            me = TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles, observ, datapoint->label);
        }else{
            observ->reg.t_channel = 0;
            me = TMatrixElement::GetMatrixElement(w, k, costheta_k, pk, particles_2, observ, datapoint->label);
        }
    }

    // real and imaginary part of an amplitude
    double amp[2];
    // now we really calculate the hadronic amplitude
    std::complex<double> amplitude = me->calculateM(L,Lp,Ly);

    amp[0] = real(amplitude);
    amp[1] = imag(amplitude);

    // Determine amplitude
    if(option == 0)
        calcpoint = amp[0]*amp[0] + amp[1]*amp[1];
    else if(option == 1)
        calcpoint = amp[0]; // real part
    else if(option == 2)
        calcpoint = amp[1]; // imag part

  }// end calculation of M amplitude

  else if(!strcmp(datapoint->observable, "diff_l"))
    {
      // We only need longitudinal responce
      respcalc.L = 1;

      // Calculate resp.L
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.L *spinaver;

      // Choose correct convention
      if(e_L)
	calcpoint *= e_el_conv_fac;

    } // end "diff_l" observable

  else if(!strcmp(datapoint->observable, "diff_t"))
    {
      // We only need the transverse response
      respcalc.T = 1;

      // Calculate resp.T
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.T *spinaver;

    } // end "diff_t" observable

  else if(!strcmp(datapoint->observable, "diff_t+l"))
    {
      // We need transverse and longitudinal response
      respcalc.L = 1;
      respcalc.T = 1;

      // Calculate resp.T and resp.L
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;
    // eps = 0 default!

      // Determine the observable
      calcpoint = ( resp.T + eps*resp.L ) *spinaver;

    } // end "diff_t+l" observable

  else if(!strcmp(datapoint->observable, "diff_r_lt"))
    {
      // We need transverse and longitudinal response
      respcalc.L = 1;
      respcalc.T = 1;

      // Calculate resp.T and resp.L
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.L / resp.T;

    } // end "diff_r_lt" observable

  else if(!strcmp(datapoint->observable, "diff_tt_unpol"))
    {
      // We need c^R_TT response
      respcalc.c_TT = 1;

      // Calculate resp.c_TT
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.c_TT *spinaver;

    } // end "diff_tt_unpol" observable

  else if(!strcmp(datapoint->observable, "diff_tl_unpol"))
    {
      // We need c^R_TT response
      respcalc.c_TL = 1;

      // Calculate resp.c_TL
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.c_TL *spinaver;

      // Choose correct convention
      if(e_L)
	calcpoint *= sqrt(e_el_conv_fac);
      if(TL_2)
	calcpoint /= sqrt(2);

    } // end "diff_tl_unpol" observable

  else if(!strcmp(datapoint->observable, "diff_tt_epol"))
    {
      // We need R_TT' response
      respcalc.TT_pol = 1;

      observ->elec.elec_pol = 1;

      // Calculate resp.TT_pol
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.TT_pol *spinaver;

    } // end "diff_tt_epol" observable

  else if(!strcmp(datapoint->observable, "diff_tl_epol"))
    {
      // We need s^R_TL' response
      respcalc.s_TL_pol = 1;

      observ->elec.elec_pol = 1;

      // Calculate resp.s_TL_pol
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine the observable
      calcpoint = resp.s_TL_pol *spinaver;

      // Choose correct convention
      if(e_L)
	calcpoint *= sqrt(e_el_conv_fac);
      if(TL_2)
	calcpoint /= sqrt(2);

    } // end "diff_tl_epol" observable

  else if(!strcmp(datapoint->observable, "diff_phi"))
    {
      // We need transverse, longitudinal, and unpolarized TT and TL response
      respcalc.L = 1;
      respcalc.T = 1;
      respcalc.c_TT = 1;
      respcalc.c_TL = 1;

      // Calculate resp.T, resp.L, resp.c_TT and resp.c_TL,
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;

      // Determine the observable
      calcpoint = spinaver*(resp.T + eps*resp.L + eps*resp.c_TT*cos(2.*phi)
			    + sqrt(eps*(1. + eps))*resp.c_TL*cos(phi));

    } // end "diff_phi" observable

  else if(!strcmp(datapoint->observable, "diff_int"))
    {
      // We need transverse, longitudinal, and unpolarized TT and TL response
      respcalc.L = 1;
      respcalc.T = 1;
      respcalc.c_TT = 1;
      respcalc.c_TL = 1;

      // Calculate resp.T, resp.L, resp.c_TT and resp.c_TL,
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;

      // Determine the observable
      calcpoint = spinaver*
	((resp.T + eps*resp.L)*(phi - phiMin)
	 + eps*resp.c_TT*(sin(2.*phi) - sin(2.*phiMin))/2.
	 + sqrt(eps*(1. + eps))*resp.c_TL*(sin(phi) - sin(phiMin)));

    } // end "diff_int" observable

  else if(!strcmp(datapoint->observable, "a_lu"))
    {
      // We need transverse, longitudinal and polarized TL response
      respcalc.L = 1;
      respcalc.T = 1;
      respcalc.s_TL_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;

      // Calculate resp.T, resp.L, resp.c_TT and resp.c_TL,
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;

      // Determine the observable
      calcpoint = sqrt(eps*(1. - eps))*resp.s_TL_pol/(resp.T + eps*resp.L);

    } // end "a_l" observable

  else if(!strcmp(datapoint->observable, "induced_pol_yp") ||
	  !strcmp(datapoint->observable, "induced_pol_yh") ||
	  !strcmp(datapoint->observable, "induced_pol_n")  )
    {
      // Formula taken from Table 1 from Carman's notes:
      //
      //          R^{y'0}_T + \epsilon R^{y'0}_L
      // P^0_y' = ------------------------------
      //          R^{00}_T  + \epsilon R^{00}_L    -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.T and resp.L (with baryon polarization along y')
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.y_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.y_pol = 0;

      // Calculate observable
      calcpoint = ( resp.T + eps*resp.L ) / sigma0;


    } // end "induced_pol_(yp/yh/n)" observable


  else if(!strcmp(datapoint->observable, "induced_pol_y"))
    {
      // Formula taken from Table 2 from Carman's notes:
      //
      // P^0_y  = 1/2 * sqrt( eps( 1 + eps ) ) *
      //
      //          s^R^{x'0}_TL costhkcm + c^R^{y'0}_TL + s^R^{z'0}_TL sinthkcm
      //          ------------------------------------------------------------
      //                         R^{00}_T  + \epsilon R^{00}_L            -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.s_TL (with hyperon polarization along x')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.s_TL = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.x_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint = resp.s_TL * costheta_k;

      // Calculate resp.s_TL (with hyperon polarization along z')
      observ->elec.bar_pol.x_pol = 0;
      observ->elec.bar_pol.z_pol = 1;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint += resp.s_TL * sqrt( 1 - costheta_k * costheta_k );

      // Calculate resp.c_TL (with hyperon polarization along y')
      respcalc.s_TL = 0;
      respcalc.c_TL = 1;
      observ->elec.bar_pol.z_pol = 0;
      observ->elec.bar_pol.y_pol = 1;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint += resp.c_TL;
      calcpoint *= sqrt(eps*(1+eps)) / sigma0 / 2.0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.y_pol = 0;


    } // end "induced_pol_y" observable


  else if(!strcmp(datapoint->observable, "transf_pol_z") ||
	  !strcmp(datapoint->observable, "transf_pol_zh") )
    {
      // Formula taken from Table 2 from Carman's notes:
      //
      // P'_z   = sqrt( 1 - eps^2 ) *
      //
      //          - R^{x'0}_TT' sinthkcm + R^{z'0}_TT' costhkcm
      //          ---------------------------------------------
      //                   R^{00}_T  + \epsilon R^{00}_L         -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.TT_pol (with hyperon polarization along x')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.TT_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.x_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint = -1.0 * resp.TT_pol * sqrt(1 - costheta_k*costheta_k);

      // Calculate resp.TT_pol (with hyperon polarization along z')
      observ->elec.bar_pol.x_pol = 0;
      observ->elec.bar_pol.z_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint += resp.TT_pol * costheta_k;
      calcpoint *= sqrt(1-eps*eps) / sigma0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.z_pol = 0;

    } // end observable "transf_pol_(z/zh)"


  else if(!strcmp(datapoint->observable, "transf_pol_zp") ||
	  !strcmp(datapoint->observable, "transf_pol_l")  )
    {
      // Formula taken from Table 1 from Carman's notes:
      //
      //          sqrt(1 - eps^2)   R^{z'0}_TT'
      // P'_z'  = -----------------------------
      //          R^{00}_T  + \epsilon R^{00}_L   -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.TT_pol (with hyperon polarization along z')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.TT_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.z_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint = sqrt(1-eps*eps) * resp.TT_pol / sigma0;
      if(!strcmp(datapoint->observable, "transf_pol_l"))
	calcpoint *= -1.0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.z_pol = 0;

    } // end observable "transf_pol_(zp/l)"


  else if(!strcmp(datapoint->observable, "transf_pol_xp") ||
	  !strcmp(datapoint->observable, "transf_pol_t")  )
    {
      // Formula taken from Table 1 from Carman's notes:
      //
      //          sqrt(1 - eps^2)   R^{x'0}_TT'
      // P'_x'  = -----------------------------
      //          R^{00}_T  + \epsilon R^{00}_L   -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.TT_pol (with hyperon polarization along x')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.TT_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.x_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint = sqrt(1-eps*eps) * resp.TT_pol / sigma0;
      if(!strcmp(datapoint->observable, "transf_pol_t"))
	calcpoint *= -1.0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.x_pol = 0;

    } // end observable "transf_pol_(xp/t)"


  else if(!strcmp(datapoint->observable, "transf_pol_x"))
    {
      // Formula taken from Table 2 from Carman's notes:
      //
      // P'_x   = 1/2 * sqrt( eps( 1 - eps ) ) *
      //
      //          c^R^{x'0}_TL' costhkcm - s^R^{y'0}_TL' + c^R^{z'0}_TL' sinthkcm
      //          ---------------------------------------------------------------
      //                             R^{00}_T  + \epsilon R^{00}_L        -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.c_TL_pol (with hyperon polarization along x')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.c_TL_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.x_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint = resp.c_TL_pol * costheta_k;

      // Calculate resp.c_TL_pol (with hyperon polarization along z')
      observ->elec.bar_pol.x_pol = 0;
      observ->elec.bar_pol.z_pol = 1;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint += resp.c_TL_pol * sqrt( 1 - costheta_k * costheta_k );

      // Calculate resp.s_TL_pol (with hyperon polarization along y')
      respcalc.c_TL_pol = 0;
      respcalc.s_TL_pol = 1;
      observ->elec.bar_pol.z_pol = 0;
      observ->elec.bar_pol.y_pol = 1;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint -= resp.s_TL_pol;
      calcpoint *= sqrt(eps*(1-eps)) / sigma0 / 2.0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.y_pol = 0;


    } // end "transf_pol_x" observable


  else if(!strcmp(datapoint->observable, "transf_pol_xh"))
    {
      // Formula taken from Table 3 from Carman's notes:
      //
      // P'_xh  = sqrt( 1 - eps^2 ) *
      //
      //          R^{x'0}_TT' costhkcm + R^{z'0}_TT' sinthkcm
      //          -------------------------------------------
      //                 R^{00}_T  + \epsilon R^{00}_L         -> sigma0

      // Determine epsilon
      if (datapoint->beam_ener_input)
	eps = calc_epsilon(wlab, klab,
			   qsquared, datapoint->e_beam_ener);
      else
	eps = datapoint->eps;


      // Calculate resp.T and resp.L (without baryon polarization)
      respcalc.L = 1;
      respcalc.T = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ,datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate sigma0
      sigma0 = ( resp.T + eps*resp.L );


      // Calculate resp.TT_pol (with hyperon polarization along x')
      respcalc.L = 0;
      respcalc.T = 0;
      respcalc.TT_pol = 1;
      observ->elec.elec_pol = 1;
      observ->elec.baryon_pol = 1;
      observ->elec.bar_pol.recpol = 1;
      observ->elec.bar_pol.x_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Temporary observable
      calcpoint = resp.TT_pol * costheta_k;

      // Calculate resp.TT_pol (with hyperon polarization along z')
      observ->elec.bar_pol.x_pol = 0;
      observ->elec.bar_pol.z_pol = 1;
      observ->quant_axis_thesis = 0;
      if (!observ->reg.t_and_u_channel)
	calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			   s, mp, particles, observ, datapoint->label);
      else {
	if (costheta_k > 0.0 || datapoint->tch) {
	  observ->reg.t_channel = 1;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles, observ, datapoint->label);
	}
	else {
	  observ->reg.t_channel = 0;
	  calc_response_func(&resp, &respcalc, w, k, costheta_k, pk,
			     s, mp, particles_2, observ, datapoint->label);
	}
      }

      // Calculate the observable
      calcpoint += resp.TT_pol * sqrt(1 - costheta_k*costheta_k);
      calcpoint *= sqrt(1-eps*eps) / sigma0;

      // Reset observ
      observ->elec.elec_pol = 0;
      observ->elec.baryon_pol = 0;
      observ->elec.bar_pol.recpol = 0;
      observ->elec.bar_pol.z_pol = 0;

    } // end observable "transf_pol_xh"


  else
    {
      error:
      throw std::logic_error("Error in calculate_electro_observable():\nObservable of type '"
                             + string(datapoint->observable) + "' not implemented!\n\n");
    }


  //  If needed, conversion from dsigma / dOmega --> dsigma / dt.
  if(datapoint->ds_dt && !string(datapoint->observable).substr(0, 4).compare("diff"))
    calcpoint *=  PI * 1e6 / ( pk * k); // dOmega / dt = 2*pi / (2*pk*k)


  return calcpoint;
}

/*! Returns the chi^2 if "datapoint" is a radiative kaon capture observable.
 */
double kaoncap_chi(Kaoncap datapoint, Class* workparticles, Class* workparticles_2,
		   Observable* observ, int* modifdatacount)
{
  int weight, ipklab;
  double w, pk, pklab;
  double ampli=0.0, tmpampli=0.0, sumampli=0.0, chisquared=0.0;
  double gridnum, realgridnum=0.0;
  double mp, mk, my;

  /* Allocate the Born masses */
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  gridnum = get_gridnum(observ->fit.narrow_grid, datapoint.observable);

  /* If there is only one kaon momentum, err_pk was put = 0.0 in the
   * getdatastructure() function. In which case
   * the for-loop construction is only run through once. */

  // Loop over pklab inside bin
  for(ipklab = 0; ipklab < gridnum; ipklab++) // #calculations = gridnum
    {
      pklab = datapoint.pk + datapoint.err_pk*((2.0*ipklab+1.0)/gridnum-1.0);

      // Conversion from lab to c.m. kaon momentum (necessary in calculations)
      lorentztrans_kaon(pklab, &pk, mp, mk);

      // Calculation of the photon energy in com. frame.
      w = construct_w(pklab, mp, mk, my);

      // Test whether pk and w are physical
      if(pk < 0.0 || w <= 0.0)
        throw std::logic_error("A datapoint of type '" + string(datapoint.observable) +
                               "' lies outside the physical region!\n\n");

      // If pk is physical, proceed...
      else
	{
	  // Calculate the appropriate observable
	  tmpampli = calculate_kaoncap_observable(&datapoint,observ,workparticles,
						  workparticles_2,w,w,datapoint.cos,pk);

	  // Add observable to bin-total
	  sumampli += tmpampli;

	  // Count gridpoints (realgridnum can be different from gridnum!)
	  realgridnum += 1;

	} // end pk = physical

	if (datapoint.label>=0) datapoint.label++; // increment it either way, because we have already taken into consideration gridnum points

    } // end loop in energy-bin

  // Average observable
  ampli = sumampli / realgridnum;


  // General check for unphysical behavior. !
  if(ampli < 0.0 && (!strcmp(datapoint.observable, "diffcs") ||
		     !strcmp(datapoint.observable, "totcs") ||
		     !strcmp(datapoint.observable, "ptcs")))
    {
      fprintparticspec(stderr, workparticles, *observ);
      std::stringstream ss;

      ss << "Error: neg. modulus !!\n(pklab = " << datapoint.pk
         << "; costheta = " << datapoint.cos << ")\n";

      throw std::logic_error(ss.str());
    }

  weight = 1;

  /* weight procedure according to
     H.Thom, Phys.Rev. 151, 1322 (1966). */

  // Multiple isospin channels
  if(observ->iso.nr_iso_channels > 1)
    weight *= observ->fit.isoWeight[datapoint.iso];

  // Calculate the chi^2:
  chisquared = weight* pow(((ampli - datapoint.ampli) / datapoint.error), 2);

  // For testing purposes, now OK (TVC20090205)
  //printf("%lf\t%lf\t%lf\t%lf\n", datapoint.pk, datapoint.cos, ampli, datapoint.ampli);
  //printf("%lf\t%lf\t%lf\n", datapoint.pk, datapoint.cos, ampli);

  /* According to the weight procedure, the total number of degrees of freedom
   * has to increase. */
  (*modifdatacount) += weight - 1;

  observ->kaoncapture = 0;

  return chisquared;
}

/*!
 * Calculates one kaon capture observable in one kinematic point
 * based on the the info in the Data-struct and the Observable-struct.
 *
 * We can calculate different observables:
 * \verbatim
 datapoint->observable = bran        branching ratio for stopped kaons
                         dcs         differential cross section
                         tcs         total cross section
                         ptcs        partial total cross section
 \endverbatim
 */
double calculate_kaoncap_observable(Kaoncap *datapoint, Observable *observ, Class *particles,
				    Class *particles_2, double w, double k, double costhk, double pk)
{
  double amplitude;
  double pseudopot;
  double mp, mk, my;
  int iang, no_intervals;
  double kincoeff, spinaver;
  double a_delta, angb, ange, anglestep, ang;

  /* kaon capture calculations require isospin 11 or 12 */
  if(datapoint->iso != 11 && datapoint->iso != 12)
    throw std::logic_error("ERROR in calculate_kaoncap_observable(..): kaon capture calculations require isospin channel 11 or 12!\n");

  observ->iso.isospin = datapoint->iso;

  /* Allocate the Born masses */
  allocate_born_masses(&mp,&mk,&my,observ->iso.isospin);

  // Set the observable specification for kaon capture !!
  observ->photoprod = 0;
  observ->electroprod = 0;
  observ->kaoncapture = 1;

  // Coeff for conversion from MeV^-2 to mubarn
  double convertcoeff = 197.32 * 197.32 * 1e4;


  /* Calculate observable
   * -------------------- */

  // Fit to branching ratio for stopped kaons
  // ----------------------------------------
  if(!strcmp(datapoint->observable, "bran"))
    {
      observ->photo.pol.nopol = 1;          // no polarization (for the photon)
      observ->kaoncapture = 1;
      observ->kaoncap.stoppedkaon = 1;
      observ->photoprod = 0;
      observ->electroprod = 0;

      spinaver = spinaveraging(observ);

      /* stopped kaons ==> pk=0 */
      pk = 0.0;

      /* photon energy (lab frame is com. frame because pk = 0) */
      w = 0.5 * ( mk+mp-my*my/(mk+mp) );

      /* (K^- p) pseudo-potential according to Burkhardt et al. */
      pseudopot = 560.0;

      /* Kinematical coefficient for ratio decay widths */
      kincoeff = PI * w / (2.0 * pseudopot * mp * mk * (mp + mk));
      amplitude = determine_M2_photo(w, w, 1, pk, particles, observ, datapoint->label);

      /* Calculate ratio of decay width in (MeV fm)^-3 */
      amplitude *= kincoeff * spinaver * pow(197.32,3);
    }

  // Fit to differential cross section data
  // --------------------------------------
  else if(!strcmp(datapoint->observable, "diffcs"))
    {
      observ->photo.pol.nopol = 1;          // no polarization
      spinaver = spinaveraging(observ);     // spinaveraging
      kincoeff = kincoeff_kaoncap(w, pk, mp, mk); // kinematic coefficient

      // Determine amplitude
      amplitude = determine_M2_photo(w, w, costhk, pk, particles, observ, datapoint->label);

      // For testing purposes, now OK (TVC20090205)
      //printf("%lf\t%lf\t%lf\t\n", datapoint->pk, datapoint->cos, amplitude);

      amplitude *= kincoeff * convertcoeff * spinaver;

      observ->photo.pol.nopol = 0;
    }

  // Fit to total cross section
  // --------------------------
  else if(!strcmp(datapoint->observable, "ptcs") ||
	  !strcmp(datapoint->observable, "totcs"))
    {
      observ->photo.pol.nopol = 1;          // no polarization
      spinaver = spinaveraging(observ);     // spinaveraging
      kincoeff = kincoeff_kaoncap(w, pk, mp, mk); // kinematic coefficient

      amplitude = 0.; // reset amplitude

      /*
       * In the case of the 2001 analysis, the available data for the
       * kaon capture process is from Crystal Ball (dissertation from
       * Phaisangittisakul). In the analysis, only partial total
       * cross sections were presented where the summation was done
       * over: -0.8 < cos(theta) < 0.8
       * The switch observ->kaoncap.partialtotcs tells us whether to
       * integrate over the interval [-0.8,0.8] for costhk, or over the
       * full interval [-1,1].
       * The interval for costhk to be integrated over is specified in the
       * datapoint->cosmin and ->cosmax variables.
       */

      /* cosmax gives upper limit for costhk, but lower limit for theta-angle */
      if(!strcmp(datapoint->observable, "ptcs"))
	{
	  angb = acos(datapoint->cosmax)*180.0/PI; /* lower angle in degrees */
	  ange = acos(datapoint->cosmin)*180.0/PI; /* upper angle in degrees */
	}
      else /* for tcs */
	{
	  angb = 0.0;
	  ange = 180.0;
	}
      no_intervals = 100;
      anglestep = (ange-angb)/no_intervals;

      // Loop over angle
      for (iang = 0; iang <= no_intervals; iang++) // #calculations *= no_intervals+1  = 101
	{
	  /* simple trapezium rule */
	  if (iang == 0 || iang == no_intervals)
	    a_delta = anglestep / 2.0;
	  else
	    a_delta = anglestep;

	  ang = PI/180.0*(angb + (double) iang * (ange - angb)/no_intervals);

	  amplitude += determine_M2_photo(w, w, cos(ang), pk,
					  particles, observ,datapoint->label)
	      * (2 * PI) * sin(ang) * a_delta * PI/180.0;
	  if (datapoint->label>=0) datapoint->label++;
	} // end loop over angle


      amplitude *= kincoeff * convertcoeff * spinaver;

      observ->photo.pol.nopol = 0;
    }

  // Error, datapoint->observable has not been implemented yet!!
  else
    throw std::logic_error("Observable " + string(datapoint->observable) + "is not defined in chifunc()!\n");

  return amplitude;
}


/*!
 * Import the dummy tot. cs. data (iso 5 and 6)
 */
int import_exp_tot_dummy_5_6(Data datapoints[], int*  datacount,
			     char* file, int iso,Observable* observ)
{
  int i, lncount, startpoint;
  double energy, error;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;
  long pos;

  startpoint = *datacount;

  if(iso == 5)
    strcpy(iso_char, "5");
  else if(iso == 6)
    strcpy(iso_char, "6");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "photo/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);


  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* dump header of the file */

  for(i=0; i<11; i++)
    fscanf(ifp, "%s", dump);

  getc(ifp);
  getc(ifp);

  pos = ftell(ifp);

  lncount = line_count(ifp, pos);

  fseek(ifp, pos, SEEK_SET);

  for(i = startpoint; i < lncount + startpoint; i++)
    {
      datapoints[i].photo.iso = iso;

      strcpy(datapoints[i].photo.observable, "totcs");
      datapoints[i].photo.ds_dt = 0;
      datapoints[i].photo_prod = 1;
      datapoints[i].electro_prod = 0;
      datapoints[i].kaoncapture = 0;

      fscanf(ifp, "%lf", &energy);

      datapoints[i].photo.emin = energy;
      datapoints[i].photo.emax = energy;

      fscanf(ifp, "%lf", &datapoints[i].photo.ampli);
      fscanf(ifp, "%s%lf", dump, &error);

      datapoints[i].photo.error = error;

      (*datacount)++;
    }

  fclose(ifp);

  return 0;
}


/*
 * Import the radiative capture data from the brookhaven experiment
 * (iso 1 and 2):
 *
 * D.A. Whitehouse et al. Phys. Rev. Lett. 63, 1352 (1989)
 */
int import_exp_brookhaven_stopped_kaoncap_11_12(Data datapoints[],
						int* datacount, char* file,
						int iso,Observable* observ)
{
  int i, start, headerlength;
  double ratio, staterr, syserr;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;


  start = *datacount;

  if(iso == 11)
    strcpy(iso_char, "11");
  else if(iso == 12)
    strcpy(iso_char, "12");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);


  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");


  headerlength = 23;

  /* dump header of the file */

  for(i=0; i<headerlength; i++)
    fscanf(ifp, "%s", dump);


  strcpy(datapoints[start].kaoncap.observable,"bran");
  datapoints[start].photo_prod = 0;
  datapoints[start].electro_prod = 0;
  datapoints[start].kaoncapture = 1;


  fscanf(ifp, "%lf%s%lf%s%lf", &ratio, dump, &staterr, dump, &syserr);

  datapoints[start].kaoncap.pk = 0.0;
  /* Used to be
     datapoints[start].kaoncap.ratio = ratio * 1e-3;
  */
  datapoints[start].kaoncap.ampli = ratio * 1e-3;
  datapoints[start].kaoncap.error =
    sqrt(staterr*staterr + syserr*syserr) * 1e-3;

  (*datacount)++;

  datapoints[start].kaoncap.iso = iso;

  fclose(ifp);

  return 0;
}

/*!
 * Import the diff. cs. data from the Crystal Ball experiment (iso 1)
 * for the kaon capture process.
 *
 * Dissertation of Nakorn Phaisangittisakul
 */
int import_exp_diff_crystalball_kaoncap_11(Data datapoints[], int* datacount,
					   char* file,Observable* observ)
{
  int i, lncount, startpoint, headerlength;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;
  long pos;

  startpoint = *datacount;

  strcpy(iso_char, "11");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);

  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  headerlength = 20;

  /* dump header of the file */

  for(i=0; i < headerlength; i++)
    fscanf(ifp, "%s", dump);

  getc(ifp);
  getc(ifp);

  pos = ftell(ifp);
  lncount = line_count(ifp, pos);
  fseek(ifp, pos, SEEK_SET);

  for(i = startpoint; i < lncount + startpoint; i++)
    {
      datapoints[i].kaoncap.iso = 11;

      strcpy(datapoints[i].kaoncap.observable,"diffcs");

      datapoints[i].photo_prod = 0;
      datapoints[i].electro_prod = 0;
      datapoints[i].kaoncapture = 1;

      fscanf(ifp, "%lf", &datapoints[i].kaoncap.pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.cos);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);

      (*datacount)++;

    }

  fclose(ifp);

  return 0;
}

/*!
 * Import the tot. cs. data from the Crystal Ball experiment (iso 1)
 * for the kaon capture process.
 *
 * Dissertation of Nakorn Phaisangittisakul
 */
int import_exp_tot_crystalball_kaoncap_11(Data datapoints[], int*  datacount,
					  char* file,Observable* observ)
{
  int i, lncount, startpoint, headerlength;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;
  long pos;

  startpoint = *datacount;

  strcpy(iso_char, "11");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);


  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* dump header of the file */

  headerlength = 18;

  for(i=0; i < headerlength; i++)
    fscanf(ifp, "%s", dump);

  getc(ifp);
  getc(ifp);

  pos = ftell(ifp);

  lncount = line_count(ifp, pos);

  fseek(ifp, pos, SEEK_SET);

  for(i = startpoint; i < lncount + startpoint; i++)
    {
      datapoints[i].kaoncap.iso = 11;

      strcpy(datapoints[i].kaoncap.observable, "ptcs");

      datapoints[i].photo_prod = 0;
      datapoints[i].electro_prod = 0;
      datapoints[i].kaoncapture = 1;
      datapoints[i].kaoncap.cosmax = 0.8;
      datapoints[i].kaoncap.cosmin = -0.8;

      fscanf(ifp, "%lf", &datapoints[i].kaoncap.pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);

      (*datacount)++;
    }

  fclose(ifp);

  return 0;
}

/*!
 * Import the diff. cs. data from the Crystal Ball experiment (iso 11)
 * for the kaon capture process. Analysis from 2009.
 *
 * Contact/paper from Sergey Prakhov.
 */
int import_exp_diff_cb2009_kaoncap_11(Data datapoints[], int* datacount,
				      char* file,Observable* observ)
{
  int i, j, k, lncount, startpoint;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;
  double pk[8], err_pk[8];

  startpoint = *datacount;

  strcpy(iso_char, "11");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);

  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* headerlength = 66; */

  /* look for first row with the values of the kaon momenta */
  strcpy(dump,"start");
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"pK"));

  /* read the eight values of pK */
  for (i=0; i<8; i++)
    fscanf(ifp, "%lf", &(pk[i]));

  /* dump the err_pk string */
  fscanf(ifp, "%s", dump);

  /* read the eight values of err_pK */
  for (i=0; i<8; i++)
    fscanf(ifp, "%lf", &(err_pk[i]));

  /* dump some more text */
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"bars."));

  /* number of datapoints in this file */
  lncount = 92;

  i = startpoint;
  for (k=0; k<11; k++)
    {
      fscanf(ifp, "%s", dump);
      for (j=0; j<8; j++)
	{
	  datapoints[i].kaoncap.iso = 11;

	  strcpy(datapoints[i].kaoncap.observable,"diffcs");

	  datapoints[i].photo_prod = 0;
	  datapoints[i].electro_prod = 0;
	  datapoints[i].kaoncapture = 1;

	  datapoints[i].kaoncap.pk = pk[j];
	  datapoints[i].kaoncap.err_pk = err_pk[j];
	  datapoints[i].kaoncap.cos = -11.0/12.0+((double)k)/6.0;
	  datapoints[i].kaoncap.err_cos = 1.0/12.0;

	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
	  i++;
	}
      i -= 8;
      fscanf(ifp, "%s", dump);
      for (j=0; j<8; j++)
	{
	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);
	  i++;
	}
      (*datacount) += 8;
    }
  for (k=11; k<12; k++)
    {
      fscanf(ifp, "%s", dump);
      for (j=0; j<4; j++)
	{
	  datapoints[i].kaoncap.iso = 11;

	  strcpy(datapoints[i].kaoncap.observable,"diffcs");

	  datapoints[i].photo_prod = 0;
	  datapoints[i].electro_prod = 0;
	  datapoints[i].kaoncapture = 1;

	  datapoints[i].kaoncap.pk = pk[j];
	  datapoints[i].kaoncap.err_pk = err_pk[j];
	  datapoints[i].kaoncap.cos = -11.0/12.0+((double)k)/6.0;
	  datapoints[i].kaoncap.err_cos = 1.0/12.0;

	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
	  i++;
	}
      i -= 4;
      for (j=0; j<5; j++)
	fscanf(ifp, "%s", dump);
      for (j=0; j<4; j++)
	{
	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);
	  i++;
	}
      (*datacount) += 4;
    }

  /*  Test  */
  if ( (startpoint+lncount) != (*datacount) )
    error_exit("Wrong counting of the data points in this file!\n");

  fclose(ifp);

  return 0;
}

/*!
 * Import the total cs. data from the Crystal Ball experiment (iso 11)
 * for the kaon capture process. Analysis from 2009.
 *
 * Contact/paper from Sergey Prakhov.
 */
int import_exp_tot_cb2009_kaoncap_11(Data datapoints[], int* datacount,
				     char* file,Observable* observ)
{
  int i, lncount, startpoint;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;

  startpoint = *datacount;

  strcpy(iso_char, "11");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);

  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* look for the ending of the text */
  strcpy(dump,"start");
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"err_cs"));

  /* number of datapoints in this file */
  lncount = 8;

  for (i=startpoint; i<startpoint+lncount; i++)
    {
      datapoints[i].kaoncap.iso = 11;

      strcpy(datapoints[i].kaoncap.observable,"totcs");

      datapoints[i].photo_prod = 0;
      datapoints[i].electro_prod = 0;
      datapoints[i].kaoncapture = 1;
      /* For the moment, we're treating these data points as total cross sections
	 and not partial total cross sections, awaiting comments from Sergey.
      datapoints[i].kaoncap.cosmax = 1.0;
      datapoints[i].kaoncap.cosmin = -1.0;
      */

      fscanf(ifp, "%lf", &datapoints[i].kaoncap.pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.err_pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);

      (*datacount)++;
    }

  /*  Test  */
  if ( (startpoint+lncount) != (*datacount) )
    error_exit("Wrong counting of the data points in this file!\n");

  fclose(ifp);

  return 0;
}

/*!
 * Import the diff. cs. data from the Crystal Ball experiment (iso 12)
 * for the kaon capture process. Analysis from 2008-2009.
 *
 * Contact/paper from Sergey Prakhov.
 */
int import_exp_diff_cb2009_kaoncap_12(Data datapoints[], int* datacount,
				      char* file,Observable* observ)
{
  int i, j, k, lncount, startpoint;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;
  double pk[8], err_pk[8];

  startpoint = *datacount;

  strcpy(iso_char, "12");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);

  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* headerlength = 66; */

  /* look for first row with the values of the kaon momenta */
  strcpy(dump,"start");
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"pK"));

  /* read the eight values of pK */
  for (i=0; i<8; i++)
    fscanf(ifp, "%lf", &(pk[i]));

  /* dump the err_pk string */
  fscanf(ifp, "%s", dump);

  /* read the eight values of err_pK */
  for (i=0; i<8; i++)
    fscanf(ifp, "%lf", &(err_pk[i]));

  /* dump some more text */
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"bars."));

  /* number of datapoints in this file */
  lncount = 96;

  i = startpoint;
  for (k=0; k<12; k++)
    {
      fscanf(ifp, "%s", dump);
      for (j=0; j<8; j++)
	{
	  datapoints[i].kaoncap.iso = 12;

	  strcpy(datapoints[i].kaoncap.observable,"diffcs");

	  datapoints[i].photo_prod = 0;
	  datapoints[i].electro_prod = 0;
	  datapoints[i].kaoncapture = 1;

	  datapoints[i].kaoncap.pk = pk[j];
	  datapoints[i].kaoncap.err_pk = err_pk[j];
	  datapoints[i].kaoncap.cos = -11.0/12.0+((double)k)/6.0;
	  datapoints[i].kaoncap.err_cos = 1.0/12.0;

	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
	  i++;
	}
      i -= 8;
      fscanf(ifp, "%s", dump);
      for (j=0; j<8; j++)
	{
	  fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);
	  i++;
	}
      (*datacount) += 8;
    }

  /*  Test  */
  if ( (startpoint+lncount) != (*datacount) )
    error_exit("Wrong counting of the data points in this file!\n");

  fclose(ifp);

  return 0;
}

/*!
 * Import the total cs. data from the Crystal Ball experiment (iso 11)
 * for the kaon capture process. Analysis from 2009.
 *
 * Contact/paper from Sergey Prakhov.
 *
 * Publication: Stanislaus et. al., PRC79, 015203 (2009), 8 pages.
 *
 */
int import_exp_tot_cb2009_kaoncap_12(Data datapoints[], int* datacount,
				     char* file,Observable* observ)
{
  int i, lncount, startpoint;
  char dump[DUMP_SIZE], filename[MAX_LOCATION_STRING], iso_char[3];
  FILE* ifp;

  startpoint = *datacount;

  strcpy(iso_char, "12");

  strcpy(filename, observ->dataFolder);
  strcat(filename, "radcap/iso.");
  strcat(filename, iso_char);
  strcat(filename, "/");
  strcat(filename, file);

  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("The data file does not exist!\n");

  /* look for the ending of the text */
  strcpy(dump,"start");
  do
    fscanf(ifp, "%s", dump);
  while(strcmp(dump,"err_cs"));

  /* number of datapoints in this file */
  lncount = 8;

  for (i=startpoint; i<startpoint+lncount; i++)
    {
      datapoints[i].kaoncap.iso = 12;

      strcpy(datapoints[i].kaoncap.observable,"totcs");

      datapoints[i].photo_prod = 0;
      datapoints[i].electro_prod = 0;
      datapoints[i].kaoncapture = 1;

      fscanf(ifp, "%lf", &datapoints[i].kaoncap.pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.err_pk);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.ampli);
      fscanf(ifp, "%lf", &datapoints[i].kaoncap.error);

      (*datacount)++;
    }

  /*  Test  */
  if ( (startpoint+lncount) != (*datacount) )
    error_exit("Wrong counting of the data points in this file!\n");

  fclose(ifp);

  return 0;
}


int line_count(FILE* ifp, long pos)
{

  /* This function counts the number of lines in a file
   * starting from a certain position. */

  int lncount;
  char c;

  /* set file position indicator to a value that
   * represents pos. bytes from the beginning
   * of the file (SEEK_SET). */
  fseek(ifp, pos, SEEK_SET);

  lncount = 0;

  while((c = getc(ifp)) != EOF)
    if(c == 10) /* 10 = value of newline character \n */
      lncount++;

  return lncount;
}

///! Determine how many parameters are free
int get_nr_par(int classindex)
{
  // Born+A channels
  if(classindex >= 0 && classindex <= 3)
    return 1;
  // B- and C-channel
  else if(classindex >= 4 && classindex <= 5)
    return 2;
  // D,E,F,G-channel
  else if(classindex >= 6 && classindex <= 9)
    return 1;
  // H,I,J,L-channel
  else if(classindex >= 10 && classindex <= 13)
    return 5;
  // M,N,O,Q-channel
  else if(classindex >= 14 && classindex <= 17)
    return 5; // no more off-shell couplings for spin 5/2 particles!
  // R-channel
  else if(classindex == 18)
    return 5;
  // V-channel
  else if(classindex == 19)
    return 2;
  // W-channel
  else if(classindex == 20)
    return 2;

  else
    {
      std::stringstream ss;
      ss << "Error in get_nr_par(int class_index): " << classindex << " has not been defined./n";
      throw std::logic_error(ss.str());
    }
}

/*!
 * Sometimes there are large error bars for the energy
 * in the datafile. In those cases, we calculate an averaged
 * value over the energy interval. "gridnum" determines how narrow
 * this averaging procedure is.
 *
 * We allow two sizes of the grid. A narrow and a wide one. This option
 * has to be specified before the fitting procedure or a chi squared
 * calculation is started.*/
double get_gridnum(int narrow_grid, char* observable)
{
  if(narrow_grid) // Narrow grid
  {
    if(!strcmp(observable, "diffcs"))
      return 9.0; /* diff. cross sec. */
    else if(!strcmp(observable, "totcs"))
      return 1.0; /* tot. cross sec. */
    else if(!strcmp(observable, "P"))
      return 9.0; /* rec. pol. */
    else if(!strcmp(observable, "S"))
      return 10.0; /* pho. pol. */

  }
  else // Wide grid
  {
    if(!strcmp(observable, "diffcs"))
      return 2.0; /* diff. cross sec. */
    else if(!strcmp(observable, "totcs"))
      return 1.0; /* tot. cross sec. */
    else if(!strcmp(observable, "P"))
      return 2.0; /* rec. pol. */
    else if(!strcmp(observable, "S"))
      return 2.0; /* pho. pol. */
  }
  // non-grid-specific
//   if(!strcmp(observable, "T"))
//     return 1.0; /* tar. pol. */
//   else if(!strcmp(observable, "C_x"))
//     return 1.0; /* beam-rec. pol. */
//   else if(!strcmp(observable, "C_z"))
//     return 1.0; /* beam-rec. pol. */
//   //kaoncapture-specific:
//   else if(!strcmp(observable, "ptcs"))
//     return 1.0; /* partial tot. cross sec. */
//   else if(!strcmp(observable, "bran"))
//     return 1.0; /* branching ratio for stopped kaon,
// 	    no grid necessary, but put here for completeness. */
  return 1.0; // default

}

int setlabels(Data** datapoints, int* datacount, Observable* observ)
{
  int ilabel= 0;
  int iso, d;
  for (iso=1; iso < ISOMAX; ++iso)
  {
    for(d=0;d<datacount[iso-1];++d)
    {
      int increment =1;
      if (datapoints[iso-1][d].photo_prod)
      {
	datapoints[iso-1][d].photo.label = ilabel;
	if (datapoints[iso-1][d].photo.emax > datapoints[iso-1][d].photo.emin) // number of E-bins
	  increment *= (get_gridnum(observ->fit.narrow_grid, datapoints[iso-1][d].photo.observable)+1);
	if (!strcmp(datapoints[iso-1][d].photo.observable,"totcs"))
	  increment *= 61; // number of angular bins
	else if (!strcmp(datapoints[iso-1][d].photo.observable,"P") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"S") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"T") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"C_xp") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"C_zp") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"O_xp") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"O_zp") )
	  increment *= 2; // number of angular bins
	else if (!strcmp(datapoints[iso-1][d].photo.observable,"C_x") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"C_z") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"O_x") ||
		!strcmp(datapoints[iso-1][d].photo.observable,"O_z") )
	  increment *= 3;
      }
      else if (datapoints[iso-1][d].electro_prod)
      {
	datapoints[iso-1][d].elec.label = ilabel;
/*	if (!strcmp(datapoints[iso-1][d].photo.observable,"P") ||
	  !strcmp(datapoints[iso-1][d].photo.observable,"S") ||
	  !strcmp(datapoints[iso-1][d].photo.observable,"T"))
	  increment *= 2; // number of angular bins*/
      }
      else if (datapoints[iso-1][d].kaoncapture)
      {
	datapoints[iso-1][d].kaoncap.label = ilabel;
	increment *= get_gridnum(observ->fit.narrow_grid, datapoints[iso-1][d].kaoncap.observable); // number of pk bins
	if (!strcmp( datapoints[iso-1][d].photo.observable,"totcs") ||
	  !strcmp(datapoints[iso-1][d].photo.observable,"ptcs"))
	  increment *= 101; // number of angular bins
      }
      ilabel += increment;
    }
  }
  set_datasize(ilabel);
  return ilabel;
}
