/*!
 * \file simanneal.cpp
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

#include "version.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include "strange_func.h"
#include "fitting.h"
#include "simanneal.h"


long idum;
double tt;



int simannealing(Class ** particles, 
		 Observable* observ, Data** datapoints, 
		 int datacount[])
{
  int i, j, iter, ndim;
  char logfile[MAX_LOCATION_STRING], tmpfile[MAX_LOCATION_STRING], couplfile[MAX_LOCATION_STRING], isostr[2];
  double chibest = 100000;
  double chiset[MAXNDIM+1] = {0}, tmpvertex[MAXNDIM] = {0}; 
  double bestvertex[MAXNDIM] = {0}; 
  double vertexset[MAXNDIM+1][MAXNDIM] = {{0}};
  double temptr, ftol = FTOL;
  Class printparticles[CLASSMAX] = {{0}};
  Varinfo varinfo[CLASSMAX+1] = {{{{{0}}}}};  /* 1 additional element for formfac. */
  Limits limits[MAXNDIM] = {{0}};


  /* initialization of the random (ran1) function */
  idum = - abs(observ->fit.seed);  

  /* 
   * Make a start vertex for the simplex method. If more than one isospin
   * channel is involved in the fitting procedure, construct a start vertex
   * based on the "iso_base" of the process. If just one isospin channel 
   * is under investigation, the choice is trivial.
   */


  // Determine isospin channel to construct initial vertices
  if(observ->iso.nr_iso_channels > 1)
    observ->iso.isospin = observ->iso.iso_base;
  else
    observ->iso.isospin = observ->iso.iso_channel[0];
  

  /* **************************************************************
   * Determine which variables are allowed to vary and which have 
   * to be static. If a variable is allow to vary, it can be a free 
   * parameter or a bound one. For the latter, determine the limits.
   * *************************************************************** */
  get_variable_info(varinfo, particles[observ->iso.isospin], observ);


  /* *********************
   * Make the startvertex.
   * ********************* */
  makestartvertex(particles[observ->iso.isospin], varinfo, *observ, 
		  vertexset, limits, &ndim);


  /* **********************************************************************
   * Make an initial set of chi-squares, corresponding to the vertex set of
   * ndim+1 dimensions.
   * And determine the vertex with best chi-squared
   * ********************************************************************** */

  // loop over ndim+1 start vertices
  for(i=0; i<ndim+1; i++)
    {
      /* loop over free parameters
       * and put vertexset[i] -> tmpvertex */
      for(j=0; j<ndim; j++)
	tmpvertex[j] = vertexset[i][j];
      
      // calculate chi-squared
      chiset[i] = chifunc(tmpvertex,particles, varinfo,
			  observ, datapoints, datacount);

      // look for best chi-squared of ndim+1 start vertices
      if(chiset[i] < chibest)
	{
	  for(j=0; j<ndim; j++)
	    bestvertex[j] = tmpvertex[j];
	  chibest = chiset[i];
	}

    } // end loop over initial vertices



  /* **********************************
   * Open log.fit.* and tmp.fit.* files
   * ********************************** */
  
  strcpy(logfile, observ->outFolder);
  strcat(logfile,"log/log.fit.");
  strcat(logfile, observ->fit.fit_nr);
  strcpy(tmpfile, observ->outFolder);
  strcat(tmpfile, "log/tmp.fit.");
  strcat(tmpfile, observ->fit.fit_nr);
  
  isostr[0] = observ->iso.iso_base + 48;
  isostr[1] = '\0';
  strcpy(couplfile, observ->outFolder);
  strcat(couplfile,"log/coupl.iso.");
  strcat(couplfile, isostr);
  strcat(couplfile, ".");
  strcat(couplfile, observ->fit.fit_nr);
  
  // print first log entry
  fprintlog(logfile, *observ, particles, varinfo, ndim);
  


  /* ***************
   * Start Annealing
   * *************** */

  iter = -1; // tracks iterations ( i think )

  // Run through different temperatures
  for(temptr = observ->fit.starttemp ; iter <= 0 ; )
    {
      iter = observ->fit.iter;

      /* Search for a "chibest" with a fixed temperature. */ 
      fixtempsearch(vertexset, chiset, ndim, bestvertex, &chibest, ftol, 
		    &iter, temptr,particles, varinfo, 
		    *observ, datapoints, datacount, limits);
      
      /* 
       * If there is more than one channel used in the chi squared 
       * calculation, print the base set, if no, print the used set.
       */
      if(observ->iso.nr_iso_channels > 1)
	insertvertex(printparticles, particles[observ->iso.iso_base], 
		     varinfo, observ, bestvertex);
      else
	insertvertex(printparticles, particles[observ->iso.isospin], 
		     varinfo, observ, bestvertex);		  

      // print log entry
      fprintlogstatus(logfile, temptr, chibest, printparticles, *observ);
      
      // print intermediate vertexset (for restart)
      fprinttmpinfo(tmpfile, vertexset, ndim);
           
      // decrease temperature
      temptr = observ->fit.tempdecl * temptr;
    }



  /* *******************************
   * Print final result to log.fit.*
   * ******************************* */
 
  /* 
   * If there is more than one channel used in the chi squared 
   * calculation, print the base set, if not, print the used set.
   */
  
  if(observ->iso.nr_iso_channels > 1) 
    insertvertex(printparticles, particles[observ->iso.iso_base], 
		 varinfo, observ, bestvertex);
  else
    insertvertex(printparticles, particles[observ->iso.isospin], 
		 varinfo, observ, bestvertex);


  // print final log entry
  fprintfinallog(logfile, couplfile, temptr, chibest, printparticles, *observ);

  return 0;
  
}

/*!
 * Creates a starting vertex for the simulated annealing
 * fitting process at the base of the particles[] array and the information
 * of the varinfo[]. 
 *
 * If observ.fit.restart = 1, the file given in observ.fit.input is read in
 * an  the calculation is restarted with this values. 
 * If observ.fit.restart = 0, N elements of the vertex are generated 
 * randomly. For the first element there are two possibilities: it is also
 * generated randomly or it is the original set of coupling constants 
 * from the "coupl.iso.*" file. 
 *
 * The limits[] array is generated in order to be used by the 
 * scalesimplex() function. The insertvertex() function uses the
 * full information in the varinfo[] structure. 
 */
int makestartvertex(Class particles[], Varinfo varinfo[], Observable observ, 
		    double vertexset[][MAXNDIM], Limits limits[], int* ndim)
{
  int i, j, k, limdim, rowdim, coldim; /* loop; */
  int born_cc_stored = 0;
  double tempvertex[MAXNDIM] = {0};
  char tmp[30];
  FILE* ifp;
  
  
  *ndim = 0;
  
  // When making a fresh start
  if(!observ.fit.restart)
    {
      /*
       * Create a new vertexset starting from the values in the coupl.iso.*
       * file; all free parameters are stored in the first vertex 
       * of this set.
       */
      
      // Loop over all diagrams
      for(i=0; i<CLASSMAX; i++)
	{
	  // Loop over all particles of same CLASS
	  for(j=0; j<particles[i].particount; j++)
	    {
	      // Born channels
	      if(i <= 2 && !born_cc_stored)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][0]);
		  /* the strong coupling constant g_KYN is stored
		   * only once in vertex[], because it is the same
		   * for all born terms. See insertvertex() to 
		   * understand */
		  born_cc_stored = 1;
	      }
	      // A-channel
	      else if(i == 3)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][0]);
		}
	      // B- and C-channel
	      else if(i == 4 || i == 5)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][0]);
		  
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].I,
				     varinfo[i].partic[j][1]);
		}
	      // D,E,F,G-channel
	      else if(i >= 6 && i <= 9)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][0]);
		}
	      // H,I,J,L-channel
	      else if(i >= 10 && i <= 13)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].G,
				     varinfo[i].partic[j][0]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][1]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].X,
				     varinfo[i].partic[j][2]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].Y,
				     varinfo[i].partic[j][3]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].Z,
				     varinfo[i].partic[j][4]);
		}
	      // M,N,O,P-channel
	      else if(i >= 14 && i <= 17)
		{
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].G,
				     varinfo[i].partic[j][0]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].H,
				     varinfo[i].partic[j][1]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].X,
				     varinfo[i].partic[j][2]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].Y,
				     varinfo[i].partic[j][3]);
		  store_vertex_value(tempvertex, limits, ndim, 
				     particles[i].partic[j].Z,
				     varinfo[i].partic[j][4]);
		}
	      if(*ndim > MAXNDIM)
		error_exit("ndim > MAXNDIM!!\n");
	    
	    } // end particle loop
	} // end diagram loop
      
      
      // Store the hadronic form factor values
      if(observ.hadronformfac)
	{
	  // form factors are stored in CLASSMAX position
	  i = CLASSMAX;
	  j = 0;
	  
	  store_vertex_value(tempvertex, limits, ndim, 
			     observ.ffac.born_cutoff,
			     varinfo[i].partic[j][0]);
	  
	  store_vertex_value(tempvertex, limits, ndim, 
			     observ.ffac.res_cutoff,
			     varinfo[i].partic[j][1]);
	}
      
      
      if(*ndim > MAXNDIM)
	error_exit("ndim > MAXNDIM!!\n");
      
      
      
      rowdim = *ndim + 1;
      coldim = *ndim;
      
      
      /* Generate the other vertices */
      
      // Loop over all ndim+1 vertices
      for(i=0; i<rowdim; i++)
	{
	  /* First set (i=0) of n cc's is the original one
	   * if a randomstart is not wanted */
	  if(i == 0 && !observ.fit.randomstart)
	    {
	      // loop over all free parameters
	      for(j=0; j<coldim; j++)
		vertexset[i][j] = tempvertex[j];
	    }

	  // Other n sets of n cc's (i=1..n) are generated randomly.
	  else
	    {
	      // loop over all free parameters
	      for(j=0; j<coldim; j++)
		{ 	   
		  // generate random number between limits in varinfo
		  if(limits[j].bound)
		    {
		      vertexset[i][j] = limits[j].low 
			+ ran1(&idum)*(limits[j].up - limits[j].low);
		      
		      // Check if it is indeed between limits
		      if(vertexset[i][j] < limits[j].low ||
			 vertexset[i][j] > limits[j].up)
			error_exit("Randomly generated cc out of boundaries!\n");
		    }
		  
		  /* generate random number between 0 and 2*tempvertex[j]
		   * when parameter is unbound */
		  else
		    vertexset[i][j] = tempvertex[j] * (1.0 + ran1(&idum) -
						       ran1(&idum)); 
		  
		} // end loop over free parameters
	    } // end random selection parameter values

	} // end loop over ndim+1 vertices
    } // end if(!observ.fit.restart)

  
  // When we restart a previous fit
  else 
    {
      /* 
       * Restart a fitting process and take the values from the file 
       * in "observ.fit.input" (tmp.fit.* file).
       */
      
      ifp = fopen(observ.fit.input, "r");
      if(ifp == NULL)
	error_exit("Restart file does not exist!\n");
      
      fscanf(ifp, "%s", tmp);

      i = 0;
      j = 0;

      while(strcmp(tmp, "@"))
	{
	  vertexset[i][j++] = atof(tmp);
	  fscanf(ifp, "%s", tmp);
	}
      
      *ndim = j;
      
      for(i=1; i<(*ndim + 1); i++)
	for(j=0; j<*ndim; j++)
	  {
	    fscanf(ifp, "%s", tmp);
	    vertexset[i][j] = atof(tmp);
	  }
       
      
      /* 
       * Generate the limits[] array and check if all the values
       * from the file are between the (new) given limits.
       */  

      limdim = 0;
      
      for(i=0; i<(CLASSMAX+1); i++)
	if(i != 1 && i != 2)
	  for(j=0; j<PARTICLEMAX; j++)
	    for(k=0; k<MAXNRFVAR; k++)
	      if(varinfo[i].partic[j][k].var)
		{
		  limits[limdim].bound = varinfo[i].partic[j][k].bound;
		  limits[limdim].up = varinfo[i].partic[j][k].up;
		  limits[limdim++].low = varinfo[i].partic[j][k].low;
		}

      if(*ndim != limdim)
	error_exit("Dimensions of \"varinfo\" and \"tmp.fit\" don't match\n"); 


      
      for(i=0; i<(*ndim+1); i++)
	for(j=0; j<*ndim; j++)
	  if(limits[j].bound && (vertexset[i][j] < limits[j].low ||
				 vertexset[i][j] > limits[j].up))
	    {
	      fprintf(stderr, "Restart value %d (%lf) of set %d ",j+1, 
		      vertexset[i][j], i+1);
	      fprintf(stderr, "out of limits (%lf <-> %lf)\n", limits[j].low,
		      limits[j].up);
	      exit(1);
	    }
    }

      
  return 0;
}

/*!
 * Looks for a vertex set that gives a best chi-square at 
 * a specific "temperature" = temptr.
 * For a full explanation of this function, see: "Numerical Recipes in C" 
 * W.H. Press et al. Cambridge (1992), p. 451.
 */
void fixtempsearch(double vertexset[][MAXNDIM], double chiset[], int ndim, 
		   double bestvertex[], double *chibest, double ftol, 
		   int *iter, double temptr,Class ** particles, Varinfo varinfo[], 
		   Observable observ, Data ** datapoints, 
		   int datacount[], Limits limits[])
{
  char statusfile[MAX_LOCATION_STRING];
  int i, ihi, ilo, j, m, n, mpts = ndim+1;
  double rtol, sum=0.0, swap, chihi, chilo, chinhi, chisave, chitry, chit;
  double vertexsetsum[MAXNDIM] = {0};
  FILE* ofp;

  
  strcpy(statusfile, observ.outFolder);
  strcat(statusfile, "log/status.");
  strcat(statusfile, observ.fit.fit_nr);
  
  tt = -temptr;
  for(n=0; n < ndim; n++) 
    {
      for(sum=0.0, m=0; m < mpts; m++) 
	sum += vertexset[m][n];
      vertexsetsum[n] = sum;
    }
  
  for(;;) 
    {
      
      ofp = fopen(statusfile,"w");
      fprintf(ofp, "chi best = %.4f;  temp = %.4f;  iter = %d\n" 
	      , *chibest, temptr, *iter);
      fprintchiset(ofp, chiset, ndim);
      fclose(ofp);
      

      /*
	printf("chi best = %.2f; temp = %.3f\n", *chibest, temptr);
	printchiset(chiset, ndim);
      */

      ilo = 0;
      ihi = 1;
      chinhi = chilo = chiset[0] + tt * log(ran1(&idum));
      chihi = chiset[1] + tt * log(ran1(&idum));
      if(chilo > chihi) 
	{
	  ihi = 0;
	  ilo = 1;
	  chinhi = chihi;
	  chihi = chilo;
	  chilo = chinhi;
	}

      for(i=2; i < mpts; i++) 
	{
	  chit = chiset[i] + tt * log(ran1(&idum));
	  if(chit <= chilo) 
	    {
	      ilo = i;
	      chilo = chit;
	    }
	  if(chit > chihi) 
	    {
	      chinhi = chihi;
	      ihi = i;
	      chihi = chit;
	    } 
	  else if(chit > chinhi) 
	    {
	      chinhi = chit;
	    }
	}

      rtol = 2.0 * fabs(chihi-chilo) / (fabs(chihi) + fabs(chilo));
     
      if(rtol < ftol || *iter < 0) 
	{
	  swap = chiset[0];
	  chiset[0] = chiset[ilo];
	  chiset[ilo] = swap;
	  for(n=0; n < ndim; n++) 
	    {
	      swap = vertexset[0][n];
	      vertexset[0][n] = vertexset[ilo][n];
	      vertexset[ilo][n] = swap;
	    }
	  break;
	}

      *iter -= 2;

      chitry = scalesimplex(vertexset, chiset, vertexsetsum, ndim, bestvertex,
			    chibest, ihi, &chihi, -1.0,particles, 
			    varinfo, observ, datapoints, datacount, limits);

      /*printchiset(chiset, ndim);*/

      if(chitry <= chilo) 
	{
	  chitry = scalesimplex(vertexset, chiset, vertexsetsum, ndim, 
				bestvertex, chibest, ihi, &chihi, 2.0,
				particles, varinfo, observ, 
				datapoints, datacount, limits);

	  /*printchiset(chiset, ndim);*/
	}
      else if(chitry >= chinhi) 
	{
	  chisave = chihi;
	  chitry = scalesimplex(vertexset, chiset, vertexsetsum, ndim, 
				bestvertex, chibest, ihi, &chihi, 0.5,
 				particles, varinfo, observ, 
				datapoints, datacount, limits);
	  
	  /*printchiset(chiset, ndim);*/
	  
	  if(chitry >= chisave) 
	    {
	      for(i=0; i < mpts; i++) 
		if(i != ilo) 
		  {
		    for(j=0; j < ndim; j++) 
		      {
			vertexsetsum[j] = 0.5 * (vertexset[i][j] + 
						 vertexset[ilo][j]);
			vertexset[i][j] = vertexsetsum[j];
		      }
		    
		    chiset[i] = chifunc(vertexsetsum,particles, 
					varinfo, &observ, datapoints, 
					datacount);
		    
		    /*printchiset(chiset, ndim);*/
		  }
	      
	      *iter -= ndim;
	      
	      for(n=0; n < ndim; n++) 
		{
		  for(sum=0.0, m=0; m < mpts; m++) 
		    sum += vertexset[m][n];
		  vertexsetsum[n] = sum;
		}
	    }
	}
      else ++(*iter);
    }

}

/*!
 * Rescales a simplex. 
 *
 * For a full explanation of this function, see: "Numerical Recipes in C" 
 * W.H. Press et al. Cambrigde (1992), p. 451.
 */
double scalesimplex(double vertexset[][MAXNDIM], double chiset[], 
		    double vertexsetsum[], int ndim, double bestvertex[], 
		    double *chibest, int ihi, double *chihi, double fac, 
		    Class ** particles, 
		    Varinfo varinfo[], Observable observ, 
		    Data ** datapoints, int datacount[], 
		    Limits limits[])
{
  int j;
  double fac1, fac2, chiflu, chitry, tryvertex[MAXNDIM] = {0};
  

  fac1 = (1.0 - fac)/ndim;
  fac2 = fac1 - fac;

  for(j=0; j < ndim; j++)
    {
      tryvertex[j] = vertexsetsum[j] * fac1 - vertexset[ihi][j] * fac2;
     
      /*
       * Restriction of some of the free parameters. If they are out of their 
       * boundaries, do not change this coordinate.
       * Maybe, some other procedures can be followed here.
       */ 
     
      if(limits[j].bound && (tryvertex[j] < limits[j].low || 
			     tryvertex[j] > limits[j].up))
	  tryvertex[j] =  vertexset[ihi][j];
    }


  
  chitry = chifunc(tryvertex,particles, varinfo, &observ, 
		   datapoints, datacount);
  
  if(chitry <= *chibest) 
    {
      for(j=0; j < ndim; j++) 
	bestvertex[j] = tryvertex[j];
      *chibest = chitry;
    }
 
  chiflu = chitry - tt * log(ran1(&idum));

  while(chiflu < 0) 
    chiflu = chitry - tt * log(ran1(&idum));

  if(chiflu < *chihi) 
    {
      chiset[ihi] = chitry;
      *chihi = chiflu;
      for(j=0; j < ndim; j++) 
	{
	  vertexsetsum[j] += tryvertex[j] - vertexset[ihi][j];
	  vertexset[ihi][j] = tryvertex[j];
	}
    }


  return chiflu;
}

/*!
 * Generates a random number between 0 and 1.
 * For a full explanation of this function, see: "Numerical Recipes in C" 
 * W.H. Press et al. Cambridge (1992), p. 280.
 */
double ran1(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if(*idum <= 0 || !iy) 
    {
      if(-(*idum) < 1) 
	*idum = 1;
      else 
	*idum = -(*idum);
	
      for(j=NTAB+7; j >= 0; j--) 
	{
	  k = (*idum)/IQ;
	  *idum = IA*(*idum-k*IQ)-IR*k;
	  if(*idum < 0) 
	    *idum += IM;
	  if(j < NTAB) 
	    iv[j] = *idum;
	}

      iy = iv[0];
    }

  k =(*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;

  if(*idum < 0) 
    *idum += IM;

  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;

  if((temp=AM*iy) > RNMX) 
    return RNMX;
  else 
    return temp;

}


int fprintlog(char* statusfile, Observable observ, 
	      Class ** particles, Varinfo varinfo[], int ndim)
{
  int i, j, k, iso;
  int born_cc_plot = 0;
  int print_em_ff = 0;
  FILE* ofp;
  time_t t;
  
  t = time(NULL);
  
  ofp = fopen(statusfile, "w");

  fprintf(ofp, "\n\n*********  Log file fit %s ************\n\n", 
	  observ.fit.fit_nr);
  fprintf(ofp, "strangecalc: version %s\n", VERSION);
  fprintf(ofp, "Calculation started at: %s\n", ctime(&t));
  fprintf(ofp, "Fit to the processes: ");

  for(i=0; i<observ.iso.nr_iso_channels; i++)
    fprintf(ofp, "(%hd)", observ.iso.iso_channel[i]);
  fprintf(ofp, "\n");
  
  fprintf(ofp, "Start temperature: %.3f\n", observ.fit.starttemp);
  fprintf(ofp, "Declination of the temp by %.2f\n", observ.fit.tempdecl);
  fprintf(ofp, "Number of iteration for each temp: %d\n\n", observ.fit.iter);
  
  if(observ.fit.polweight != 0)
    fprintf(ofp, "Weight factor for polarization observables: %d\n", 
	    observ.fit.polweight);

  if(observ.photo.kin.anglestep != 0)
    fprintf(ofp, "Angle step in angular integration: %d\n", 
	    observ.photo.kin.anglestep);
    
  if(observ.fit.narrow_grid)
    fprintf(ofp, "Grid size: narrow\n");
  else 
    fprintf(ofp, "Grid size: wide\n");

  if(observ.cgln)
    fprintf(ofp, "Calculation method: CGLN decomposition\n");
  else
    fprintf(ofp, "Calculation method: Explicit diagrams\n");

  fprintf(ofp, "Model type: %s \n", observ.modelType);
  
  fprintf(ofp, "Background width in 'u' and 't' channel: ");

  if(observ.backgr_width)
    fprintf(ofp, "yes\n");
  else
    fprintf(ofp, "no\n");


  fprintf(ofp, "\nUsed data base:\n");
 

  for(i=0; i<observ.iso.nr_iso_channels; i++)
    {
      iso = observ.iso.iso_channel[i];
      
      fprintf(ofp, "\tiso channel (%hd):\n", iso);
      
      for(j=0 ; j<observ.fit.nr_experiments[iso]; j++)
	fprintf(ofp, "\t-->\t%s\n", observ.fit.database_info[iso][j]); 

      if(observ.fit.elec_diffcs[iso] == 1)
	print_em_ff = 1;
    }

  
  if(print_em_ff)
    {
      fprintf(ofp, "\nUsed e.m. form factors:\n");

      if(observ.elec.emff.nuc_gk)
	fprintf(ofp, "\tNucleon form fac: Gari-Krumpelmann\n");
      else if(observ.elec.emff.nuc_lo)
	fprintf(ofp, "\tNucleon form fac: Lomon\n");
      else if(observ.elec.emff.nuc_gm)
	fprintf(ofp, "\tNucleon, K & K* form fac: Guidal monopole\n");


      if(observ.elec.emff.kplus_david)
	fprintf(ofp, "\tK+ form fac: David\n");
      else if(observ.elec.emff.kplus_monopole)
	fprintf(ofp, "\tK+ form fac: Monopole (Lambda = %.3f)\n", 
		observ.elec.emff.kplus_monopole_cutoff);

      if(observ.elec.emff.kstar_david)
	fprintf(ofp, "\tK* form fac: David\n");
      else if(observ.elec.emff.kstar_munz)
	fprintf(ofp, "\tK* form fac: Munz\n");
      else if(observ.elec.emff.kstar_monopole)
	fprintf(ofp, "\tK+ form fac: Monopole (Lambda = %.3f)\n", 
		observ.elec.emff.kstar_monopole_cutoff);

      if(observ.elec.emff.gauge_gr)
	fprintf(ofp, "Gauge restoration by Gross-Riska\n");
      else if(observ.elec.emff.gauge_modif_ff)
	fprintf(ofp, "Gauge restoration by modified f.f.\n");



    }




  fprintf(ofp, "\nParameter info:\n");
  fprintf(ofp,   "***************\n");


  if(observ.iso.nr_iso_channels > 1)
    {
      iso = observ.iso.iso_base;
    }
  else
    {
      iso = observ.iso.iso_channel[0];
    }

 
  k = 1;
  for(i=0; i<CLASSMAX; i++)
    for(j=0; j<particles[iso][i].particount; j++)
      {
	if((i <= 2 && !born_cc_plot) || i == 3)
	  {
	    give_param_info(ofp, &k, particles[iso][i].partic[j].H,
			    varinfo[i].partic[j][0]);
	    if(i <= 2)
	      born_cc_plot = 1;
	  }
	else if(i == 4 || i == 5)
	  {
	    give_param_info(ofp, &k, particles[iso][i].partic[j].H,
			    varinfo[i].partic[j][0]);

	    give_param_info(ofp, &k, particles[iso][i].partic[j].I,
			    varinfo[i].partic[j][1]);
	  }
	else if(i >= 6 && i <= 9)
	  {
	    give_param_info(ofp, &k, particles[iso][i].partic[j].H,
			    varinfo[i].partic[j][0]);
	  }
	else if(i >= 10 && i <= 13)
	  {
	    give_param_info(ofp, &k, particles[iso][i].partic[j].G,
			    varinfo[i].partic[j][0]);

	    give_param_info(ofp, &k, particles[iso][i].partic[j].H,
			    varinfo[i].partic[j][1]);

	    give_param_info(ofp, &k, particles[iso][i].partic[j].X,
			    varinfo[i].partic[j][2]);

	    give_param_info(ofp, &k, particles[iso][i].partic[j].Y,
			    varinfo[i].partic[j][3]);

	    give_param_info(ofp, &k, particles[iso][i].partic[j].Z,
			    varinfo[i].partic[j][4]);
	  }
	
      }	//FIXME -> extend to spin 5/2 particles, and gauge invariant couplings!!!

  i = CLASSMAX;
  j = 0;

  give_param_info(ofp, &k, observ.ffac.born_cutoff, varinfo[i].partic[j][0]);

  give_param_info(ofp, &k, observ.ffac.res_cutoff, varinfo[i].partic[j][1]);


  fprintf(ofp, "Number of free parameters:\t%d\n", ndim);

  fclose(ofp);

  return 0;
  
}

/*!
 * Plots the value and some information concerning the 
 * freedom of the parameter in the fit.
 */
int give_param_info(FILE* ofp, int* count, double value, Celinfo partic)
{
  fprintf(ofp, "Param.%3d:  %12.4e   ", (*count)++, value);

  if(partic.var == 0)
    fprintf(ofp, "\t(fixed)\n");
  else if(partic.bound)
    fprintf(ofp, "\t(%.2f  <->  %.2f)\n", partic.low, partic.up);
  else
    fprintf(ofp, "\t(free)\n");
  

  return 0;
}





int printchiset(double chiset[], int ndim)
{
  int i;
  
  printf("chiset[] = ");
  
  for(i=0; i<ndim+1; i++)
    printf("%.3f  ", chiset[i]);
  
  printf("\n");

  return 0;
}


int fprintchiset(FILE* ofp, double chiset[], int ndim)
{
  int i;
  
  fprintf(ofp, "chiset[] = ");
  
  for(i=0; i<ndim+1; i++)
    fprintf(ofp, "%.3f  ", chiset[i]);
  
  fprintf(ofp, "\n");
 
  return 0;
  
}

/*!
 * Prints the status of the fitting process to the logfile.
 */
int fprintlogstatus(char* logfile, double temptr, double chibest, Class 
		    printparticles[], Observable observ)
{
  FILE* ofp;
  

  ofp = fopen(logfile, "a");
  
  fprintf(ofp,"\nThe chi squared by temperature %6.2f is:\n", temptr);
  fprintf(ofp,"%f\n", chibest);
  fprintf(ofp,"The coupling constants producing this chi squared are:\n");
  

 
  fprintparticspec(ofp, printparticles, observ);
  
  fclose(ofp);      
  
  
  return 0;
}

/*!
 * Prints the FINAL status of the fitting process to the logfile.
 */
int fprintfinallog(char* logfile, char* couplfile, double temptr, double chibest, Class 
		   printparticles[], Observable observ)
{
  FILE* ofp;
  FILE* ccf;
  
  ofp = fopen(logfile, "a");
  ccf = fopen(couplfile, "w");
  
  fprintf(ofp,"\nThe FINAL chi squared by temperature %6.2f is:\n", temptr);
  fprintf(ofp,"%f\n", chibest);
  fprintf(ofp,"The coupling constants producing this chi squared are:\n");
  fprintparticspec(ofp, printparticles, observ);

  
  fprintf(ccf, "Class  NickName   G*           H*          I*          ");
  fprintf(ccf, "X*          Y*          Z*\n");
  fprintf(ccf, "-----  --------   --           --          --          ");
  fprintf(ccf, "--          --          --\n");
  fprintparticspec(ccf, printparticles, observ);
  
  return 0;
}





int fprinttmpinfo(char* tmpfile, double vertexset[][MAXNDIM], int ndim)
{

  int i, j;
  FILE* ofp;
  
  
  ofp = fopen(tmpfile, "w");
  

  for(i=0; i<ndim+1; i++)
    {
      for(j=0; j<ndim; j++)
	fprintf(ofp, "%15.6e", vertexset[i][j]);
      
      if(i == 0)
	fprintf(ofp, "  @\n");
      else
	fprintf(ofp, "   \n");
    }
  

  fclose(ofp);
  
  return 0;
}





