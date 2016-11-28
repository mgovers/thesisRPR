/*!
 * \file io_specific.cpp
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

#include "version.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include "strange_func.h"
#include "io_specific.h"


/*!
 * \brief Initializes the "observ" structure.
 * The structure "observ" contains all the information of the 
 * process and the observables that are calculated.
 */
int observablespecification(Observable* observ, FILE* ifp, FILE* ofp1, 
			    FILE* ofp2)
{
  char ans[DUMP_SIZE];
     
  

  fprintf(ofp1,"\n\n");
  fprintf(ofp1,"\t::::::::::::::::::::::::::: \n\n");
  fprintf(ofp1,"\t STRANGECALC (version %s)\n\n",VERSION);
  fprintf(ofp1,"\t::::::::::::::::::::::::::: \n\n");

  /*
   * Specify the input and output folder
   */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Location of the input folder? \t\t\t: ");
  fscanf(ifp,"%s", observ->inFolder);
  fprintf(ofp2,"%s\n", observ->inFolder);

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Location of the output folder? \t\t\t: ");
  fscanf(ifp,"%s", observ->outFolder);
  fprintf(ofp2,"%s\n", observ->outFolder);

  /*
   * Specify the data folder
   */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Location of the data folder? \t\t\t: ");
  fscanf(ifp,"%s", observ->dataFolder);
  fprintf(ofp2,"%s\n", observ->dataFolder);

  /* 
   * Calculation type
   */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Specific Calculation, Fitting procedure \n");
  fprintf(ofp1,"or chi Squared calc.? (c/f/s)\t\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);


  if(!strcmp(ans, "f"))
     {
       /* fitting */

       observ->fitting = 1;
       observ->chisquare = 0;
      
       isospinspecification(observ, ifp, ofp1, ofp2);
       fittingspecification(observ, ifp, ofp1, ofp2);
     } 


  else if(!strcmp(ans, "s"))
    {
      /* chi-squared calculation */

      observ->chisquare = 1;
      observ->fitting = 0;
      
      isospinspecification(observ, ifp, ofp1, ofp2);
      chicalcspecification(observ, ifp, ofp1, ofp2);
    }


  else if(!strcmp(ans, "c"))
    {
      error_exit("This option is depricated.\n");
    }
  else 
    error_exit("Answer c(alc.), f(itting) procedure or (chi) s(quared)\n");


  /* 
   * Regge specifications 
   */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Isobar theory or Regge theory? (i/r)\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);
  

  if(!strcmp(ans, "i"))
    {
      /*-------------
       * Isobar model 
       *-------------*/

      observ->regge = 0;

      /* 
       * Width in u- and t-channel?
       */
     
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Zero or finite width in u- and t- channel? (0/f): ");
      fscanf(ifp,"%s", ans);
      fprintf(ofp2,"%s\n", ans);

      if(!strcmp(ans, "0"))
	observ->backgr_width = 0;
      else if(!strcmp(ans, "f"))
	observ->backgr_width = 1;
      else
	error_exit("Answer 0 or f(inite)!!\n");
    }


  else if(!strcmp(ans, "r"))
    {
      /*------------
       * Regge model
       *------------*/

      observ->regge = 1;

      reggespecification(observ, ifp, ofp1, ofp2);

    }
  else
    error_exit("Answer i(sobar) or r(egge)\n");


  /*
   * Form factor specifications 
   */


  /* In electroproduction, there are EM ff's */
  // HAS CHANGED !!!!
  /* EM form factors used to be specified here (in Stijn's and Tamara's code)
   * This has changed. They are now specified in an input file "formfactors".
   * This file is read in and stored in the particles[] array later on in numeric()
   */
  
  /* Unless we work in the pure Regge model, hadronic ffs can be added */
  if(!observ->regge || observ->reg.res_on_bg)
    hadronformfacspecification(observ, ifp, ofp1, ofp2);
  



  /* 
   * PS or PV coupling for spin 1/2 resonances
   */  

  observ->pvcoupling = 0; 

  /*
    fprintf(ofp1,"---\n");
    fprintf(ofp1,"PV-coupling for the spin 1/2 resonances? (y/n)\n");
    fscanf(ifp,"%s", ans);
    fprintf(ofp2,"%s", ans);
    
    if(!strcmp(ans, "y"))
    observ->pvcoupling = 1;  
    else if(!strcmp(ans, "n"))
    observ->pvcoupling = 0;  
    else 
    error_exit("Answer y(es) or n(o)!\n"); 
  */


    
  /* 
   * Gauge-invariant or offshell couplings
   */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,
	  "Enter model type (for old model choose: (n) offshell  model (y) consistent couplings) \t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp1,"%s\n",ans);
  
  // Full backward compatibility.
  observ->cgln =0;
  if(!strcmp(ans, "n"))
    strcpy(observ->modelType,"offshell");
  else if(!strcmp(ans, "y"))
    strcpy(observ->modelType,"consistent");
  else 
  {
    observ->cgln=1;
    strcpy(observ->modelType,ans);
  }
  
  
  fprintf(ofp1,"\n");
  
  return 0;

}




/*!
 * Determines the isospin specifications from "ifp"
 * and prints the messages to "ofp1". 
 */
int isospinspecification(Observable* observ, FILE* ifp, FILE* ofp1, FILE* ofp2)
{
  int i=0,j,k;
  char ans[DUMP_SIZE];
  char anstmp[DUMP_SIZE];

  if(observ->fitting || observ->chisquare)
    {	


      /*-----------------------
       * Fitting or chi squared
       *-----------------------*/


      /*
       * Determine isospin channels to be included 
       */

      fprintf(ofp1,"---\n");
      fprintf(ofp1,"\t(1)   ===>   g + p   ->  K+ + L0 \n");
      fprintf(ofp1,"\t(2)   ===>   g + p   ->  K+ + S0 \n");
      fprintf(ofp1,"\t(3)   ===>   g + p   ->  K0 + S+ \n");
      fprintf(ofp1,"\t(4)   ===>   g + n   ->  K0 + L0 \n");
      fprintf(ofp1,"\t(5)   ===>   g + n   ->  K0 + S0 \n");
      fprintf(ofp1,"\t(6)   ===>   g + n   ->  K+ + S- \n");
      fprintf(ofp1,"\t(7)   ===>   g + p   ->  P0 + p  \n");
      fprintf(ofp1,"\t(8)   ===>   g + p   ->  P+ + n  \n");
      fprintf(ofp1,"\t(9)   ===>   g + n   ->  P0 + n  \n");
      fprintf(ofp1,"\t(10)  ===>   g + n   ->  P- + p  \n");
      fprintf(ofp1,"\t(11)  ===>   K- + p  ->  g  + L0 \n");
      fprintf(ofp1,"\t(12)  ===>   K- + p  ->  g  + S0 \n");
      fprintf(ofp1,"Give the numbers of the processes in the fit\n");
      fprintf(ofp1,"in increasing order.\n");
      fprintf(ofp1,"Combined fit only possible for processes derived from same base channel!\n");
      fprintf(ofp1,"(multiple nrs as one line separated by \",\")\t: ");


      fscanf(ifp,"%s", ans);
      fprintf(ofp2,"%s\n", ans);
  


      /* observ->iso.iso_channel[] array (comp's=0,1,...iso.nr_iso_channels-1)
       * contains isospin channel numbers, from 1 to 12 (numbers
       * must be given in the right order!!!) */

       observ->iso.nr_iso_channels = 1;
      for(i=0; ans[i] != '\0'; i++)
	if (ans[i] == ',')
	  observ->iso.nr_iso_channels++;
      k=0;
      j=0;
      i=0;
      do{
	if(isdigit(ans[i])) 
	  {
	    anstmp[j++] = ans[i];
	  }
	else
	  {
	    anstmp[j] = '\0';
	    j = 0;
	    observ->iso.iso_channel[k++] = (short) atoi(anstmp);
	  }
      }while(ans[i++] != '\0');

      if (k != observ->iso.nr_iso_channels)
	error_exit("Wrong counting of number of isospin channels!\n");


      /* observ->iso.iso_base is either 1 (g p -> K+ L) 
       *                             or 2 (g p -> K+ S0)
       *			     or 7 (g p -> P0 p)
       *			     or 8 (g p -> P+ n)
       * observ->iso.isospin can be 1...12, but is set equal to iso_base unless
       * only a non-base-channel is selected */

      observ->iso.isospin = observ->iso.iso_channel[0];
      observ->iso.iso_base = getisospinbasechannel(observ->iso.iso_channel[0]);

  
      /* WARNING: only for processes derived from the same base channel 
       * (ie same cc's), a combined fit can be made */
      
      for(i=1; i<observ->iso.nr_iso_channels; ++i) {
	if( getisospinbasechannel(observ->iso.iso_channel[i]) != observ->iso.iso_base) {
	  fprintf(stderr,"ERROR in isospinspecification(..): in fit or chi squared calc. all isospin channels should be of the same base channel!\n");
	  exit(1);
	}
      }

    }
  else
    {

      /*------------------------------------------
       * Calculation with given coupling constants
       *------------------------------------------*/

      /* 
       * Select (only one!) isospin channel 
       */

      fprintf(ofp1,"---\n");
      fprintf(ofp1,"\t(1)   ===>   g + p   ->  K+ + L0 \n");
      fprintf(ofp1,"\t(2)   ===>   g + p   ->  K+ + S0 \n");
      fprintf(ofp1,"\t(3)   ===>   g + p   ->  K0 + S+ \n");
      fprintf(ofp1,"\t(4)   ===>   g + n   ->  K0 + L0 \n");
      fprintf(ofp1,"\t(5)   ===>   g + n   ->  K0 + S0 \n");
      fprintf(ofp1,"\t(6)   ===>   g + n   ->  K+ + S- \n");
      fprintf(ofp1,"\t(7)   ===>   g + p   ->  P0 + p  \n");
      fprintf(ofp1,"\t(8)   ===>   g + p   ->  P+ + n  \n");
      fprintf(ofp1,"\t(9)   ===>   g + n   ->  P0 + n  \n");
      fprintf(ofp1,"\t(10)  ===>   g + n   ->  P- + p  \n");
      fprintf(ofp1,"\t(11)  ===>   K- + p  ->  g + L0  \n");
      fprintf(ofp1,"\t(12)  ===>   K- + p  ->  g + S0  \n");
      fprintf(ofp1,"Give a process number\t\t\t\t: ");

      fscanf(ifp,"%hd", &observ->iso.isospin);
      fprintf(ofp2,"%hd\n", observ->iso.isospin);

      if(observ->iso.isospin < 1 || observ->iso.isospin > 12)
	error_exit("Answer with number between 1 and 12!\n");
     

      observ->iso.nr_iso_channels = 1;


      /* determine base isospin channel (1,2,7,8) */
      observ->iso.iso_base = getisospinbasechannel(observ->iso.isospin);
	
      if(observ->iso.isospin == 11 || observ->iso.isospin == 12 ) {
	// observ->iso.isospin = observ->iso.iso_base; if have no idea why this line was here?? PVC
	observ->kaoncapture = 1;
      }
    }
  
  return 0;
}

/*!
 * Determines the fitting specifications.
 */
int fittingspecification(Observable* observ, FILE* ifp, FILE* ofp1, FILE* ofp2)
{
  char ans[DUMP_SIZE];

  /* will be changed in chicalcspecification() */
  observ->photoprod = 1;


  /* Start from *.log file produced by a previous (ev. aborted) fit? */
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Restart a previous process? (y/n)\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);


  if(!strcmp(ans, "n"))
    {
      observ->fit.restart = 0;

      /* Start with random cc's (following varinfo.*), or use coupl.iso.* ? */
      fprintf(ofp1,"---\n");      
      fprintf(ofp1,"Complete random start? (y/n)\t\t\t: ");
      fscanf(ifp,"%s", ans);
      fprintf(ofp2,"%s\n", ans);
      
      
      if(!strcmp(ans, "y"))
	observ->fit.randomstart = 1;
      else if(!strcmp(ans, "n"))
	observ->fit.randomstart = 0;
      else 
	error_exit("Answer yes or no!\n");

    }
  else if(!strcmp(ans, "y"))
    {
      observ->fit.restart = 1;
      
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Give the name of the input file (full path)\t: ");
      fscanf(ifp,"%s", observ->fit.input);
      fprintf(ofp2,"%s\n", observ->fit.input);

      fprintf(ofp1,"Specify the restarted process...\n");
    }
  else
    error_exit("Answer y(es) or n(o)!!\n");
  
            
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Fit number? (no longer than 9 characters!) \t\t\t\t\t: ");
  fscanf(ifp,"%s", observ->fit.fit_nr);
  fprintf(ofp2,"%s\n", observ->fit.fit_nr);
  
  
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Start temperature?\t\t\t\t: ");
  fscanf(ifp,"%lf", &observ->fit.starttemp);
  fprintf(ofp2,"%lf\n", observ->fit.starttemp);
  
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Declination of temperature?\t\t\t: ");
  fscanf(ifp,"%lf", &observ->fit.tempdecl);
  fprintf(ofp2,"%lf\n", observ->fit.tempdecl);
  
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Number of iterations by every temp?\t\t: ");
  fscanf(ifp,"%d", &observ->fit.iter);
  fprintf(ofp2,"%d\n", observ->fit.iter);
  
  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Integer random seed?\t\t\t\t: ");
  fscanf(ifp,"%d", &observ->fit.seed);  
  fprintf(ofp2,"%d\n", observ->fit.seed);  



  chicalcspecification(observ, ifp, ofp1, ofp2);

 
  return 0;
}

  
/*!
 * Determines the specifications to calculate a chi squared.
 */  
int chicalcspecification(Observable* observ, FILE* ifp, FILE* ofp1, FILE* ofp2)
{
  int i;
  short tot_angle = 0, pol_data = 0, e_grid = 0;
  char ans[4];

  databasespecification(observ, &tot_angle, &pol_data, &e_grid, 
			ifp, ofp1, ofp2);

  if (observ->iso.nr_iso_channels > 1)
    {
      for(i=0; i<observ->iso.nr_iso_channels; ++i) {
	  fprintf(ofp1,"---\n");
	  fprintf(ofp1,"Weight factor for iso-channel %d?\t\t\t: ",
		  observ->iso.iso_channel[i]);
	  fscanf(ifp,"%hd", &observ->fit.isoWeight[observ->iso.iso_channel[i]]);
	  fprintf(ofp2,"%hd\n", observ->fit.isoWeight[observ->iso.iso_channel[i]]);

	  if( observ->fit.isoWeight[observ->iso.iso_channel[i]] <= 0.0 ) 
	    error_exit("Weight factor should be strictly positive!\n");
      }
    }
  else
    {
      observ->fit.isoWeight[observ->iso.iso_channel[0]] = 1.0;
    }  
    
  if(tot_angle)
    {
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Give theta integration step (in deg.)\t\t: ");
      fscanf(ifp,"%hd", &observ->photo.kin.anglestep);
      fprintf(ofp2,"%hd\n", observ->photo.kin.anglestep);
      /* in case of kaon capture, we use the same anglestep */
      observ->kaoncap.kin.anglestep = observ->photo.kin.anglestep;
    }
  
  if(pol_data)
    {
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Weight factor for polarization?\t\t\t: ");
      fscanf(ifp,"%hd", &observ->fit.polweight);
      fprintf(ofp2,"%hd\n", observ->fit.polweight);
    }

  if(e_grid)
    {
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Narrow or wide grid to average? (n/w)\t\t: ");
      //fscanf(ifp,"%s", ans);
      strcpy(ans,"n");
      fprintf(ofp2,"%s\n", ans);

      if(!strcmp(ans, "n"))
	observ->fit.narrow_grid = 1;
      else if(!strcmp(ans, "w"))
	observ->fit.narrow_grid = 0;
      else 
	error_exit("Answer n(arrow) or w(ide)!\n");
    }

  return 0;
}

  
/*!
 * Determines the databases used in the chi squared calculation. For
 * each isospin channel included in the fit, several datasets can be 
 * selected. The names of these sets are stored in database_info[iso][exp].
 */
int databasespecification(Observable* observ, short* tot_angle, 
			  short* pol_data, short* e_grid, 
			  FILE* ifp, FILE* ofp1, FILE* ofp2)
{
  char exp[5], ans[400];
  int i, j, k, l, iso, experiment;


  /* why ??? */
  observ->photoprod = 1;



  /* For each isospin channel, determine the datasets to which the results
   * should be fitted */
  	  
  for(i=0; i<observ->iso.nr_iso_channels; i++)
    {      
      iso = observ->iso.iso_channel[i];

      
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Give data points to be used in chi^2 calculation ");
      fprintf(ofp1,"for isospin (%d):\n", iso); 

      if(iso == 1)
	{
	  /*-------------------------------
	   * Isopin channel 1 (g p -> K+ L)
	   *-------------------------------*/

	  fprintf(ofp1,
		  "Medium-energy photoproduction data\n");
	  		  
	  fprintf(ofp1,
		  "\t(101)\t-->\tsaphir '98\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(102)\t-->\tsaphir '98\t  -- photo --\ttot. cross sec. \n");

	  fprintf(ofp1,
		  "\t(103)\t-->\tsaphir '98\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(104)\t-->\tsaphir '03\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(105)\t-->\tforw. saphir '03\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(106)\t-->\tsaphir '03\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(107)\t-->\tforw. saphir '03\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(108)\t-->\tclas-cmu '03\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(109)\t-->\tforward clas '03  -- photo ");
	  fprintf(ofp1,"--\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(110)\t-->\tclas-cmu '03\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(111)\t-->\tforward clas '03 \n");
	  fprintf(ofp1,"\t\t\t(cos > -0.15)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(112)\t-->\tforward clas '03\n");
	  fprintf(ofp1, "\t\t\t(cos > 0.35)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(113)\t-->\tclas '05\t  -- photo ");
	  fprintf(ofp1,"--\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(114)\t-->\tforward clas '05\n");
	  fprintf(ofp1,"\t\t\t(cos > -0.15)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(115)\t-->\tforward clas '05\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.35)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(116)\t-->\tclas '07\t  -- photo ");
	  fprintf(ofp1,"--\tbeam-rec. asym (c)\n");

	  fprintf(ofp1,
		  "\t(117)\t-->\tspring-8 '03\t  -- photo --\tpho. asym.\n");

	  fprintf(ofp1,
		  "\t(118)\t-->\tspring-8 '04\t  -- photo --\tdiff.cross sec.\n");

	  fprintf(ofp1,
		  "\t(119)\t-->\tgraal '06\t  -- photo --\tpho. asym.\n");
	  fprintf(ofp1,
		  "\t(120)\t-->\tforward graal '06\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\tpho. asym.\n");

	  fprintf(ofp1,
		  "\t(121)\t-->\tgraal '06\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(122)\t-->\tforward graal '06\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(123)\t-->\tbonn '78\t  -- photo --\ttar. asym. \n");

	  fprintf(ofp1,
		  "\t(124)\t-->\tLEPS '07\t  -- photo --\tdiff.cross sec. (u)\n");

	  fprintf(ofp1,
		  "\t(125)\t-->\tLEPS '07\t  -- photo --\tdiff.cross sec. (costhkcm)\n");	  

	  fprintf(ofp1,
		  "\t(126)\t-->\tLEPS '07\t  -- photo --\tpho. asym. \n");

	  fprintf(ofp1,
		  "\t(127)\t-->\tGRAAL '08\t  -- photo --\ttar. asym. \n");

	  fprintf(ofp1,
		  "\t(128)\t-->\tGRAAL '08\t  -- photo --\tbeam-rec. asym. (o) \n");
// positive costhkcm data 	  
	  fprintf(ofp1,
		  "\t(129)\t-->\tforward clas '03 \n \t\t\t (cos > 0.0) \t -- photo -- \t rec. asym.\n");
	  fprintf(ofp1,
		  "\t(130)\t-->\tforward clas '05 \n \t\t\t (cos > 0.0) \t -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(131)\t-->\tforward clas '07 \n \t\t\t (cos > 0.0) \t -- photo -- \t beam-rec. asym (c)\n");
	  fprintf(ofp1,
		  "\t(132)\t-->\tGRAAL '08 \n \t\t\t (cos > 0.0) \t -- photo -- \t tar. asym.\n");
	  fprintf(ofp1,
		  "\t(133)\t-->\tGRAAL '08 \n \t\t\t (cos > 0.0) \t -- photo -- \t beam-rec. asym. (o)\n");
	  fprintf(ofp1,
		  "\t(134)\t-->\tCLAS '09 \n \t\t\t \t -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(135)\t-->\tCLAS '09 \n \t\t\t \t -- photo -- \t rec. asym. (o)\n");
	  fprintf(ofp1,
		  "\t(136)\t-->\tforward CLAS '09 \n \t\t\t (cos > 0.0) \t -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(137)\t-->\tforward CLAS '09 \n \t\t\t (cos > 0.0) \t -- photo -- \t rec. asym. (o)\n");
	  fprintf(ofp1,
		  "\t(138)\t-->\tex fwd CLAS '09 \n \t\t\t (cos > 0.35) \t -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(139)\t-->\tex fwd CLAS '09 \n \t\t\t (cos > 0.35) \t -- photo -- \t rec. asym. (o)\n");


	  fprintf(ofp1,
		  "Medium-energy electroproduction data\n");

	  fprintf(ofp1,
		  "\t(201)\t-->\tcea '72 \t  -- elec -- \tT+L cross sec.\n");

	  fprintf(ofp1,
		  "\t(202)\t-->\tharv.-corn. '74\t  -- elec -- \tT+L cross ");
	  fprintf(ofp1,"sec.\n");

	  fprintf(ofp1,
		  "\t(203)\t-->\tharv.-corn. '77\t  -- elec -- \tT+L cross ");
	  fprintf(ofp1,"sec.\n");

 	  fprintf(ofp1,
		  "\t(204)\t-->\tjlab '03\t  -- elec -- \tT & L cross sec.\n");

 	  fprintf(ofp1,
		  "\t(205)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (2.567 GeV, narrow bins)\n");

 	  fprintf(ofp1,
		  "\t(206)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (2.567 GeV, wide bins)\n");

 	  fprintf(ofp1,
		  "\t(207)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (4.056 GeV)\n");

 	  fprintf(ofp1,
		  "\t(208)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT, L, TT, LT cross sec.\n");

	  fprintf(ofp1,
		  "\t(209)\t-->\tjlab '07\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tLT' cross sec.\n");

	  fprintf(ofp1,
		  "\t(210)\t-->\tCLAS '09\t  -- elec -- \tT+L cross ");
	  fprintf(ofp1,"sec.\n");

	  fprintf(ofp1,
		  "\t(211)\t-->\tCLAS '09\t  -- elec -- \tT cross ");
	  fprintf(ofp1,"sec.\n");

	  fprintf(ofp1,
		  "\t(212)\t-->\tCLAS '09\t  -- elec -- \tL cross ");
	  fprintf(ofp1,"sec.\n");

	  fprintf(ofp1,
		  "\t(213)\t-->\tjlab '09\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\ttransf. pol (unprimed ref. frame, 4.261 GeV)\n");

	  fprintf(ofp1,
		  "\t(214)\t-->\tjlab '09\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\ttransf. pol (unprimed ref. frame, 5.754 GeV)\n");

	  fprintf(ofp1,
		  "\t(215)\t-->\tjlab '09\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tratio longitudinal/transverse\n");

	  fprintf(ofp1,
		  "\t(216)\t-->\tjlab '00\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L cross sec (Ph.D. Cha, unpublished)\n");
	  fprintf(ofp1,
		  "\t(217)\t-->\tCLAS '13\t  -- elec -- \tU, LT, TT cross sections.\n");
	  fprintf(ofp1,
		  "\t(218)\t-->\tCLAS '13\t  -- elec -- \tLT' cross sections.\n");
	  fprintf(ofp1,
		  "\t(219)\t-->\tCLAS '14\t  -- elec -- \tinduced Lambda polarization.\n");


	  fprintf(ofp1,"High-energy photoproduction data\n");

	  fprintf(ofp1,
		  "\t(401)\t-->\tslac '69\t  -- photo --\tdiff. ");
	  fprintf(ofp1,"cross sec. (t)\n");

	  fprintf(ofp1,
		  "\t(402)\t-->\tslac '69\t  -- photo --\tdiff. ");
	  fprintf(ofp1,"cross sec. (u)\n");

	  fprintf(ofp1,
		  "\t(403)\t-->\tslac '79\t  -- photo --\tpho.");
	  fprintf(ofp1,"asym. (t)\n");

	  fprintf(ofp1,
		  "\t(404)\t-->\tdesy '72\t  -- photo --\trec.");
	  fprintf(ofp1,"asym. (t)\n");
 
	  fprintf(ofp1,
		  "\t(405)\t-->\tex fwd CLAS '09 (W > 2.6 GeV, cos > 0.35)\n \t\t\t\t\t  -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(406)\t-->\tex fwd CLAS '09 (W > 2.6 GeV, cos > 0.35)\n \t\t\t\t\t  -- photo -- \t rec. asym. (o)\n");
	  fprintf(ofp1,"High-energy electroproduction data (W > 2.5 GeV, -t < 1 GeV^2)\n");
	  fprintf(ofp1,
		  "\t(501)\t-->\tCornell '74\t  -- elec -- \tunseparated cross sections.\n");
	  fprintf(ofp1,
		  "\t(502)\t-->\tCornell '77\t  -- elec -- \tunseparated cross sections.\n");
	  fprintf(ofp1,
		  "\t(503)\t-->\tCLAS '13\t  -- elec -- \tU, LT, TT cross sections.\n");
	  fprintf(ofp1,
		  "\t(504)\t-->\tCLAS '13\t  -- elec -- \tLT' cross sections.\n");
	  fprintf(ofp1,
		  "\t(505)\t-->\tCLAS '14\t  -- elec -- \tinduced Lambda polarization.\n");
	}	          					     


      else if(iso == 2)
	{
	 
	  /*--------------------------------
	   * Isopin channel 2 (g p -> K+ S0)
	   *--------------------------------*/

	  fprintf(ofp1,"Medium-energy photoproduction data\n");	

	  fprintf(ofp1,
		  "\t(101)\t-->\tsaphir '98\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(102)\t-->\tsaphir '98\t  -- photo --\ttot. cross sec. \n");

	  fprintf(ofp1,
		  "\t(103)\t-->\tsaphir '98\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(104)\t-->\tsaphir '03\t  -- photo --\tdiff. cross sec. \n");

	  fprintf(ofp1,
		  "\t(105)\t-->\tforw. saphir '03\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\tdiff. cross sec. \n");

	  fprintf(ofp1,
		  "\t(106)\t-->\tsaphir '03\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(107)\t-->\tforw. saphir '03\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\trec. asym.\n");
	  
	  fprintf(ofp1,
		  "\t(108)\t-->\tclas-cmu '03\t  -- photo --\tdiff. cross sec. \n");

	  fprintf(ofp1,
		  "\t(109)\t-->\tforward clas '03  -- photo ");
	  fprintf(ofp1,	"--\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(110)\t-->\tclas-cmu '03\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(111)\t-->\tforward clas '03\n");
	  fprintf(ofp1,"\t\t\t(cos > -0.15)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(112)\t-->\tforward clas '03\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.35)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(113)\t-->\tclas '05\t  -- photo ");
	  fprintf(ofp1,	"--\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(114)\t-->\tforward clas '05\n");
	  fprintf(ofp1,	"\t\t\t(cos > -0.15)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(115)\t-->\tforward clas '05\n");
	  fprintf(ofp1,	"\t\t\t(cos > 0.35)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(116)\t-->\tclas '07\t  -- photo ");
	  fprintf(ofp1,"--\tbeam-rec. asym (c)\n");

	  fprintf(ofp1,
		  "\t(117)\t-->\tspring-8 '03\t  -- photo --\tpho. asym.\n");

	  fprintf(ofp1,
		  "\t(118)\t-->\tspring-8 '04\t  -- photo --\tdiff.cross sec.\n");

	  fprintf(ofp1,
		  "\t(119)\t-->\tgraal '06\t  -- photo --\tpho. asym.\n");

	  fprintf(ofp1,
		  "\t(120)\t-->\tforward graal '06\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\tpho. asym.\n");


	  fprintf(ofp1,
		  "\t(121)\t-->\tgraal '06\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(122)\t-->\tforward graal '06\n");
	  fprintf(ofp1,"\t\t\t(cos > 0.0)\t  -- photo --\trec. asym.\n");

	  fprintf(ofp1,
		  "\t(123)\t-->\tleps '06\t  -- photo --\tdiff.cross sec.\n");

	  fprintf(ofp1,
		  "\t(124)\t-->\tleps '06\t  -- photo --\tpho. asym.\n");
// positive costhkcm data 	  
	  fprintf(ofp1,
		  "\t(129)\t-->\tforward clas '03 \n \t\t\t (cos > 0.0) \t  -- photo -- \trec. asym.\n");
	  fprintf(ofp1,
		  "\t(130)\t-->\tforward clas '05 \n \t\t\t (cos > 0.0) \t  -- photo -- \tdiff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(131)\t-->\tforward clas '07 \n \t\t\t (cos > 0.0) \t  -- photo -- \tbeam-rec. asym (c)\n");
	  fprintf(ofp1,
		  "\t(132)\t-->\tclas '10\t  -- photo --\tdiff.cross sec.\n");
	  fprintf(ofp1,
		  "\t(133)\t-->\tforward clas '10 \n \t\t\t (cos > 0.0) \t  -- photo -- \tdiff.cross sec.\n");
	  fprintf(ofp1,
		  "\t(134)\t-->\tclas '10\t  -- photo --\trec. asym.\n");
	  fprintf(ofp1,
		  "\t(135)\t-->\tforward clas '10 \n \t\t\t (cos > 0.0) \t  -- photo -- \trec. asym.\n");

	  fprintf(ofp1,"Medium-energy electroproduction data\n");	

	  fprintf(ofp1,
		  "\t(201)\t-->\tcea '72   \t  -- elec -- \tT+L cross sec. \n");

	  fprintf(ofp1,
		  "\t(202)\t-->\tharv.-corn. '74\t  -- elec -- \tT+L cross ");
	  fprintf(ofp1,"sec.\n");

	  fprintf(ofp1,
		  "\t(203)\t-->\tharv.-corn. '77\t  -- elec -- \tT+L cross ");
	  fprintf(ofp1,"sec.\n");

 	  fprintf(ofp1,
		  "\t(204)\t-->\tjlab '03\t  -- elec -- \tT & L cross sec.\n");

 	  fprintf(ofp1,
		  "\t(205)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (2.567 GeV, narrow bins)\n");

 	  fprintf(ofp1,
		  "\t(206)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (2.567 GeV, wide bins)\n");

 	  fprintf(ofp1,
		  "\t(207)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L, TT, LT cross sec. (4.056 GeV)\n");

 	  fprintf(ofp1,
		  "\t(208)\t-->\tjlab '06\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT, L, TT, LT cross sec.\n");

	  fprintf(ofp1,
		  "\t(209)\t-->\tjlab '09\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\ttransf. pol (unprimed ref. frame)\n");

	  fprintf(ofp1,
		  "\t(210)\t-->\tjlab '00\t  -- elec --\n");
	  fprintf(ofp1,"\t\t\tT+L cross sec (Ph.D. Cha, unpublished)\n");


	  fprintf(ofp1,"High-energy data\n");

	  fprintf(ofp1,
		  "\t(401)\t-->\tslac '69\t  -- photo --\tdiff. cross ");
	  fprintf(ofp1,"sec. (t)\n");

	  fprintf(ofp1,
		  "\t(402)\t-->\tslac '69\t  -- photo --\tdiff. cross ");
	  fprintf(ofp1,"sec. (u)\n");

	  fprintf(ofp1,
		  "\t(403)\t-->\tslac '79\t  -- photo --\tpho. ");
	  fprintf(ofp1,"asym. (t)\n");
	  fprintf(ofp1,
		  "\t(404)\t-->\tex fwd clas '10 (W > 2.6 GeV)\n \t\t\t (cos > 0.35) \t  -- photo -- \t diff. cross sec.\n");
	  fprintf(ofp1,
		  "\t(405)\t-->\tex fwd clas '10 (W > 2.6 GeV)\n \t\t\t (cos > 0.35) \t  -- photo -- \t rec. asym. (o)\n");
	}	


      else if(iso == 3)					     
	{	  
	  /*--------------------------------
	   * Isopin channel 3 (g p -> K0 S+)
	   *--------------------------------*/        					     
	  fprintf(ofp1,"Medium-energy photoproduction data\n");	

	  fprintf(ofp1,
		  "\t(101)\t-->\tsaphir '99\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(102)\t-->\tsaphir '99\t  -- photo --\ttot. cross sec. \t!! Not for fitting !!\n");

	  fprintf(ofp1,
		  "\t(103)\t-->\tsaphir '99\t  -- photo --\trec. asym. \t\t!! Not for fitting !!\n");

	  fprintf(ofp1,
		  "\t(104)\t-->\tsaphir '05\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(105)\t-->\tforw. saphir '05\n");
	  fprintf(ofp1,	"\t\t\t(cos > 0.0)\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(106)\t-->\tsaphir '05\t  -- photo --\trec. asym. \t\t!! Not for fitting !!\n");

	  fprintf(ofp1,
		  "\t(107)\t-->\tforw. saphir '05\n");
	  fprintf(ofp1,	"\t\t\t(cos > 0.2)\t  -- photo --\trec. asym.\t\t!! Not for fitting !!\n");

	  fprintf(ofp1,
		  "\t(108)\t-->\tsaphir '05\t  -- photo --\ttot. cross sec. \t!! Not for fitting !!\n");

	  fprintf(ofp1,
		  "\t(109)\t-->\tCLAS '03\t  -- photo --\tdiff. cross sec.\n");

	  fprintf(ofp1,
		  "\t(110)\t-->\tCB-ELSA/TAPS '07  -- photo --\tdiff. cross sec.\n");
	  
	  fprintf(ofp1,
		  "\t(111)\t-->\tCB-ELSA/TAPS '07  -- photo --\ttot. cross sec.\n");

	  fprintf(ofp1,
		  "\t(112)\t-->\tCB-ELSA/TAPS '07  -- photo --\trec. asym.\n");	  

	  fprintf(ofp1,
		  "\t(113)\t-->\tMAMI-C '12        -- photo --\tdiff. cross sec.\t!! Preliminary data !!\n");	  

	  fprintf(ofp1,
		  "\t(114)\t-->\tMAMI-C '12        -- photo --\trec. asym.\t\t!! Preliminary data !!\n");	  
	}

      else if(iso == 4)
	{
	  /*--------------------------------
	   * Isopin channel 4 (g n -> K0 L0)
	   *--------------------------------*/	
	  fprintf(ofp1,
		  "\t(101)\t-->\tCLAS '09\t  -- photo --\tpho. asym. \t!! For our eyes only !!\n");
	}

      else if(iso == 5)
	{
	  /*--------------------------------
	   * Isopin channel 5 (g n -> K0 S0)
	   *--------------------------------*/	
	  fprintf(ofp1,
		  "\t(101)\t-->\tCLAS '09\t  -- photo --\tpho. asym. \t!! For our eyes only !!\n");
	}


      else if(iso == 6)
	{
	  /*--------------------------------
	   * Isopin channel 6 (g n -> K+ S-)
	   *--------------------------------*/	

	  fprintf(ofp1,
		  "\t(101)\t-->\tleps '06\t  -- photo --\tdiff.cross sec.\n");

	  fprintf(ofp1,
		  "\t(102)\t-->\tleps '06\t  -- photo --\tpho. asym.\n");

	  fprintf(ofp1,
		  "\t(103)\t-->\tCLAS 2010\t  -- photo --\tdiff.cross sec.\n");
	}

      else if(iso == 8)
	{
	  fprintf(ofp1,
		  "\t(801)\t-->\tFpi-1 (JLab) '97\t  -- elec -- \tL, T, TT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(802)\t-->\tFpi-2 (JLab) '03\t  -- elec -- \tvirtual-photoproduction cross section (W = 1911 - 2001 MeV).\n");
	  fprintf(ofp1,
		  "\t(803)\t-->\tFpi-2 (JLab) '03\t  -- elec -- \tvirtual-photoproduction cross section (W = 2127 - 2308 MeV).\n");
	  fprintf(ofp1,
		  "\t(804)\t-->\tFpi-2 (JLab) '03\t  -- elec -- \tL, T, TT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(805)\t-->\tpi-CT (JLab) '07\t  -- elec -- \tunseparated cross section.\n");
	  fprintf(ofp1,
		  "\t(806)\t-->\tDESY-114 (Brauel) '76\t  -- elec -- \tunseparated and TT, TL separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(807)\t-->\tCLAS (JLab) '04 \t  -- elec -- \tbeam single spin asymmetry.\n");
	  fprintf(ofp1,
		  "\t(808)\t-->\tDESY-126 (Ackermann) '78  -- elec -- \tL, T, TT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(809)\t-->\tDESY (Brauel) '79 \t  -- elec -- \tL, T, TT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(810)\t-->\tHERMES (DESY) '08 \t  -- elec -- \tunseparated cross section.\n");
	  fprintf(ofp1,
		  "\t(811)\t-->\tpi-CT (JLab) '07 \t  -- elec -- \tL, T, TT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(812)\t-->\tCEA (Brown) '73 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	  fprintf(ofp1,
		  "\t(813)\t-->\tCEA (Brown) '73 \t  -- elec -- \tLT separated cross section.\n");
	  fprintf(ofp1,
		  "\t(814)\t-->\tCEA (Brown) '73 \t  -- elec -- \tTT, LT separated cross sections.\n");
	  fprintf(ofp1,
		  "\t(815)\t-->\tCornell (Bebek) '74 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	  fprintf(ofp1,
		  "\t(816)\t-->\tCornell (Bebek) '74 \t  -- elec -- \tunseparated and TT, TL separated cross sections.\n");
	  fprintf(ofp1,
	  	  "\t(817)\t-->\tCornell (Bebek) '76 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	  fprintf(ofp1,
	  	  "\t(818)\t-->\tCornell (Bebek) '76 \t  -- elec -- \tL, T separated cross sections.\n");
	  fprintf(ofp1, 
		  "\t(819)\t-->\tCornell (Bebek) '78 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	  fprintf(ofp1,
		  "\t(820)\t-->\tCLAS (JLab) '13 \t  -- elec -- \tunseparated cross section.\n");
	}

      else if(iso == 10)
	{
	  fprintf(ofp1,
		  "\t(1001)\t-->\tDESY-114 (DAS) '76\t  -- elec -- \tunseparated and TT, TL separated cross sections.\n");
	  fprintf(ofp1,
	  	  "\t(1002)\t-->\tCornell (Bebek) '76 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	  fprintf(ofp1,
	  	  "\t(1003)\t-->\tCornell (Bebek) '78 \t  -- elec -- \tvirtual-photoproduction cross section.\n");
	}

      else if(iso == 11)
	{
	  /*---------------------------------
	   * Isopin channel 11 (K- p -> g L0)
	   *---------------------------------*/	

	  fprintf(ofp1,
		  "Medium-energy kaon capture data\n");

	  fprintf(ofp1,
		  "\t(301)\t-->\tbrookhaven '89\t  -- radcap -- \tratio\n");

	  fprintf(ofp1,
		  "\t(302)\t-->\tcrystal ball 2001 -- radcap -- \tdiff. ");
	  fprintf(ofp1,"cross sec.\n");

	  fprintf(ofp1,
		  "\t(303)\t-->\tcrystal ball 2001 -- radcap -- \ttot. ");
	  fprintf(ofp1,"cross sec.\n");

	  fprintf(ofp1,
		  "\t(304)\t-->\tcrystal ball 2009 -- radcap -- \tdiff. ");
	  fprintf(ofp1,"cross sec.\n");

	  fprintf(ofp1,
		  "\t(305)\t-->\tcrystal ball 2009 -- radcap -- \ttot. ");
	  fprintf(ofp1,"cross sec.\n");
	}

      else if(iso == 12)
	{
	  /*---------------------------------
	   * Isopin channel 12 (K- p -> g S0)
	   *---------------------------------*/	

	  fprintf(ofp1,"Medium-energy kaon capture data\n");

	  fprintf(ofp1,
		  "\t(301)\t-->\tbrookhaven '89\t  -- radcap -- \tratio\n");

	  fprintf(ofp1,
		  "\t(304)\t-->\tcrystal ball 2009 -- radcap -- \tdiff. ");
	  fprintf(ofp1,"cross sec.\n");

	  fprintf(ofp1,
		  "\t(305)\t-->\tcrystal ball 2009 -- radcap -- \ttot. ");
	  fprintf(ofp1,"cross sec.\n");

	  fprintf(ofp1,
		  "\t(306)\t-->\tcrystal ball 2009 (Stanislaus et al.) -- radcap -- \ttot. ");
	  fprintf(ofp1,"cross sec.\n");
	}


      else
	error_exit("No data available for this experiment\n");


      fprintf(ofp1,"Give nr. separated with \",\"\t\t\t: ");

      fscanf(ifp,"%s", ans);
      fprintf(ofp2,"%s\n", ans);


      k = 0;
      l = 0;
      while(ans[l] != '\0')
	{

	  j = 0;
	  while(ans[l] != ',' && ans[l] != '\0')
	    exp[j++] = ans[l++];
	  exp[j] = '\0';

	  if(ans[l] == ',')
	    l++;


	  experiment = atoi(exp);


	  /* 
	   * For each dataset in this isospin channel, determine which kind of 
	   * variables should be computed (photo/electro/kaoncap, diff/tot/...)
	   */

	  if(iso == 1)
	    {
	      /*-------------------------------
	       * Isopin channel 1 (g p -> K+ L)
	       *-------------------------------*/

	      switch (experiment){
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.saphir.98");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 102:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.tot.saphir.98");
		observ->fit.photo_totcs[iso] = 1;
		observ->photoprod = 1;
		*tot_angle = 1;
		*e_grid = 1;
		break;
	      case 103:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.saphir.98");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 104:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.saphir.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 105:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.saphir.fwang.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 106:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.saphir.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 107:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.saphir.fwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 108:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas-cmu.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 109:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas-cmu.fwang.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 110:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas-cmu.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 111:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas-cmu.fwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 112:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas-cmu.exfwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 113:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 114:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.fwang.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 115:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.exfwang.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 116:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.cxz_clas");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 117:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.pho.spring8.03");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 118:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.spring8.04");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 119:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.pho.graal.05");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 120:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.pho.fwang.graal.05");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 121:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.graal.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 122:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.fwang.graal.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 123:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.tar.bonn.78");
		observ->fit.photo_tarpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 124:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.leps.07");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 125:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.leps.07_cos");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 126:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.pho.leps.07");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 127:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.tar.graal.08");
		observ->fit.photo_tarpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 128:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.olt.graal.08");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;  
	      case 129:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas-cmu.poscos.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
		/* no E-grid BUT cos-grid! */
	      case 130:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.poscos.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 131:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.cxz_clas.poscos");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 132:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.tar.graal.poscos.08");
		observ->fit.photo_tarpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 133:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.olt.graal.poscos.08");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;  
	      case 134:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.09");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;  
              case 135:
                strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas.09");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;  
	      case 136:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.poscos.09");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;  
              case 137:
                strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas.poscos.09");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;  
	      case 138:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.exfwang.09");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
              case 139:
                strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas.exfwang.09");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;    		
		
	      case 201:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.cea_72");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 202:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.harvard_cornell_74");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 203:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.harvard_cornell_77");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 204:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.jlab_e93_018");
		observ->fit.elec_diffcs[iso] = 1; /* this indicates diff cs 
						   * without polarized baryons! 
						   * fit to polarized electro-data
						   * not yet implemented! */
		observ->electroprod = 1;
		break;
	      case 205:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.lambda-3str-2.567.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 206:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.lambda-3str-2.567-lt.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 207:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.lambda-3str-4.056.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 208:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.lambda-4str-epsphi.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 209:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.clas_E-00-112_iso-1_sigma-ltprime");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 210:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.coman_2009_arxiv0911-3943_sigma_l+t.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;	
	      case 211:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.coman_2009_arxiv0911-3943_sigma_t.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;	
	      case 212:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.coman_2009_arxiv0911-3943_sigma_l.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 213:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.carman_2009_PRC79-065205_transPol_unprimed_ebeam4gev.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 214:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.carman_2009_PRC79-065205_transPol_unprimed_ebeam6gev.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 215:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.carman_2009_PRC79-065205_lt_ratio.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 216:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.cha_2000_phd_exp-E91-016_sigma_l+t.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 217:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_U-LT-TT_E=5.499");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 218:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_LT'_E=5.499");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 219:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_induced_Lambda_pol");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 401:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.slac.t.69");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 402:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.slac.u.69");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 403:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.pho.slac.t.79");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 404:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.desy.t.72");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 405:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.diff.clas.exfwang.HE.09");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 406:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.1.rec.clas.exfwang.HE.09");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 501:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.Cornell_unseparated_1974-kaonHighE");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 502:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.Cornell_unseparated_1977-kaonHighE");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 503:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_U-LT-TT_E=5.499-kaonHighE");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 504:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_LT'_E=5.499-kaonHighE");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 505:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.1.CLAS_induced_Lambda_pol-kaonHighE");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  else if(iso == 2)
	    {
	      /*--------------------------------
	       * Isopin channel 2 (g p -> K+ S0)
	       *--------------------------------*/

	      switch (experiment) {
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.saphir.98");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 102:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.tot.saphir.98");
		observ->fit.photo_totcs[iso] = 1;
		observ->photoprod = 1;
		*tot_angle = 1;
		*e_grid = 1;
		break;
	      case 103:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.saphir.98");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 104:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.saphir.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 105:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.saphir.fwang.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 106:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.saphir.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 107:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.saphir.fwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 108:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas-cmu.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 109:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas-cmu.fwang.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 110:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas-cmu.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 111:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas-cmu.fwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 112:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas-cmu.exfwang.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		/* no E-grid BUT cos-grid! */
		break;
	      case 113:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 114:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.fwang.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 115:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.exfwang.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 116:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.cxz_clas");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 117:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.pho.spring8.03");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 118:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.spring8.04");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 119:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.pho.graal.05");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 120:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.pho.fwang.graal.05");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 121:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.graal.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 122:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.fwang.graal.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 123:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.leps.06");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 124:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.pho.leps.06");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 129:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas-cmu.poscos.03");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
		/* no E-grid BUT cos-grid! */
              case 130:
                strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.poscos.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 131:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.cxz_clas.poscos");
		observ->fit.photo_brpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 132:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.10");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 133:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.poscos.10");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 134:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas.10");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 135:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas.poscos.10");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
		
	      case 201:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.cea_72");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 202:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.harvard_cornell_74");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 203:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.harvard_cornell_77");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 204:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.jlab_e93_018");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 205:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.sigma0-3str-2.567.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 206:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.sigma0-3str-2.567-lt.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 207:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.sigma0-3str-4.056.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 208:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.sigma0-4str-epsphi.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 209:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.carman_2009_PRC79-065205_transPol_unprimed.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 210:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.2.cha_2000_phd_exp-E91-016_sigma_l+t.dat");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 401:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.slac.t.69");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 402:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.slac.u.69");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 403:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.pho.slac.t.79");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 404:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.diff.clas.exfwang.HE.10");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 405:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.2.rec.clas.exfwang.HE.10");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  else if(iso == 3)
	    {
	      /*--------------------------------
	       * Isopin channel 3 (g p -> K0 S+)
	       *--------------------------------*/  

	      switch (experiment) {
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.saphir.99");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 102:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.tot.saphir.99");
		observ->fit.photo_totcs[iso] = 1;
		observ->photoprod = 1;
		*tot_angle = 1;
		*e_grid = 1;
		break;
	      case 103:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.rec.saphir.99");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 104:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.saphir.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 105:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.fwangles.saphir.05");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		*e_grid = 1;
		break;
	      case 106:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.rec.saphir.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 107:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.rec.fwangles.saphir.05");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		*e_grid = 1;
		break;
	      case 108:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.tot.saphir.05");
		observ->fit.photo_totcs[iso] = 1;
		observ->photoprod = 1;
		*tot_angle = 1;
		*e_grid = 1;
		break;
	      case 109:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.clas.03");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 110:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.cbelsa-taps.07");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 111:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.tot.cbelsa-taps.07");
		observ->fit.photo_totcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 112:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.rec.cbelsa-taps.07");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 113:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.diff.mami-c.2012");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 114:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.3.rec.mami-c.2012");
		observ->fit.photo_recpol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }

	  else if(iso == 4)
	    {
	      /*--------------------------------
	       * Isopin channel 4 (g n -> K0 L0)
	       *--------------------------------*/	          					     
	      switch (experiment) {
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.4.pho.clas.09");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }

	  else if(iso == 5)
	    {
	      /*--------------------------------
	       * Isopin channel 5 (g n -> K0 S0)
	       *--------------------------------*/	          					     
	      switch (experiment) {
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.5.pho.clas.09");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  else if(iso == 6)
	    {
	      /*--------------------------------
	       * Isopin channel 6 (g n -> K+ S-)
	       *--------------------------------*/

	      switch (experiment) {
	      case 101:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.6.diff.leps.06");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      case 102:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.6.pho.leps.06");
		observ->fit.photo_phopol[iso] = 1;
		observ->photoprod = 1;
		*pol_data = 1;
		break;
	      case 103:
		strcpy(observ->fit.database_info[iso][k++], "photo.iso.6.diff.clas.2010");
		observ->fit.photo_diffcs[iso] = 1;
		observ->photoprod = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }
	  

	  else if(iso == 8)
	    {
	      /*--------------------------------
	       * Isopin channel 8 (g p -> pi+ n)
	       *--------------------------------*/

	      switch (experiment) {
	      case 801:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.JLab_Fpi-1");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 802:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Fpi-2(Tadevosyan)-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 803:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Fpi-2(Horn)-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 804:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.JLab_Fpi-2");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 805:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.JLab_pi-CT_unseparated");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 806:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.DESY_proton-U_TT_LT");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 807:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.JLab-a_lu");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 808:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.DESY_separated_Q2=0.35_W=2.10");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 809:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.DESY_separated_Q2=0.70_W=2.19");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 810:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.HERMES-diff_U");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 811:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.JLab_pi-CT_separated");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 812:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-CEA_1973-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 813:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-CEA_1973-diff_LT");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 814:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-CEA_1973-diff_TT-LT");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 815:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Cornell_1974-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 816:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Cornell_1974-diff_U-TT-LT");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 817:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Cornell_1976-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 818:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Cornell_1976-diff_L-T");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 819:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.SAID_electro-Cornell_1978-diff_phi");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 820:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.8.CLAS2013_unseparated");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  else if(iso == 10)
	    {
	      /*---------------------------------
	       * Isopin channel 10 (g n -> pi- p)
	       *---------------------------------*/

	      switch (experiment) {
	      case 1001:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.10.DESY_neutron-U_TT_LT");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 1002:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.10.SAID_electro-Cornell_1976-diff_phi-neutron");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      case 1003:
		strcpy(observ->fit.database_info[iso][k++], "electro.iso.10.SAID_electro-Cornell_1978-diff_phi-neutron");
		observ->fit.elec_diffcs[iso] = 1;
		observ->electroprod = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  else if(iso == 11)
	    {
	      /*---------------------------------
	       * Isopin channel 11 (K- p -> g L0)
	       *---------------------------------*/
	      
	      switch (experiment) {
	      case 301:
		strcpy(observ->fit.database_info[iso][k++], "brookhaven_89");
		observ->fit.kaoncap_branrat[iso] = 1;
		observ->kaoncapture = 1;
		break;
	      case 302:
		strcpy(observ->fit.database_info[iso][k++], "diff.crystalball_2001");
	        observ->fit.kaoncap_diffcs[iso] = 1;
	        observ->kaoncapture = 1;
		break;
	      case 303: 
		strcpy(observ->fit.database_info[iso][k++], "tot.crystalball_2001");
	        observ->fit.kaoncap_totcs[iso] = 1;
	        observ->kaoncapture = 1;
		*tot_angle = 1;
		break;
	      case 304:
		strcpy(observ->fit.database_info[iso][k++], "diff.cbc_ags_2009");
	        observ->fit.kaoncap_diffcs[iso] = 1;
	        observ->kaoncapture = 1;
		break;
	      case 305: 
		strcpy(observ->fit.database_info[iso][k++], "tot.cbc_ags_2009");
	        observ->fit.kaoncap_totcs[iso] = 1;
	        observ->kaoncapture = 1;
		*tot_angle = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }

	  else if(iso == 12)
	    {
	      /*---------------------------------
	       * Isopin channel 12 (K- p -> g S0)
	       *---------------------------------*/
	      
	      switch (experiment) {
	      case 301:
		strcpy(observ->fit.database_info[iso][k++], "brookhaven_89");
		observ->fit.kaoncap_branrat[iso] = 1;
		observ->kaoncapture = 1;
		break;
	      case 304:
		strcpy(observ->fit.database_info[iso][k++], "diff.cbc_ags_2009");
	        observ->fit.kaoncap_diffcs[iso] = 1;
	        observ->kaoncapture = 1;
		break;
	      case 305: 
		strcpy(observ->fit.database_info[iso][k++], "tot.cbc_ags_2009");
	        observ->fit.kaoncap_totcs[iso] = 1;
	        observ->kaoncapture = 1;
		*tot_angle = 1;
		break;
	      case 306: 
		strcpy(observ->fit.database_info[iso][k++], "tot.cbc_ags_stanislaus_2009");
	        observ->fit.kaoncap_totcs[iso] = 1;
	        observ->kaoncapture = 1;
		*tot_angle = 1;
		break;
	      default:
		error_exit("Not a valid experiment!\n");
	      }
	    }


	  if(k > MAXEXP)
	    error_exit("More experiments than MAXEXP\n");
	  
	} 

      /* 
       * Determine how many experiments there are in each isospin channel
       */
      
      observ->fit.nr_experiments[iso] = k;
      
    }


  return 0;
}

/*!
 * Reads in the hadronic formfactor specifications.
 * Remark that the cutoff values are always read in from the coupl.iso.# file.
 */
int hadronformfacspecification(Observable* observ, FILE* ifp, FILE* ofp1, 
			       FILE* ofp2)
{
  char ans[3];

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Hadronic form factors (y/n)\t\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);
  
  if(!strcmp(ans, "y"))
    {
      
      observ->hadronformfac = 1;
      
      fprintf(ofp1,"---\n");
      fprintf(ofp1,"Dipole, Gaussian or Lorentz hadronic form factors? (d/g/l)\t\t: ");
      fscanf(ifp,"%s", ans);
      fprintf(ofp2,"%s\n", ans);
      if (!strcmp(ans, "d"))
	error_exit("Dipole is no longer supported in version 10.0\n");
      else if(strcmp(ans, "l") && strcmp(ans, "g"))
  	error_exit("Answer Dipole, Gaussian or Lorentz!\n");
      observ->ffac.gauss_strongff = 1; // by default, set to Gaussian

      /* Isobar model: only here, a gauge-restoration procedure is
       * necessary; in the (mixed) Regge model, only the resonances 
       * can have hadronic ff's */	  
      if(!observ->regge)
	{	      
	  /* select gauge-restoration Fhat (here DW) */
	  observ->ffac.davidson_fhat = 1;
	  observ->ffac.haberzettl_fhat = 0;
	  observ->ffac.ohta_fhat = 0;
	  
	  if(observ->ffac.haberzettl_fhat)
	    {
	      fprintf(ofp1,"---\n");
	      fprintf(ofp1,"Fraction\t\t\t\t\t: ");
	      fscanf(ifp,"%lf", &observ->ffac.fraction);
	      fprintf(ofp2,"%lf\n", observ->ffac.fraction);
	    }    
	  
	  if(observ->ffac.haberzettl_fhat || observ->ffac.ohta_fhat)
	    {
	      printf("\n\n!!!! ATTENTION: Davdison-Workman \"F_hat\" ");
	      printf("is NOT used !!!!\n\n");
	    }
	}
    }
  
  else if(!strcmp(ans, "n"))
    observ->hadronformfac = 0;
  
  else 
    error_exit("Answer y(es) or n(o)!\n");
  
  return 0;
}

/*!
 * Reads in the specifications for the Regge propagator.
 */
int reggespecification(Observable* observ, FILE* ifp, FILE* ofp1, FILE* ofp2)
{
  char ans[60], val[20];
  int i, j, l, value;

  

  /* T- or U-channel reggeization, or both? */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"T or U channel reggeization, or both? (t/u/b)\t\t: ");
  //fscanf(ifp,"%s", ans);
  strcpy(ans,"t");
  fprintf(ofp2,"%s\n", ans);
  fprintf(ofp1,"Only T channel reggeization is implemented\n");
  
  if(!strcmp(ans, "t"))
    {
      /* 
       * T-channel Reggeization 
       */

      observ->reg.t_channel = 1;
    }
  else if(!strcmp(ans, "u"))
    {
      /* 
       * U-channel Reggeization 
       */
      
      observ->reg.t_channel = 0;

    }

  else if(!strcmp(ans, "b"))
    {
      /* 
       * T- and U-channel Reggeization, depending on cos_theta
       */

      observ->reg.t_and_u_channel = 1;
            
    }
  else
    error_exit("Answer T or U channel or both!\n");


  /* Exact legendre form, or high-energy, low-t/u limit? */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Exact or asymptotic behavior? (x/a)\t\t: ");
  //fscanf(ifp,"%s", ans);
  strcpy(ans,"a");
  fprintf(ofp2,"%s\n", ans);
  fprintf(ofp1,"Only asymptotic behaviour is implemented\n");

  if(!strcmp(ans, "x"))
    observ->reg.asymp = 0;
  else if(!strcmp(ans, "a"))
    observ->reg.asymp = 1;
  else
    error_exit("Answer eXact or Asymptotic\n");
  

  /* Replacing s by (s-u)/2 */

  /* observ->reg.s_modif = 0; */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"'s' or '(s-u)/2' dependence? (s/u)\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);
  
  if(!strcmp(ans, "s"))
    observ->reg.s_modif = 0;
  else if(!strcmp(ans, "u"))
    observ->reg.s_modif = 1;
  else
    error_exit("Answer 's' or 'u' -> (s-u)/2 dependence\n");



  /* Linear trajectories or saturating ones? */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Linear or Saturating trajectories? (l/s)\t: ");
  //fscanf(ifp,"%s", ans);
  strcpy(ans,"l");
  fprintf(ofp2,"%s\n", ans);
  fprintf(ofp1,"Only linear trajectories are implemented\n");

  if(!strcmp(ans, "l"))
    observ->reg.sat_traj = 0;
  else if(!strcmp(ans, "s"))
    observ->reg.sat_traj = 1;
  else
    error_exit("Answer l(inear) or s(aturating) trajectories!!\n");


  /* S. Janssen's or M. Guidal's recipe for Regge propagator of 
   * spin>0 particles? */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Spin-shift recipe (alpha->alpha-n)\n");
  fprintf(ofp1,"of Guidal or Janssen? (j/g)\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);

  if(!strcmp(ans, "g"))
    observ->reg.spinshift_guidal = 1;
  else if(!strcmp(ans, "j"))
    observ->reg.spinshift_guidal = 0;
  else
    error_exit("Answer g(uidal) or j(anssen) recipe!!\n");  


  /* Pure Regge (eventually with duality corrections), 
   * or reggeized background plus isobar-model resonances? */

  fprintf(ofp1,"---\n");
  fprintf(ofp1,"Pure regge theory, or resonances Superimposed\n");
  fprintf(ofp1,"on reggeized background? (p/s)\t\t: ");
  fscanf(ifp,"%s", ans);
  fprintf(ofp2,"%s\n", ans);
  
  
  if(!strcmp(ans, "p"))
	{
	  /* Pure Regge: reggeized background only */

	  observ->reg.res_on_bg = 0;


	  /* Duality corrections (Veneziano Model)? (for now, I'm only sure this works
	   * with Stijn's spin-shift recipe, not with Guidal's.*/

	  if(observ->reg.rot)
	    {
	      fprintf(ofp1,"---\n");
	      fprintf(ofp1,"Duality corrections (caution: for now, only correct\n");
	      fprintf(ofp1,"with Janssen spin-shift recipe)? (y/n)\t\t\t: ");
	      fscanf(ifp,"%s", ans);
	      fprintf(ofp2,"%s\n", ans);
	      
	      if(!strcmp(ans, "y"))
		observ->reg.dual_corr = 1;
	      else if(!strcmp(ans, "n"))
		observ->reg.dual_corr = 0;
	      else
		error_exit("Answer Yes or No duality corrections!\n");
	    }
	  
	  if(observ->reg.dual_corr)
	    {
	      fprintf(ofp1,"---\n");
	      fprintf(ofp1,"Indicate the s-channel trajectories:\n ");
	      fprintf(ofp1,"\t(1)  --> N^* (1/2)^+ (a)\n");
	      fprintf(ofp1,"\t(2)  --> N^* (3/2)^- (a)\n");
	      fprintf(ofp1,"\t(3)  --> N^* (1/2)^+ (b)\n");
	      fprintf(ofp1,"\t(4)  --> N^* (3/2)^- (b)\n");
	      fprintf(ofp1,"\t(5)  --> N^* (1/2)^- (a)\n");
	      fprintf(ofp1,"\t(6)  --> N^* (3/2)^+ (a)\n");
	      fprintf(ofp1,"\t(7)  --> N^* (1/2)^- (b)\n");
	      fprintf(ofp1,"\t(8)  --> N^* (3/2)^+ (b)\n");
	      fprintf(ofp1,"\t(9)  --> N^* (1/2)^- (c)\n");
	      fprintf(ofp1,"\t(10) --> N^* (3/2)^+ (c)\n");
	      
	      if(observ->iso.iso_base == 2)
		{
		  fprintf(ofp1,"\t(11) --> D^* (1/2)^+\n");
		  fprintf(ofp1,"\t(12) --> D^* (3/2)^-\n");
		  fprintf(ofp1,"\t(13) --> D^* (1/2)^-\n");
		  fprintf(ofp1,"\t(14) --> D^* (3/2)^+\n");
		}
	      
	      fprintf(ofp1,"(give numbers separated with \",\" or type \"all\")\t: ");
	      
	      fscanf(ifp,"%s", ans);
	      fprintf(ofp2,"%s\n", ans);
	      
	      if(!strcmp(ans, "all"))
		{
		  if(observ->iso.iso_base == 1)
		    strcpy(ans, "1,2,3,4,5,6,7,8,9,10");
		  else if(observ->iso.iso_base == 2)
		    strcpy(ans, "1,2,3,4,5,6,7,8,9,10,11,12,13,14");
		}
	      
	      
	      i = 0;
	      l = 0;
	      while(ans[l] != '\0')
		{
		  j = 0;
		  while(ans[l] != ',' && ans[l] != '\0')
		    val[j++] = ans[l++];
		  
		  if(ans[l] == ',')
		    l++;
		  
		  val[j] = '\0';
		  value = atoi(val);
		  
		  observ->reg.s_traj_type[i++] = value;
		  
		}
	    }
	}
  
  else if(!strcmp(ans, "s"))
    {
      /* Resonances superimposed on Reggeized background */
      observ->reg.res_on_bg = 1;
      
      /* No duality corrections */
      observ->reg.dual_corr = 0;
    }
  
  else
    error_exit("Answer p(ure regge) or s(uperimposed resonances)\n");



  return 0;
}




