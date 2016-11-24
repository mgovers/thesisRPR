/*!
 * \file strange_func.cpp
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
#include <stdexcept>
#include <Structures.h>
#include "strange_func.h"
#include "fitting.h"
#include "calcMatrixElement.h"
#include "FormFactorSpecification.h"

/*!
 * Returns the isospin base channel that goes with isospin channel 'iso'
 */
int getisospinbasechannel(int iso)
{
  switch(iso) {
  case 1: return 1;
    break;
  case 2: return 2;
    break;
  case 3: return 2;
    break;
  case 4: return 1;
    break;
  case 5: return 2;
    break;
  case 6: return 2;
    break;
  case 7: return 7;
    break;
  case 8: return 8;
    break;
  case 9: return 7;
    break;
  case 10: return 8;
    break;
  case 11: return 1;
    break;
  case 12: return 2;
    break;
  default:
    fprintf(stderr,"ERROR in getisospinbasechannel(int): isospin channel %d unkown.\n",
	    iso);
    exit(1);    
  }

  return 0;
}

/*!
 * In this function, the particle specifications (mass, width, 
 * coupling constants, ...) are read in from the files 
 * "./input/numinput/mass.in" and "./input/numinput/coupl.iso.*"
 * (where * is the isospin channel), and stored in an easily accessible way  
 * in the array of Class, particles[]. Every element of the particles[] 
 * array contains a class of particles. For each class, the different 
 * particles (with their properties) are stored in the "partic" array 
 * of Properties, and "particount" gives the number of particles of one 
 * class. Remark that in the file "mass.in" all the particles are stored,
 * the particle selection for a specific situation is made in the
 * "coupl.iso.*" file. At the end off the coupl.iso.* file the hadronic 
 * form factor cutoff values are read in and stored in the "observ" 
 * structure.
 */
int particlespecify(Class particles[], Observable* observ, FILE* ifp2, FILE* ofp)
{
  int i;
  
  
  /*-----------------
   * Read coupl.iso.*
   *-----------------*/

  /* 
   * Read the cc's and cutoffs from the "coupl.iso.*" file, and store them 
   * in the particles[] array and the observ structure (which has already 
   * been partly filled by observable-specification()). 
   * Input is taken from the base isospin channel 1 or 2.
   */

  /* Initialize no of particles in each class to 0 */
  for(i=0; i<CLASSMAX; i++)
    particles[i].particount = 0;


  /* Store the coupl.iso.* info in particles[] and observ. */
  import_baseclass_particles(particles, observ, ifp2, ofp);
	
  //fclose(ifp);
 
  return 0;
}

/*!
 * - Read the NICKNAME and the G,H,X,Y,Z values for each particle j 
 *   in each class i from the coupl.iso.* file,
 *   and store them in the PARTICLES[] array 
 *   (particles[i].partic[j].nickname, particles[i].partic[j].G/H/X/Y/Z).
 * - Read the Born and res. CUTOFFS, and store them in the OBSERV structure 
 *   (observ->ffac.born/res_cutoff).
 * - Determine the NUMBER of particles in each CLASS and store 
 *   in the PARTICLES[] array (particles[i].particount) 
 * - For all Regge trajectories, specify whether the phase is constant or rotating
 *   by prompting the user. (This could be specified in an input file?)
 * - If particount for one of the Born terms is zero, the nickname is added.
 *   Without changing particount!
 * 
 * This function has as a parameter the base isospin channel (K+L or K+S0): 
 * ifp = input/numinput/coupl.iso.1 or 2; isospin = 1 or 2
 */
int import_baseclass_particles(Class particles[], Observable* observ,
			       FILE* ifp2, FILE* ofp) 
{
  int i,j;
  char dump[DUMP_SIZE], classlabel[7], filename[MAX_LOCATION_STRING];



  /* ----------------------------
   * Determine and open inputfile
   * ----------------------------*/
  int isospin = observ->iso.iso_base;

  sprintf(filename,"%s%s%d",
	  observ->inFolder,
	  "numinput/coupl.iso.",
	  isospin);

  FILE* ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("error in opening \"coupl.iso-file\"\n");


  /* Dump the header of the coupl.iso.*-file (7 column titles and 7 lines) */
  
  for(i=0; i<16; i++)
    fscanf(ifp, "%s", dump);

  

  /*----------------------------- 
   * Store the coupling constants 
   *-----------------------------*/
  
  fscanf(ifp, "%s", classlabel);

  while(strcmp(classlabel, "cutoff"))
    {

      /* Find the index i in the particles[] array of Classes 
	 corresponding to the classlabel. */

      i = findclassindex(classlabel[0]);


      /* Scan the nickname and G,H,I;X,Y,Z,X2,Y2,Z2 values, and store them in the 
       * particles[i].partic[j] structure corresp. to partic. j in class i */

      fscanf(ifp, "%s", particles[i].partic[particles[i].particount].nickname);

      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].G);
      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].H);
      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].I);
      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].X);
      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].Y);
      fscanf(ifp, "%le", &particles[i].partic[particles[i].particount].Z);


      /* Specify wether gauge-invariant couplings are used or not*/
      if( (!strcmp(observ->modelType,"offshell")) && i < 14) // Always use gauge invariant couplings for spin 5/2 particles (i =,14,15,16,17)
	  particles[i].partic[particles[i].particount].gic = 0;
      else 
	particles[i].partic[particles[i].particount].gic = 1;


      /* Count the particles in the class i while storing their properties */  

      particles[i].particount++;

      
      if(particles[i].particount > PARTICLEMAX)
	error_exit("More particles in a class than allowed by PARTICLEMAX\n");
      
      
      fscanf(ifp, "%s", classlabel);
      
    }



  /*--------------
   * Store cutoffs
   *--------------*/

  /* 
   * Store the form factor values, specified at the end of the file.
   */
 
  if(!strcmp(classlabel, "cutoff"))
    {
      fscanf(ifp, "%s%lf", dump, &observ->ffac.born_cutoff);
      fscanf(ifp, "%s%s%lf", dump, dump, &observ->ffac.res_cutoff);
    }
  else
    {
      fprintf(stderr, 
	      "Error: coupl.iso.%d file is not cleanly terminated!\n",
	      isospin);
      exit(1);
    }

  fclose(ifp);

  

  /*-----------------------------
   * Phase for Regge trajectories
   *-----------------------------*/
  if(observ->regge)
    {
      printf("\n---");
      /* T-channel reggeization */
      if(observ->reg.t_channel || observ->reg.t_and_u_channel)
	{
	  int Nclasses = 6;
	  short phase1, phase2, Nphases;
	  char phases[3];
	  
	  int classindex[] = 
	    {1,   // t-channel Born exchange
	     4,   // t-channel vector meson exhange
	     5,   // t-channel axial vector meson exchange
	     18,  // t-channel Born exchange (pionHighE)
	     19,  // t-channel vector meson exchange (pionHighE)
	     20}; // t-channel axial vector meson exchange (pionHighE)

	  for(i = 0; i < Nclasses; i++)
	    {
	      for(j = 0; j < particles[classindex[i]].particount; j++)
		{
		  Nphases = 0;
		  printf("\n%s exchange in t-channel (diagram %c): %s, %s or %s phase? ",
			 particles[classindex[i]].partic[j].nickname,
			 findclasslabel(classindex[i]),
			 "constant (0)",
			 "rotating (1)",
			 "non-degenerate (2)");
		  fscanf(ifp2,"%s",phases);
		  fprintf(ofp,"%s\n",phases);
		  
		  Nphases = (phases[1] == ',' ? 2 : 1);

		  phase1 = phases[0] - 48;
		  phase2 = (Nphases == 2 ? phases[2] - 48 : phase1);
			    
		  if(Nphases == 1)
		    {
		      particles[classindex[i]].partic[j].regge_phase = phase1;
		      particles[classindex[i]].partic[j].regge_phase_nonbase = phase1;
		    }
		  else if(Nphases == observ->iso.nr_iso_channels)
		    {
		      short firstIsBase = (observ->iso.iso_channel[0] == observ->iso.iso_base);
		      particles[classindex[i]].partic[j].regge_phase = 
			(firstIsBase ? phase1 : phase2);
		      particles[classindex[i]].partic[j].regge_phase_nonbase =
			(firstIsBase ? phase2 : phase1);
		    }
		  else
		    error_exit("Number of Regge phases exceeds number of isospin channels.\n");
		}
	    }
	}

      /* U-channel reggeization */
      if(!observ->reg.t_channel || observ->reg.t_and_u_channel)
	{
	  // Loop over all U-channel born terms
	  for(i=0; i<particles[2].particount; i++)
	    {
	      printf("\n%s exchange in u-channel born term: constant (0) or rotating (1) phase? ",
		     particles[2].partic[i].nickname);
	      fscanf(ifp2,"%hd",&particles[2].partic[i].regge_phase);
	      fprintf(ofp,"%d\n",particles[2].partic[i].regge_phase);
	    }
	  // Loop over all U-channel extended born terms
	  for(i=0; i<particles[3].particount; i++)
	    {
	      printf("\n%s exchange in u-channel extended born term: constant (0) or rotating (1) phase? ",
		     particles[3].partic[i].nickname);
	      fscanf(ifp2,"%hd",&particles[3].partic[i].regge_phase);
	      fprintf(ofp,"%d\n",particles[3].partic[i].regge_phase);
	    }
	}
    }


  /*------------------------------------------------
   * Longitudinal couplings for spin-1/2 resonances?
   *------------------------------------------------*/
  
  /* In case of electroproduction, longitudinal couplings
   * can be included for spin-1/2 resonances.
   * The user has specified this in the input and this was
   * stored in observ.
   * Now this info is stored in the relevant places of the
   * particles[] array.
   */

  if(observ->electroprod && observ->elec.long_coupl)
    {
      // Loop over all resonant diagrams and all possible resonances
      for(i=6; i<CLASSMAX; i++)
	{
	  for(j=0; j<particles[i].particount; j++)
	    {
	      // and store that longitudinal couplings are needed
	      particles[i].partic[j].long_coupling = 1;
	    }
	}
    }


  /*----------------------
   * Store Born properties
   *----------------------*/

  /* 
   * If one of the Born terms is not included in coupl.iso, the nicknames 
   * are added (though particount remains 0) so the Born-term masses, 
   * needed later, can be read in. 
   */
  

  if(particles[0].particount == 0)
    {
      /* Which particle is in the s-channel... */
      
      if(observ->iso.isospin <= 3 || observ->iso.isospin == 7 || 
	 observ->iso.isospin == 8 || observ->iso.isospin >= 11 ) 
	strcpy(particles[0].partic[0].nickname, "P");
      else 
	strcpy(particles[0].partic[0].nickname, "N");
    }
  
  if(particles[1].particount == 0)
    {
      /* Which particle is in the t-channel... */
      
      if(observ->iso.isospin <= 2 || observ->iso.isospin == 6 )
	strcpy(particles[1].partic[0].nickname, "K+");
      else if(observ->iso.isospin <= 5 )
	strcpy(particles[1].partic[0].nickname, "K0");
      else if(observ->iso.isospin == 7 ||
	      observ->iso.isospin == 9 )
	strcpy(particles[1].partic[0].nickname, "P0");
      else if(observ->iso.isospin == 8 )
	strcpy(particles[1].partic[0].nickname, "P+");
      else if(observ->iso.isospin == 10 )
	strcpy(particles[1].partic[0].nickname, "P-");
      else
	strcpy(particles[1].partic[0].nickname, "K+");
    }
  
  if(particles[2].particount == 0)
    {
      /* Which particle is in the u-channel... */
      
      if(observ->iso.isospin == 1 || 
	 observ->iso.isospin == 4 ||
	 observ->iso.isospin == 11 )
	strcpy(particles[2].partic[0].nickname, "L");
      else if(observ->iso.isospin == 2 || 
	      observ->iso.isospin == 5 ||
	      observ->iso.isospin == 12 )
	strcpy(particles[2].partic[0].nickname, "S0");
      else if(observ->iso.isospin == 3)
	strcpy(particles[2].partic[0].nickname, "S+");
      else if(observ->iso.isospin == 6)
	strcpy(particles[2].partic[0].nickname, "S-");
      else if(observ->iso.isospin == 7 ||
	      observ->iso.isospin == 10 )
	strcpy(particles[2].partic[0].nickname, "P");
      else if(observ->iso.isospin == 8 ||
	      observ->iso.isospin == 9 )
	strcpy(particles[2].partic[0].nickname, "N");
    }
  
  return 0;

}	

/*!
 * If other isospin channels are calculated than the "main" ones 
 * (1 or 2), the base set of coupling constant can be used as well.
 * A number of changed need to be applied however:
 * - The nicknames of certain exchanged particles in the 
 *   (extended) Born terms have to change.
 * - Next the properties (mass,width,charge) of the particles
 *   are read in from the "mass.in" file
 * - The coupling constants change using isospin symmetry and
 *   helicity amplitudes.
 */
int change_isospin_channel(Class particles[], Observable* observ)
{
  /* ----------------
   * modify nicknames 
   * ---------------- */
  
  modify_nicknames(particles,observ);


  /* --------------------------------------------- *
   * Store the particle masses, widths and charge, *
   * read from the mass.in file.                   *
   * --------------------------------------------- */

  import_particle_properties(particles,observ);


  /*-----------------------------------------
   * Magnetic transition moment ratios p/n
   *-----------------------------------------*/

  /*
   * In the case that production off the neutron needs to be calculated 
   * starting from one of the two main sets of coupling constants, 
   * determine the magnetic transition 
   * moment ratios r_kappa_n_p, r_kappa_1_n_p, r_kappa_2_n_p (in Properties 
   * structure) for each particle, to convert the EM vertices: 
   * kappa_{n,N*}/kappa_{p,N*} = photohel amp. ratio.
   */

  if( (observ->iso.isospin >= 4 && observ->iso.isospin <= 6) ||
      (observ->iso.isospin >= 9 && observ->iso.isospin <= 10) ) 
    magnetic_transition_ratio(particles,observ);


  /* ----------------------------- *
   * Change coupling constants for *
   * different isospin channels    *
   * ----------------------------- */

  change_coupling_constants(particles,observ,observ->iso.isospin);


  return 0;
}

/*!
 * Receives a list of particles (particleclass).
 * We perform a loop over all particles and when a particle has
 * an old nickname (oldNick) it is replaced by the new one (newNick).
 */
int change_nickname(Class* particleclass, const char* oldNick, const char* newNick)
{
  int particle;

  for(particle=0; particle<particleclass->particount; particle++)
    {
      if( !strcmp(particleclass->partic[particle].nickname, oldNick) )
	strcpy(particleclass->partic[particle].nickname, newNick);
    }


  return 0;
}

/*!
 * Modifies the nicknames of the exchanged particles in the Born terms
 * (S,T,U) and the extended Born terms (A,B). For these 5 classes, the coupl.iso.*
 * files contain the nicknames of the particles needed in the base isospin
 * channel. In the derived channels, however, other particles are needed 
 * (eg n instead of p). This makes a difference, because in the "mass.in" 
 * file the masses and widths of the different isospin partners are NOT 
 * the same (eg p vs. n); they ARE the same for the other resonances.
 * So, we must make sure these other masses and widths are effectively 
 * taken, when later on the mass.in info is read in! 
 */
int modify_nicknames(Class particles[], Observable* observ)
{
  // All none base classes
  // g + p -> K^0 + S^+  <=  g + p -> K+ S0
  if(observ->iso.isospin == 3) 
    {
      change_nickname(&particles[1],"K+","K0");
      change_nickname(&particles[2],"S0","S+");
      change_nickname(&particles[4],"B1","B2");
      change_nickname(&particles[4],"B4","B3");
      change_nickname(&particles[5],"K1","K10");
      change_nickname(&particles[5],"C1","C2");
    }
  // g + n -> K0 + L  <=  g + p -> K+ L0
  else if(observ->iso.isospin == 4) 
    {
      change_nickname(&particles[0],"P","N");
      change_nickname(&particles[1],"K+","K0");
      change_nickname(&particles[4],"B1","B2");
      change_nickname(&particles[4],"B4","B3");
      change_nickname(&particles[5],"K1","K10");
      change_nickname(&particles[5],"C1","C2");
    }
  // g + n -> K0 + S0  <=  g + p -> K+ S0
  else if(observ->iso.isospin == 5) 
    {
      change_nickname(&particles[0],"P","N");
      change_nickname(&particles[1],"K+","K0");
      change_nickname(&particles[4],"B1","B2");
      change_nickname(&particles[4],"B4","B3");
      change_nickname(&particles[5],"K1","K10");
      change_nickname(&particles[5],"C1","C2");
    }
  // g + n -> K+ + S-  <=  g + p -> K+ S0
  else if(observ->iso.isospin == 6) 
    {
      change_nickname(&particles[0],"P","N");
      change_nickname(&particles[2],"S0","S-");
    }
  // g + n -> P0 + n  <=  g + p -> P0 + p
  else if(observ->iso.isospin == 9) 
    {
      change_nickname(&particles[0],"P","N");
      change_nickname(&particles[2],"P","N");
    }
  // g + n -> P- + p  <=  g + p -> P+ + n
  else if(observ->iso.isospin == 10)
    {
      change_nickname(&particles[0],"P","N");
      change_nickname(&particles[1],"P+","P-");
      change_nickname(&particles[2],"N","P");
    }

  /*
   * We need to make sure the names of the s-,t- and u-channel diagrams are
   * correct for the current isospin channel, even when particount==0. 
   */

  // Double-check the Born channels
  if(particles[0].particount == 0)
    {
      // S-channel      
      if(observ->iso.isospin <= 3 || observ->iso.isospin == 7 || 
	 observ->iso.isospin == 8 || observ->iso.isospin >= 11 ) 
	strcpy(particles[0].partic[0].nickname, "P");
      else 
	strcpy(particles[0].partic[0].nickname, "N");
    }
  if(particles[1].particount == 0)
    {
      // T-channel      
      if(observ->iso.isospin <= 2 || observ->iso.isospin == 6 )
	strcpy(particles[1].partic[0].nickname, "K+");
      else if(observ->iso.isospin <= 5 )
	strcpy(particles[1].partic[0].nickname, "K0");
      else if(observ->iso.isospin == 7 ||
	      observ->iso.isospin == 9 )
	strcpy(particles[1].partic[0].nickname, "P0");
      else if(observ->iso.isospin == 8 )
	strcpy(particles[1].partic[0].nickname, "P+");
      else if(observ->iso.isospin == 10 )
	strcpy(particles[1].partic[0].nickname, "P-");
      else
	strcpy(particles[1].partic[0].nickname, "K+");
    }
  if(particles[2].particount == 0)
    {
      // U-channel      
      if(observ->iso.isospin == 1 || 
	 observ->iso.isospin == 4 ||
	 observ->iso.isospin == 11 )
	strcpy(particles[2].partic[0].nickname, "L");
      else if(observ->iso.isospin == 2 || 
	      observ->iso.isospin == 5 ||
	      observ->iso.isospin == 12 )
	strcpy(particles[2].partic[0].nickname, "S0");
      else if(observ->iso.isospin == 3)
	strcpy(particles[2].partic[0].nickname, "S+");
      else if(observ->iso.isospin == 6)
	strcpy(particles[2].partic[0].nickname, "S-");
      else if(observ->iso.isospin == 7 ||
	      observ->iso.isospin == 10 )
	strcpy(particles[2].partic[0].nickname, "P");
      else if(observ->iso.isospin == 8 ||
	      observ->iso.isospin == 9 )
	strcpy(particles[2].partic[0].nickname, "N");
    }

  return 0;
}

/*!
 * Stores the particle masses, widths and charge read from the mass.in file
 * and sets the particles' spin.
 */
int import_particle_properties(Class particles[], Observable* observ)
{
  int i, k, rol, match;
  char classlabel[2], filename[MAX_LOCATION_STRING];
  char tmpnickname[4];
  char dump[DUMP_SIZE];
  
  FILE* ifp;

  strcpy(filename, observ->inFolder);
  strcat(filename, "numinput/mass.in");
  ifp = fopen(filename, "r");
  if (ifp==NULL)
    error_exit("Cannot open numinput/mass.in");
  

  /* Dump the file header (6 column titles, 6 lines) */
  
  for(i=0; i<14; i++)
    fscanf(ifp, "%s", dump);
  


  fscanf(ifp, "%s", classlabel);

  while(strcmp(classlabel, "X")) /* EOF character */
    {

      /* Find index in particles[] array of Class */
      i = findclassindex(classlabel[0]); 
     
      match = 0;
	
      fscanf(ifp, "%s", tmpnickname);

      
      if(particles[i].particount != 0) 
	{
	  /* Loop over part's in class partic[i], and compare their nicknames 
	   * with tmpnickname: check if the mass.in particle is also 
	   * in coupl.iso */

	  for(rol = 0; rol < particles[i].particount && !match; rol++)


	    /* if match is found, store mass and width */

	    if(!strcmp(particles[i].partic[rol].nickname,tmpnickname)) 
	      {
		fscanf(ifp, "%s", dump);
		fscanf(ifp, "%s", dump);
		fscanf(ifp, "%lf", 
		       &particles[i].partic[rol].mass);
		fscanf(ifp, "%lf", 
		       &particles[i].partic[rol].width);
		fscanf(ifp, "%lf",
		       &particles[i].partic[rol].E);
		match = 1;
	   
	      }
	}

      else if(i >= 0 && i <= 2)  
	{
	  if(!strcmp(particles[i].partic[0].nickname,tmpnickname)) 
	    {
	      /* 
	       * Although this Born term is not included in coupl.iso, store 
	       * its mass and width (necessary for conversions, normalizations)
	       */
	      
	      fscanf(ifp, "%s", dump);
	      fscanf(ifp, "%s", dump);
	      fscanf(ifp, "%lf", 
		     &particles[i].partic[0].mass);
	      fscanf(ifp, "%lf", 
		     &particles[i].partic[0].width);
	      fscanf(ifp, "%lf",
		     &particles[i].partic[0].E);
	      match = 1;
	    }
	}

      
      if(!match) 
	{
	  /* Not in coupl.iso, no Born term: skip */
	  for(k=0; k<5; k++)
	    fscanf(ifp, "%s", dump);
	}


      fscanf(ifp, "%s", classlabel);
    }

  // Finally set the spin of all particles
  for ( i=0; i<CLASSMAX; i++ ) {
    if(i == 1 || i == 18)
      for ( k=0; k<PARTICLEMAX; k++ ) particles[i].partic[k].spin = 0.0;
    else if(i == 4 || i == 5 || i == 19 || i == 20)
      for ( k=0; k<PARTICLEMAX; k++ ) particles[i].partic[k].spin = 1.0;
    else if(i == 0 || i == 2 || i == 3 || (i >= 6 && i < 10))
      for ( k=0; k<PARTICLEMAX; k++ ) particles[i].partic[k].spin = 1.0/2.0;
    else if(i >= 10 && i < 14)
      for ( k=0; k<PARTICLEMAX; k++ ) particles[i].partic[k].spin = 3.0/2.0;
    else if(i >= 14 && i < 18)
      for ( k=0; k<PARTICLEMAX; k++ ) particles[i].partic[k].spin = 5.0/2.0;
    else error_exit("Spin for particles exceeding CLASSMAX is not defined!");
  }

  fclose(ifp);

  return 0;
}

/*! 
 * Reads the photocoupling helicity amplitudes 
 * for some transitions "N* -> p" and "N* -> n" (N* = N1,N2,...N8) from the 
 * file "./input/numinput/phohelapm.in". Starting from these helicity
 * amplitudes, the ratio of the magnetic transition moments for this N* to 
 * the proton and the neutron are calculated according to:
 *   
 * Spin 1/2 N*:
 *\verbatim
     Kappa_{n,N*}    A^n_{1/2} 
     ------------ =  --------
     Kappa_{p,N*}    A^p_{1/2}
 \endverbatim
 *
 * Spin 3/2 N*:
 *\verbatim
     Kappa^1_{n,N*}    sqrt(3) * A^n_{1/2} +/- A^n_{3/2} 
     -------------  =  ---------------------------------
     Kappa^1_{p,N*}    sqrt(3) * A^p_{1/2} +/- A^p_{3/2} 
 
     Kappa^2_{n,N*}    sqrt(3) * A^n_{1/2} - M_p/M_N* * A^n_{3/2} 
     -------------  =  -----------------------------------------
     Kappa^2_{p,N*}    sqrt(3) * A^p_{1/2} - M_p/M_N* * A^p_{3/2} 
 \endverbatim
 *
 * For Delta resonances, we have Kappa_{n,D*} = Kappa_{p,D*} etc,
 * so all ratios for Deltas are set to equal to 1.
 *
 * The function also checks if the "./input/numinput/phohelapm.in" file is
 * correct, i.e. if all N*'s of class D,E,H,I have a conversion 
 * coefficient != 0, so the magn. moment in the base isospin channel
 * can be converted to another isospin channel.
 */
int magnetic_transition_ratio(Class particles[], Observable* observ)
{
  FILE* ifp;
  int i, match, diagram, particle;
  fpos_t startposition;
  char dump[DUMP_SIZE], tmpnickname[4], filename[MAX_LOCATION_STRING];
  double a_n_1, a_p_1, a_n_3, a_p_3, mp_mNstar, sign;
  
  /* Open the input file with the helicity amplitudes */
  strcpy(filename, observ->inFolder);
  strcat(filename, "numinput/phohelamp.in");
  ifp = fopen(filename, "r");
  if(ifp == NULL)
    error_exit("File \"phohelamp.in\" is not found!!\n");


  /* dump the header of the phohelamp.in-file (6 column titles, 6 lines) */
  for(i=0; i<12; i++)
    fscanf(ifp, "%s", dump);

  /* store this starting position */
  fgetpos( ifp, &startposition );
  
  /* ------------------------------------------------------------------ *
   * N* resonances                                                      *
   *                                                                    *
   * Loop over all spin-1/2 N* particles and find the corresponding     *
   * helicity amplitudes in the phohelamp.in file.                      *
   * ------------------------------------------------------------------ */
  
  for(diagram=6; diagram<CLASSMAX; ++diagram) {
    for(particle = 0; particle < particles[diagram].particount; ++particle) 
    {
      if( particles[diagram].partic[particle].nickname[0] == 'D' ) {
	/* Delta resonances have isospin independent EM couplings */
	if(diagram == 6 || diagram == 7) {
	  particles[diagram].partic[particle].r_kappa_n_p = 1.;
	} else if(diagram == 10 || diagram == 11 
			 || diagram ==14 || diagram ==15 ) //FIXME!!
	{
	  particles[diagram].partic[particle].r_kappa_1_n_p = 1.;
	  particles[diagram].partic[particle].r_kappa_2_n_p = 1.;
	}
	else if(!(diagram >= 12)) {
	  error_exit("Unknown type of delta resonance.\n");
	}
      }
      
      else {
	/* reposition to the first relevant line */
	fsetpos( ifp, &startposition );

	/* Read nicknames untill we find the relevant one */
	match = 0;
	fscanf(ifp, "%s", dump);
      
	while( strcmp( dump, "X" ) && !match ) {
	  fscanf(ifp, "%s%lf%lf%lf%lf", 
		 tmpnickname, &a_n_1, &a_p_1, &a_n_3, &a_p_3);
	  
	  if(!strcmp(particles[diagram].partic[particle].nickname, tmpnickname))
	    match=1;
	  else
	    fscanf(ifp, "%s", dump);
	}
         
	/* Calculate the correct conversion coefficient.
	   This depends on the spin of the resonance. */
	if(diagram == 6 || diagram == 7) {
	  /* D and E diagrams: spin 1/2 */

	  if(!match) error_exit("Some helicity amplitudes are missing!\n");
	
	  particles[diagram].partic[particle].r_kappa_n_p = a_n_1 / a_p_1;
	} // D and E diagrams

	else if(diagram == 10 || diagram == 11) {
	  /* H and I: spin 3/2 */

	  if(!match) error_exit("Some helicity amplitudes are missing!\n");
	  
	  /* m_p/m_{N*} */
	  mp_mNstar = M_P / particles[diagram].partic[particle].mass; 
	  
	  if(diagram == 10) sign = 1; /* H: 3/2^+ */
	  else sign = -1;             /* I: 3/2^- */
	  
	  particles[diagram].partic[particle].r_kappa_1_n_p = 
	    (sqrt(3) * a_n_1 + sign * a_n_3) / 
	    (sqrt(3) * a_p_1 + sign * a_p_3);
	  
	  particles[diagram].partic[particle].r_kappa_2_n_p = 
	    (sqrt(3) * a_n_1 - mp_mNstar * a_n_3) / 
	    (sqrt(3) * a_p_1 - mp_mNstar * a_p_3);
	  
	} // H and I diagrams

	else if( diagram == 8 || diagram == 9 || diagram == 12 || diagram == 13 ) {
	  /* diagrams F, G, J and L */
	  printf("INFO in magnetic_transition_ratio(Class*,Observable*): %s %s",
		 "we assume the hyperon's EM coupling constants to be",
		 "isospin independent.\n");
	  if( diagram == 8 || diagram == 9 )
	    particles[diagram].partic[particle].r_kappa_n_p = 1.;
	  else {
	    particles[diagram].partic[particle].r_kappa_1_n_p = 1.;
	    particles[diagram].partic[particle].r_kappa_2_n_p = 1.;
	  }
	}
	else {
	  fprintf(stderr,
		  "INFO in magnetic_transition_ratio(Class*,Observable*): %s %s.\n",
		  "No conversion rules for EM coupling constants of diagram type",
		  dump);
	}
      } // end nucleon resonances
    } // loop over all particles
  } // loop over all non-Born diagrams
  
  fclose(ifp);

  return 0; 
}

/*!
 * There are two independent sets of coupling constants, one for the
 * (g + p -> L + K+) process (read from the coupl.iso.1 files), and one 
 * for the (g + p -> S0 + K+) process (read from coupl.iso.2).
 * From these two sets, one can construct the sets of coupling constants 
 * for the four other processes:
 * - Supposing isospin symmetry at the hadronic vertices, the coefficients 
 *   for conversion p to n, and S0 to S+ and S-, are just Clebsch-Gordans. 
 * - For the EM vertices, we use the ratios of measured photocoupling 
 *   helicity amplitudes (for baryon magn. moments) and of measured decay 
 *   widths (for mesons) to determine conversion factors.
 *
 * Remark that in some channels (3 and 6: charged S production), particle 
 * A is thrown away, because there is no extended Born term!
 *
 * In principle the coupl.iso.* files contain all exchanged particles for the two
 * base isospin classes. And all other isospin cases can be derived using that
 * information.
 * This does not work in practice, e.g. the K0S+ channel needs a K*(1410) in the
 * t channel, this diagram is not needed in it's isospin brother. This implies that
 * in the coupl.iso.* files the user can add diagrams that are not needed in that base
 * class. 
 * This could result in errors. Therefore we will first try to check that all 
 * exchanged particles are compatible with the base class. This check is not exhaustive
 * and the user should always double check the input files!!!
 *
 * For the moment we only check one thing:
 * - some vector mesons are needed in  K0 production channels,
 *   but are not needed in the base classes. Therefore we remove these particles from
 *   the particles[] structure for the K+ production channels.
 */
int change_coupling_constants(Class particles[], Observable* observ, int iso)
{
  double kappa_p, kappa_s0, kappa_s, kappa_u, r_kappa_sig_res, 
    clebsch, clebsch_n, clebsch_d;
  int i, j;

  // No A-channel for charged sigma production
  if(observ->iso.isospin == 3 || observ->iso.isospin == 6) 
      particles[3].particount = 0;

  // Remove neutral vector mesons from charged-kaon production channels
  if(observ->iso.isospin == 1 || 
     observ->iso.isospin == 2 ||
     observ->iso.isospin == 6 )
    {
      /* Remove all B3,K10 and C2 particles from the B or C channel
       * by putting all coupling constants to zero
       * This could be done by actualy removing the 
       * Properties from the partic[] array. */
      int particle;

      for(particle=0; particle<particles[4].particount; particle++)
	{
	  if( !strcmp(particles[4].partic[particle].nickname,"B3") )
	    {
	      particles[4].partic[particle].G = 0.0;
	      particles[4].partic[particle].H = 0.0;
	      particles[4].partic[particle].I = 0.0;
	    }
	}
      for(particle=0; particle<particles[5].particount; particle++)
	{
	  if( !strcmp(particles[5].partic[particle].nickname,"K10") ||
	      !strcmp(particles[5].partic[particle].nickname,"C2") )
	    {
	      particles[5].partic[particle].G = 0.0;
	      particles[5].partic[particle].H = 0.0;
	      particles[5].partic[particle].I = 0.0;
	    }
	}
    }


  /* Make sure the proton and S^+ magnetic moment are
   * available, even when particles[0,2].partic[0].G
   * are 0 (e.g. if the S or U channel does not occur
   * in coupl.iso, or G is just 0 there)
   * (PS: particles[0].partic[0].G = kappa_p now; it is only 
   * changed into kappa_N for iso>3 later in this function) */

  if (particles[0].particount == 0 || 
      particles[0].partic[0].G == 0) 
    kappa_p = KAPPA_P;
  else
    kappa_p = particles[0].partic[0].G;


  if (particles[2].particount == 0 || 
      particles[2].partic[0].G == 0) 
    kappa_s0 = KAPPA_S0;
  else
    kappa_s0 = particles[2].partic[0].G;




  /*---------------------------------------- *
   * L production on the neutron             *
   *                                         *
   * Transform the set of coupling constants *
   * for the g + p -> K+ + L0 process        *
   * into a set for:                         *
   *                                         *
   * g + n -> K0 + L0                        *
   *                                         *
   *---------------------------------------- */

  if(iso == 4) 
    {
      /* 
       * Born terms: 
       * - Main strong coupling constant g_{K,L,n} 
       *   remains unchanged since g_{K+,L,p} = g_{K0,L,n}.
       * - Elastic magnetic moments: Kappa_{L0} remains unchanged, 
       *                             Kappa_{p} changes into Kappa_{n}
       * => only s-channel is modified.
       */

      // S channel
      for(j=0; j<particles[0].particount; j++) 
	{
	  /* S channel:
	   * - G: Kappa_n
	   * - H: g_{K0,L,n} */
	  particles[0].partic[j].G = KAPPA_N;
	}

      // A channel
      for(j=0; j<particles[3].particount; j++)
	{
	  /* A-particle: Y' = S0 exch -> cc's =
	   * - G: Kappa_{L,S0}
	   * - H: g_{K0,S0,n}
	   * Hadr:  g(K0,S0,n) = - g_{K+,S0,p} 
	   * EM:    unchanged */ 
	  particles[3].partic[j].H *= -1.0;
	}

      // B channel
      for(j=0; j<particles[4].particount; j++)
	{

	  if (!strcmp(particles[4].partic[j].nickname, "B2"))
	    {
	      /* B-particle: K*0(892) exch -> cc's =
	       * - G: Kappa_{K0,K*0}
	       * - H: g_{K*0,L,n}^v
	       * - I: g_{K*0,L,n}^t
	       * Hadr:  g_{K*0,L,n}^{v,t} = g{K*+,L,p)^{v,t}
	       * EM:    Kappa_{K*0,K0}= Kappa_{K*+,K+} * KSP2KS0 */
	      
	      particles[4].partic[j].G *= KSP2KS0;
	    }
	  else if (!strcmp(particles[4].partic[j].nickname, "B3"))
	    {
	      /* B-particle: K*0(1410) exch -> cc's =
	       * - G: Kappa_{K0,K*0}
	       * - H: g_{K*0,L,n}^v
	       * - I: g_{K*0,L,n}^t
	       * Hadr:  g_{K*0,L,n}^{v,t} = g{K*+,L,p)^{v,t}
	       * EM:    Kappa_{K*0,K0}= Kappa_{K*+,K+} * B3P2B30 */
	      
	      particles[4].partic[j].G *= B4toB3;
	    }
	  else /* No other particles are implemented */
	    {
	      printf("%s%s","Conversion factor for EM coupling constant ",
		     "of particles in B channel are needed!\n");
	      exit(1);
	    }
	}

      // C channel
      for(j=0; j<particles[5].particount; j++)
	{
	  if (!strcmp(particles[5].partic[j].nickname, "K10"))
	    {
	      /* C-particle: K1(1270) exch -> like K* */
 
	      particles[5].partic[j].G *= K1P2K10; 
	      
	    }
	  else if (!strcmp(particles[5].partic[j].nickname, "C2"))
	    {
	      /* C-particle: K1(1400) exch -> like K* */

	      particles[5].partic[j].G *= K1400P2K14000; 
	    
	    }
	  else /* No other particles are implemented */
	    {
	      printf("%s%s","Conversion factor for EM coupling constant ",
		     "of particles in C channel are needed!\n");
	      exit(1);
	    }
	}


      for(i=6; i<=13; i++) /* resonances: N*,L*,S* 
			    * (in Lambda production, no Delta resonances 
			    * due to isospin conservation at the (D*,K,L) 
			    * vertex!) */

	for(j=0; j<particles[i].particount; j++)

	  if(i == 6 || i == 7) /* D, E */
	    {
	      /* D,E: spin-1/2 N* resonances -> cc's = 
	       * - G: kappa_{n,N*0}
	       * - H: g_{K0,L,N*0}
	       * Hadr:  g_{K0,L,N*0} = g_{K+,L,N*+}
	       * EM:    Kappa_{n,N*0} = Kappa{p,N*+} * r_kappa_n_p 
	       */
	      
	      particles[i].partic[j].G *= particles[i].partic[j].r_kappa_n_p;
	    }
	  else if(i == 8 || i == 9) /* F,G */
	    {
	      /* F,G: spin-1/2 Y* resonances -> cc's =
	       * - G: kappa_{L,Y*0)}
	       * - H: g_{K0,Y*0,n}
	       * Hadr:  Y* = Lambda*: g_{K0,L*,n} = g_{K+,L*,p}
	       *        Y* = Sigma*0: g_{K0,S*,n} = - g_{K+,S*,p}
	       * EM: unchanged */

	      if(particles[i].partic[j].nickname[0] == 'S')
		particles[i].partic[j].H *= -1;
	    }
      	  else if(i == 10 || i == 11) /* H,I*/
	    {
	      /* H,I: spin-3/2 N* resonances -> cc's =
	       * - G: kappa^(1)_{n,N*0}
	       * - H: kappa^(2)_{n,N*0}
	       * - I: f_{K0,L,N*0}
	       * same as above, but with f<->g and (1/2) with the Kappa's */
	
	      particles[i].partic[j].G *= 
		particles[i].partic[j].r_kappa_1_n_p;
	      particles[i].partic[j].H *= 
		particles[i].partic[j].r_kappa_2_n_p;
	    }
	  else if(i == 12 || i == 13) /* J,L */
	    {
	      if(particles[i].partic[j].nickname[0] == 'S')
		{
		  /* J,L: spin-3/2 Y* resonances
		   * same as above, but with f<->g and (1/2) for Kappa's */
		  
		  particles[i].partic[j].I *= -1;
		}
	    }
    }


  /*---------------------------------------- *
   * Sigma production                        *
   *                                         *
   * Transform the set of coupling constants *
   * for the g + p -> K+ + S0 process        *
   * into a set for:                         *
   *                                         *
   * 3: g + p -> K0 + S+                     *
   * 5: g + n -> K0 + S0                     *
   * 6: g + n -> K+ + S-                     *
   *                                         *
   *---------------------------------------- */

  else if(iso == 3 || iso == 5 || iso == 6) 
    {

      /* 
       * "Switch" statement: fix the values of conversion coefficients 
       * for hadronic and EM vertex, for each of the 3 channels.
       *
       *
       * Hadronic: for each isospin channel we need only two conversion 
       * relations, i.e. 
       *
       * g_{K?(*),S?(*),N?(*)} = clebsch_n * g_{K+(*),S0(*),N+(*)}
       * g_{K?(*),S?(*),D?(*)} = clebsch_d * g_{K+(*),S0(*),D+(*)}
       * where clebsch_{n,p} depend only on the isospin (I,I3) of the different
       * particles, i.e. on the channel (3, 5 or 6)
       * 
       *
       * EM: 
       * - Elastic magnetic moments change:  
       *   Kappa_{p} changes into kappa_s = Kappa_{n} for channel 5,6
       *   Kappa_{S0} changes into kappa_u = Kappa_{S+/-} for channel 3/6.
       *
       * - Transition magnetic moments change:
       *   the ratios kappa_{N*n}/kappa{N*p} are determined by helicity 
       *   amplitudes of the transitions N*->{p/n}+g. For the ratios 
       *   kappa_{Y*S+/-}/kappa_{Y*S0}, those amplitudes are not well known,
       *   so there we use the recipe:
       *
       *                    mu_{S0}                  mu_{S0}
       *   kappa_{S0*,S0} = ------- kappa_{S+*,S+} = ------- kappa_{S-*,S-}
       *                    mu_{S+}                  mu_{S-}
       */

      switch(iso) {

      
	/* ---------------- *
	 * g + p -> K0 + S+ *
	 * ---------------- */
      case 3 : 

	/* strong */
	clebsch_n = sqrt(2);
	clebsch_d = -1/sqrt(2);

	/* EM */
	kappa_s = kappa_p; /* Kappa_p */
	kappa_u = KAPPA_SP; /* Kappa_{S+} */

	r_kappa_sig_res = (KAPPA_SP + 1.0) / kappa_s0;

	/* HAS BEEN MODIFIED! The previous line used to be:
	 * "r_kappa_sig_res = particles[0].partic[0].H / (KAPPA_SP + 1.0);"
	 * i.e. kappa_{S*,S+} / kappa_{S*,S0} = kappa_p / mu_{S+}
	 *      ??? why kappa_p instead of kappa_{S0} ??? 
	 * NB: This change_isospin function is called before the cc's are
	 * normalized, so H is just kappa_{S0}. */

	break;
	
	/* ---------------- *
	 * g + n -> K0 + S0 *
	 * ---------------- */
      case 5 :
	
	/* strong (only for S, not for L resonances!!!) */
	clebsch_n = -1;
	clebsch_d = 1; 

	/* EM */
	kappa_s = KAPPA_N; /* Kappa_n */
	kappa_u = kappa_s0; /* Kappa_{S0} */
	r_kappa_sig_res = 1.0; 

	break;

	/* ---------------- *
	 * g + n -> K+ + S- *
	 * ---------------- */
      case 6 :

	/* strong */
	clebsch_n = sqrt(2);
	clebsch_d = 1/sqrt(2);

	/* EM */
	kappa_s = KAPPA_N; /* Kappa_n */
	kappa_u = KAPPA_SM; /* Kappa_{S-} */


	r_kappa_sig_res = (KAPPA_SM - 1.0) / kappa_s0;

	/* HAS BEEN MODIFIED! The previous line used to be:
	 * "r_kappa_sig_res = particles[0].partic[0].H / (KAPPA_SM - 1.0);"
	 * kappa_{Y*,S+/-} / kappa_{Y*,S0} = kappa_n / mu_{S-} (>0) */
	/* ??? why kappa_n instead of kappa_{S0} ??? */

	break;
      }
  

      /* ---------------------------------- *
       * Recalculate the coupling constants *
       * ---------------------------------- */

      for(j=0; j<particles[0].particount; j++)
	{
	  /* S channel: N = p/n exch -> cc's = 
	   * - G: Kappa_{N?}
	   * - H: g_{K?,S?,N?}
	   */	  

	  particles[0].partic[j].H *= clebsch_n; 
	  /* HAS BEEN MODIFIED */ /* Correct??? */
	  if(particles[0].partic[j].G != 0)
	    particles[0].partic[j].G = kappa_s;
	}

	
      for(j=0; j<particles[1].particount; j++)
	{
	  /* T channel: K=K+/K0 exch -> cc's =
	   * - H: g_{K?,S?,N?}
	   */

	  particles[1].partic[j].H *= clebsch_n;
	}

	 
      for(j=0; j<particles[2].particount; j++)
	{
	  /* U channel: Y=S0,S+,S- exch -> cc's =
	   * - G: Kappa_{S?}
	   * - H: g_{K?,S?,N?}
	   */

	  particles[2].partic[j].H *= clebsch_n;
	  /* HAS BEEN MODIFIED */
	  if(particles[2].partic[j].G != 0)
	    particles[2].partic[j].G = kappa_u;
	}
    

      if(iso == 3 || iso == 6)
	{
	  /* A-particle: Y'= L exch; 
	   * remains in channel 5 (S0 prod), vanishes in 3,6 -> cc's: 
	   * - G: kappa_{L,S0}
	   * - H: g_{K?,L,N?}
	   * Hadr:  g_{K0,L0,n} = g{K+,L0,p} 
	   * EM:    unchanged */
      
	  particles[3].particount = 0;
	}
      

      for(j=0; j<particles[4].particount; j++)
	{

	  /* B-particle: K* exchange -> cc's = 
	   * - G: Kappa_{K?,K*?}
	   * - H: g_{K*?,S?,N?}^v
	   * - I: g_{K*?,S?,N?}^t
	   * Hadr:  in all cases g_{K*?,S?,N?} = g_{K*+,S0,p} * clebsch_n.
	   * EM:    in channel 3 and 5 -> K0 production instead of K+ 
	   *        so Kappa_{K0,K*0} = Kappa_{K+,K*+} * KSP2KS0.
	   *
	   *        in channel 3 a K*(1410) is possible. This exchange
	   *        does not occur in the base class and the given
	   *        coupling constant is Kappa_(KO,K*0(1410)) already
	   */

	  particles[4].partic[j].H *= clebsch_n;
	  particles[4].partic[j].I *= clebsch_n;

	  /* K0 production */
	  if(iso == 3 || iso == 5)
	    {
	      /* K*(892) */
	      if (!strcmp(particles[4].partic[j].nickname,"B2") )
		{
		  particles[4].partic[j].G *= KSP2KS0;
		}
	      
	      /* K*(1410) */
	      else if (!strcmp(particles[4].partic[j].nickname,"B3") )
		{
		  particles[4].partic[j].G *= B4toB3;
		}
	      
	      /* No other particles are implemented */
	      else 
		{
		  printf("%s%s","Conversion factor for EM coupling constant ",
			 "of particles in B channel are needed!\n");
		  exit(1);
		}
	    }
	}

      for(j=0; j<particles[5].particount; j++)
	{

	  /* C-particle: ISOBAR: K1^+(1270) = K1^0(1270) exch
	   *               -> analogue to particles B 
	   */

	  particles[5].partic[j].H *= clebsch_n;
	  particles[5].partic[j].I *= clebsch_n;

	  /* K0 production */
	  if(iso == 3 || iso == 5)
	    {
	      /* K1(1270), only in isobar model! */
	      if (!strcmp(particles[5].partic[j].nickname,"K10"))
		particles[5].partic[j].G *= K1P2K10;
 
	      /* K1(1400) exch */
	      else if (!strcmp(particles[5].partic[j].nickname, "C2"))
		particles[5].partic[j].G *= K1400P2K14000; 
	      
	      /* No other particles are implemented */
	      else  
		{
		  printf("%s%s","Conversion factor for EM coupling constant ",
			 "of particles in C channel are needed!\n");
		  exit(1);
		}
	    }
	}
	
          
      
      for(i=6; i<CLASSMAX; i++)  /* D,E,F,G,H,I,J,L */
	for(j=0; j<particles[i].particount; j++)
	  {
	    /* Determine if the resonance is a delta. */
	    if(particles[i].partic[j].nickname[0] == 'D')
	      clebsch = clebsch_d;
	    else
	      clebsch = clebsch_n;
	    
	    
	    if(i == 6 || i == 7)
	      {
		/* D,E: spin-1/2 N*,D* resonances  -> cc's = 
		 * - G: kappa_{N?,(N*,D*)?}
		 * - H: g_{K?,S?,(N*,D*)?}
		 */
		
		/* Hadr:  in all cases g_{K?,S?,(N*,D*)?} 
		 *        = g_{K+,S0,(N*,D*)+} * clebsch_{n,d} */
		particles[i].partic[j].H *= clebsch;
		
		/* EM:    in channel 5 and 6 (prod. off a n), Kappa_{N*0,n} 
		 * = Kappa{N*+,p} * r_kappa_n_p (with A_{1/2}, A_{3/2}). */
		if(iso == 5 || iso == 6)
		  particles[i].partic[j].G *= particles[i].partic[j].r_kappa_n_p;

	      }
	    else if(i == 10 || i == 11)
	      {
		/* H,I: spin-3/2 (N*,D*) resonances -> cc's =
		 * - G: kappa^(1)_{N?,(N*,D*)}
		 * - H: kappa^(2)_{N?,(N*,D*)}
		 * - I: f_{K?,S?,(N*,D*)?}
		 * same as above, but with f<->g and (1/2) with the Kappa's */

		particles[i].partic[j].I *= clebsch;

		if(iso == 5 || iso == 6)
		  {
		    particles[i].partic[j].G *= particles[i].partic[j].r_kappa_1_n_p;
		    particles[i].partic[j].H *= particles[i].partic[j].r_kappa_2_n_p;
		  }
	      }

	    /* Careful for hyperon resonances!
	     * 
	     * This used to be implemented in Stijns's and Tamara's code. Since
	     * we do not need hyperon exchange for t-channel Regge I've removed 
	     * this part of the code.
	     */

	    else if(i == 8 || i == 9)
	      {
		/* Spin 1/2 Y^* */
		
		printf("%s%s","Changing isospin channels with hyperon resonances ",
		       "is not implemented!\n");
		exit(1);
	      }

	    else if(i == 12 || i == 13)
	      {		
		/* Spin 3/2 Y^* */

		printf("%s%s","Changing isospin channels with hyperon resonances ",
		       "is not implemented!\n");
		exit(1);
	      }
	  }
    }

  /*---------------------------------------- *
   * pi^- production                         *
   *                                         *
   * Transform the set of coupling constants *
   * for the g + p -> pi+ + n process        *
   * into a set for the g + n -> pi- + p     *
   * process.                                *
   *                                         *
   *---------------------------------------- */

  if(iso == 10)
    {
      for(j = 0; j < particles[18].particount; j++)
	particles[18].partic[j].E *= -1.;

      for(j = 0; j < particles[20].particount; j++)
	particles[20].partic[j].I *= -1.;
    }

  return 0;
}

/*!
 * Determines the index i (0...13) in the particles[] array of 
 * Classes from the "classlabel" character (S,T,U,A,B,...,J,L).
 */ 
int findclassindex(char classlabel)
{
  int i;
  
  /* index = 0,1,2 = Born terms */
  if(classlabel == 'S' || classlabel == 's' || classlabel == 'p' )
    i = 0;
  else if(classlabel == 'T' || classlabel == 't' || classlabel == 'k')
    i = 1;
  else if(classlabel == 'U' || classlabel == 'u' || classlabel == 'y')
    i = 2;
  
  /* there is a 'L' class, but no 'K' class */
  else if(classlabel == 'L' || classlabel == 'l') 
    i = 13; 
  else if(classlabel == 'M' || classlabel == 'm') 
    i = 14; 
  else if(classlabel == 'N' || classlabel == 'n') 
    i = 15; 
  else if(classlabel == 'O' || classlabel == 'o') 
    i = 16; 
  else if(classlabel == 'Q' || classlabel == 'q') 
    i = 17; 
  else if(classlabel == 'R') // pi/b1 t-channel 
    i = 18; 
  else if(classlabel == 'V') // rho/a2 t-channel
    i = 19; 
  else if(classlabel == 'W') // a1 t-channel
    i = 20; 
  
  /* index = 3...12 */
  else if(isupper(classlabel)) 
    {
      /* classlabel = 'A'...'J'; 'A' = 65 */
      i = classlabel - 65 + 3; 
    }
  else 
    {
      /* classlabel = 'a'...'j'; 'a'= 97 */
      i = classlabel - 97 + 3; 
    }
  
  return i;
}

/*!
 * Determines the "classlabel" character from the index i 
 * in the particles[] array of Classes.
 */ 
char findclasslabel(int index)
{
  char label;
  
  if(index == 0)
    label = 'S';
  else if(index == 1)
    label = 'T';
  else if(index == 2)
    label = 'U';
  else if(index == 13)
    label = 'L';
  else if(index == 14)
    label = 'M';
  else if(index == 15)
    label = 'N';
  else if(index == 16)
    label = 'O';
  else if(index == 17)
    label = 'Q';
  else if(index == 18) // pi/b1 t-channel 
    label = 'R';
  else if(index == 19) // rho/a2 t-channel
    label = 'V';
  else if(index == 20) // a1 t-channel
    label = 'W';
  else /* if index = 3 -> 'A', etc */
    label = 'A' + index - 3; 
  
  return label;
}

/*!
 * Print the stored particle information to file (ofp),
 * as read from the coupl.iso.* and mass.new files: class, 
 * nickname, cc's, offshell par's, born and resonance cutoffs.
 */
int fprintparticspec(FILE* ofp, Class particles[], Observable observ)
{
  int i,j;
//   printf("\n%s\n","--------------------------------------------------");
//   printf("%s\n\n","Exchanged particles and their coupling constants");

  for(i=0; i<CLASSMAX; i++)
    if(particles[i].particount != 0)
      for(j=0; j<particles[i].particount; j++)
	{
	  fprintf(ofp, "%-5c%-4s", findclasslabel(i), 
		  particles[i].partic[j].nickname);
	  
	  /*
	  printf("%.4e  %.4e", particles[i].partic[j].mass, 
	  particles[i].partic[j].width);
	  */
	  fprintf(ofp, "%15.6e%15.6e%15.6e", particles[i].partic[j].G, 
		  particles[i].partic[j].H,particles[i].partic[j].I);
	  fprintf(ofp, "%15.6e%15.6e%15.6e", particles[i].partic[j].X, 
		  particles[i].partic[j].Y, particles[i].partic[j].Z);
	  fprintf(ofp, "\n");
  	  /*
	  printf("%12.3e%12.3e%12.3e", particles[i].partic[j].r_kappa_n_p, 
	  particles[i].partic[j].r_kappa_1_n_p, 
	  particles[i].partic[j].r_kappa_2_n_p);
	  */
	}

  if(observ.hadronformfac)
    fprintf(ofp,"cutoff born:\t%.6f\ncutoff res:\t%.6f\n",  
	    observ.ffac.born_cutoff, observ.ffac.res_cutoff);
  else	
    fprintf(ofp,"cutoff born:\tX\ncutoff res:\tX\n");
 
//   printf("%s\n","--------------------------------------------------");
   
  return 0;
}

/*! 
 * Checks if the Born coupling constants g_{KYN} are equal for the three Born
 * channels.
 */
int check_born_cc(Class particles[])
{
  int s_ch, t_ch, u_ch;

  s_ch = particles[0].particount;
  t_ch = particles[1].particount;
  u_ch = particles[2].particount;

  if(s_ch + t_ch + u_ch > 1) /* at least 2 of the 3 Born-channels exist */
    {
      if(s_ch && t_ch) /* s and t equal? */
	{
	  if(particles[0].partic[0].H != particles[1].partic[0].H)
	    {
	      printf("s-channel cc: %f <--> t-channel cc: %f\n", 
		     particles[0].partic[0].H, particles[1].partic[0].H);
	      error_exit("Differences in Born cc !!\n");
	    }
	}
      
      if(s_ch && u_ch) /* s and u equal? */
	{
	  if(particles[0].partic[0].H != particles[2].partic[0].H)
	    {
	      printf("s-channel cc: %f <--> u-channel cc: %f\n", 
		     particles[0].partic[0].H, particles[2].partic[0].H);     
	      error_exit("Differences in Born cc !!\n");
	    }
	}    
      
      if(t_ch && u_ch) /* t and u equal? */
	{
	  if(particles[1].partic[0].H != particles[2].partic[0].H)
	    {
	      printf("t-channel cc: %f <--> u-channel cc: %f\n", 
		     particles[1].partic[0].H, particles[2].partic[0].H);
	      error_exit("Differences in Born cc !!\n");
	    }
	}
    }
  
  return 0;

}

/*!
 * We want to store the mass of the struck nucleon and the produced
 * kaon and hyperon in mp,mk and my respectively. This can be done
 * with:
 * mp = particles[0].partic[0].mass;
 * mk = particles[1].partic[0].mass;
 * my = particles[2].partic[0].mass;
 *
 * If the born masses are needed however before "mass.in" has been
 * read in, this function gets called
 */
int allocate_born_masses(double* mp, double* mk, double* my,int iso)
{
  // g p -> K+ L
  if (iso==1)
    {
      *mp = M_P;
      *mk = M_KP;
      *my = M_L;
    }
  // g p -> K+ S0
  else if (iso==2)
    {
      *mp = M_P;
      *mk = M_KP;
      *my = M_S0;
    }
  // g p -> K0 S+
  else if (iso==3)
    {
      *mp = M_P;
      *mk = M_K0;
      *my = M_SP;
    }
  // g n -> K0 L
  else if (iso==4)
    {
      *mp = M_N;
      *mk = M_K0;
      *my = M_L;
    }
  // g n -> K0 S0
  else if (iso==5)
    {
      *mp = M_N;
      *mk = M_K0;
      *my = M_S0;
    }
  // g n -> K+ S-
  else if (iso==6)
    {
      *mp = M_P;
      *mk = M_KP;
      *my = M_SM;
    }
  // g + p   ->  P0 + p
  else if (iso==7)
    {
      *mp = M_P;
      *mk = M_PI0;
      *my = M_P;
    }
  // g + p   ->  P+ + n
  else if (iso==8)
    {
      *mp = M_P;
      *mk = M_PIP;
      *my = M_N;
    }  
  // g + n   ->  P0 + n
  else if (iso==9)
    {
      *mp = M_N;
      *mk = M_PI0;
      *my = M_N;
    }
  // g + n   ->  P- + p
  else if (iso==10)
    {
      *mp = M_N;
      *mk = M_PIM;
      *my = M_P;
    }

  // K- p -> g L0
  else if (iso==11)
    {
      *mp = M_P;
      *mk = M_KP;
      *my = M_L;
    }
  // K- p -> g S0
  else if (iso==12)
    {
      *mp = M_P;
      *mk = M_KP;
      *my = M_S0;
    }
  // other isospin channels
  else
    {
      printf("\n%s\n","Isospin channel 9, 10 and >12 and not implemented in allocate_born_masses()");
      exit(1);
    }

  return 0;
}

/*!
 * Calculates epsilon, the transverse polarization of
 * the virtual photon, for a certain electron beam (lab) energy eps_1
 */
double calc_epsilon(double wlab, double klab, double qq, 
		    double e_beam_ener)
{
  double eps=1., costheta_e;
  
  if( qq>0. ) {
    /* theta_e: stijn notes (131) */
    costheta_e = 1.0 - qq / ( 2. * e_beam_ener * (e_beam_ener - wlab));

    if( fabs(costheta_e)-1. > 0. ) {
      fprintf(stderr,"ERROR in calc_epsilon(..): wrong kinematics!\n");
      exit(1);
    }
    
    /* epsilon: stijn notes (99) */
    eps = 1.0 / (1.0 + 2.*klab*klab/qq * (1.-costheta_e) / (1.+costheta_e));
  }
  
  return eps;

}

/*!
 * \return Kinematic coefficient in the photoproduction differential cross section.
 */
double kincoeff_photo(double w, double pk, double mp)
{
  return (pk / w) * 1./4. / ( pow((sqrt(mp*mp + w*w) + w), 2) *16.*PI*PI) ;

}

/*!
 * \return Kinematic coefficient in the radiative kaon capture diff. cross section.
 */
double kincoeff_kaoncap(double w, double pk, double mp, double mk)
{
  return (w / pk) * 1./4. / ( pow((sqrt(mp*mp+pk*pk) + sqrt(mk*mk+pk*pk)), 2) *16.*PI*PI);
}

/*! \return The spin averaging factor */
double spinaveraging(const Observable* observ)
{
  
  /* 
   * The spin averaging factor is determined by the type and number of 
   * particles in the process of which the polarization is measured. 
   *
   *
   * Some considerations about the spin averaging:
   *
   *        *******************
   *        * Photoproduction *
   *        *******************
   *
   * -> Unpolarized photoproduction (diff or total):
   *          averaging over 2 initial particle spins: photon and proton -> 1/4
   *          determine_M2_photo() computes <s>, code returns <S>
   *         
   *          so:            <S> = <s>/4 
   * 
   *
   * -> Polarized beam photoproduction:
   *          polarized cross section is defined as:
   *                         <S^+/-> = <S>_unpol +/- <G_5>
   *
   *          with  <S>_unpol = (<S^+> + <S^->)/2
   *          and <G_5> = (<S^+> - <S^->)/2
   *
   *          so the asymmetry is
   *                         <A> = (<S^+> - <S^->)/(<S^+> + <S^->)
   *                         <A> = (<S^+> - <S^->)/(2 * <S>)
   *                         <A> = (2 * <G_5>)/(2 * <S>_unpol)
   *                         <A> = <G_5>/<S>_unpol
   *
   *          determine_M2_photo() computes <g_5> = (<s^+> - <s^->),
   *          the code returns <G_5>
   *
   *          because of
   *          the averaging over initial proton spin ( <S^+/-> = <s^+/->/2 )
   *          and the factor 1/2 in the definition of <G_5> we get:
   *         
   *          so:            <G_5> = <g_5>/4

   *
   * -> Polarized target photoproduction:
   *          identical reasoning as for polarized beam photoproduction
   *
   * -> Polarized recoil photoproduction:
   *          polarized cross section is defined as:
   *                         <S^+/-> = 1/2*<S>_unpol +/- <G_5>
   *
   *          with  <S>_unpol = (<S^+> + <S^->)
   *          and <G_5> = (<S^+> - <S^->)/2
   *
   *          so the asymmetry is
   *                         <A> = (<S^+> - <S^->)/(<S^+> + <S^->)
   *                         <A> = (2 * <G_5>)/<S>_unpol
   *
   *          determine_M2_photo() computes <g_5> = (<s^+> - <s^->),
   *          the code returns <G_5>
   *
   *          because of
   *          the averaging over initial proton and photon spin ( <S^+/-> = <s^+/->/4 )
   *          and the factor 1/2 in the definition of <G_5> we get:
   *         
   *          so:            <G_5> = <g_5>/8
   *
   *
   * -> Polarized beam-target photoproduction:
   *          double polarized cross section is defined as:
   *                         <S^+/-|+> = <S^+/-> + <G_5^+/->
   *                         <S^+/-|-> = <S^+/-> - <G_5^+/->
   *          
   *          with <G_5^+/-> = ( <S^+/-|+> - <S^+/-|-> )/2
   *          and <S>_unpol = (<S^++> + <S^+-> + <S^-+> + <S^-->)/4
   *                        = <S^+> = <S^->
   *
   *          so the asymmetry is
   *                         <A> = (<S^++> - <S^-+>)/(<S^++> + <S^-+>)
   *                         <A> = (<S^+> + <G_5^+> - <S^+> + <G_5^+>)
   *                                 /(2*<S>_unpol)
   *                         <A> = (2 * <G_5^+>)/(2 * <S>_unpol)
   *                         <A> = <G_5^+>/<S>_unpol
   *
   *          determine_M2_photo() computes <g_5^+> = ( <s^+|+> - <s^+|-> )
   *          the code returns <G_5^+>
   *
   *          there's no averaging over inital spins
   *          but there is an extra factor 1/2 in the definition of <G_5^+>
   *
   *          so:            <G_5> = <g_5>/2
   *
   *
   * -> Polarized beam-recoil photoproduction:
   *          double polarized cross section is defined as:
   *                         <S^+/-|+> = 1/2*<S^+/-> + <G_5^+/->
   *                         <S^+/-|-> = 1/2*<S^+/-> - <G_5^+/->
   *          
   *          with <G_5^+/-> = ( <S^+/-|+> - <S^+/-|-> )/2
   *          and <S>_unpol = (<S^++> + <S^+-> + <S^-+> + <S^-->)/2
   *                        = <S^+> = <S^->
   *
   *          so the asymmetry is
   *                         <A> = (<S^++> - <S^-+>)/(<S^++> + <S^-+>)
   *                         <A> = (<S^+>/2 + <G_5^+> - <S^+>/2 + <G_5^+>)
   *                                 /<S>_unpol
   *                         <A> = (2 * <G_5^+>)/<S>_unpol
   *
   *          determine_M2_photo() computes <g_5^+> = ( <s^+|+> - <s^+|-> )
   *          the code returns <G_5^+>
   *
   *          because of
   *          the averaging over inital proton spin (factor 1/2)
   *          and an extra factor 1/2 in the definition of <G_5^+>
   *
   *          so:            <G_5> = <g_5>/4
   *                        
   *
   *
   * -> Polarized target-recoil photoproduction: not yet implemented !!!
   *
   *
   *        *********************
   *        * Electroproduction *
   *        *********************
   *
   *  -> Unpolarized electroproduction:
   *          averaging over initial electron and proton spin -> 1/4.
   *          snlink/ file contains <h>_ll'
   *
   *
   * -> Polarized electron beam electroproduction: 
   *          averaging over initial proton spin + extra fact. 1/2 from 
   *          (1+g_5)/2 in trace -> 1/4. 
   *          snlink/ file contains <h>_ll', code uses <h>_ll'/4 : 
   *          when chi (stijn notes 95) is calculated in the code, 
   *          a factor 1/4 is included instead of 1/16, so spinaverage 
   *          is not included in chi!!! 
   *          
   *          so: <ds/domega>_ll' ~ <h>_ll'/4   
   *
   *
   * -> Polarized target electroproduction: ???
   *          averaging over initial electron spin + factor 1/2 from 
   *          (1+g5)/2 in trace -> 1/4.
   *          snlink/ file contains <g5>_ll' (<h>_ll' with 
   *          [1 -> 1 +/- g_5(n.g)]), 
   *          code uses <G5>_ll'= <g5>_ll'/4:
   *          when chi (stijn notes 95) is calculated, a factor 1/4 is 
   *          included instead of 1/16, so spinaverage is not included 
   *          in chi !!! 
   *
   *          so: R^(0,*) = <G5>
   *
   *          ??? NOTE: polarized response functions (ie which do not occur in
   *          the unpolarized cross-section) are asymmetries!!
   *              This is in agreement with the relation:
   *          
   *                    R_T^(0,y) = T . R_T^(0,0)
   *           
   *              since:
   *           
   *                    R_T^(0,y) = 2 <G5> / (2 <S>) * <S>
   *                    R_T^(0,y) = <G5>
   * 
   *
   * -> Polarized recoil electroproduction: ???
   *          averaging over target and electron spin + 1/2 from 
   *          (1+g5)/2 -> 1/8
   *          (factor 1/4 is ALWAYS included in xi !!)
   *          snlink/ file contains <g5>_ll' (<h>_ll' with 
   *          [1 +/- g_5(n.g)]???),
   *          code uses <G5>_ll'= <g5>_ll'/4 to compute the cross-section
   *
   *          so:       R^(*,0) =  1/2 * <G5> ??? (is what I=Tamara thought)
   *
   *                    R^(*,0) =  2 * <G5>  with: -> <G5> = <g5>/8 
   *                    (what Stijn wrote; but maybe it's wrong??? see below)
   *                            
   *
   * -> Radiative capture: 
   *          Kaon is pseudoscalar -> no spin, only averaging over initial 
   *          proton spin -> 1/2 
   *
   *
   */
  
  if(observ->photoprod)
    {
      if(observ->photo.pol.nopol)
	return 1/4.;
      else if(observ->photo.pol.sinpol.phopol)
	return 1/4.;
      else if(observ->photo.pol.sinpol.tarpol)
	return 1/4.;
      else if(observ->photo.pol.sinpol.recpol)
	return 1/8.;
      else if(observ->photo.pol.doubpol.beamtar)
	return 1/2.;
      else if(observ->photo.pol.doubpol.beamrec)
	return 1/4.;
    }
  else if(observ->electroprod)
    {
      /* This could be wrong.
       * This is what my (=Pieter's) intuition tells me */
      return 1/4. ;
    }
  else if(observ->kaoncapture)
    return 1/2.;
  else
    error_exit("Error in \"spinaveraging\"!");

  return 0;
  
}

/*!
 * Lorentz transformation for photon lab-energy to C.M.-energy for the
 * photoproduction process.
 */
int lorentztrans_photo(double wlab, double* w, double mp)
{
  *w = wlab * mp / sqrt(pow((wlab + mp),2) - wlab * wlab);

  return 0;

}

/*!
 * Lorentz transformation for lab energy and momentum
 * to c.m. energy and momentum for the electroproduction process.
 */
int lorentztrans_electro(double wlab, double* w, double klab, 
			 double* k, double mp)
{
  *w = (wlab*(wlab+mp) - klab*klab) / sqrt(pow((wlab + mp),2) - klab*klab);
  
  *k = klab * mp / sqrt(pow((wlab + mp),2) - klab*klab);
  
  return 0;
 
}

/*!
 * Lorentz transformation for kaon-momentum 
 * to C.M.-momentum for the radiative kaon capture process.
 */
int lorentztrans_kaon(double pklab, double* pk, double mp, double mk)
{
  /* This can be done faster! */
  /*
    double s;*/
 /*Mandelstam (lorentzinvariant)*/

  /*
  s = pow( ( sqrt( pklab*pklab + mk*mk ) + mp ) , 2 ) - pklab*pklab;
  
  *pk = sqrt( (s*s + mp*mp*mp*mp + mk*mk*mk*mk 
	       - 2*mk*mk*s - 2*mp*mp*s - 2*mk*mk*mp*mp) / (4*s) );
  */

  /* This is faster code! */
  double EKlab = sqrt(pklab*pklab+mk*mk);
  *pk = mp*pklab/sqrt(mk*mk+mp*mp+2*mp*EKlab);
  return 0;

}

/*!
 * The absolute value |pk| of the momentum of the two escaping particles
 * is determined through the energy conservation (mass) relation.
 */
double construct_pk(double w, double k, double mp, double mk, double my)
{
  /* Check if a solution for pk exists!If massrelation(pk=0) > 0, this would
   * mean that the mass relation, which is a rising function of |pk|^2,
   * still gives a value larger than zero for the smallest physical
   * value of |pk|^2 (nl 0), and will only become 0 at negative |pk|^2
   * values -> a physical |pk|^2 cannot be found!!! */
  if(massrelation(w, k, 0., mp, mk, my) > 0.)
    return -1.;
  
  // the invariant mass
  double sqrts = sqrt(mp*mp + k*k) + w;

  // CM energy of the kaon
  double ek = (sqrts + (mk*mk-my*my)/sqrts)/2.;
  
  // Check for close to threshold
  if( (ek-mk)<STRANGEUFLOW ) return 0.;

  return sqrt(ek*ek-mk*mk);
}

/*!
 * Mass relation for photo- and electroproduction process in the c.m. system
 * (RHS - LHS of expression (74) in Stijn's notes.
 * For photoproduction w == k.
 */
double massrelation(double w, double k, double pk, double mp,
			   double mk, double my)
{
  double temp;

  
  temp = sqrt(mk*mk + pk*pk) + sqrt(my*my + pk*pk)
    - sqrt(mp*mp + k*k) - w;

  
  return temp;
}

/*!
 * The momentum of the two escaping particles in the case of kaon
 * capture is calculated
 */
double construct_w(double pklab, double mp, double mk, double my)
{
  /* This can be done faster */
  /*  double s; */
  /*Mandelstam (lorentzinvariant)*/
  /*double w;
  
  s = pow( ( sqrt( pklab*pklab + mk*mk ) + mp ) , 2 ) - pklab*pklab;

  w = (s - my*my)/(2*sqrt(s));

  return w;

  */

  /* This is faster code */
  double sqrt_s = sqrt(mp*mp+mk*mk+2.0*mp*sqrt(mk*mk+pklab*pklab));
  return 0.5*(sqrt_s - my*my/sqrt_s);

}

/*!
 * Conversion: dsigma / dt = domega_to_dt_conversion dsigma / dOmega.
 *
 * Jacobian | dOmega / dt | = 2*pi / (2*pk*k).
 *
 * Factor 1e-6 since we need dsigma / dt in (mubarn * GeV^-2), 
 * while pk,k are in MeV
 */
double domega_to_dt_conversion(double pk, double k)
{
  return PI / ( pk * k * 1e-6 );
}

/*!
 * Conversion:  dsigma / du = domega_to_du_conversion . dsigma / dOmega.
 *
 * Jacobian | dOmega / du | = - 2*pi / (2*pk*k). 
 *
 * Remark that we omit the sign, since the cross section 
 * has to be positive !!
 *
 * Factor 1e-6 since we need dsigma / du in (mubarn * GeV^-2), 
 * while pk,k are in MeV  
 */
double domega_to_du_conversion(double pk, double k)
{
   return PI / ( pk * k * 1e-6 );
}

/*! Error exit after printing message */
void error_exit(const std::string& message)
{
  throw std::logic_error("\nERROR: " + message + "\n===> BYE!\n#########\n\n\n");
}

