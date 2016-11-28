/*!
 * \file FormFactorSpecification.cpp
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

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <typeinfo>
#include <FormFactor.h>
#include "FormFactorParametrization.h"
#include "FormFactorSpecification.h"

using namespace std;



/*! ff_specification contains information about a formfactor parametrization 
 * and it's parameters extracted from 'input/numinput/formfactors'.
 *
 * \return FormFactor based on this information
 * and the type of coupling constant (cc_type) with which it will be used.
 */
FormFactor* specify_ff(string ff_specification, CC cc_type, int gic)
{
  /* Extract formfactor parametrization code and parameters from ff_specification
   * **************************************************************************** */
  
  /* Array of doubles (size=15) to store all form factor parametrization parameters.
   * Initialized to zero */
  double parameters[15];
  for(int i=0; i<15; i++)
    parameters[i] = 0.0;
  
  // ff_specification is stored in a C-string ff_spec for tokenization
  char ff_spec[ff_specification.size()+1];
  strcpy( ff_spec , ff_specification.c_str() );
  
  /* Tokenize ff_spec:
   * - the first and second token are the name and J^P of the particle
   * - the third token is the parametrization code and is stored in parametrization_type
   * - all others are parameters and are stored in the parameters[] array.
   *   The number of parameters given in ff_specification is usually less than 15.
   *
   * When only two tokens are present, the user didn't provide sufficient formfactor information.
   */
  char* token; /* C-string to temporarely store parts of ff_spec before conversion
		* to double or int */
  
  // Dump first and second token
  token = strtok (ff_spec," \t"); // first tokenization
  token = strtok(NULL," \t");
  
  // Use third token to check whether user gave enough info
  token = strtok(NULL," \t");
  if( token == NULL ) {
    cerr << "Form factor information in 'input/numinput/formfactors' is incomplete."
	 << endl;
    exit(1);
  }
  int parametrization_type = atoi(token); // parametrization code
  
  int parameter_nr = 0; // keep track of the number of imported parameters
  
  token = strtok (NULL," \t"); // next tokenization  
  while ( token != NULL)
    {
      parameters[parameter_nr++] = atof(token); // store the parameter
      token = strtok (NULL," \t"); // next tokenization
    }
  
  
  /* Choose the correct form factor parameterization
   * *********************************************** */
  ffCalculator ff_parametrization;
  /* The parametrizations are stored in a 2dim. array defined in 
   * 'FormFactorParametrization.h' */
    if(gic==0) 
      ff_parametrization = ff_parametrizations[parametrization_type-1][cc_type];
    else 
      ff_parametrization = gic_ff_parametrizations[parametrization_type-1][cc_type];

  // Allocate some dynamic memory and create the correct FormFactor object
  FormFactor* ptr_to_formfactor = new FormFactor(ff_parametrization,
						 parameters[0],
						 parameters[1],
						 parameters[2],
						 parameters[3],
						 parameters[4],
						 parameters[5],
						 parameters[6],
						 parameters[7],
						 parameters[8],
						 parameters[9],
						 parameters[10],
						 parameters[11],
						 parameters[12],
						 parameters[13],
						 parameters[14]);
  
  // Notify when something went wrong
  if (ptr_to_formfactor==NULL)
    {
      cout << "Error in allocation of formfactor!" << endl << endl;
      exit(1);
    }

  // return the pointer to the correct FormFactor
  return ptr_to_formfactor;
}

/* ******************************************************************************** */



/*! The EM formfactors are specified in the 'formfactors' file in the 'input/numinput'
 * folder.
 * Using that info, all exchanged particles in the particles[] array are given the
 * correct EM formfactor.
 */
int em_formfactor_specification(Observable* observ, Class particles[])
{
  
  if (observ->iso.isospin==11 || observ->iso.isospin==12) // kaon capture: no FF needed!
  {
    // we will set all form factors to the null pointer.
    // loop over all types of diagrams
    for(int diagram=0; diagram<CLASSMAX; diagram++)
    {
      // loop over all exchanged particles for a specific diagram
      for(int particle=0; particle<particles[diagram].particount; particle++)
      {
	particles[diagram].partic[particle].formfactorE = NULL;
	particles[diagram].partic[particle].formfactorG = NULL;
	particles[diagram].partic[particle].formfactorH = NULL;
	particles[diagram].partic[particle].formfactorI = NULL;
      }
    }
	      
    
    return 0;
  }
  
  cout << endl << "Storing EM formfactors in particles[] array..." << endl << endl;	
  // In case of electroproduction, the EM formfactors need to be specified
  if(observ->electroprod)
    {

      /* Open 'formfactors' file
       * *********************** */
      
      // Contruct name of file containing the form factor specifications
      char ff_input_file[MAX_LOCATION_STRING];
      strcpy(ff_input_file, observ->inFolder);
      strcat(ff_input_file, "numinput/formfactors.");
      if( (observ->iso.isospin>=1 && observ->iso.isospin<=3) // kaon production from the proton
	 || (observ->iso.isospin==7) || (observ->iso.isospin==8) ) // pion production from the proton
 	strcat(ff_input_file, "proton");
      else if( (observ->iso.isospin>=4 && observ->iso.isospin<=6)// kaon production from the neutron
	 || (observ->iso.isospin==9) || (observ->iso.isospin==10) ) // pion production from the neutron
	strcat(ff_input_file, "neutron");
      else
      {
	cerr << "ERROR in em_formfactor_specification(..): "
	     << "Isospin channel " << observ->iso.isospin
	     << " unknown!\n";
	exit(1);
      }
	
      
      // Open the file
      ifstream ff_input( ff_input_file );
      
      // Exit if error opening file
      if( !ff_input )
	{
	  cout << "Error opening formfactor specification file ("
	       << ff_input_file << ")" << endl << endl;
	  exit(1);
	}
      
      string ff_specification; /* this C++ string temporarely holds info on the type
				* of formfactor parametrization to use and it's
				* parameters */
      
      // Skip first 3 lines
      for(int count=0; count<3; count++)
	getline( ff_input , ff_specification );
      
      // Store in begin_position where the info about formfactors begins
      int begin_position = ff_input.tellg();
      
      
      
      /* Loop over all diagrams and exchanged particles
       * ********************************************** */
      
      // loop over all types of diagrams
      for(int diagram=0; diagram<CLASSMAX; diagram++)
	{
	  // loop over all exchanged particles for a specific diagram
	  for(int particle=0; particle<particles[diagram].particount; particle++)
	    {
	      
	      /* Look for formfactor parametrization for specific particle
	       * ********************************************************* */
	      
	      // position the stream to the beginning of the input file
	      ff_input.seekg( begin_position );
	      
	      // Run through the file untill we find the relevant info
	      ff_input >> ff_specification;
	      while( ff_specification != particles[diagram].partic[particle].nickname &&
		     ff_specification != "X" ) 
		{
		  getline( ff_input , ff_specification );
		  ff_input >> ff_specification;
		}
	      
	      // When EOF was reached before info was found -> exit
	      if ( ff_specification == "X" )
		{
		  cout << "No info about EM form factor for "
		       << particles[diagram].partic[particle].nickname 
		       << " was found!" << endl << endl;
		  exit(1);
		}
	      else
		{
		  // Store remaining line in ff_specifaction
		  getline( ff_input , ff_specification );
		}
	      
	      
	      /* for each specific type of diagram, certain coupling constants
	       * get form factors as specified in 'input/numinput/formfactors' 
	       * ************************************************************* */
	      
	      int gic = particles[diagram].partic[particle].gic;

	      switch(diagram)
		{
		  // S-channel born diagram
		case 0:
		  particles[diagram].partic[particle].formfactorE = 
		    specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;
		  
		  // T-channel born diagram
		case 1:
		  particles[diagram].partic[particle].formfactorE = 
		    specify_ff(ff_specification,E,gic);
		  break;
		  
		  // U-channel born diagram
		case 2:
		  particles[diagram].partic[particle].formfactorE = 
		    specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // A-channel extended born diagram
		case 3:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;
		  
		  // B-diagram
		case 4:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // C-diagram
		case 5:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // D-diagram
		case 6:
		  if(particles[diagram].partic[particle].long_coupling)
		    particles[diagram].partic[particle].formfactorE = 
		      specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // E-diagram
		case 7:
		  if(particles[diagram].partic[particle].long_coupling)
		    particles[diagram].partic[particle].formfactorE = 
		      specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // F-diagram
		case 8:
		  if(particles[diagram].partic[particle].long_coupling)
		    particles[diagram].partic[particle].formfactorE = 
		      specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // G-diagram
		case 9:
		  if(particles[diagram].partic[particle].long_coupling)
		    particles[diagram].partic[particle].formfactorE = 
		      specify_ff(ff_specification,E,gic);
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // H-diagram
		case 10:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // I-diagram
		case 11:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // J-diagram
		case 12:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // l-diagram
		case 13:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // M-diagram
		case 14:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // N-diagram
		case 15:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // O-diagram
		case 16:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // Q-diagram
		case 17:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  particles[diagram].partic[particle].formfactorH = 
		    specify_ff(ff_specification,H,gic);
		  break;

		  // R-diagram
		case 18:
		  particles[diagram].partic[particle].formfactorE = 
		    specify_ff(ff_specification,E,gic);
		  break;

		  // V-diagram
		case 19:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		  // W-diagram
		case 20:
		  particles[diagram].partic[particle].formfactorG = 
		    specify_ff(ff_specification,G,gic);
		  break;

		default:
		  cout << "EM form factors not yet implemented for diagram nr."
		       << diagram << "." << endl << endl;
		  exit(1);
		  
		} // end of switch(diagram)
	      
	    } // end particle loop
	} // end diagram loop

} // end of electroproduction case

  // No EM formfactors needed in photoproduction
  else
    {
      cout << "No EM form factors needed in photoproduction" << endl;
      exit(1);
    }

  return 0;
}

/* ************************************************************************************** */



/*! The user has specified his choice of strong form factors
 * in io_specific and those were stored in observ.
 * In this function those specifications are implemented in the
 * particles array
 */
int strong_formfactor_specification(Observable* observ, Class particles[])
{
  cout << endl << "Storing hadronic formfactors in particles[] array..." << endl << endl;


  /* When doing pure Regge exchange without added resonances
   * strong formfactors are not needed */
  if(observ->regge && !observ->reg.res_on_bg)
    {
      cout << "In pure Regge-exchange no hadronic form factors needed!"
	   << endl << endl;
      exit(1);
    }


  /* Determine type of hadronic formfactors
   * ************************************** */

  // strongFF is obviously the strong formfactor parametrization
  ffCalculator strongFF;

  // Make some changes when the lorentz modeltype was chosen
  if(!strcmp(observ->modelType,"lorentz"))
    strongFF = strong_lorentz;
  else if (!strcmp(observ->modelType,"lorentz_gaussian"))
    strongFF = strong_lorentz_gaussian;
    // Make sure the strong formfactor is set to Lorentz in case of the Lorentz modelType 
  else
    strongFF = strong_gaussian;


  /* Loop over all diagrams and exchanged particles
   * ********************************************** */

  /* It is not always necessary to run through all types of diagrams.
   * In the RPR case, the (extended) born terms do not need
   * strong formfactors, so we begin the loop at diagram nr.6.
   */
  int startDiagram = 6;


  /* !! ATTENTION !!
   * In the isobar case we provide hadronic formfactors for both
   * resonant as background terms. One should beware however that
   * strong formfactors break gauge invariance in the born terms!
   * Several gauge restauration procedures exist (see Stijn's thesis).
   * None of them are implemented in this code, so the results are 
   * NOT correct.
   * For isobar calculations use the old code.
   */
  if(!observ->regge)
    {
      startDiagram = 0;
      cout << "Using hadronic form factors in the isobar model. "
	   << "No gauge restauration procedure implemented, "
	   << "so results are probably incorrect!" << endl
	   << " Use Stijn's old code for isobar model calculations."
	   << endl << endl;
    }
  
  
  // loop over all types of diagrams
  for(int diagram=startDiagram; diagram<CLASSMAX; diagram++)
    {

      // loop over all exchanged particles for a specific diagram
      for(int particle=0; particle<particles[diagram].particount; particle++)
	{
	  
	  /* for each specific type of diagram, certain coupling constants
	   * get form factors.
	   * The cutoff mass is stored in the observ structure.
	   * ************************************************************* */
	  
	  switch(diagram)
	    {
	      // S-channel born diagram
	    case 0:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	      // T-channel born diagram
	    case 1:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	      // U-channel born diagram
	    case 2:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;

	      // A-diagram
	    case 3:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	      // B-diagram
	    case 4:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;

	      // C-diagram
	    case 5:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;

	      // D-diagram
	    case 6:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // E-diagram
	    case 7:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // F-diagram
	    case 8:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // G-diagram
	    case 9:
	      particles[diagram].partic[particle].formfactorH =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // H-diagram
	    case 10:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // I-diagram
	    case 11:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // J-diagram
	    case 12:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // L-diagram
	    case 13:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // M-diagram
	    case 14:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // N-diagram
	    case 15:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // O-diagram
	    case 16:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // Q-diagram
	    case 17:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.res_cutoff);
	      break;

	      // R-diagram
	      case 18:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	      // V-diagram
	      case 19:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	      // W-diagram
	      case 20:
	      particles[diagram].partic[particle].formfactorI =
		new FormFactor(strongFF,
			       particles[diagram].partic[particle].mass,
			       observ->ffac.born_cutoff);
	      break;
	      
	    default:
	      cout << "Hadronic form factors not yet implemented for diagram nr."
		   << diagram << "." << endl << endl;
	      exit(1);
	      
	    } // end of switch(diagram)
	  
	} // end particle loop
    } // end diagram loop
  
  return 0;
}


/* ************************************************************************************** */


/*! Store the cutoff masses for the strong form factors in all relevant formfactors
 * in the particles[] array.
 *
 * This function looks similar to strong_formfactor_specification(). See that function
 * for more clarifying comments
 */
int update_strong_formfactor(Observable* observ, Class particles[])
{
  

  /* Loop over all diagrams and exchanged particles
   * ********************************************** */

  /* In the RPR case, the (extended) born terms do not have
   * strong formfactors, so we begin the loop at diagram nr.6.
   */
  int startDiagram = 6;


  /* !! ATTENTION !!
   * In the isobar case we provide hadronic formfactors for both
   * resonant as background terms.
   */
  if(!observ->regge)
    startDiagram = 0;

  
  
  // loop over all types of diagrams
  for(int diagram=startDiagram; diagram<CLASSMAX; diagram++)
    {

      // loop over all exchanged particles for a specific diagram
      for(int particle=0; particle<particles[diagram].particount; particle++)
	{
	  
	  /* We update the cutoff mass in the relevant FormFactor objects
	   * The cutoff mass is stored in the observ structure.
	   * ************************************************************* */
	  
	  switch(diagram)
	    {
	      // S-channel born diagram
	    case 0:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;
	      
	      // T-channel born diagram
	    case 1:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;
	      
	      // U-channel born diagram
	    case 2:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;

	      // A-diagram
	    case 3:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;
	      
	      // B-diagram
	    case 4:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;

	      // C-diagram
	    case 5:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;

	      // D-diagram
	    case 6:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // E-diagram
	    case 7:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // F-diagram
	    case 8:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // G-diagram
	    case 9:
	      particles[diagram].partic[particle].formfactorH->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // H-diagram
	    case 10:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // I-diagram
	    case 11:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // J-diagram
	    case 12:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // L-diagram
	    case 13:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // M-diagram
	    case 14:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // N-diagram
	    case 15:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // O-diagram
	    case 16:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // Q-diagram
	    case 17:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.res_cutoff);
	      break;

	      // R-channel
	    case 18:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;
	      // V-channel
	    case 19:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;
	      // W-channel
	    case 20:
	      particles[diagram].partic[particle].formfactorI->
		setParameters(particles[diagram].partic[particle].mass,
			      observ->ffac.born_cutoff);
	      break;	      
	      
	    default:
	      cout << "Hadronic form factors not yet implemented for diagram nr."
		   << diagram << "." << endl << endl;
	      exit(1);
	      
	    } // end of switch(diagram)
	  
	} // end particle loop
    } // end diagram loop

  return 0;
}

/* ************************************************************************************** */

/*! In specify_ff() memory is dynamically allocated for FormFactor objects.
 * When the program terminates this memory should be freed
 */
int release_formfactors(Class particles[])
{
  // loop over all types of diagrams
  for(int diagram=0; diagram<CLASSMAX; diagram++)
    {
      // loop over all exchanged particles for a specific diagram
      for(int particle=0; particle<particles[diagram].particount; particle++)
      {
	    delete particles[diagram].partic[particle].formfactorE;
	    particles[diagram].partic[particle].formfactorE = NULL;
	    delete particles[diagram].partic[particle].formfactorG;
	    particles[diagram].partic[particle].formfactorG = NULL;
	    delete particles[diagram].partic[particle].formfactorH;
	    particles[diagram].partic[particle].formfactorH = NULL;
	    delete particles[diagram].partic[particle].formfactorI;
	    particles[diagram].partic[particle].formfactorI = NULL;
	}
    }

  return 0;
}
