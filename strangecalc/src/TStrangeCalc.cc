/*!
 * \file TStrangeCalc.cc
 * \ingroup wrapper
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Dave Ireland <d.ireland@physics.gla.ac.uk>
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

#define TStrangeCalc_cxx

#include "TStrangeCalc.h"
#include <FormFactorSpecification.h>
#include <TMatrixElement.h>
#include <cstring>
#include <cfloat>
#include <cassert>
#include <vector>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::make_pair;

/*!
 * \class TStrangeCalc
 * Class to wrap around the functions of strangecalc.
 * 
 * \date started Thu May 10 16:07:43 BST 2001
 * 
 * \author Dave Ireland <d.ireland@physics.gla.ac.uk>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 */

/*! \brief Constructor for the strangecalc wrapper class TStrangeCalc
 * \param log file pointer for output
 * \param fit_specification char array with the name of the fit_specification file used to initialise an object
 * \param datafolder optionally an alternative datafolder can be specified. In this case the one in fit_specification is ignored.
 *
 */
TStrangeCalc::TStrangeCalc ( FILE *log, const char *fit_specification, const char* datafolder ) :
  datapoints ( new ( Data* [ISOMAX] ) ),
  particles ( new ( Class* [ISOMAX] ) ),
  printparticles ( new Class[CLASSMAX] ),
  varinfo ( new Varinfo[CLASSMAX+1] ) /* 1 additional element for formfac*/
{
  // Define well-initialized structs only once
  static Observable emptyobserv = {0};
  static Data emptydata = {0};
  static Class emptyclass = {0};
  static Varinfo emptyvarinfo = {{{{0}}}};
  
  // initialize observ to all 'zeros'
  observ = emptyobserv;
  observ.reg.t_and_u_channel = 0;

  // Weight for fitting should be initialized to 1 by default
  observ.fit.polweight = 1;
  for(int i=0; i<ISOMAX; ++i)
    observ.fit.isoWeight[i] = 1;


  // Constructor
  // Also performs initialisation for strangeness calculations.
  // Dynamically allocate large arrays.
  int i, j, k;

  for ( i=0; i<CLASSMAX; i++ )
    {
      // Initialize printparticles[] to zero
      printparticles[i] = emptyclass;
      // Initialize varinfo[] to zero
      varinfo[i]=emptyvarinfo;
    }
  varinfo[CLASSMAX]=emptyvarinfo;

  for ( i=0; i<ISOMAX; i++ )
    {
      datacount[i]=0;
      datapoints[i] = new Data[DATAMAX];
      for ( j=0; j<DATAMAX; j++ )
	{
	  // Initialize datapoints[][] to zero
	  datapoints[i][j] = emptydata;
	  datapoints[i][j].photo.label = -1;
	  datapoints[i][j].elec.label = -1;
	  datapoints[i][j].kaoncap.label = -1;
	}
      particles[i] = new Class[CLASSMAX];
      for ( j=0; j<CLASSMAX; j++ )
	{
	  // Initialize particles[][] to zero
	  particles[i][j] = emptyclass;
	  for ( k=0; k<PARTICLEMAX; k++ )
	    {
	      strcpy ( particles[i][j].partic[k].nickname,"" );
	      particles[i][j].partic[k].mass = 0.0;
	      particles[i][j].partic[k].width = 0.0;
	      particles[i][j].partic[k].E = 0.0;
	      particles[i][j].partic[k].G = 0.0;
	      particles[i][j].partic[k].H = 0.0;
	      particles[i][j].partic[k].I = 0.0;
	      particles[i][j].partic[k].J = 0.0;
	      particles[i][j].partic[k].X = 0.0;
	      particles[i][j].partic[k].Y = 0.0;
	      particles[i][j].partic[k].Z = 0.0;
	      particles[i][j].partic[k].r_kappa_n_p = 0.0;
	      particles[i][j].partic[k].r_kappa_1_n_p = 0.0;
	      particles[i][j].partic[k].r_kappa_2_n_p = 0.0;
	      particles[i][j].partic[k].regge_phase = 0;
	      particles[i][j].partic[k].regge_phase_nonbase = 0;
	      particles[i][j].partic[k].long_coupling = 0;
	      particles[i][j].partic[k].formfactorE = NULL;
	      particles[i][j].partic[k].formfactorG = NULL;
	      particles[i][j].partic[k].formfactorH = NULL;
	      particles[i][j].partic[k].formfactorI = NULL;

	      if(j == 1 || j == 18)
		particles[i][j].partic[k].spin = 0.0;
	      else if(j == 4 || j == 5 || j == 19 || j == 20)
		particles[i][j].partic[k].spin = 1.0;
	      else if(j == 0 || j == 2 || j == 3 || (j >= 6 && j < 10))
		particles[i][j].partic[k].spin = 1.0/2.0;
	      else if(j >= 10 && j < 14)
		particles[i][j].partic[k].spin = 3.0/2.0;
	      else if(j >= 14 && j < 18)
		particles[i][j].partic[k].spin = 5.0/2.0;
	      else error_exit("Spin for particles exceeding CLASSMAX is not defined!");
	    }
	}
    }

  // initialize fVertex to zero
  for ( i=0; i<MAXNDIM; ++i )
    fVertex[i] = 0.0;

  //  The fitting specifications are read in from the
  //  <fit_specifications> file.
  //
  //  When the fit specifications are read in, the setup is done
  //  of all the structures and variables needed by the calculations.

  for ( i=0 ; i< ( CLASSMAX+1 ) ; i++ )
    {
      for ( j=0 ; j<PARTICLEMAX ; j++ )
	{
	  for ( k=0 ; k<MAXNRFVAR ; k++ )
	    {
	      varinfo[i].partic[j][k].var   = short ( 0 );
	      varinfo[i].partic[j][k].bound = short ( 0 );
	      varinfo[i].partic[j][k].low   = double ( 0 );
	      varinfo[i].partic[j][k].up    = double ( 0 );
	    }
	}
    }

  pho_diff_setup = pho_rec_setup = pho_tar_setup = pho_pho_setup = 0;
  elec_diff_setup = 0;
  isospin = 0;

  // Get the fit specifications
  /* Open configuation file */
  FILE *ifp, *ofp2;
  ifp = fopen ( fit_specification, "r" );
  if ( ifp == NULL )
    error_exit ( "TStrangeCalc::TStrangeCalc: Could not open file:'fit_specification'\n" );

  ofp2 = log;

  int istat = observablespecification ( &observ, ifp, log, ofp2 );
     
  fprintf(log, "\nTStrangeCalc::TStrangeCalc: \nGetting fit specifications...\t\t");
  fflush(log);
  TestStatus ( istat );

  // Check whether the user specified an alternative data folder
  if( datafolder ) strcpy( observ.dataFolder, datafolder );
  
  // Import and construct the particle properties and  data points
  fprintf(log, "Importing data...\t");
  fflush(log);
  istat = SetupStrangeStruc ( ifp );
  TestStatus ( istat );
  fclose ( ifp );

  // Determine which isospin channels are needed and what is the
  // "isospin_base".
  // If more than one isospin channel is involved, construct a 
  // start vertex based on the "iso_base" of the process.
  // If just one isospin channel is under investigation, the
  // choice is trivial.
  if(observ.iso.nr_iso_channels > 1) {
    isospin = observ.iso.isospin = observ.iso.iso_base;
  }
  else {
    isospin = observ.iso.isospin = observ.iso.iso_channel[0];
  }

  // Determine which variables are allowed to vary and which have to
  // be static. If a variable is allowed to vary, it can be a free
  // parameter or a bounded one. In the latter case, determine the
  // limits between it can vary. Note that in the GA, all parameters
  // are bounded
  fprintf(log,"Getting variable info...\t\t");
  fflush(log);
  istat = get_variable_info ( varinfo, particles[observ.iso.iso_base], &observ );
  TestStatus ( istat );


  //     Plot the "particles[]" and "varinfo[]" information in the
  //     log-file.
  PlotParameterInfo( log, particles[observ.iso.iso_base] );


  // Make a start vertex determined by the "particles" structure. "fDim" is
  // the dimension of the vertex.
  fprintf(log, "Making a start vertex...\t\t");
  fflush(log);
  istat = MakeStartvertex();
  TestStatus ( istat );

  // Give a label to each parameter
  AssignVertexName();

  // The limits for each parameter need to be output to the log
  // file...
  fprintf ( log, "\n<limits>\n" );
  for ( int i=0 ; i<fDim ; i++ )
    {
      fprintf ( log, "%d\t%f\t%f\n", ( i+1 ), fLimits[i].low, fLimits[i].up );
    }
  fprintf ( log, "</limits>\n\n" );

  fprintf ( log, "##########################################\n\n" );

  // A bit more pseudo-XML to log...
  fprintf ( log, "<vertexdim>\n" );
  fprintf ( log, "\t%d\n", fDim );
  fprintf ( log, "</vertexdim>\n" );


  // A first calculation of the chi squared based on the coupling
  // constants in the "./input/numinput/coupl.iso.*" file.
  fprintf(log, "Initial chi-squared evaluation...\t");
  fflush(log);
  SetChiSquared();
  fprintf(log,"done\n");
  fflush(log);
  copyparticles ( printparticles,particles[observ.iso.iso_base],&observ );
   
  // set the dataset weight
  fDataWeight = setlabels(datapoints, datacount, &observ);
   
}


/*! \brief Copy constructor
 * \param toCopy TStrangeCalc object used to initialise this instance
 */
TStrangeCalc::TStrangeCalc ( const TStrangeCalc& toCopy ) :
  fDim ( toCopy.fDim ),
  fChiSquared ( toCopy.fChiSquared ),
  observ ( toCopy.observ ),
  datapoints ( new ( Data* [ISOMAX] ) ),
  particles ( new ( Class* [ISOMAX] ) ),
  printparticles ( new Class[CLASSMAX] ),
  varinfo ( new Varinfo[CLASSMAX+1] ),    /* 1 additional element for formfac.*/
  pho_diff_setup ( toCopy.pho_diff_setup ),
  pho_rec_setup ( toCopy.pho_rec_setup ),
  pho_tar_setup ( toCopy.pho_tar_setup ),
  pho_pho_setup ( toCopy.pho_pho_setup ),
  elec_diff_setup ( toCopy.elec_diff_setup ),
  isospin ( toCopy.isospin )
{
  // The explicit copy constructor

  for ( int i=0; i<MAXNDIM; ++i )
    {
      fVertex[i] = toCopy.fVertex[i];
      fLimits[i] = toCopy.fLimits[i];
      for ( int j=0; j<kNameLength; ++j )
	fName[i][j] = toCopy.fName[i][j];
    }

  for ( int i=0; i<ISOMAX; ++i )
    {
      datacount[i] = toCopy.datacount[i];
      datapoints[i] = new Data[DATAMAX];
      particles[i] = new Class[CLASSMAX];
      for ( int j=0; j<DATAMAX; ++j )
	// copy all the Data structs
	datapoints[i][j] = toCopy.datapoints[i][j];
      for ( int j=0; j<CLASSMAX; ++j )
	{
	  // copy all the Class structs
	  // The following code is almost identical to copyparticles.
	  // Important difference: pointers to Formfactors are dereferenced!
	  particles[i][j]= toCopy.particles[i][j];

	  for ( int k=0; k< PARTICLEMAX; k++ )
            // the Properties struct
            //(in array called partic) contains some pointers to formfactors.
	    {
	      // formfactor info (different from copyparticles(...): dereferenced!)
	      if ( toCopy.particles[i][j].partic[k].formfactorE ==NULL )
		particles[i][j].partic[k].formfactorE = NULL;
	      else  particles[i][j].partic[k].formfactorE = new
		      FormFactor ( * ( toCopy.particles[i][j].partic[k].formfactorE ) );
	      if ( toCopy.particles[i][j].partic[k].formfactorG ==NULL )
		particles[i][j].partic[k].formfactorG = NULL;
	      else  particles[i][j].partic[k].formfactorG = new
		      FormFactor ( * ( toCopy.particles[i][j].partic[k].formfactorG ) );
	      if ( toCopy.particles[i][j].partic[k].formfactorH ==NULL )
		particles[i][j].partic[k].formfactorH = NULL;
	      else  particles[i][j].partic[k].formfactorH = new
		      FormFactor ( * ( toCopy.particles[i][j].partic[k].formfactorH ) );
	      if ( toCopy.particles[i][j].partic[k].formfactorI ==NULL )
		particles[i][j].partic[k].formfactorI = NULL;
	      else  particles[i][j].partic[k].formfactorI = new
		      FormFactor ( * ( toCopy.particles[i][j].partic[k].formfactorI ) );

	    }
	}
      // the above completes particles[i][j] = toCopy.particles[i][j];
    }

  copyparticles ( printparticles,particles[observ.iso.iso_base],&observ );
   
  for ( int i=0; i<CLASSMAX; ++i )
    {
      varinfo[i] = toCopy.varinfo[i];
    }
  varinfo[CLASSMAX] = toCopy.varinfo[CLASSMAX];
}


/*! \brief Assignment operator
 * \param \toCopy TStrangeCalc object on the right hand side
 */
TStrangeCalc& TStrangeCalc::operator= ( const TStrangeCalc& toCopy )
{
  // The assignment operator

  if ( this != &toCopy )
    {
      fDim = toCopy.fDim;
      fChiSquared =toCopy.fChiSquared;
      observ =toCopy.observ;
      pho_diff_setup =toCopy.pho_diff_setup;
      pho_rec_setup =toCopy.pho_rec_setup;
      pho_tar_setup =toCopy.pho_tar_setup;
      pho_pho_setup =toCopy.pho_pho_setup;
      elec_diff_setup =toCopy.elec_diff_setup;
      isospin = toCopy.isospin;

      // Note that arrays (if allocated on the stack) are copied
      // automagically when you copy a struct that contains them. They do not
      // however like to be copied explicitly, that's why we need this loop:
      for ( int i=0; i<MAXNDIM; ++i )
	{
	  fVertex[i] = toCopy.fVertex[i];
	  fLimits[i] = toCopy.fLimits[i];
	  for ( int j=0; j<kNameLength; ++j )
            fName[i][j] = toCopy.fName[i][j];
	}

      for ( int i=0; i<ISOMAX; ++i )
	{
	  datacount[i] = toCopy.datacount[i];
	  for ( int j=0; j<DATAMAX; ++j )
            // copy all the Data structs
            datapoints[i][j] = toCopy.datapoints[i][j];
	  for ( int j=0; j<ISOMAX; ++j )
	    {
	      // copy all the Class structs
	      // The following code is almost identical to copyparticles.
	      // Important difference: pointers to Formfactors are dereferenced!
	      particles[i][j]= toCopy.particles[i][j];

	      for ( int k=0; k< PARTICLEMAX; k++ )
		// the Properties struct (in array called partic)
		// contains some pointers to formfactors.
		{
		  // formfactor info (different from copyparticles(...): dereferenced!)
		  if ( toCopy.particles[i][j].partic[k].formfactorE ==NULL )
		    {
		      delete particles[i][j].partic[k].formfactorE;
		      particles[i][j].partic[k].formfactorE = NULL;
		    }
		  else * ( particles[i][j].partic[k].formfactorE ) =
			 * ( toCopy.particles[i][j].partic[k].formfactorE );

		  if ( toCopy.particles[i][j].partic[k].formfactorG ==NULL )
		    {
		      delete particles[i][j].partic[k].formfactorG;
		      particles[i][j].partic[k].formfactorG = NULL;
		    }
		  else * ( particles[i][j].partic[k].formfactorG ) =
			 * ( toCopy.particles[i][j].partic[k].formfactorG );

		  if ( toCopy.particles[i][j].partic[k].formfactorH ==NULL )
		    {
		      delete particles[i][j].partic[k].formfactorH;
		      particles[i][j].partic[k].formfactorH = NULL;
		    }
		  else * ( particles[i][j].partic[k].formfactorH ) =
			 * ( toCopy.particles[i][j].partic[k].formfactorH );

		  if ( toCopy.particles[i][j].partic[k].formfactorI ==NULL )
		    {
		      delete particles[i][j].partic[k].formfactorI;
		      particles[i][j].partic[k].formfactorI = NULL;
		    }
		  else * ( particles[i][j].partic[k].formfactorI ) =
			 * ( toCopy.particles[i][j].partic[k].formfactorI );
		}
	    }
	  // the above completes particles[i][j] = toCopy.particles[i][j];
	}
      for ( int i=0; i<CLASSMAX; ++i )
	{
	  printparticles[i] = toCopy.printparticles[i];
	  for ( int k=0; k< PARTICLEMAX; k++ )
            // the Properties struct (in array called partic)
            // contains some pointers to formfactors.
	    {
	      * ( printparticles[i].partic[k].formfactorE ) =
		* ( toCopy.printparticles[i].partic[k].formfactorE );
	      * ( printparticles[i].partic[k].formfactorG ) =
		* ( toCopy.printparticles[i].partic[k].formfactorG );
	      * ( printparticles[i].partic[k].formfactorH ) =
		* ( toCopy.printparticles[i].partic[k].formfactorH );
	      * ( printparticles[i].partic[k].formfactorI ) =
		* ( toCopy.printparticles[i].partic[k].formfactorI );
	    }
	  varinfo[i] = toCopy.varinfo[i];
	}
      varinfo[CLASSMAX] = toCopy.varinfo[CLASSMAX];
    }
  return *this;
}

/*! \brief Prints the current values of the variable coupling constants:
 * index and value is printed to stdout.
 */
void TStrangeCalc::PrintVertex() const
{
  cout << "\nPrintVertex()\n";
  for ( int i=0 ; i<fDim ; i++ )
    cout << i << "\t" << fVertex[i] << endl;
}

/*! \brief Prints the limits for the fitting process:
 * prints upper and lower bounds to stdout.
 */
void TStrangeCalc::PrintLimits() const
{
  cout << "\nPrintLimits()\n";
  for ( int i=0 ; i<fDim ; i++ )
    cout << i << "\t" << fLimits[i].low
	 << "\t" << fLimits[i].up << endl;
}

/*! \brief Prints particles information.
 * The particle classes, nicknames and coupling constants are 
 * printed in the coupl.iso.x format.
 * \param stream file pointer to write to
 * 
 */

void TStrangeCalc::PrintParticles ( FILE* stream )
{
  //  cout << "\nPrintParticles():\n";
  insertvertex ( printparticles, particles[observ.iso.iso_base],
		 varinfo, &observ, fVertex );
  fprintparticspec ( stream, printparticles, observ );
}

//----------------------------------------------------------------------

/*! \brief Changes the variable parameters with an array
 * \param array array of new parameters
 * 
 */
void TStrangeCalc::SetVertex ( const double array[] )
{
  for ( int i=0 ; i<fDim ; i++ )
    {
      fVertex[i] = array[i];
    }
}


/*! \brief Changes the variable parameters with a vector
 * \param ind vector of new parameters
 * 
 */
void TStrangeCalc::SetVertex ( const individual ind )
{
  for ( int i=0 ; i<fDim ; i++ )
    {
      fVertex[i] = ind[i];
    }
}


/*! \brief Calculates the internal chisquared 
 */
void TStrangeCalc::SetChiSquared()
{
  fChiSquared = chifunc ( fVertex, particles, varinfo,
			  &observ, datapoints, datacount );

}


//----------------------------------------------------------------------
/*! \brief Determines the weighted number of degrees of freedom of the fit.
 * When all observables are given the same weight, we have
 * trivially noDegreesOfFreedom = no. of datapoints 
 * This is not always the case however.
 * \return number of degrees of freedom (weights included)
 */
int TStrangeCalc::GetNDF() const
{

  int noDegreesOfFreedom = 0; // the final result
  int weight = 1;             // intermediate result

  // Loop over all relevant isospin channels
  for ( int isoChannels=0; isoChannels<observ.iso.nr_iso_channels; isoChannels++ )
    {
      int iso = observ.iso.iso_channel[isoChannels];  // set isospin
      noDegreesOfFreedom += datacount[iso]; // add no. of datapoints

      // Loop over datapoints
      for ( int data=0; data<datacount[iso]; data++ )
	{
	  weight = 1;

	  // weight procedure according to H.Thom, Phys.Rev. 151, 1322 (1966).
	  if ( datapoints[iso][data].photo_prod &&
               ( !strcmp ( datapoints[iso][data].photo.observable, "rec" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "pho" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "tar" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "c_x" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "c_xp") ||
                 !strcmp ( datapoints[iso][data].photo.observable, "c_z" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "c_zp") ||
                 !strcmp ( datapoints[iso][data].photo.observable, "o_x" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "o_xp") ||
                 !strcmp ( datapoints[iso][data].photo.observable, "o_z" ) ||
                 !strcmp ( datapoints[iso][data].photo.observable, "o_zp") )
	       )
	    {
	      if ( observ.fit.polweight != 0 )
		weight *= observ.fit.polweight;
	      else
		error_exit ( "error in polarization weight!!\n" );
	    }

	  if ( observ.iso.nr_iso_channels > 1 )
            weight *= observ.fit.isoWeight[datapoints[iso][data].iso];

	  noDegreesOfFreedom += weight -1;

	} // end loop datapoints
    } // end loop isospin channels

  return noDegreesOfFreedom;
}

//----------------------------------------------------------------------

double TStrangeCalc::GetCalcpoint(Data* datapoint,double w, double k,
				  double cos, double pk, double s ,
				  double phi, double phiMin)
{
  double calcpoint = -DBL_MAX;

  // check if object was initialized for requested isospin channel
  if ( getisospinbasechannel ( datapoint->iso ) != observ.iso.iso_base )
    {
      cerr << "ERROR in TStrangeCalc::GetCalcpoint(...): "
	   << "TStrangeCalc object has not been initialized for "
	   << "isospin channel no." << datapoint->iso << ".\n";
      return calcpoint;
    }

  // initialize particles structure
  observ.iso.isospin = datapoint->iso;
  insertvertex ( printparticles,particles[datapoint->iso], varinfo, &observ,
		 fVertex );

  // Transform coupling constants to correct isospin channel
  change_coupling_constants ( printparticles,&observ,observ.iso.isospin );

  // calculate datapoint
  if ( datapoint->electro_prod )
    // no need to check this as form factors are initialized anyway in
    //TStrangeModel!
    // maybe we should display warning in case GetCalcpoint is used outside
    //TStrangeViewer.
    // if(observ.fit.elec_diffcs[datapoint->iso])
    calcpoint =
      calculate_electro_observable ( &datapoint->elec,&observ,
				     printparticles,0,w,k,cos,pk,s,phi,phiMin);

  //     else
  //       cerr << "ERROR in TStrangeCalc::GetCalcpoint(...): "
  //     << "TStrangeCalc object has not been initialized for "
  //     << "electroproduction calculations.\n";

  else if ( datapoint->photo_prod )
    calcpoint =
      calculate_photo_observable ( &datapoint->photo,&observ,
				   printparticles,0,w,k,cos,pk );

  else if ( datapoint->kaoncapture )
    calcpoint =
      calculate_kaoncap_observable ( &datapoint->kaoncap,&observ,
				     printparticles,0,w,k,cos,pk );

  else
    {
      cerr << "ERROR in TStrangeCalc::GetCalcpoint(...): "
	   << "Only photo and electro production are implemented!\n";
      exit ( 1 );
    }

  return calcpoint;
}


//----------------------------------------------------------------------
/*! \brief Sets the matrix element singleton's kinematics and return it with a copy of observ and particles.
 * Also correctly set the iso.isospin property of the observ member, insert the fVertex 
 * and transform the cc's for this isospin channel.
 * \param iso isospin channel for which it will return a matrix element.
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \return Pointer to the TMatrixElement singleton. The class TMatrixElement owns the returned TMatrixElement*
 */
TMatrixElement* TStrangeCalc::GetMatrixElement ( int iso,double w, double k,
						 double cos, double pk )
{
  // check if object was initialized for requested isospin channel
  if ( getisospinbasechannel ( iso ) != observ.iso.iso_base )
    {
      cerr << "ERROR in TStrangeCalc::GetMatrixElement(...): "
	   << "TStrangeCalc object has not been initialized for "
	   << "isospin channel no." << iso << ".\n";
      exit ( 1 );
    }

  // Make duplicates of some data member
  double vertex[MAXNDIM];
  for(int i=0; i<MAXNDIM; ++i) vertex[i] = fVertex[i];

  // initialize particles structure
  observ.iso.isospin = iso;
  insertvertex ( printparticles,particles[iso], varinfo, &observ, vertex );

  // Transform coupling constants to correct isospin channel
  change_coupling_constants ( printparticles,&observ,observ.iso.isospin );

  return TMatrixElement::GetMatrixElement( w,k,cos,pk,printparticles,&observ);
}

//----------------------------------------------------------------------
/*! \brief Sets the matrix element singleton's kinematics and return it with a copy of observ and particles. The kinematics are allowed to have off-mass-shell particles.
 * Also correctly set the iso.isospin property of the observ member, insert the fVertex 
 * and transform the cc's for this isospin channel.
 * \param iso isospin channel
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param mn nucleon mass (in MeV)
 * \param mk kaon mass (in MeV)
 * \param my hyperon mass (in MeV)
 * \return Pointer to the TMatrixElement singleton. The class TMatrixElement owns the returned TMatrixElement*
 */
TMatrixElement* TStrangeCalc::GetMatrixElement(int iso,double w, double k,
					       double cos, double pk, double mn,
					       double mk, double my)
{
  // check if object was initialized for requested isospin channel
  if ( getisospinbasechannel ( iso ) != observ.iso.iso_base )
    {
      cerr << "ERROR in TStrangeCalc::GetMatrixElement(...): "
	   << "TStrangeCalc object has not been initialized for "
	   << "isospin channel no." << iso << ".\n";
      exit ( 1 );
    }

  // Make duplicates of some data member
  double vertex[MAXNDIM];
  for(int i=0; i<MAXNDIM; ++i) vertex[i] = fVertex[i];

  // initialize particles structure
  observ.iso.isospin = iso;
  insertvertex ( printparticles,particles[iso], varinfo, &observ, vertex );

  // Transform coupling constants to correct isospin channel
  change_coupling_constants ( printparticles,&observ,observ.iso.isospin );

  return TMatrixElement::GetMatrixElement( w,k,cos,pk,printparticles,&observ,
					   mn, mk, my);
}

//----------------------------------------------------------------------
/*! \brief Returns the hadronic current J^\mu.
 * \param iso isospin channel
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 */
FourVector< GammaStructure > TStrangeCalc::GetCurrent(int iso, double w, double k, double cos, double pk) const
{
  // check if object was initialized for requested isospin channel
  if ( getisospinbasechannel ( iso ) != observ.iso.iso_base )
    {
      cerr << "ERROR in TStrangeCalc::GetCurrent(...): "
	   << "TStrangeCalc object has not been initialized for "
	   << "isospin channel no." << iso << ".\n";
      exit ( 1 );
    }
  
  // Make duplicates of some data member
  Observable tempObserv = observ;
  double vertex[MAXNDIM];
  for(int i=0; i<MAXNDIM; ++i) vertex[i] = fVertex[i];
  
  // initialize particles structure
  tempObserv.iso.isospin = iso;
  insertvertex ( printparticles,particles[iso], varinfo, &tempObserv, vertex );
  
  // Transform coupling constants to correct isospin channel
  change_coupling_constants ( printparticles,&tempObserv,tempObserv.iso.isospin );
  
  // return the current
  return TMatrixElement::GetMatrixElement( w,k,cos,pk,printparticles,&tempObserv)->GetCurrent();
}

//----------------------------------------------------------------------
/*! \brief Returns the hadronic current J^\mu. The kinematics are allowed to have off-mass-shell particles.
 * \param iso isospin channel
 * \param w photon energy in the CMF
 * \param k photon momentum in the CMF
 * \param costhkcm cos(theta_K) in the CMF
 * \param pK Kaon momentum in the CMF
 * \param mn nucleon mass (in MeV)
 * \param mk kaon mass (in MeV)
 * \param my hyperon mass (in MeV)
 */
FourVector< GammaStructure > TStrangeCalc::GetCurrent(int iso, double w, double k, double cos, double pk, double mn, double mk, double my) const
{
  // check if object was initialized for requested isospin channel
  if ( getisospinbasechannel ( iso ) != observ.iso.iso_base )
    {
      cerr << "ERROR in TStrangeCalc::GetCurrent(...): "
	   << "TStrangeCalc object has not been initialized for "
	   << "isospin channel no." << iso << ".\n";
      exit ( 1 );
    }
  
  // Make duplicates of some data member
  Observable tempObserv = observ;
  double vertex[MAXNDIM];
  for(int i=0; i<MAXNDIM; ++i) vertex[i] = fVertex[i];
  
  // initialize particles structure
  tempObserv.iso.isospin = iso;
  insertvertex ( printparticles,particles[iso], varinfo, &tempObserv, vertex );
  
  // Transform coupling constants to correct isospin channel
  change_coupling_constants ( printparticles,&tempObserv,tempObserv.iso.isospin );
  return TMatrixElement::GetMatrixElement( w,k,cos,pk,printparticles,&tempObserv,mn, mk, my)->GetCurrent();
}

//----------------------------------------------------------------------
/*! \brief Returns the varinfo-vector for the parameters of fVertex
 * This is necessary e.g. for Minuit when the error bounds are printed.
 * \return vector of the form vector< pair < pair< string,string>, vector < string> > >,
 * i.e. particlelist < pair < pair< class,nickname>, parameterlist < "fixed"|"limits"|"nolimits" > > >
 */
varinfovector TStrangeCalc::GetVarinfoVector() const
{
  varinfovector vvect;
  // loop over diagrams
  for(int i=0; i< CLASSMAX; i++) // CLASSMAX+1 -> formfact
    {
      // loop over particles of same CLASS
      for(int j=0; j<printparticles[i].particount; j++)
	{
	  std::vector < string> particleinfo;
	  for (int k=0 ; k<get_nr_par(i); k++ )// loop over possible free variables (<MAXNRFVAR)
	    particleinfo.push_back( (varinfo[i].partic[j][k].var==0) ? "fixed": ((varinfo[i].partic[j][k].bound==0)?"nolimits": "limits"));	
	  vvect.push_back(make_pair( make_pair( string( 1, findclasslabel(i) ), printparticles[i].partic[j].nickname), particleinfo));
	}
    }
  
  // cutoff info is situated under the index CLASSMAX of the varinfo array
  if (observ.hadronformfac)
    {
      for (int i=0;i<2;i++)
	{ 
	  std::vector < string > particleinfo;
	  particleinfo.push_back( (varinfo[CLASSMAX].partic[0][i].var==0) ? "fixed": ((varinfo[CLASSMAX].partic[0][i].bound==0)? "nolimits": "limits"));
	  vvect.push_back( make_pair( make_pair( "cutoff",((i==0)?"born:":"res:")), particleinfo));
	}
    }
  
  return vvect;
}

/*! \brief Assigns a new datapoint array to the datapoints pointer.
 * Delete original datapoints arrays and replace by new newdatapoints and 
 * newdatacount. Note that these should be dynamically allocated (with new) 
 * and not be deleted after this function is called!!
 */
void TStrangeCalc::SetData ( Data** newdatapoints, int* newdatacount )
{
  for ( int iso=ISOMAX; iso>0; --iso )
    {
      datacount[iso-1]=newdatacount[iso-1];
      delete[] datapoints[iso-1];
    }
  delete[] datapoints; // discard memory allocated for old datapoints.
  datapoints = newdatapoints; // let datapoints point to newdatapoints.
  
  fDataWeight = setlabels(datapoints, datacount, &observ);
  
  //TODO Clear cache arrays...
}

/* \brief Gets the weight of the dataset attributed to this node.
 * Returns the number of matrix elements that have to be evaluated to get the
 * chi squared of the current dataset (NOT counting e.g. different components 
 * L/Lp/Ly in the case of an unpolarised variable)
 */



//----------------------------------------------------------------------

void TStrangeCalc::TestStatus ( int istat )
{
  if ( istat )
    {
      cout << "failed\n";
      exit ( 1 );
    }
  else
    {
      cout << "OK\n";
    }
}

/*! \brief Destructor for TStrangeCalc
 * releases dynamically allocated memory
 */

TStrangeCalc::~TStrangeCalc()
{
  for ( int i=0; i<ISOMAX; i++ )
    release_formfactors ( particles[i] );

  for ( int i=ISOMAX; i>0; --i )
    {
      delete[] datapoints[i-1];
      delete[] particles[i-1];
    }
  delete[] datapoints;
  delete[] particles;
  delete[] printparticles;
  delete[] varinfo;
}

/*! \brief  Get the value of a particle's property. 
 * Known properties are (case sensitive):
 *\verbatim
 - mass:          particle's mass
 - width:         particle's total decay width
 (should be 0 for all excepts-channel resonances)
 - G:             EM c.c. (kappa) for diagrams S,U,A,B,C,D,E,F,G
 EM c.c. (kappa_1) for diagrams H,I,J,L
 - H:             Strong c.c. for diagrams S,T,U,A,D,E,F,G
 Strong vector c.c. for diagrams B,C
 EM c.c. (kappa_2) for diagrams H,I,J,L
 - I:             Strong tensor c.c. for diagrams B,C
 Strong c.c. for diagrams H,I,J,L
 - J:             (none)
 - X:             off-shell parameter for diagrams H,I,J,L
 - Y:             off-shell parameter for diagrams H,I,J,L
 - Z:             off-shell parameter for diagrams H,I,J,L
 - r_kappa_n_p:   ratio of EM c.c. (kappa) neutron/proton for diagrams D,E,F,G
 - r_kappa_1_n_p: ratio of EM c.c. (kappa_1) neutron/proton for diagrams H,I,J,L
 - r_kappa_2_n_p: ratio of EM c.c. (kappa_2) neutron/proton for diagrams H,I,J,L
 - cutoff:        the hadronic formfactor's cutoff value. 
 Note that particlename in this case* is either "res" (resonances)
 or "born" (born diagrams)
 \endverbatim
 * \param particlename nickname of the particle you wish to modify
 * \param property parameter to vary, see above list
 * \param value value to which the parameter is set
 * 
 */
double TStrangeCalc::GetParticleProperty(const char* particlename,
					 const char* property)
{
  if(!strcmp(property, "cutoff"))
    {
      if(!strcmp(particlename, "res"))
	return observ.ffac.res_cutoff;
      else if(!strcmp(particlename, "born"))
	return observ.ffac.born_cutoff;
      else
	{
	  cerr << "WARNING in TStrangeModel::GetParticleProperty"
	       << "(const TString&,const TString&): "
	       << "unknown cutoff type.\n";
	  exit(1);
	}
    }
  else
    for(int iso = 0; iso < ISOMAX; iso++)
      {
	Class *particles = GetParticles(iso);

	for(int diagram = 0; diagram < CLASSMAX; diagram++)
	  for(int particle = 0; particle < particles[diagram].particount; particle++)
	    if(!strcmp(particlename, particles[diagram].partic[particle].nickname))
	      {
		if(!strcmp(property, "mass"))
		  return particles[diagram].partic[particle].mass;
		else if(!strcmp(property, "width"))
		  return particles[diagram].partic[particle].width;
		else if(!strcmp(property, "G"))
		  return particles[diagram].partic[particle].G;
		else if(!strcmp(property, "H"))
		  return particles[diagram].partic[particle].H;
		else if(!strcmp(property, "I"))
		  return particles[diagram].partic[particle].I;
		else if(!strcmp(property, "J"))
		  return particles[diagram].partic[particle].J;
		else if(!strcmp(property, "X"))
		  return particles[diagram].partic[particle].X;
		else if(!strcmp(property, "Y"))
		  return particles[diagram].partic[particle].Y;
		else if(!strcmp(property, "Z"))
		  return particles[diagram].partic[particle].Z;
		else if(!strcmp(property, "r_kappa_n_p"))
		  return particles[diagram].partic[particle].r_kappa_n_p;
		else if(!strcmp(property, "r_kappa_1_n_p"))
		  return particles[diagram].partic[particle].r_kappa_1_n_p;
		else if(!strcmp(property, "r_kappa_2_n_p"))
		  return particles[diagram].partic[particle].r_kappa_2_n_p;
		else
		  {
		    cerr << "WARNING in TStrangeModel::GetParticleProperty"
			 << "(const TString&,const TString&): "
			 << "unknown particle property.\n";
		    exit(1);
		  }
	      }
      }

  return 0; // dummy return
}

/*! \brief  Sets the value of a particle's property. 
 * Known properties are (case sensitive):
 *\verbatim
 - mass:          particle's mass
 - width:         particle's total decay width
 (should be 0 for all excepts-channel resonances)
 - G:             EM c.c. (kappa) for diagrams S,U,A,B,C,D,E,F,G
 EM c.c. (kappa_1) for diagrams H,I,J,L
 - H:             Strong c.c. for diagrams S,T,U,A,D,E,F,G
 Strong vector c.c. for diagrams B,C
 EM c.c. (kappa_2) for diagrams H,I,J,L
 - I:             Strong tensor c.c. for diagrams B,C
 Strong c.c. for diagrams H,I,J,L
 - J:             (none)
 - X:             off-shell parameter for diagrams H,I,J,L
 - Y:             off-shell parameter for diagrams H,I,J,L
 - Z:             off-shell parameter for diagrams H,I,J,L
 - r_kappa_n_p:   ratio of EM c.c. (kappa) neutron/proton for diagrams D,E,F,G
 - r_kappa_1_n_p: ratio of EM c.c. (kappa_1) neutron/proton for diagrams H,I,J,L
 - r_kappa_2_n_p: ratio of EM c.c. (kappa_2) neutron/proton for diagrams H,I,J,L
 - cutoff:        the hadronic formfactor's cutoff value. 
 Note that particlename in this case* is either "res" (resonances)
 or "born" (born diagrams)
 \endverbatim
 * \param particlename nickname of the particle you wish to modify
 * \param property parameter to vary, see above list
 * \param value value to which the parameter is set
 * 
 */
void TStrangeCalc::SetParticleProperty ( const char* particlename,
					 const char* property,
					 double value )
{
  if (!strcmp(property,"cutoff"))
    {
      if(!strcmp(particlename,"res"))
	{
	  observ.ffac.res_cutoff = value;
	  for ( int iso=0; iso<ISOMAX; ++iso )
	    {
	      update_strong_formfactor(&observ, GetParticles(iso));
	    }
	}
      else if (!strcmp(particlename,"born"))
	{
	  observ.ffac.born_cutoff = value;
	  for ( int iso=0; iso<ISOMAX; ++iso )
	    {
	      update_strong_formfactor(&observ, GetParticles(iso));
	    }
	}
      else
	{
	  cerr << "WARNING in TStrangeModel::SetParticleProperty"
	       << "(const TString&,const TString&,double): "
	       << "unknown cutoff type.\n";
	}
    }
  else
    {
      for ( int iso=0; iso<ISOMAX; ++iso )
	{
	  Class *particles = GetParticles ( iso );

	  for ( int diagram=0; diagram<CLASSMAX; ++diagram )
	    {
	      for ( int particle=0; particle<particles[diagram].particount; ++particle )
		{

		  if ( !strcmp ( particlename,particles[diagram].partic[particle].nickname ) )
		    {

		      if ( !strcmp ( property,"mass" ) )
			particles[diagram].partic[particle].mass = value;
		      else if ( !strcmp ( property,"width" ) )
			particles[diagram].partic[particle].width = value;
		      else if ( !strcmp ( property,"G" ) )
			particles[diagram].partic[particle].G = value;
		      else if ( !strcmp ( property,"H" ) )
			particles[diagram].partic[particle].H = value;
		      else if ( !strcmp ( property,"I" ) )
			particles[diagram].partic[particle].I = value;
		      else if ( !strcmp ( property,"J" ) )
			particles[diagram].partic[particle].J = value;
		      else if ( !strcmp ( property,"X" ) )
			particles[diagram].partic[particle].X = value;
		      else if ( !strcmp ( property,"Y" ) )
			particles[diagram].partic[particle].Y = value;
		      else if ( !strcmp ( property,"Z" ) )
			particles[diagram].partic[particle].Z = value;
		      else if ( !strcmp ( property,"r_kappa_n_p" ) )
			particles[diagram].partic[particle].r_kappa_n_p = value;
		      else if ( !strcmp ( property,"r_kappa_1_n_p" ) )
			particles[diagram].partic[particle].r_kappa_1_n_p = value;
		      else if ( !strcmp ( property,"r_kappa_2_n_p" ) )
			particles[diagram].partic[particle].r_kappa_2_n_p = value;
		      else
			cerr << "WARNING in TStrangeModel::SetParticleProperty"
			     << "(const TString&,const TString&,double): "
			     << "unknown particle property.\n";

		    }
		  //if(!strcmp(particlename,particles[diagram].partic[particle].nickname))

		} // for(int particle=0; particle<particles[diagram].particount;
	      //++particle)
	    } // for(int diagram=0; diagrams<CLASSMAX; ++diagrams)
	} // for(int iso=1; i<ISOMAX; ++i)
    }
}

/*! \brief Sets the model and implementation type of a TStrangeCalc.
 *  \param modelType theory, e.g. offshell or consistent
 *  \param cgln implementation type:
 *		cgln decomposition (cgln=true) 
 *  		or the explicit diagram summation (cgln=false).
 * All isospin channels initialised by the fit_specification are changed
 * The electroproduction variables will only be initialized in the case of electroproduction!
 */
void TStrangeCalc::SetModelType(const char* modelType, const bool cgln) 
{
  SetModelImplementation(cgln);
  
  strcpy(observ.modelType,modelType);
  

  
  // Loop over all isospins
  for(int i =0; i<observ.iso.nr_iso_channels; i++) 
    {
      int iso = observ.iso.iso_channel[i];
      observ.iso.isospin = iso;
      // Loop over classes
      for(int diag = 10; diag<14; diag++) // couplings for spin 5/2 are always gauge invariant. hence diag<14
	{
	  // Loop over all particles per class
	  for (int p = 0; p < particles[iso][diag].particount; p++)
	    {
	      //Specify wether gauge-invariant couplings are used or not, 
	      if (!strcmp(observ.modelType,"offshell"))
		particles[iso][diag].partic[p].gic = 0;
	      else
		particles[iso][diag].partic[p].gic = 1;
	    }  
	}
      release_formfactors(particles[iso]);
      if ( observ.fit.elec_diffcs[iso] )
	{
	  observ.electroprod = 1;
	  em_formfactor_specification(&observ,particles[iso]);
	  observ.electroprod = 0;
	}
      // set correct strong form factors
      if( observ.hadronformfac )
	strong_formfactor_specification(&observ,particles[iso]);

    }
}

/*! \brief Sets the implementation type of a TStrangeCalc.
 *  \param cgln implementation type:
 *		cgln decomposition (cgln=true) 
 *  		or the explicit diagram summation (cgln=false).
 */
void TStrangeCalc::SetModelImplementation(const bool cgln) 
{
  observ.cgln = cgln?1:0;
}

/*!
 * \brief Initializes fVertex[] and fLimits[]
 * based on the information in the particles[observ.iso.iso_base]
 * and the varinfo[] arrays.
 *
 * The elements of fVertex[] are generated on the basis of the set
 * of coupling constants from the "coupl.iso.*" file.
 *
 * The fLimits[] array is generated in order to be used by the fitter.
 * The information from the varinfo[] structure is used.
 */
int TStrangeCalc::MakeStartvertex()
{
  int born_cc_stored = 0;
  
  fDim = 0;
  
  /*
   * Create a vertex starting from the values in the coupl.iso.*
   * file which are stored in the "particles[]" structure.
   */
  
  // Loop over all diagrams
  for(int i=0; i<CLASSMAX; i++)
    {
      // Loop over all particles of same CLASS
      for(int j=0; j<particles[observ.iso.iso_base][i].particount; j++)
	{
	  // Born channels
	  if(i <= 2 && !born_cc_stored)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
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
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][0]);
	    }
	  // B- and C-channel
	  else if(i == 4 || i == 5)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][0]);
	      
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].I,
				 varinfo[i].partic[j][1]);
	    }
	  // D,E,F,G-channel
	  else if(i >= 6 && i <= 9)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][0]);
	    }
	  // H,I,J,L-channel
	  else if(i >= 10 && i <= 13)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].G,
				 varinfo[i].partic[j][0]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][1]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].X,
				 varinfo[i].partic[j][2]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Y,
				 varinfo[i].partic[j][3]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Z,
				 varinfo[i].partic[j][4]);
	    }
	  // M,N,O,P-channel
	  else if(i >= 14 && i <= 17)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].G,
				 varinfo[i].partic[j][0]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][1]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].X,
				 varinfo[i].partic[j][2]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Y,
				 varinfo[i].partic[j][3]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Z,
				 varinfo[i].partic[j][4]);
            }
	  // R-channel
	  else if(i == 18)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][0]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].I,
				 varinfo[i].partic[j][1]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].X,
				 varinfo[i].partic[j][2]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Y,
				 varinfo[i].partic[j][3]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].Z,
				 varinfo[i].partic[j][4]);

	    }
	  // V-channel
	  else if(i == 19)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].G,
				 varinfo[i].partic[j][0]);
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].H,
				 varinfo[i].partic[j][1]);
	    }
  	  // W-channel
	  else if(i == 20)
	    {
	      store_vertex_value(fVertex, fLimits, &fDim, 
				 particles[observ.iso.iso_base][i].partic[j].G,
				 varinfo[i].partic[j][0]);
	    }

	  
	  if(fDim > MAXNDIM)
	    error_exit("ndim > MAXNDIM!!\n");
	  
	} // end particle loop
    } // end diagram loop
  
  
  // Store the hadronic form factor values
  if(observ.hadronformfac)
    {
      // form factors are stored in CLASSMAX position
      int i = CLASSMAX;
      int j = 0;
      
      store_vertex_value(fVertex, fLimits, &fDim, 
			 observ.ffac.born_cutoff,
			 varinfo[i].partic[j][0]);
      
      store_vertex_value(fVertex, fLimits, &fDim, 
			 observ.ffac.res_cutoff,
			 varinfo[i].partic[j][1]);
    }
  
  
  if(fDim > MAXNDIM)
    error_exit("ndim > MAXNDIM!!\n");

  return 0;
}


/*!
 * \brief Plots a particles[] structure
 * in the form of  the "coupl.iso.*" file and the varinfo in the form of the
 * "varinfo.*" file.
 *
 * \param particles a particles structure for one isospin channel
 * \param ofp output file stream						\
 */
int TStrangeCalc::PlotParameterInfo(FILE* ofp, Class particles[])
{

  FILE *ifp;
  char filename[MAX_LOCATION_STRING], string[200], isospin[2];


  fprintf(ofp, "\n###############################################\n\n");
  fprintf(ofp, "<coupl.iso>\n");
  fprintf(ofp, "Class  NickName   G*           H*          I*          ");
  fprintf(ofp, "X*          Y*          Z*\n");
  fprintf(ofp, "-----  --------   --           --          --          ");
  fprintf(ofp, "--          --          --\n");
  fprintparticspec(ofp, particles, observ);
  fprintf(ofp, "</coupl.iso>\n\n\n");
    

  isospin[0] = observ.iso.iso_base + 48;
  isospin[1] = '\0';

  strcpy(filename, observ.inFolder );
  strcat(filename, "numinput/varinfo.");
  strcat(filename, isospin); 
  ifp = fopen(filename, "r");
      
  fprintf(ofp, "\n###############################################\n\n");
  fprintf(ofp, "<varinfo>\n");
    
  if(fgets(string, 200, ifp) == NULL)
    error_exit("Error in reading \"varinfo\" while making log header");
    
  while(string[0] != 'c')
    {	
      fprintf(ofp, "%s", string);
      fgets(string, 200, ifp);
    }
    
  fprintf(ofp, "%s", string);
  fgets(string, 199, ifp);
  fprintf(ofp, "%s", string);
  fprintf(ofp, "</varinfo>\n\n");

  fprintf(ofp, "##########################################\n\n");	

  fclose(ifp);
  return 0;

}


/*!
 * \brief Sets the name of a label
 * The name will look like "(particle name).(suffix)" and ndim will be incremented if the
 * property is variable.
 *
 * \param ndim Count of variable particles (will be incremented if particle is variable)
 * \param suffix Most often the coupling constant name
 * \param partic_prop particle coupling constant
 * \param partic_info tells whether the  coupling constant is variable 
 */
void TStrangeCalc::StoreName(int& ndim, char suffix[], 
			     Properties partic_prop, Celinfo partic_info)	
{
  if(partic_info.var)
    {
      strcpy(fName[ndim], partic_prop.nickname);
      strcat(fName[ndim], ".");
      strcat(fName[ndim], suffix);
	
      ndim++;
    }
}


/*!
 * \brief Initializes fName[][]
 * on the basis of the varinfo[] structure, particular name is
 * assigned to every vertex point, in order to keep track of its
 * variation while fitting.
 */
int TStrangeCalc::AssignVertexName()
{
  int ndim = 0;
  int born_cc_stored = 0;
  static char label[5];
    
  for(int i=0; i<CLASSMAX; i++)
    for(int j=0; j<particles[observ.iso.iso_base][i].particount; j++)
      if(i <= 2 && !born_cc_stored)
	{
	  strcpy(label, "H");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	    
	  born_cc_stored = 1;
	}	
      else if(i == 3)
	{
	  strcpy(label, "H");
	  StoreName(ndim, label, 
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	}
      else if(i == 4 || i == 5)
	{
	  strcpy(label, "H");
	  StoreName(ndim,label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	    
	    
	  strcpy(label, "I");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	}
      else if(i >= 6 && i <= 9)
	{
	  strcpy(label, "H");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	}
      else if(i >= 10 && i <= 13)
	{
	  strcpy(label, "G");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	  strcpy(label, "H");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	  strcpy(label, "X");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][2]);
	  strcpy(label, "Y");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][3]);
	  strcpy(label, "Z");
	  StoreName(ndim, label, 
		    particles[observ.iso.iso_base][i].partic[j],	
		    varinfo[i].partic[j][4]);
	}
      else if(i >= 14 && i <= 17)
	{
	  strcpy(label, "G");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	  strcpy(label, "H");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	  strcpy(label, "X");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][2]);
	  strcpy(label, "Y");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][3]);
	  strcpy(label, "Z");
	  StoreName(ndim, label, 
		    particles[observ.iso.iso_base][i].partic[j],	
		    varinfo[i].partic[j][4]);
	}
      else if(i == 18)
	{
	  strcpy(label, "G");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	  strcpy(label, "I");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	    
	}
      else if(i == 19)
	{
	  strcpy(label, "G");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	  strcpy(label, "H");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	  strcpy(label, "I");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][2]);
	}
      else if(i == 20)
	{
	  strcpy(label, "G");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][0]);
	  strcpy(label, "I");
	  StoreName(ndim, label,
		    particles[observ.iso.iso_base][i].partic[j],
		    varinfo[i].partic[j][1]);
	}

    
    
  if(observ.hadronformfac)
    {
      /* Store the hadronic form factor values */

      int i = CLASSMAX;
      int j = 0;

      if(varinfo[i].partic[j][0].var)
	strcpy(fName[ndim++], "B.cutoff");
	    
      if(varinfo[i].partic[j][1].var)
	strcpy(fName[ndim++], "R.cutoff");
	    
    }	
      
  
  if(ndim > fDim) {
    cerr << "ERROR in TStrangeCalc::AssignVertexName(): "
	 << "dimension mismatch." << endl;
    assert( ndim==fDim );
  }
      
  return 0;	
}


/*!
 * \brief Sets up the datapoints structures and particle properties
 */
int TStrangeCalc::SetupStrangeStruc(FILE* ifp)
{
  /* First we construct the
   * particle specifications in the "particles[iso][]" array
   * and import the data points for the base isospin channel.
   */
  SetupStrangeStruc(ifp,observ.iso.iso_base);
  
  for(int i=0; i<observ.iso.nr_iso_channels; i++) 
    {
      /*
       * Next, for all the channels under investigation, construct the
       * particle specifications in the "particles[iso][]" array
       * and import the data points.
       */

      int iso = observ.iso.iso_channel[i];
      
      if( iso!=observ.iso.iso_base )
	SetupStrangeStruc(ifp,iso);
    }

  setlabels(datapoints,datacount,&observ);

  return 0;
}

/*!
 * \brief Sets up the datapoints structures and particle properties for one isospin channel
 * Calling this function for a non-base isospin channel *before* the base isospin channel
 * has been initialized will result in an error.
 */
void TStrangeCalc::SetupStrangeStruc(FILE* ifp, int iso)
{
  /*
   * This function does the "setup" of the datapoints structures
   * and particle properties
   */
  
  observ.photoprod = 0;
  observ.electroprod = 0;
  observ.kaoncapture = 0;

  /*
   * For the channel under investigation, construct the
   * particle specifications in the "particles[iso][]" array
   * and import the data points.
   */
 
  observ.iso.isospin = iso;
  
  /* 
   * The particle properties are read in from the 
   * ./input/numinput/coupl.iso.* file. 
   */
  if( iso==observ.iso.iso_base )
    particlespecify(particles[iso], &observ, ifp, stdout);
  else
    copyparticles(particles[iso],particles[observ.iso.iso_base],&observ);
  
  
  /*
   * We need to change isospin channel
   * change_isospin_channel() can not be used however,
   * because the coupling constants may not be changed yet!
   * (this would disturb makestartvertex()).
   */
  
  modify_nicknames(particles[iso],&observ);
  import_particle_properties(particles[iso],&observ);
  if(iso > 3) 
    magnetic_transition_ratio(particles[iso],&observ);
  
  /*
   * Store Form factor info in particles[]
   */
  if( observ.hadronformfac )
    strong_formfactor_specification(&observ,particles[iso]);
  
  if ( observ.fit.elec_diffcs[iso] )
    {
      observ.electroprod = 1;
      em_formfactor_specification(&observ,particles[iso]);
      observ.electroprod = 0;
    }
      
      
  /* Store the data in an easy accessible structure */
  observ.iso.isospin = iso;

  getdatastructure(datapoints[iso], &datacount[iso], &observ, 
		   particles[iso]);
}
