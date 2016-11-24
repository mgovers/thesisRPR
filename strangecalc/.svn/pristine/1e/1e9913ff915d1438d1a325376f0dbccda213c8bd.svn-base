/*!
 * \file TStrangeCalc.h
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

#ifndef TStrangeCalc_h
#define TStrangeCalc_h

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <strange_func.h>
#include <fitting.h>
#include <io_specific.h>
#include <Structures.h>

#include <FourVector.h>
#include <GammaStructure.h>
#include <iostream>
#include <vector>
#include <string>
#include <utility>

typedef std::vector<double> individual;
typedef std::pair < std::string , std::string > stringpair ;
//typedef std::vector < std::pair < stringpair , std::vector < double > > > couplisovect;
typedef std::vector < std::pair < stringpair , std::vector < std::string> > > varinfovector;

// Forward declarations
class TMatrixElement;

class TStrangeCalc {
private:
  friend class TStrangeModel;
  friend class TMultiModel;

private:
    int    fDim;       		//!< the number of used elements in fVertex[]
    int    fDataWeight;		//!< computational load of the current dataset (no of ME's to evaluate)
    double fChiSquared;	       	//!< the current chi-squared value 
    double fVertex[MAXNDIM];	//!< array of parameters 
    Limits fLimits[MAXNDIM];	//!< limits between which parameters can vary
  static const int kNameLength = 20;  //!< length of parameter labels
    char fName[MAXNDIM][kNameLength]; //!< parameter labels
    void   TestStatus(int);     //!< Not needed outwith class

   // Interface to strangecalc...
    Observable observ;
    Data** datapoints;          // [ISOMAX][DATAMAX]
    int datacount[ISOMAX];   
    Class** particles;          //[ISOMAX][CLASSMAX]
    Class* printparticles;      //[CLASSMAX]
    Varinfo* varinfo;           //[CLASSMAX+1]
    int pho_diff_setup, pho_rec_setup, pho_tar_setup, pho_pho_setup;
    int elec_diff_setup;
    short isospin;


public:
  // Class operations
  TStrangeCalc() {};                        // default constructor
  TStrangeCalc(FILE *log, const char* fit_spec, const char* datafolder=0);     // constructor, calls Init(...)
public:
  virtual ~TStrangeCalc();	            // default destructor
  TStrangeCalc(const TStrangeCalc& toCopy); // explicit copy constructor
  TStrangeCalc& operator= (const TStrangeCalc& toCopy); // explicit assignment operator
  
  // Class arrays and others are printed directly from class
  // methods. Meant for debugging...
  void   	PrintVertex() const;
  void   	PrintLimits() const;
  void   	PrintParticles(FILE* stream=stdout);
    
  // Class "numbers" are Set and Got as per convention...
  Observable      *GetObserv() {return &observ;}	//!< get the observ struct
  Data           **GetDatapoints() {return datapoints;} //!< get the array of data points
  int             *GetDatacount() {return datacount;}	//!< get the number of datapoints
  Class           *GetParticles(int iso) const {return particles[iso];} //!< get the particles of a certain isospin channel
  double           GetParticleProperty(const char* particlename, const char* property);
  double          *GetVertex() {return fVertex;};		//!< get array of current values of the variable parameters
  double 	   GetChiSquared() const {return fChiSquared;};	//!< get the current chi squared value
  Limits 	   GetLimit(int i) const {return fLimits[i];};	//!< get the limits of the i'th variable parameter
  int    	   GetDim() const {return fDim;};		//!< get the number of variable parameters
  const char      *GetName(int i) const {return fName[i];};	//!< get the i'th parameter label
  int		   GetDataWeight() const {return fDataWeight;}  //!< get the number of matrix elements to evaluate to get the chi squared for this dataset
  
  int              GetNDF() const;	//!< get the weighted amount of degrees of freedom
  double           GetCalcpoint(Data*,double,double,double,double,double,double=0.0,double=0.0); //!< calculate the model prediction for a datapoint
  TMatrixElement*  GetMatrixElement(int iso, double w, double k, double cos, double pk); // do not delete the MatrixElement*
  TMatrixElement*  GetMatrixElement(int iso, double w, double k, double cos, double pk, double mn, double mk, double my); // do not delete the MatrixElement*
  FourVector<GammaStructure> GetCurrent(int iso, double w, double k, double cos, double pk) const ;
  FourVector<GammaStructure> GetCurrent(int iso, double w, double k, double cos, double pk, double mn, double mk, double my) const;
  varinfovector    GetVarinfoVector() const; //!< get varinfo vector< pair < pair< string,string>, vector < string> > >

  
  void   	   SetChiSquared();	// wrapper around chifunc
  void             SetData ( Data** newdatapoints, int newdatacount[] ); //!< assign new block of datapoints
  void   	   SetVertex( const double [] );
  void             SetVertex( const individual);
  void             SetParticleProperty(const char* particlename,const char* property,double value);
  void 		   SetModelType(const char* modelType, const bool cgln);//!< set the theory and implementation type only for the assigned isospins and reaction types!
void 		   SetModelImplementation(const bool cgln);

private:
  int              MakeStartvertex();
  int              PlotParameterInfo(FILE* ofp, Class particles[]);
  int              AssignVertexName();
  void             StoreName(int& ndim, char suffix[], Properties partic_prop, Celinfo partic_info);
  int              SetupStrangeStruc(FILE* ifp);
  void             SetupStrangeStruc(FILE* ifp, int iso);
};

//----------------------------------------------------------------------

#endif // #ifdef TStrangeCalc_h

