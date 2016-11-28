/*!
 * \file CCBoundaries.cpp
 *
 * Check the "physical" boundaries for coupling constants.
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 
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

#include <TStrangeModel.h>
#include <TStrangeCalc.h>
#include <TKinematics.h>
#include <TCalcInfo.h>
#include <TMatrixElement.h>
#include <TCGLNDiagnostic.h>
#include <TString.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::ios_base;
int main ( int argc, char *argv[] )
{
  if ( argc>1 )
  {
    TStrangeModel model ( "model","testing sensible boundaries for different particles." );
    model.SetStrangeModel ( argv[1] );
    TString modeltype = argc>2 ? argv[2] : "consistent";
    int iso = argc>3 ? atoi ( argv[3] ) : 1;
    model.SetModelType ( iso, modeltype, 1 );
    TKinematics kinematics ( "tk","tk",iso,"costhkcm:wlab:qsquared",1.,1500.,0. );
    kinematics.SetVarRange ( 2,1500.,3000.,100 );
    TCalcInfo tc(TCalcInfo::kPhoto, iso,"tcs" );

    TStrangeCalc* calc = &(model.GetTStrangeCalc ( iso ));
    int ndim = calc->GetDim() - 1;
    vector<double> ind ( ndim+1,1.0 );
    vector<double> max_cc (ndim,0.0);
    vector<double> min_cc (ndim,0.0);
    double tcs_limit = 25.; // 5 microbarn is the order of magnitude

    int iterations=9; // accuracy of 0.25
    for ( int sign=1; sign>=-1; sign-=2)
    {
      cout << "Sign: " <<  sign << ".\n";
      for ( int cc=0;cc<ndim;cc++ )
      {
	cout << "CC:" <<  cc << ".\n";
	for (int cutoff=1000; cutoff<=2500 ; cutoff+=300)
	{
	  cout << "Cutoff:" <<  cutoff << ".\n";
	  ind[ndim] = cutoff; // set the correct cutoff value

	  double range = sign*64.;// actually this is the halfway point (range is actually 128)
	  double val = range;
	  double cclimitvalue = range*2; //lowest value at which the tcs is too high.
	  cout << "Setting value to...  ";
	  // now perform a binary search to determine the maximum total cross section as a function of energy.
	  for (int step = 0; step < iterations; step++)
	  {
	    ind[cc]=val;
	    cout << val<< '\t';
	    
	    double xmax=0.;// maximum for this cutoff and parameter cc.
	    double ymax=0.;  
	    
	    calc->SetVertex(ind);// set the vertex.
	    
	    double *x = kinematics.GetPhysicalVarArray ( "wlab" );
	    double *y = model.GetPhysicalCalcpoints ( kinematics,&tc );
	    
	    for ( int i = 0; i < kinematics.GetNumberOfPhysicalSteps(); i++ )
	    {
	      
	      if ( y[i]> ymax )
	      {
		// A new maximum has been found!
		xmax=x[i];
		ymax=y[i];
		cout << "; max("<< xmax << ","<< ymax<< ","<< cutoff << ")";
	      }
	    }
	    delete [] x;
	    delete [] y;
	    range *=0.5; // halve the range in which were are searching
	    
	    if (ymax >tcs_limit)
	    {
	      // set the new upper bound
	      cclimitvalue=val;
	      val = val-range;
	    }
	    else
	      val = val+range; 
	  }
	  // Now we check whether this new value is less restrictive than previously found values; in that case we have to broaden our boundaries.
	  if (sign>0) 
	  {
	    cout << "\nMax value for cc "<< cc<< " at cutoff "<< cutoff << ": "<< cclimitvalue << "\n";
	    if (cclimitvalue > max_cc[cc])
	    {
	      max_cc[cc] = cclimitvalue;
	      cout << "updated!\n";
	    }
	  }
	  
	  else 
	  {
	    cout << "\nMin value for cc "<< cc<< " at cutoff "<< cutoff << ": "<< cclimitvalue<< "\n";
	    if(cclimitvalue < min_cc[cc]) 
	    {
	      min_cc[cc]=cclimitvalue;
	      cout << "updated!\n";
	    }
	  }
	  //set back to 1.0
	  ind[cc]=1.0;
	}
	
      }
    }

    TString basename = "limits_"+ modeltype;
    ofstream limits;
    
    Class * particles = calc->GetParticles(iso);
    for(int i=0; i< CLASSMAX; i++) // CLASSMAX+1 -> formfact
    {
      // loop over particles of same CLASS
      if (particles[i].particount==1)
      {
	TString fname = basename + "_" + particles[i].partic[0].nickname + ".txt";
	limits.open(fname, ios_base::out | ios_base::app ); // open file
	
	limits << findclasslabel(i) << '\t' << particles[i].partic[0].nickname;
	for (int cc=0; cc<ndim; cc++)
	{
	  limits << '\t'<< "f,l "<< min_cc[cc] << " " << max_cc[cc];
	}
	if(get_nr_par(i) > 2)  limits << "\t x \t x \t x";
	limits << "\n";
	limits.close(); 
      }
    }
   
   
  }
  else
  {
    cerr << "Please provide a valid fit_specification file! "<< endl;
  }
  return 0;
}


