/*!
 * \file GaugeInvariance_test.cpp
 *
 * Relatively low-level test for a MatrixElement object.
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 
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
#include <TDataset.h>
#include <TMatrixElement.h>
#include <TCGLNDiagnostic.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;
using std::fixed;
using std::scientific;
using std::ios_base;
using std::complex;

void printMatrixElement(TMatrixElement* me)
{
  ios_base::fmtflags originalFormat = cout.flags(); // use cout flags function to save original format

  // Calculate average
  double average = 0.;
  double reals[12];
  double imags[12];
  for (int lp=0;lp<2;++lp) {
    for ( int ly=0;ly<2;++ly) {
      for (int l=0;l<3;++l) {
	const complex<double> element = me->calculateM(l,lp,ly);
	reals[lp+2*ly+4*l] = element.real();
	imags[lp+2*ly+4*l] = element.imag();
	average += (fabs(reals[lp+2*ly+4*l]) + fabs(imags[lp+2*ly+4*l]));
      }
    }
  }
  const int precision = abs((int)(floor(log10(average/24.e11))));
 
  // Print
  for (int lp=0;lp<2;++lp) {
    for ( int ly=0;ly<2;++ly) {
      for (int l=0;l<3;++l) {
	cout << "MatrixElement("<<lp <<","<<ly<<") = ("
	     << fixed << setprecision(precision) << reals[lp+2*ly+4*l]
	     << ","
	     << fixed << setprecision(precision) << imags[lp+2*ly+4*l]
	     << ")"
	     << endl;
      }
    }
  }

  cout.flags( originalFormat ); // restore format
}

int main(int argc, char *argv[])
{
  if (argc>1)
  {
    TDataset::SetDataFolder("../../share/data/");
    TStrangeModel model("model","testing Gauge Invariance");
    model.SetStrangeModel(argv[1]);
    int iso = argc>2 ? atoi(argv[2]): 1;
    if (argc>3) 
    {
      if (!strcmp(argv[3],"n"))
	model.SetModelType(iso,"offshell",0);
      else if (!strcmp(argv[3],"y"))
	model.SetModelType(iso,"consistent",0);
      else
	model.SetModelType(iso,argv[3],1);
    }  
    
    Class * particles = model.GetTStrangeCalc(iso).GetParticles(iso); 
    for ( int diagram=3; diagram<CLASSMAX; diagram++ )
    {
      for ( int particle=0; particle<particles[diagram].particount; particle++ )
      {
	particles[diagram].partic[particle].mass=1800.;
	cout << "Setting mass of particle " <<  particle << "(diagram "<< diagram << ") to 1800 MeV for testing purposes\n";
      } // particle loop
    } // diagram loop
    TKinematics tk("tk","tk",iso,"costhkcm:w:qsquared",0.4,1900.,0.);   
    
    for (int w=1900;w<=2000;w+=100)
    {
      tk.SetVar(2,w);
      
      cout << "-------------------------------------\n W = "<<w<< " MeV \n";
      TMatrixElement* me = model.GetTStrangeCalc(iso).GetMatrixElement(iso, tk.GetWcm(), tk.GetKcm(), tk.GetCosthkcm(), tk.GetPkcm());
      printMatrixElement(me);
      
      cout << setprecision (4);
      cout << "\nCurrent: \n";
      me->GetCurrent().print();	
      me->Print();
    }
    
    //-------------------------
    if (argc>4)
    {
      if (!strcmp(argv[4],"full"))
      {
	cout <<"_____________________________________________________________________________\n"
	    <<"Checking whether the difference between the old and new currents is proportional to k4vect...\n";
	
	// Diagrammatic consistent couplings
	model.SetModelType(iso,"consistent",false);
        FourVector<GammaStructure> oldcurrent = model.GetCurrent(iso,tk); 
	// CGLN consistent couplings
	model.SetModelType(iso,"consistent",true);
        FourVector<GammaStructure> newcurrent = model.GetCurrent(iso,tk);
	
	cout << "The old current is:\n";
	oldcurrent.print();
	cout<<"The new current is:\n";
	newcurrent.print();

	newcurrent -= oldcurrent;
	FourVector< double> k4vect  = FourVector<double> ( tk.GetWcm(),0.,0.,tk.GetKcm());             // photon 4vector
	bool proportional=true;
	double uflow = 1.e-7;
	
	cout << "Photon four-momentum is :\n";
	for (int mu=0; mu<4;mu++)
	  cout << k4vect[mu] << '\n';
	
	cout << "The difference between the two currents is:\n";
	for (int mu=0; mu<4;mu++)
	{
	  newcurrent[mu].value().print();
	  if (abs(k4vect[mu])>uflow)
	  {
	    newcurrent[mu] *= 1./k4vect[mu];
	    if (mu>0)
	    {
	      GammaStructure diff  = newcurrent[mu]-newcurrent[mu-1];
	      for(int i=0;i<4;i++)
	      {
		for(int j=0;j<4;j++)
		{
		  if (abs(diff.value()(i,j)) > uflow)
		  {
		    proportional=false;
		    break;		
		  }
		}
	      }
	    }
	  }
	  else
	    for(int i=0;i<4;i++)
	    {
	      for(int j=0;j<4;j++)
	      {
		if (abs(newcurrent[mu].value()(i,j))>uflow)
		{
		  proportional=false;
		  break;
		}
	      }
	      if (!proportional) 
		break;
	    }
	}
	cout << "Difference between old and new currents "<< (proportional?"is ": "is not ") << "proportional to k4vect!\n";

	//-------------------------
	
	cout <<"_____________________________________________________________________________\n"
	    <<"Checking whether the CGLN basis is zero for k=epsilon (photoproduction only)\n.";
	TCGLNDiagnostic diagme( tk.GetWcm(), tk.GetKcm(), tk.GetCosthkcm(), tk.GetPkcm(), model.GetTStrangeCalc(iso).GetParticles(iso),model.GetTStrangeCalc(iso).GetObserv());
	
	for (int i=0;i<6;i++)
	{
	  diagme.CalcMwithkvect(i).print();
      //  tk.SetVarRange(1,-.8,1.,5);
      //  tk.SetVarRange(3,.1e6,1.e6,5);
	}
      }

    }

  }
  else
  {
    cerr << "Usage: " << argv[0] << " fit-specification [iso [modeltype [full]]] "<< endl;
  }
  return 0;
}


