/*!
 * \file Speed_test.cpp
 *
 * Test the relative speed of old and new calculation methods.
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
#include <TMatrixElement.h>
#include <iostream>
#include <ctime>
#include <assert.h>
#include <fstream>
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
int main(int argc, char *argv[])
{  
  int cgln_start=0;
  int cgln_stop=1;
  if (argc>1)
  {
    TStrangeModel model("model","testing implementation speeds");
    model.SetStrangeModel(argv[1]);
    int iso = argc>2 ? atoi(argv[2]): 1;
    TKinematics tk("tk","tk",iso,"costhkcm:wlab:qsquared",1.,1500.,0.);
    if (argc<4 || !strcmp(argv[3],"test") || !strcmp(argv[3],"stress"))// quick test
    {
      int numCalcs=10;
      if (argc>3 && !strcmp(argv[3],"stress"))
	numCalcs = 1000;
	
      for (int cgln=cgln_start;cgln<=cgln_stop;cgln++)
      {
	model.SetModelType(iso,"consistent",cgln);
	
	clock_t start = clock();
	for (int k=0;k<numCalcs;k++)
	{
	  TMatrixElement* me = model.GetTStrangeCalc(iso).GetMatrixElement(iso, tk.GetWcm(), tk.GetKcm(), tk.GetCosthkcm(), tk.GetPkcm()); // do not delete
	  for (int lp=0;lp<2;lp++)
	  {
	    for ( int ly=0;ly<2;ly++)
	    {
	      for (int l=0;l<3;l++)
	      {
		me->calculateM(l,lp,ly);
		if (k==numCalcs-1)
		  cout << "("<< l << "," <<lp <<","<<ly<<") = "<< me->calculateM(l,lp,ly)<< "; ";
	      }
	      if (k==numCalcs-1) cout << "\n";
	    }
	  }
	}
	
	clock_t stop = clock();
	cout << "\n" << numCalcs << " * 12 Matrix Element calculations for cgln "<< cgln <<" took " << ((double) (stop-start) )/ (1.0*CLOCKS_PER_SEC) <<"s.\n\n";
      }
      
      cout << "\n-------------------------------------------------------------------\n\n"
	      << "Performing "<< numCalcs << " chi^2 calculations...." << endl;
      double old_chisquared = 0.0;
      for (int cgln=cgln_start;cgln<=cgln_stop;cgln++)
      {
	model.SetModelType(iso,"consistent",cgln);
	
	clock_t start = clock();
	for (int k=0;k<numCalcs;k++)
	{
	  model.GetTStrangeCalc(iso).SetChiSquared();
	}
	clock_t stop = clock();
	cout << "Chi squared = " << model.GetTStrangeCalc(iso).GetChiSquared() << "\n";

	cout  << numCalcs << " chi^2 calculations (" 
	      << model.GetTStrangeCalc(iso).GetDataWeight() 
	      << " evaluations) for cgln "<< cgln << " took "
	      << ((double) (stop-start) )/ (1.0*CLOCKS_PER_SEC)<< "s."
	      << endl;
	
	if (!cgln)
	  old_chisquared = model.GetTStrangeCalc(iso).GetChiSquared();
	else if (old_chisquared != 0.0) 
	  assert(abs(old_chisquared-model.GetTStrangeCalc(iso).GetChiSquared())<1e-5);
      }
    }
    else if(!strcmp(argv[3],"full")) // check the behaviour for more chi^2 values
    {
      ofstream statfile;
      statfile.open (argc>4?argv[4]:"Speed_test_output.dat");
      statfile << "# Number of ME evaluations: " << model.GetTStrangeCalc(iso).GetDataWeight() << "\n"
	       << "# numCalc "<< (cgln_start?"":"time_cgln_0 ")
	       << (cgln_stop? "time_cgln_1":"") << endl;
      for (int numCalcs=1; numCalcs <=2048; numCalcs*=2)
      {
	statfile << numCalcs;
	for (int cgln=cgln_start;cgln<=cgln_stop;cgln++)
	{
	  model.SetModelType(iso,"consistent",cgln);
	  clock_t start = clock();
	  for (int k=0;k<numCalcs;k++)
	  {
	    model.GetTStrangeCalc(iso).SetChiSquared();
	  }
	  clock_t stop = clock();
	  
	  statfile<< "\t" << ((double) (stop-start) )/ (1.0*CLOCKS_PER_SEC);
	}
	statfile << endl;
      }
      statfile.close();
    }
  }
  else
  {
    cerr << "Usage: "<< argv[0] << "fit_specification_path [isopin [test|stress|full [full_outputfile]] ] \n Default mode: test "<< endl;
  }
  return 0;
}


