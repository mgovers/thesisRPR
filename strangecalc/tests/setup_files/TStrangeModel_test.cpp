/*!
 * \file TStrangeModel_test.cpp
 *
 * Verify the functionality of the TStrangeModel class.
 *
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
#include <TKinematics.h>
#include <TCalcInfo.h>
#include <TFile.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
  // two ways to call this: either with no arguments, or with 3 arguments:
  // root-binary modelType cgln (0 or 1).
  TFile *file = new TFile(argc>3 ? argv[1] : "test_strangecalc-wrapper.root");

  // Create new models
  TStrangeModel model_iso1("model_iso1","testing TStrangeModel");
  TStrangeModel model_iso2("model_iso2","testing TStrangeModel");
  TStrangeModel model_iso11("model_iso11","testing TStrangeModel");
  
  model_iso1.SetStrangeModel("../setup_files/fit_specification_test_iso1");
  model_iso2.SetStrangeModel("../setup_files/fit_specification_test_iso2+3");
  model_iso11.SetStrangeModel("../setup_files/fit_specification_test_iso11");
  
  if (argc>3)
  {
    cerr << "Setting model types...\n"; 
     model_iso1.SetModelType(1,argv[2],atoi(argv[3]));
     model_iso2.SetModelType(2,argv[2],atoi(argv[3]));
     model_iso11.SetModelType(11,argv[2],atoi(argv[3]));
  }
     
  // Load old models
  TStrangeModel *oldmodel_iso1
    = (TStrangeModel*) file->Get("model_iso1");
  TStrangeModel *oldmodel_iso2
    = (TStrangeModel*) file->Get("model_iso2");
  TStrangeModel *oldmodel_iso11
    = (TStrangeModel*) file->Get("model_iso11");

  TCalcInfo ci(TCalcInfo::kElectro,1,"transf_pol_xp");
  ci.SetElectroBeamEnergy(2567.);
  TKinematics tk("tk","tk",1,"costhkcm:wlab:qsquared",1.,1500.,0.);
  tk.SetVarRange(1,-.8,1.,5);
  tk.SetVarRange(3,.1e6,1.e6,5);

  do {
    double calcpoint = model_iso1.GetCalcpoint(tk,&ci);
    double oldcalcpoint = oldmodel_iso1->GetCalcpoint(tk,&ci);
    if( fabs(calcpoint) > 1e-8 ) 
    cout << calcpoint << " ?= " << oldcalcpoint << endl;
  } while( tk.Next() );

  tk.SetIsospin(5);
  tk.GoTo(0);
  ci.SetIsospin(5);
  ci.SetElectroBeamEnergy(5567.);
  tk.SetVarRange(1,-.8,1.,5);
  tk.SetVarRange(3,.1,.2e6,5);
  do {
    double calcpoint = model_iso2.GetCalcpoint(tk,&ci);
    double oldcalcpoint = oldmodel_iso2->GetCalcpoint(tk,&ci);
    if( fabs(calcpoint) > 1e-8 ) 
    cout << calcpoint << " ?= " << oldcalcpoint << endl;
  } while( tk.Next() );

  TKinematics tk2("tk2","tk2",11,"costhkcm:pklab:qsquared",1.,300.,0.);
  tk2.SetVarRange(1,-1.,1.,5);
  tk2.SetVarRange(2,200.,600.,5);
  TCalcInfo ci2(TCalcInfo::kCapture,11,"dcs");
  do {
    double calcpoint = model_iso11.GetCalcpoint(tk2,&ci2);
    double oldcalcpoint = oldmodel_iso11->GetCalcpoint(tk2,&ci2);
    if( fabs(calcpoint) > 1e-8 ) 
      cout << calcpoint << " ?= " <<  oldcalcpoint << endl;
  } while( tk2.Next() );
  

  delete oldmodel_iso1;
  delete oldmodel_iso2;
  delete oldmodel_iso11;
  file->Close();
  delete file;
  
  return 0;
}
