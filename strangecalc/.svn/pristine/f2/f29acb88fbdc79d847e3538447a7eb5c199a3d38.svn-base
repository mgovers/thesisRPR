/*!
 * \file CreateTStrangeModel.cpp
 *
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

#include <TStrangeModel.h>
#include <TKinematics.h>
#include <TCalcInfo.h>
#include <TFile.h>
#include <iostream>
using namespace std;

int main()
{
  TFile *file = new TFile("test_strangecalc-wrapper_offshell_noCGLN.root","RECREATE");

  // Create new models
  TStrangeModel model_iso1("model_iso1","testing TStrangeModel");
  TStrangeModel model_iso2("model_iso2","testing TStrangeModel");
  TStrangeModel model_iso11("model_iso11","testing TStrangeModel");
  
  model_iso1.SetStrangeModel("../setup_files/fit_specification_test_iso1");
  model_iso2.SetStrangeModel("../setup_files/fit_specification_test_iso2+3");
  model_iso11.SetStrangeModel("../setup_files/fit_specification_test_iso11");
  
  model_iso1.Write();
  model_iso2.Write();
  model_iso11.Write();
  
  file->Close();
  delete file;
  
  //**************************************************************************/  

  file = new TFile("test_strangecalc-wrapper_consistent_noCGLN.root","RECREATE");
  
  model_iso1.SetModelType(1,"consistent",false);
  model_iso2.SetModelType(2,"consistent",false);
  model_iso11.SetModelType(11,"consistent",false);
  
  model_iso1.Write();
  model_iso2.Write();
  model_iso11.Write();
  
  file->Close();
  delete file;  

  //**************************************************************************/
  
  file = new TFile("test_strangecalc-wrapper_consistent_CGLN.root","RECREATE");
  
  model_iso1.SetModelType(1,"consistent",true);
  model_iso2.SetModelType(2,"consistent",true);
  model_iso11.SetModelType(11,"consistent",true);
  
  model_iso1.Write();
  model_iso2.Write();
  model_iso11.Write();
  
  file->Close();
  delete file;  
  
  //**************************************************************************/
  return 0;
}
