/*!
 * \file makerootfile.cpp
 *
 * Use this file to generate "test_strangecalc-wrapper.root"
 * with the command [] strangeViewer -q makerootfile.cpp
 *
 * "test_strangecalc-wrapper.root" is used by the other test apps.
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

{
  TFile *file = new TFile("test_strangecalc-wrapper.root","recreate");

  TStrangeModel model_iso1("model_iso1","testing TStrangeModel");
  TStrangeModel model_iso2("model_iso2","testing TStrangeModel");
  TStrangeModel model_iso11("model_iso11","testing TStrangeModel");
  
  model_iso1.SetStrangeModel("fit_specification_test_iso1");
  model_iso2.SetStrangeModel("fit_specification_test_iso2+3");
  model_iso11.SetStrangeModel("fit_specification_test_iso11");

  model_iso1.Write();
  model_iso2.Write();
  model_iso11.Write();

 file->Close();
 delete file;
  
}
