/*!
 * \file Observable_test.cpp
 *
 * Test program to verify all observables implemented in strangecalc.
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

#include <vector>
#include <TString.h>
#include <TStrangeModel.h>
#include <TKinematics.h>
#include <TCalcInfo.h>
#include <TDataset.h>
#include <TFile.h>
#include <strange_func.h>
#include <iostream>
#include <cassert>
using namespace std;


int main(int argc, char* argv[])
{
  // Set isospin channel
  const int isospin = 2;
  const int kaoncap_isospin = 11;

  // Create new models
  TDataset::SetDataFolder("../../share/data/");
  TStrangeModel model("model","testing");
  TStrangeModel kaoncap_model("kaoncap_model","testing");

  if( getisospinbasechannel(isospin) == 1 ) {
    model.SetStrangeModel("../setup_files/fit_specification_test_iso1",
			  TStrangeModel::kSilent);
  } else if( getisospinbasechannel(isospin) == 2 ) {
    model.SetStrangeModel("../setup_files/fit_specification_test_iso2+3",
			  TStrangeModel::kSilent);
  } else {
    cerr << "ERROR: unexpected isospin channel" << endl;
    assert( 1==0 );
  }

  if( getisospinbasechannel(kaoncap_isospin) == 1 ) {
    kaoncap_model.SetStrangeModel("../setup_files/fit_specification_test_iso11",
				  TStrangeModel::kSilent);
  } else {
    cerr << "ERROR: unexpected isospin channel" << endl;
    assert( 1==0 );
  }

  // Prepare kinematics
  TKinematics tk("tk","tk",isospin,"costhkcm:wlab",1.,1500.);
  tk.SetVarRange(1,-1.,1.,3);
  tk.SetVarRange(2,1500.,7500.,3);

  // Check all photoproduction observables
  TCalcInfo photo_obs(TCalcInfo::kPhoto,isospin,"");
  vector<TString> photo_obs_list;
  photo_obs_list.push_back("dcs");
  photo_obs_list.push_back("tcs");
  photo_obs_list.push_back("rec");
  photo_obs_list.push_back("pho");
  photo_obs_list.push_back("tar");
  photo_obs_list.push_back("c_x");
  photo_obs_list.push_back("c_xp");
  photo_obs_list.push_back("c_z");
  photo_obs_list.push_back("c_zp");
  photo_obs_list.push_back("o_x");
  photo_obs_list.push_back("o_xp");
  photo_obs_list.push_back("o_z");
  photo_obs_list.push_back("o_zp");
  photo_obs_list.push_back("b1^2");
  photo_obs_list.push_back("b2^2");
  photo_obs_list.push_back("b3^2");
  photo_obs_list.push_back("b4^2");
  photo_obs_list.push_back("b1_R");
  photo_obs_list.push_back("b2_R");
  photo_obs_list.push_back("b3_R");
  photo_obs_list.push_back("b4_R");
  photo_obs_list.push_back("b1_I");
  photo_obs_list.push_back("b2_I");
  photo_obs_list.push_back("b3_I");
  photo_obs_list.push_back("b4_I");
  photo_obs_list.push_back("H1^2");
  photo_obs_list.push_back("H2^2");
  photo_obs_list.push_back("H3^2");
  photo_obs_list.push_back("H4^2");
  photo_obs_list.push_back("H1_R");
  photo_obs_list.push_back("H2_R");
  photo_obs_list.push_back("H3_R");
  photo_obs_list.push_back("H4_R");
  photo_obs_list.push_back("H1_I");
  photo_obs_list.push_back("H2_I");
  photo_obs_list.push_back("H3_I");
  photo_obs_list.push_back("H4_I");
  photo_obs_list.push_back("diffcs");
  photo_obs_list.push_back("totcs");
  photo_obs_list.push_back("S");
  photo_obs_list.push_back("T");
  photo_obs_list.push_back("P");
  photo_obs_list.push_back("C_x");
  photo_obs_list.push_back("C_xp");
  photo_obs_list.push_back("C_z");
  photo_obs_list.push_back("C_zp");
  photo_obs_list.push_back("O_x");
  photo_obs_list.push_back("O_xp");
  photo_obs_list.push_back("O_z");
  photo_obs_list.push_back("O_zp");
  photo_obs_list.push_back("E");
  photo_obs_list.push_back("F");
  photo_obs_list.push_back("G");
  photo_obs_list.push_back("H");
  photo_obs_list.push_back("T_x");
  photo_obs_list.push_back("T_xp");
  photo_obs_list.push_back("T_z");
  photo_obs_list.push_back("T_zp");
  photo_obs_list.push_back("L_x");
  photo_obs_list.push_back("L_xp");
  photo_obs_list.push_back("L_z");
  photo_obs_list.push_back("L_zp");

  for(int i=0; i<photo_obs_list.size(); ++i) {

    photo_obs.SetObservable( photo_obs_list[i] );
    cout << photo_obs_list[i] << endl;

    tk.GoTo(0);
    do {
      double calcpoint = model.GetCalcpoint(tk,&photo_obs);
      if( fabs(calcpoint) > 1e-8 ) 
	cout << calcpoint << endl;
    } while( tk.Next() );

    cout << endl;
  }

  // Check all electroproduction observables
  TCalcInfo elec_obs(TCalcInfo::kElectro,isospin,"");
  elec_obs.SetElectroBeamEnergy(5000.);
  vector<TString> elec_obs_list;
  elec_obs_list.push_back("diff_phi");
  elec_obs_list.push_back("diff_l");
  elec_obs_list.push_back("diff_t");
  elec_obs_list.push_back("diff_t+l");
  elec_obs_list.push_back("diff_r_lt");
  elec_obs_list.push_back("diff_tt_unpol");
  elec_obs_list.push_back("diff_tl_unpol");
  elec_obs_list.push_back("diff_tl_epol");
  elec_obs_list.push_back("diff_tt_epol");
  elec_obs_list.push_back("a_lu");
  elec_obs_list.push_back("induced_pol_y");
  elec_obs_list.push_back("induced_pol_yp");
  elec_obs_list.push_back("induced_pol_n");
  elec_obs_list.push_back("induced_pol_yh");
  elec_obs_list.push_back("transf_pol_x");
  elec_obs_list.push_back("transf_pol_xp");
  elec_obs_list.push_back("transf_pol_t");
  elec_obs_list.push_back("transf_pol_xh");
  elec_obs_list.push_back("transf_pol_z");
  elec_obs_list.push_back("transf_pol_zp");
  elec_obs_list.push_back("transf_pol_l");
  elec_obs_list.push_back("transf_pol_zh");

  tk.FixVariables();
  tk.SetFormat("costhkcm:wlab:qsquared:phi");
  tk.SetVarRange(1,-1.,1.,3);
  tk.SetVar(2,2000.);
  tk.SetVarRange(3,1.e5,1.e6,3);
  tk.SetVar(4,1.);
  tk.SetPhiLimits(1.,2.);

  for(int i=0; i<elec_obs_list.size(); ++i) {

    elec_obs.SetObservable( elec_obs_list[i] );
    cout << elec_obs_list[i] << endl;

    tk.GoTo(0);
    do {
      double calcpoint = model.GetCalcpoint(tk,&elec_obs);
      if( fabs(calcpoint) > 1e-8 ) 
	cout << calcpoint << endl;
    } while( tk.Next() );

    cout << endl;
  }

  // Check all kaon capture observables
  TCalcInfo kaoncap_obs(TCalcInfo::kCapture,kaoncap_isospin,"");
  vector<TString> kaoncap_obs_list;
  kaoncap_obs_list.push_back("dcs");
  kaoncap_obs_list.push_back("tcs");
  kaoncap_obs_list.push_back("ptcs");
  kaoncap_obs_list.push_back("bran");

  TKinematics kaoncap_tk("tk","tk",kaoncap_isospin,"costhkcm:pklab",1.,1500.);
  kaoncap_tk.SetVarRange(1,-1.,1.,3);
  kaoncap_tk.SetVarRange(2,50.,250.,3);

  for(int i=0; i<kaoncap_obs_list.size(); ++i) {

    kaoncap_obs.SetObservable( kaoncap_obs_list[i] );
    cout << kaoncap_obs_list[i] << endl;

    kaoncap_tk.GoTo(0);
    do {
      double calcpoint = kaoncap_model.GetCalcpoint(kaoncap_tk,&kaoncap_obs);
      if( fabs(calcpoint) > 1e-8 ) 
	cout << calcpoint << endl;
    } while( kaoncap_tk.Next() );

    cout << endl;
  }
  
  return 0;
}
