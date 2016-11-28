/*!
 * \file pionHighE_test.cpp
 *
 * Very low-level test for the charged pion electroproduction observables
 * above the resonance region.
 *
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

#include <TStrangeModel.h>
#include <TDataset.h>
#include <FormFactorParametrization.h>
#include <Lagrangian.h>
#include <complex>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
  // Initialize TStrangeModel (and don't write to standard out)
  TDataset::SetDataFolder("../../share/data/");
  TStrangeModel model("","");
  model.SetStrangeModel("../setup_files/fit_specification_test_iso8", TStrangeModel::kSilent);

  // Initialize TKinematics and TCalcInfo
  TKinematics kin("","",8,"w:t:qsquared:phi",0.,0.,0.,0.);
  TCalcInfo cs(TCalcInfo::kElectro,8,"a_lu");
  cs.SetElectroConvention(1);
  cs.SetDsDt();

  int NSteps = 6;
  double 
    mtMin,
    mtMax,
    phiMin,
    phiMax,
    beta,
    betaMin,
    betaMax,
    xi,
    xiMin,
    xiMax,
    Q2,
    Q2Min,
    Q2Max,
    W,
    WMin,
    WMax,
    pi = 4.*atan(1.);
  complex<double> F;

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_l, diff_t | JLab-DESY (W = 2.19 GeV, Q^2 = 0.7 GeV^2)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,2.19e3);
  kin.SetVar(3,0.7e6);

  mtMin = -kin.GetTMin();
  mtMax = 0.4e6;
  
  cout << setprecision(2) << fixed
       << "======================================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "======================================" << endl
       << "       t         diff_l         diff_t" << endl
       << " (GeV^2)     (ub/GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------     ----------" << endl;

  cout << setprecision(5) << fixed;

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15);
      cs.SetObservable("diff_l");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_t");
      cout << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_l, diff_t, diff_tt, diff_lt | Fpi-2 (W = 2.274 GeV, Q^2 = 1.416 GeV^2)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,2.274e3);
  kin.SetVar(3,1.416e6);

  mtMin = -kin.GetTMin();
  mtMax = 0.25e6;
  
  cout << setprecision(3) << fixed
       << "====================================================================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "====================================================================" << endl
       << "       t         diff_l         diff_t        diff_tt        diff_lt" << endl
       << " (GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------     ----------     ----------     ----------" << endl;

  cout << setprecision(5) << fixed;

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15);
      cs.SetObservable("diff_l");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_t");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tt_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tl_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_l, diff_t, diff_tt, diff_lt | pi-CT (W = 2.22 GeV, Q^2 = 3.91 GeV^2)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,2.22e3);
  kin.SetVar(3,3.91e6);

  mtMin = -kin.GetTMin();
  mtMax = 0.84e6;
  
  cout << setprecision(2) << fixed
       << "====================================================================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "====================================================================" << endl
       << "       t         diff_l         diff_t        diff_tt        diff_lt" << endl
       << " (GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------     ----------     ----------     ----------" << endl;

  cout << setprecision(5) << fixed;

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15);
      cs.SetObservable("diff_l");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_t");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tt_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tl_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_l, diff_t, diff_tt, diff_lt | Fpi-1 (W = 1.983 GeV, Q^2 = 0.526 GeV^2)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,1.983e3);
  kin.SetVar(3,0.526e6);

  mtMin = -kin.GetTMin();
  mtMax = 0.09e6;
  
  cout << setprecision(3) << fixed
       << "====================================================================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "====================================================================" << endl
       << "       t         diff_l         diff_t        diff_tt        diff_lt" << endl
       << " (GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------     ----------     ----------     ----------" << endl;

  cout << setprecision(5) << fixed;

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15);
      cs.SetObservable("diff_l");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_t");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tt_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tl_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_phi | Fpi-1 (W = 1.942 GeV, Q^2 = 0.612 GeV^2, t = -0.05 GeV^2, epsilon = 0.74)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,1.942e3);
  kin.SetVar(2,-0.05e6);
  kin.SetVar(3,0.612e6);
  cs.SetElectroEpsilon(0.74);

  phiMin = 0.;
  phiMax = pi;;
  
  cout << setprecision(3) << fixed
       << "=======================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "t = " << kin.GetVar(2)*1e-6 << " GeV^2" << endl
       << "epsilon = " << cs.GetElectroEpsilon() << endl
       << "=======================" << endl
       << "    phi       diff_phi" << endl
       << "  (rad)     (ub/GeV^2)" << endl
       << "  -----     ----------" << endl;

  cout << setprecision(5) << fixed;
  cs.SetObservable("diff_phi");

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(4,phiMin + (phiMax - phiMin)*((double)i)/((double)(NSteps - 1)));
      cout << kin.GetVar(4)
	   << setw(15)
	   << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_U | pi-CT (W = 2.25 GeV, Q^2 = 3.9 GeV^2, epsilon = 0.39)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,2.25e3);
  kin.SetVar(3,3.9e6);
  cs.SetElectroEpsilon(0.39);

  mtMin = -kin.GetTMin();
  mtMax = 0.75e6;
  
  cout << setprecision(2) << fixed
       << "=======================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "epsilon = " << cs.GetElectroEpsilon() << endl
       << "=======================" << endl
       << "       t         diff_U" << endl
       << " (GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------" << endl;

  cout << setprecision(5) << fixed;
  cs.SetObservable("diff_t+l");

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15)
	   << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_U (neutron) | DESY (W = 2.19 GeV, Q^2 = 1.35 GeV^2, epsilon = 0.84)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetIsospin(10);
  cs.SetIsospin(10);
  kin.SetVar(1,2.19e3);
  kin.SetVar(3,1.35e6);
  cs.SetElectroEpsilon(0.84);

  mtMin = -kin.GetTMin();
  mtMax = 0.55e6;
  
  cout << setprecision(2) << fixed
       << "=====================================================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2      (neutron)" << endl
       << "epsilon = " << cs.GetElectroEpsilon() << endl
       << "=====================================================" << endl
       << "       t         diff_U        diff_tt        diff_lt" << endl
       << " (GeV^2)     (ub/GeV^2)     (ub/GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------     ----------     ----------" << endl;

  cout << setprecision(5) << fixed;

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15);
      cs.SetObservable("diff_t+l");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tt_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << setw(15);
      cs.SetObservable("diff_tl_unpol");
      cout << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  kin.SetIsospin(8);
  cs.SetIsospin(8);

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * diff_U | HERMES (x_B = 0.35, Q^2 = 5.0 GeV^2, E = 27.6 GeV)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetFormat("xb:t:qsquared:phi");
  kin.SetVar(1,0.35);
  kin.SetVar(3,5.0e6);
  cs.SetElectroBeamEnergy(27.6e3);

  mtMin = -kin.GetTMin();
  mtMax = 3.e6;
  
  cout << setprecision(2) << fixed
       << "=======================" << endl
       << "x_B = " << kin.GetVar(1) << endl
       << setprecision(1) << fixed
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << "E = " << cs.GetElectroBeamEnergy()*1e-3 << " GeV" << endl
       << "=======================" << endl
       << "       t         diff_U" << endl
       << " (GeV^2)     (ub/GeV^2)" << endl
       << " -------     ----------" << endl;

  cout << setprecision(5) << fixed;
  cs.SetObservable("diff_t+l");

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15)
	   << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  kin.SetFormat("w:t:qsquared:phi");

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * a_lu | JLab (W = 2.0 GeV, Q^2 = 1.5 GeV^2, E = 5.77 GeV)
   *---------------------------------------------------------------------------------------------
   */
  kin.SetVar(1,2.0e3);
  kin.SetVar(3,1.5e6);
  cs.SetElectroBeamEnergy(5.77e3);
  cs.SetDsDomega();

  mtMin = -kin.GetTMin();
  mtMax = 2.e6;
  
  cout << setprecision(1) << fixed
       << "=======================" << endl
       << "W = " << kin.GetVar(1)*1e-3 << " GeV" << endl
       << "Q^2 = " << kin.GetVar(3)*1e-6 << " GeV^2" << endl
       << setprecision(2) << fixed
       << "E = " << cs.GetElectroBeamEnergy()*1e-3 << " GeV" << endl
       << "=======================" << endl
       << "       t           a_lu" << endl
       << " (GeV^2)               " << endl
       << " -------           ----" << endl;

  cout << setprecision(5) << fixed;
  cs.SetObservable("a_lu");

  for(int i = 1; i < NSteps; i++)
    {
      kin.SetVar(2,-(mtMin + (mtMax - mtMin)*((double)i)/((double)(NSteps - 1))));
      cout << kin.GetVar(2)*1e-6
	   << setw(15)
	   << model.GetCalcpoint(kin,&cs)
	   << endl;
    }

  cs.SetDsDt();

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * F_s | (W = 2.2 GeV, Q^2 = 1.3 GeV^2)
   *---------------------------------------------------------------------------------------------
   */
  W = 2.2e3;
  Q2 = 1.3e6;

  betaMin = 1.;
  betaMax = 4.;
  xiMin = 0.1;
  xiMax = 0.9;
  
  cout << setprecision(1) << fixed
       << "================================================" << endl
       << "W = " << W*1e-3 << " GeV" << endl
       << "Q^2 = " << Q2*1e-6 << " GeV^2" << endl
       << "================================================" << endl
       << "beta             xi                          F_s" << endl
       << "----             --                          ---" << endl;

  for(int i = 0; i < 2*NSteps - 2; i++)
    {
      beta = betaMin + (betaMax - betaMin)*((double)i)/((double)(2*NSteps - 3));
      for(int j = 0; j < NSteps - 1; j++)
	{
	  xi = xiMin + (xiMax - xiMin)*((double)j)/((double)(NSteps - 2));
	  cout << setprecision(2) << fixed
	       << beta
	       << setw(15)
	       << setprecision(1) << fixed
	       << xi
	       << setw(29)
	       << setprecision(5) << fixed;
	  F = nucleonEMFF(Q2,W*W,beta,xi,1);
	  cout << F
	       << endl;
	}
    }

  cout << endl << endl;

  /*
   *---------------------------------------------------------------------------------------------
   * F_{g pi pi}
   *---------------------------------------------------------------------------------------------
   */
  WMin = 1.9e3;
  WMax = 2.2e3;
  Q2Min = 0.35e6;
  Q2Max = 2.0e6;
  
  cout << "=======================================" << endl
       << endl
       << "=======================================" << endl
       << "      W            Q^2      F_{g pi pi}" << endl
       << "(GeV^2)        (GeV^2)" << endl
       << "-------        -------      -----------" << endl;
       
  cout << setprecision(5) << fixed;

  for(int i = 0; i < NSteps - 1; i++)
    {
      W = WMin + (WMax - WMin)*((double)i)/((double)(NSteps - 2));
      for(int j = 0; j < NSteps - 1; j++)
	{
	  Q2 = Q2Min + (Q2Max - Q2Min)*((double)j)/((double)(NSteps - 2));
	  cout << setprecision(3) << fixed
	       << "  " << W*1e-3
	       << setw(15)
	       << setprecision(4) << fixed
	       << Q2*1e-6
	       << setw(17)
	       << setprecision(5) << fixed
	       << pionFF(Q2,0,0,W)
	       << endl;
	}
    }

  return 0;
}
