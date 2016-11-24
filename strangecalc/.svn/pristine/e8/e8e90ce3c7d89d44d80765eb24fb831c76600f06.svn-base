//TCalcInfo is actually a wrapper for Data struct (see fitting.h). A datapoint knows its TCalcInfo, but not its TKinematics.
//TKinematics determines the model graph.
//Run this script with strangeViewer and enter the following datasets: 200 201 202 203 204 205 207 208 209.
{
  TFile * storeFile = new TFile("electrodata_iso1.root","RECREATE");
  TDataset::SetDataFolder("./");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","tk",1,"w:costhkcm:qsquared",2.16,1.0,290000.0);
  int offset = 0;
  double costhkcm,qsquared,w;
//*********************************************************************************************
    
  d = new TDataset("CEA_1972__Brown__PRL_28_1086","T+L cross sec. by CEA"); //201
  // 8 datapoints
  d->ImportDataset(1,"201");
  d->Scan("*");
  
  //pick data with qsquared = 290000, variable W
  d->AddSelection("qsquared == 290000");
    //pick data with <W> = 2.16 (!= 2.15 Tamara) and costhkcm = 0, variable Q^2
  d->AddSelection("w >= 2050 && w <= 2220");
  
  d->SetSelection(1);
  tk->SetVarRange(1,1500.0,3000.0,50);
  d->SetKinematics(tk);
  

  d->SetSelection(2);
  tk->FixVariable(1);
  tk->SetVar(1,2160);
  tk->SetVarRange(3,1,4500000.0,50);
  d->SetKinematics(tk);
  
  d->Write();
  
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("Cornell_1974__Bebek__PRL_32_21","T+L cross sec. by Harvard & Cornell");//202
  // 4 datapoints
  d->ImportDataset(1,"202"); 
  d->Scan("*");
  
  // 3 datapoints with W approx 2663 (take 2.7)
  d->AddSelection("W >=2650 && W <=2700");//1
  // 1 datapoint with W approx 2200 MeV, matches previous CEA set.
  d->AddSelection("w == 2200");//2
  
  d->SetSelection(1);
  tk->SetVar(1,2663);
  tk->SetVarRange(3,1,4500000.0,50);
  d->SetKinematics(tk);
  
  d->SetSelection(2);
  tk->SetVar(1,2160);
  tk->SetVarRange(3,1,4500000.0,50);
  d->SetKinematics(tk);
  d->Write();
  
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("Cornell_1977__Bebek__PRD_15_594","T+L cross sec. by Harvard & Cornell");//203
  // 14 datapoints
  d->ImportDataset(1,"203"); 
  
  // 7 datapoints with an average W of 2.16 GeV (actually (2 170 + 2 150 + 2 200 + 2 170 + 2 100 + 2 210 + 2 220) / 7 = 2 174.28571)
  d->AddSelection("w <= 2220"); //1
  // now choose data with a Q² that is compatible with CLAS values: 
  // approx 1 GeV²
  d->AddSelection("qsquared >= 750000 && qsquared < 1275000" );//2
  // approx 1.55 GeV²
  d->AddSelection("qsquared >= 1275000 && qsquared < 1800000" );//3
  // approx 2.05 GeV²
  d->AddSelection("qsquared >= 1800000 && qsquared < 2300000" );//4
      
  // approx 2.16 GeV  (actually 2.17)
  d->SetSelection(1);
  d->Scan("*");
  tk->SetVar(1,2160);
  tk->SetVarRange(3,1,4500000.0,50);
  d->SetKinematics(tk);
  
  
  //think about these 2: they are _very_ far off... approx 1.2 GeV²
  d->SetSelection(2);
  d->Scan("*");
  tk->SetVarRange(1,1500,3000,50);
  tk->FixVariable(3);
  tk->SetVar(3,1000000.0);
  d->SetKinematics(tk);
  
  //think about these 2: they are _very_ far off... approx 1.3 GeV²
  d->SetSelection(3);
  d->Scan("*");
  tk->SetVar(3,1550000.0);
  d->SetKinematics(tk);
  
  //acceptable! approx 2 GeV²
  d->SetSelection(4);
  d->Scan("*");
  tk->SetVar(3,2050000.0);
  d->SetKinematics(tk);  

  d->Write();
  delete d;

//*********************************************************************************************

  d = new TDataset("JLAB_2003__Mohring__PRC_67_055205","T&L cross sec. by JLAB");//204
  //Q^2, t? ((W = 1.83 (!= 1.84 Tamara) and thkcm=8deg or costhkcm approx 1)
  d->ImportDataset(1,"204"); 
  
  d->MakeSelections("observable");
  // approx 1 GeV²
  d->AddSelection("observable==\"diff_t\" && qsquared == 1000000");
  d->AddSelection("observable==\"diff_l\" && qsquared == 1000000");
  // approx 2.05 GeV²
  d->AddSelection("observable==\"diff_l\" && qsquared == 2000000");
  d->AddSelection("observable==\"diff_t\" && qsquared == 2000000");
  
  
  d->SetSelection(1);
  d->Scan("*");
  tk->FixVariable(1);
  tk->SetVar(1,1830);
  tk->SetVarRange(3,1,4500000.0,50);
  d->SetKinematics(tk);
  
  d->SetSelection(2);
  d->Scan("*");
  d->SetKinematics(tk);
  
  
  //not immediately useful...
  d->SetSelection(3);
  d->Scan("*");
  tk->FixVariable(3);
  tk->SetVar(3,1000000.0);
  tk->SetVarRange(1,1500,3000,50);
  d->SetKinematics(tk);
  
  d->SetSelection(4);
  d->Scan("*");
  d->SetKinematics(tk);

  d->SetSelection(5);
  d->Scan("*");
  tk->FixVariable(3);
  tk->SetVar(3,2000000.0);
  tk->SetVarRange(1,1500,3000,50);
  d->SetKinematics(tk);
  
  d->SetSelection(6);
  d->Scan("*");
  d->SetKinematics(tk);
  
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("CLAS_2007__Ambrozewicz__PRC_75_045203_2.567","T+L, TT, LT cross sec. by CLAS (2.567 GeV)");//205
  //Q^2, W, costhkcm
  // Q² is either 0.65 or 1 GeV²!
  d->ImportDataset(1,"205"); 
  d->MakeSelections("observable:qsquared:w");
  d->MakeSelections("observable:qsquared:costhkcm");
  d->MakeSelections("observable:costhkcm:w");
  
  tk->FixVariable(1);
  tk->FixVariable(3);
  tk->SetVarRange(2,-1.0,1.0,50);
  
  double wvalues_205[8] = {1650.0, 1725.0, 1775.0, 1825.0, 1875.0, 1925.0, 1975.0, 2050.0};
  double cvalues_205[6] = {-0.6, -0.25, 0.05, 0.35, 0.65, 0.9};
  double qvalues_205[2] = {.65E6, 1E6};
    
  for(int i=0; i<8;++i)
  {
    for(int obs=0;obs<3;++obs)
    {
      for(int qq=0;qq<2;++qq)
      {
	if(!(wvalues_205[i] == 2050.0 && qq==1))
	{
	  d->SetSelection(6*i+3*qq+obs+1);
	  tk->SetVar(1,wvalues_205[i]);
	  tk->SetVar(3,qvalues_205[qq]);
	  d->SetKinematics(tk);
	}
      }
    }
  }
  offset = d->GetSelection();
  
  tk->FixVariable(2);
  tk->SetVarRange(1,1500.0,3000.0,50);
  
  for(int i=0; i<6;++i)
  {
    for(int obs=0;obs<3;++obs)
    {
      for(int qq=0;qq<2;++qq)
      {
	  d->SetSelection(offset+6*i+3*qq+obs+1);
	  tk->SetVar(2,cvalues_205[i]);
	  tk->SetVar(3,qvalues_205[qq]);
	  d->SetKinematics(tk);
      }
    }
  }
  offset = d->GetSelection();
  
  tk->FixVariable(1);
  tk->SetVarRange(3,1,4500000.0,50);
  for(int i=0; i<8;++i)
  {
    for(int obs=0;obs<3;++obs)
    {
      for(int cos=0;cos<6;++cos)
      {
	  d->SetSelection(offset+18*i+3*cos+obs+1);
	  tk->SetVar(1,wvalues_205[i]);
	  tk->SetVar(2,cvalues_205[cos]);
	  d->SetKinematics(tk);
      }
    }
  }
 
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("CLAS_2007__Ambrozewicz__PRC_75_045203_4.056","T+L, TT, LT cross sec. by CLAS (4.056 GeV)");//207
  //Q^2, W, costhkcm
  //Q² is 1, 1.55, 2.05 or 2.55 GeV²
  //d->ImportDataset(1); d->MakeSelections("qsquared:wlab");
  //d->MakeSelections("costhkcm:w");
  
  d->ImportDataset(1,"207"); 
  d->MakeSelections("observable:qsquared:w");
  d->MakeSelections("observable:qsquared:costhkcm");
  d->MakeSelections("observable:costhkcm:w");
  
  tk->FixVariable(1);
  tk->FixVariable(3);
  tk->SetVarRange(2,-1.0,1.0,50);
  
  double wvalues_207[8] = {1650.0, 1750.0, 1850.0, 1950.0,2050.0,2150.0,2250.0,2350.0};
  double cvalues_207[6] = {-0.6, -0.25, 0.05, 0.35, 0.65, 0.9};
  double qvalues_207[4] = {1.0E6, 1.55E6, 2.05E6, 2.55E6};
    
  offset=0;
  for(int i=0; i<8;++i)
  {
    for(int qq=0;qq<4;++qq)
    {
      for(int obs=0;obs<3;++obs)
      {
	if(wvalues_207[i] < 2250.0)
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_207[i]);
	  tk->SetVar(3,qvalues_207[qq]);
	  d->SetKinematics(tk);
	}
	//NB no qq=2.55 data for w=2250 and 2350!!
	else if(wvalues_207[i] == 2250.0 && qq !=3)
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_207[i]);
	  tk->SetVar(3,qvalues_207[qq]);
	  d->SetKinematics(tk);
	}
	else if(wvalues_207[i] == 2350.0 && qq!=3)
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_207[i]);
	  tk->SetVar(3,qvalues_207[qq]);
	  d->SetKinematics(tk);
	}
      }
    }
  }
  tk->FixVariable(2);
  tk->SetVarRange(1,1500.0,3000.0,50);
  
  for(int i=0; i<6;++i)
  {
    for(int qq=0;qq<4;++qq)
    {
      for(int obs=0;obs<3;++obs)
      {
	d->SetSelection(++offset);
	tk->SetVar(2,cvalues_207[i]);
	tk->SetVar(3,qvalues_207[qq]);
	d->SetKinematics(tk);
      }
    }
  }
  
  
  //NB no costhkcm=-0.6 data for w=2050 and 2150!!
  tk->FixVariable(1);
  tk->SetVarRange(3,1,4500000.0,50);//FIXME
  for(int i=0; i<8;++i)
  { 
    for(int cos=0;cos<6;++cos)
    {
      for(int obs=0;obs<3;++obs)
      {
	if(!((wvalues_207[i]==2050.0 || wvalues_207[i]==2150.0) && cvalues_207[cos]==-0.6))
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_207[i]);
	  tk->SetVar(2,cvalues_207[cos]);
	  d->SetKinematics(tk);
// 	  cout  << "Selection: " << d->GetSelection() 
// 	      << "\nW: " << wvalues[i]
// 	      << "\nCos: " << cvalues[cos]<< "\n";
	  
	}
      }
    }
  }
  cout  << "Selection: " << d->GetSelection() << "\n";
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("CLAS_2007__Ambrozewicz__PRC_75_045203_sep","T, L, TT, LT cross sec. by CLAS");//208
  //W, costhkcm
  //Q^2 = 1 gev^2
  d->ImportDataset(1,"208"); 
  d->MakeSelections("observable:w");
  d->MakeSelections("observable:costhkcm");
  
  //TKinematics* tk = new TKinematics("tk","tk",1,"w:costhkcm:qsquared",2.16,1.0,290000.0);
  
  tk->FixVariable(1);
  tk->SetVarRange(2,-1.0,1.0,50);
  tk->FixVariable(3);
  tk->SetVar(3,1E6);
  
  double wvalues_208[4] = {1650.0, 1750.0, 1850.0, 1950.0};
  double cvalues_208[6] = {-0.6, -0.25, 0.05, 0.35, 0.65, 0.9};
  
  offset=0;
  for(int i=0; i<4;++i)
  {
    for(int obs=0;obs<4;++obs)
      {
	d->SetSelection(++offset);
	tk->SetVar(1,wvalues_208[i]);
	d->SetKinematics(tk);
      }
  }
  
  tk->FixVariable(2);
  tk->SetVarRange(1,1500.0,3000.0,50);
  
  for(int i=0; i<6;++i)
  {
      for(int obs=0;obs<4;++obs)
      {
	d->SetSelection(++offset);
	tk->SetVar(2,cvalues_208[i]);
	d->SetKinematics(tk);
      }
  }
  
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("CLAS_2008__Nasseripour__PRC_77_065208","LT' cross sec. by CLAS");//209 
  //Q^2, W, costhkcm
  //Q^2 is either 1 or 0.65 GeV²!
  
  d->ImportDataset(1,"209"); 
  d->MakeSelections("qsquared:w");// all combinations of angles are present
  d->MakeSelections("qsquared:costhkcm");
  d->MakeSelections("costhkcm:w");

  tk->FixVariable(1);
  tk->FixVariable(3);
  tk->SetVarRange(2,-1.0,1.0,50);
  
  double wvalues_209[8] = {1650.0, 1725.0, 1775.0, 1825.0, 1875.0, 1925.0, 1975.0, 2025.0};
  double cvalues_209[6] = {-0.6, -0.25, 0.05, 0.35, 0.65, 0.9}; 
  double qvalues_209[2] = {.65E6, 1E6};
    
  offset = 0;
  
  for(int i=0; i<8;++i)
  {
      for(int qq=0;qq<2;++qq)
      {
	if(!(wvalues_209[i] == 2025.0 && qvalues_209[qq]==1e6))
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_209[i]);
	  tk->SetVar(3,qvalues_209[qq]);
	  d->SetKinematics(tk);
	}
      }
  }
  
  tk->FixVariable(2);
  tk->SetVarRange(1,1500.0,2050.0,50);
  
  for(int cos=0; cos<6;++cos)
  {
      for(int qq=0;qq<2;++qq)
      {
	d->SetSelection(++offset);
	tk->SetVar(2,cvalues_209[cos]);
	tk->SetVar(3,qvalues_209[qq]);
	d->SetKinematics(tk);
      }
  }
  
  tk->FixVariable(1);
  tk->SetVarRange(3,1,4500000.0,50);//FIXME
  for(int i=0; i<8;++i)
    {
      for(int cos=0;cos<6;++cos)
	{
	  d->SetSelection(++offset);
	  tk->SetVar(1,wvalues_209[i]);
	  tk->SetVar(2,cvalues_209[cos]);
	  d->SetKinematics(tk);
	}
    }
  
  d->Write();
  delete d;
  
  //*********************************************************************************************
  
  d = new TDataset("CLAS_2009__Coman__arXiv_0911-3943__sigmaT","T cross sec. by CLAS");
  // 8 datapoints
  d->ImportDataset(1,"211");
  
  d->MakeSelections("qsquared");
  
  d->SetSelection(1);
  tk->FixVariable(3);
  tk->SetVar(2,1.);
  tk->SetVar(3,1.9e+06);
  tk->SetVarRange(1,1600.,2500.,100);
  d->SetKinematics(tk);
  
  d->SetSelection(2);
  tk->SetVar(3,2.35e+06);
  d->SetKinematics(tk);
  
  d->Write();
  delete d;

  //*********************************************************************************************
  
  d = new TDataset("CLAS_2009__Coman__arXiv_0911-3943__sigmaL","L cross sec. by CLAS");
  // 8 datapoints
  d->ImportDataset(1,"212");
  
  d->MakeSelections("qsquared");
  
  d->SetSelection(1);
  tk->SetVar(2,1.);
  tk->SetVar(3,1.9e+06);
  tk->SetVarRange(1,1600.,2500.,100);
  d->SetKinematics(tk);
  
  d->SetSelection(2);
  tk->SetVar(3,2.35e+06);
  d->SetKinematics(tk);
  
  d->Write();
  delete d;

  //*********************************************************************************************
  
  d = new TDataset("CLAS_2009__Carman__PRC_79_065205__transfPol_ebeam4GeV","Transferred polarizations by CLAS");
  d->ImportDataset(1,"213");

  // Make selections for plots as a function of costhkcm
  //----------------------------------------------------
  d->MakeSelections("qsquared:w:observable");
  
  // we need to remove the surplus selections
  for(int i=26; i>0; --i) { 
    d->SetSelection(i); 
    if(d->GetNrOfEntries()<2) d->RemoveSelection(i);
  }
  
  // set the kinematics
  tk->FixVariables();
  tk->SetVarRange(2,-1.,1.,50);
  
  d->SetBranchAddress("w",&w);
  d->SetBranchAddress("qsquared",&qsquared);
  
  for(int i=1; i<=6; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk->SetVar(1,w);
    tk->SetVar(3,qsquared);
    d->SetKinematics(tk);
  }  
  d->ResetBranchAddresses();

  // Make selections for plots as a function of W
  //---------------------------------------------
  d->SetSelection(0);
  d->MakeSelections("observable");
  d->SetSelection(7);
  d->AddSelection("costhkcm>-4.10+3.92e-3*w-0.84e-6*w*w-.01 && costhkcm<-4.10+3.92e-3*w-0.84e-6*w*w+.01");
  d->SetSelection(9);
  d->SetDescription("observable==\"transf_pol_z\" && costhkcm==-4.10+3.92w-0.84w^2");
  d->SetSelection(8);
  d->AddSelection("costhkcm>-4.10+3.92e-3*w-0.84e-6*w*w-.01 && costhkcm<-4.10+3.92e-3*w-0.84e-6*w*w+.01");
  d->SetSelection(10);
  d->SetDescription("observable==\"transf_pol_x\" && costhkcm==-4.10+3.92w-0.84w^2");
  d->RemoveSelection(7);
  d->RemoveSelection(7);
  
  d->Write();
  delete d;

  //*********************************************************************************************
  
  d = new TDataset("CLAS_2009__Carman__PRC_79_065205__transfPol_ebeam6GeV","Transferred polarizations by CLAS");
  d->ImportDataset(1,"214");

  // Make selections for plots as a function of costhkcm
  //----------------------------------------------------
  d->MakeSelections("qsquared:w:observable");
  
  // we need to remove the surplus selections
  for(int i=86; i>0; --i) { 
    d->SetSelection(i); 
    if(d->GetNrOfEntries()<2) d->RemoveSelection(i);
  }
  
  // set the kinematics
  tk->FixVariables();
  tk->SetVarRange(2,-1.,1.,50);
  
  d->SetBranchAddress("w",&w);
  d->SetBranchAddress("qsquared",&qsquared);
  
  for(int i=1; i<=6; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk->SetVar(1,w);
    tk->SetVar(3,qsquared);
    d->SetKinematics(tk);
  }  
  d->ResetBranchAddresses();
  
  // Make selections for plots as a function of W
  //---------------------------------------------
  d->SetSelection(0);
  d->MakeSelections("observable");
  d->SetSelection(7);
  d->AddSelection("costhkcm>-3.74+3.48e-3*w-0.69e-6*w*w-.01 && costhkcm<-3.74+3.48e-3*w-0.69e-6*w*w+.01");
  d->SetSelection(9);
  d->SetDescription("observable==\"transf_pol_z\" && costhkcm==-3.74+3.48w-0.69w^2");
  d->SetSelection(8);
  d->AddSelection("costhkcm>-3.74+3.48e-3*w-0.69e-6*w*w-.01 && costhkcm<-3.74+3.48e-3*w-0.69e-6*w*w+.01");
  d->SetSelection(10);
  d->SetDescription("observable==\"transf_pol_x\" && costhkcm==-3.74+3.48w-0.69w^2");
  d->RemoveSelection(7);
  d->RemoveSelection(7);

  // Make selections for plots as a function of qsquared
  //----------------------------------------------------
  d->SetSelection(0);
  d->MakeSelections("costhkcm:w:observable");

  // we need to remove the surplus selections
  for(int i=122; i>0; --i) { 
    d->SetSelection(i); 
    if(d->GetNrOfEntries()<2) d->RemoveSelection(i);
  }

  // set the kinematics
  tk->FixVariables();
  tk->SetVarRange(3,1.5e6,5.e6,100);
  
  d->SetBranchAddress("w",&w);
  d->SetBranchAddress("costhkcm",&costhkcm);
  
  for(int i=9; i<=10; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk->SetVar(1,w);
    tk->SetVar(2,costhkcm);
    d->SetKinematics(tk);
  }  
  d->ResetBranchAddresses();
  
  d->Write();
  delete d;

  //*********************************************************************************************
  
  d = new TDataset("CLAS_2000__Cha__PhD-unpublished__sigmaL+T","unpolarized cross sec. by CLAS");
  // 11 datapoints
  d->ImportDataset(1,"216");

  tk->FixVariables();
  tk->SetVarRange(2,-1.,1.,50);
  
  tk->SetVar(1,1782.);
  tk->SetVar(3,.51e6);
  d->AddSelection("qsquared>400000",tk);
  
  tk->SetVar(1,1893.);
  tk->SetVar(3,.38e6);
  d->AddSelection("qsquared<400000",tk);

  d->Write();
  delete d;


  //*********************************************************************************************
  delete tk;
  storeFile->Close();
}
