//$\gamma p \rightarrow Pi^{+} n$ results

{
  TFile * file = new TFile("/home/sam/code/data/data_iso8.root","RECREATE");
  TDataset::SetDataFolder("/home/sam/code/data/");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","iso.8 t kinematics",8,"t:wlab:qsquared",0,5000,0);
  TKinematics* tk2 = new TKinematics("tk2","iso.8 kinematics",8,"costhkcm:wlab:qsquared",1,500,0);
  double cos, wlab;
  
 //*********************************************************************************************
  d = new TDataset("BONN_1983__Althoff__ZPC_18_199","Diff.cross sections by Bonn (Althoff)");

  d->ImportDataset(8,"101"); //184 datapoints
  d->MakeSelections("costhkcm"); // 2 selections
  
  tk2->SetVarRange(2,470,1400,1000);
  d->SetBranchAddress("costhkcm",&cos);
    for(int i=1; i<=2; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk2.SetVar(1,cos);
    d->SetKinematics(tk2);
  }
  d->Scan("*");
  d->Write();
  delete d;
  
  //*********************************************************************************************
  d = new TDataset("1977__Fujii__NPB_120_395","Tot. cross sections by Fujii");
  d->ImportDataset(8,"102");  
  tk2->SetVarRange(2,250,800,1000);
  d->SetKinematics(tk2);
//   d->SetSelection(0);
//   d->GoTo(0);
//   tk.SetVar(2,wlab);
//   d->SetKinematics(tk2);
  d->Scan("*");
  d->Write();
  delete d;
  

//*********************************************************************************************
  d = new TDataset("CEA_1967__BarYam__PRL_19_40","Diff.cross sections (t) by CEA (Bar-Yam)");

  d->ImportDataset(8,"401"); //5 datapoints
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.4e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("CEA_1967__Dowd__PRL_18_414","Diff.cross sections (t) by CEA (Dowd)");
  
  d->ImportDataset(8,"402"); //17 datapoints
  d->MakeSelections("wlab"); // 3 selections
  tk->SetVarRange(1,0.1e6,-0.6e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=3; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************

  d = new TDataset("CEA_1967__Joseph__PRL_19_1206","Diff.cross sections (t) by CEA (Joseph)");
  
  d->ImportDataset(8,"403"); // 6 datapoints
  d->MakeSelections("wlab"); // 2 selections
  tk->SetVarRange(1,0.1e6,-1.5e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=2; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************

  d = new TDataset("DESY_1968__Heide__PRL_21_248","Diff.cross sections (t) by DESY");
  
  d->ImportDataset(8,"404"); //18 datapoints
  d->MakeSelections("wlab"); // 3 selections
  tk->SetVarRange(1,0.1e6,-0.61e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=3; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("SLAC_1968__Boyarski__PRL_20_300", "Diff.cross sections (t) by SLAC");

  d->ImportDataset(8,"405"); //62 datapoints
  d->MakeSelections("wlab"); // 4 selections
  tk->SetVarRange(1,0.1e6,-2.13e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=4; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************

  d = new TDataset("CEA_1970__Sherden__PRL_19_40","Diff.cross sections (t) by CEA (Sherden)");

  d->ImportDataset(8,"406"); //9 datapoints
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.5e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
    
//*********************************************************************************************

  d = new TDataset("CEA_1970__BarYam__PRL_25_1053","Photon asymmetry (t) by CEA");

  d->ImportDataset(8,"407"); //6 datapoints
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.2e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
    
//*********************************************************************************************

  d = new TDataset("SLAC_1973__Burfeindt__PRL_33_509","Photon asymmetry (t) by Burfeindt");

  d->ImportDataset(8,"408"); //3 datapoints
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.0e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
  
//*********************************************************************************************

  d = new TDataset("SLAC_1973__Sherden__PRL_30_24","Photon asymmetry (t) by SLAC");

  d->ImportDataset(8,"409"); //12 datapoints
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.5e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************
  	
	
  d = new TDataset("SLAC_1970__Morehouse__PRL_25_835","Target asymmetry (t) by SLAC");
	
  d->ImportDataset(8,"410"); // 11 datapoints
  d->MakeSelections("wlab"); // 2 selections
  tk->SetVarRange(1,0.1e6,-1.1e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=2; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;

  //*********************************************************************************************
	
	
  d = new TDataset("DESY_1975__Genzel__NPB_9_196","Target asymmetry (t) by DESY");
	
  d->ImportDataset(8,"411"); // 27 datapoints
  d->MakeSelections("wlab"); // 3 selections
  tk->SetVarRange(1,0.1e6,-1.3e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=3; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
	
  //*********************************************************************************************
	
  d = new TDataset("SLAC_1976__Anderson__PRD_14_679","Diff. cross section (t) by SLAC");
	
  d->ImportDataset(8,"412"); // 43 datapoints
  d->MakeSelections("wlab"); // 3 selections
  tk->SetVarRange(1,0.1e6,-12e6,1000);
  d->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=3; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(2,wlab);
    d->SetKinematics(tk);
  }
  d->Scan("*");
  d->Write();
  delete d;
	
//*********************************************************************************************
  
 file->Close();
 delete tk;
  
}
