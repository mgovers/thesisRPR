{
  TFile * file = new TFile("data_iso1_HE.root","RECREATE");
  TDataset::SetDataFolder("./");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","iso.1 HE kinematics",1,"t:wlab:qsquared",0,5000,0);
  double wlab, w,costhkcm;
  
//*********************************************************************************************
  
  d = new TDataset("SLAC_1969__Boyarski__PRL_22_1131", "HE Diff.cross sections by SLAC");

  d->ImportDataset(1,"401"); //56 datapoints
  d->MakeSelections("wlab"); // 4 selections
  tk->SetVarRange(1,0.1e6,-2.1e6,100);
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
  d = new TDataset("SLAC_1979__Quinn__PRD_20_1553","HE Photon asymmetry by SLAC");

  d->ImportDataset(1,"403"); //9 datapoints (NOTE: check this! trailing newlines cause phantom data points)
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.1e6,100);
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  
//*********************************************************************************************
  
  d = new TDataset("DESY_1972__Vogel__PLB_40_513","HE Recoil polarization by DESY");
  d->ImportDataset(1,"404"); //7 datapoints (NOTE: check this! trailing newlines cause phantom data points)
  d->SetSelection(0);
  tk->SetVarRange(1,0.1e6,-1.1e6,100);
  
  d->SetBranchAddress("wlab",&wlab);
  d->GoTo(0);
  tk.SetVar(2,wlab);
  d->SetKinematics(tk);
  d->Scan("*");
  d->Write();
  delete d;
  delete tk;
  
  //*********************************************************************************************
  /*
   * HIGH-ENERGY part of the McCracken results -- also cos(theta) > 3.5
   * NB:No differential cross section results are given for the 1.955, 2.735, or 2.745 GeV bins.
   */

  TKinematics* tk = new TKinematics("tk","iso.1 kinematics",1,"w:costhkcm:qsquared",2600,1,0);
  d = new TDataset("CLAS_2009__McCracken__PRC_81_025201__dcs_HE", "HE Diff.cross sections by CLAS");
  
  d->ImportDataset(1,"405"); // 132 datapoints
  
  
  int nsel= d->MakeSelections("w"); // 22 selections
  d->SetBranchAddress("w",&w);
  tk->SetVarRange(2,0.35,1.0,100);
  for(int i=1; i<=nsel; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(1,w);
    d->SetKinematics(tk);
  }
  tk->FixVariable(2);
  tk->SetVarRange(1,2600,3000,100);
   d->SetBranchAddress("costhkcm",&costhkcm);
  // cos binning: 6 bins (centred 0.4 to 0.9 in steps of 0.1)
  for (int i= 0; i<6 ;i++) 
  {
    d->SetSelection(0);
    double c = 0.4+(i*0.1);
    char query[128];
    sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
    d->SetSelection(d->AddSelection(query)); 
    d->GoTo(0);
    //cout << "Setting costhckm of "<<query << " to "<< costhkcm << endl;
    tk.SetVar(2,costhkcm); // We WILL use the costhkcm from the data (!= 0.9 for fw bin)
    d->SetKinematics(tk);
  }
  d->ViewSelections();
  d->Write();
  delete d;

/*********************************************************************************************/
 d = new TDataset("CLAS_2009__McCracken__PRC_81_025201__rec_HE", "HE Recoil polarization by CLAS");
 
 d->ImportDataset(1,"406"); // 132 datapoints
 
 nsel= d->MakeSelections("w"); // 22 selections
 d->SetBranchAddress("w",&w);
 tk->FixVariable(1);
 tk->SetVarRange(2,0.35,1.0,100);
 for(int i=1; i<=nsel; ++i) {
   d->SetSelection(i);
   d->GoTo(0);
   tk.SetVar(1,w);
   d->SetKinematics(tk);
 }
 
 tk->FixVariable(2);
 tk->SetVarRange(1,2600,3000,100);
 d->SetBranchAddress("costhkcm",&costhkcm);
 // cos binning: 6 bins (centred 0.4 to 0.9 in steps of 0.1)
 for (int i= 0; i<6 ;i++) 
 {
     d->SetSelection(0);
     double c = 0.4+(i*0.1);
     char query[128];
     sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
     d->SetSelection(d->AddSelection(query));
     d->GoTo(0);
     //cout << "Setting costhckm of "<<query << " to "<< costhkcm << endl;
     tk.SetVar(2,costhkcm);
     d->SetKinematics(tk);
 }
 d->ViewSelections();
 d->Write();
 delete d;
 //*********************************************************************************************/
 delete tk;
  file->Close();
}
