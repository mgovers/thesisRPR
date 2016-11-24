{
  TFile * file = new TFile("data_iso2_HE.root","RECREATE");
  TDataset::SetDataFolder("./");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","iso.1 HE kinematics",2,"t:wlab:qsquared",0,5000,0);
  double wlab,w,costhkcm;
  
//*********************************************************************************************
  
  d = new TDataset("SLAC_1969__Boyarski__PRL_22_1131", "HE Diff.cross sections by SLAC");

  d->ImportDataset(2,"401"); //48 datapoints
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

  d->ImportDataset(2,"403"); //9 datapoints (NOTE: check this! trailing newlines cause phantom data points)
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
  
  //________________________________________________________________________
  d = new TDataset("CLAS_2010__Dey__PRC_82_025202__dcs_HE", "CLAS-PRC82(2010)025202-HE");

  d->ImportDataset(2,"404"); // 132 datapoints
  
  int nsel= d->MakeSelections("w"); // 22 selections
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
  
  // cos binning: 6 bins (centred 0.4 to 0.9 in steps of 0.1)
  d->SetBranchAddress("costhkcm",&costhkcm);
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

  //________________________________________________________________________
  d = new TDataset("CLAS_2010__Dey__PRC_82_025202__rec_HE", "CLAS-PRC82(2010)025202-HE");

  d->ImportDataset(2,"405"); // 45 datapoints
 nsel= d->MakeSelections("w"); // 16 selections
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
 
 // cos binning: 6 bins (centred 0.4 to 0.9 in steps of 0.1)
 d->SetBranchAddress("costhkcm",&costhkcm);
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

 //*********************************************************************************************/

 file->Close();
 delete tk;
  
}
