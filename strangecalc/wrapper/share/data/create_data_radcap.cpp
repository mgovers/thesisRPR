{
  /*
    Use this script to generate data_radcap.root with the command
    > strangeViewer -q create_data_radcap.cpp 
    from within your local working copy of the strangecalc-data folder.
   */

  TDataset::SetDataFolder("./");

  //
  // K- g --> p Lambda (new analysis)
  //=================================
  TDataset *iso11new = new TDataset("CB_2009__Prakhov__priv-comm__dcs11",
				    "Diff. cross sections for Kg->pLamda");
  
  iso11new->ImportDataset(11,"304");
  iso11new->MakeSelections("pklab");    // selection 1 -> 8
  iso11new->MakeSelections("costhkcm"); // selection 9 -> 20

  TKinematics kin("kin","Kinematics for $K^- p #rightarrow #gamma #Lambda", 
		  11,"pklab:costhkcm:qsquared", 
		  100, -0.916666666, 0.0);
  double cos,pklab;

  // pklab selections
  kin.SetVarRange(2,-1.,1.,30);
  iso11new->SetBranchAddress("pklab",&pklab);
  for(int i=1; i<=8; ++i) {
    iso11new->SetSelection(i);
    iso11new->GoTo(0);
    kin.SetVar(1,pklab);
    iso11new->SetKinematics(&kin);
  }
  iso11new->ResetBranchAddresses();

  // cos selections
  kin.FixVariable(2);
  kin.SetVarRange(1,10.,1000.,50);
  iso11new->SetBranchAddress("costhkcm",&cos);
  for(int i=9; i<=20; ++i) {
    iso11new->SetSelection(i);
    iso11new->GoTo(0);
    kin.SetVar(2,cos);
    iso11new->SetKinematics(&kin);
  }
  iso11new->ResetBranchAddresses();


  //
  // Total cross section
  // K- g --> p Lambda (new analysis)
  //=================================
  TDataset *total11new = new TDataset("CB_2009__Prakhov__priv-comm__tcs11",
				    "Diff. cross sections for Kg->pLamda");
  
  total11new->ImportDataset(11,"305");
  kin.SetVarRange(1,10.,1000.,50);
  total11new->SetKinematics(&kin);


  //
  // K- g --> p Sigma (new analysis)
  //================================
  TDataset *iso12new = new TDataset("CB_2009__Prakhov__priv-comm__dcs12",
				    "Diff. cross sections for Kg->pSigma");
  
  iso12new->ImportDataset(12,"304");
  iso12new->MakeSelections("pklab");    // selection 1 -> 8
  iso12new->MakeSelections("costhkcm"); // selection 9 -> 20

  kin.FixVariable(1);
  kin.SetIsospin(12);
  kin.SetTitle("Kinematics for $K^- p #rightarrow #gamma #Sigma");
 
  // pklab selections
  kin.SetVarRange(2,-1.,1.,30);
  iso12new->SetBranchAddress("pklab",&pklab);
  for(int i=1; i<=8; ++i) {
    iso12new->SetSelection(i);
    iso12new->GoTo(0);
    kin.SetVar(1,pklab);
    iso12new->SetKinematics(&kin);
  }
  iso12new->ResetBranchAddresses();

  // cos selections
  kin.FixVariable(2);
  kin.SetVarRange(1,10.,1000.,50);
  iso12new->SetBranchAddress("costhkcm",&cos);
  for(int i=9; i<=20; ++i) {
    iso12new->SetSelection(i);
    iso12new->GoTo(0);
    kin.SetVar(2,cos);
    iso12new->SetKinematics(&kin);
  }
  iso12new->ResetBranchAddresses();


  //
  // Total cross section
  // K- g --> p Sigma (new analysis)
  //=================================
  TDataset *total12new = new TDataset("CB_2009__Prakhov__priv-comm__tcs12",
				      "Total cross sections for Kg->pSigma");
  
  total12new->ImportDataset(12,"305");
  kin.SetVarRange(1,10.,1000.,50);
  total12new->SetKinematics(&kin);


  //
  // Total cross section
  // K- g --> p Sigma (old analysis)
  //=================================
  TDataset *total12old = new TDataset("CB_2009__Stanislaus__PRC79-015203",
				      "Total cross sections for Kg->pSigma");
  
  total12old->ImportDataset(12,"306");
  kin.SetVarRange(1,10.,1000.,50);
  total12old->SetKinematics(&kin);
  

  TFile *file = new TFile("data_radcap.root","RECREATE");
  iso11new->Write();
  iso12new->Write();
  total11new->Write();
  total12new->Write();
  total12old->Write();
  file->Close();
  delete iso11new;
  delete iso12new;
  delete total11new;
  delete total12new;
  delete total12old;
  delete file;

}

