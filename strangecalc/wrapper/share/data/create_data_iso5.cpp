{
  // Create data_iso5.root file with command strangeViewer -q create_data_iso5.cpp

  TDataset::SetDataFolder("./");
  TKinematics tk("tk","iso.5 kinematics",5,"costhkcm:wlab:qsquared",1,1000,0);
  double wlab,cos;

  //________________________________________________________
  // CLAS 2010  -  N.Hassall, preliminary data
  // Beam asymmetry

  TDataset *claspho = new TDataset("CLAS_2010__Hassall__prelim__pho",
				   "CLAS preliminary (Hassall)");

  claspho->ImportDataset(5,"101");
  claspho->MakeSelections("wlab"); // 6 selections
  claspho->MakeSelections("costhkcm"); // 6 selections

  // wlab selections
  tk.SetVarRange(1,-1,1,50);
  claspho->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=6; ++i) {
    claspho->SetSelection(i);
    claspho->GoTo(0);
    tk.SetVar(2,wlab);
    claspho->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,915.,2400.,100);
  claspho->ResetBranchAddresses();
  claspho->SetBranchAddress("costhkcm",&cos);
  for(int i=7; i<=12; ++i) {
    claspho->SetSelection(i);
    claspho->GoTo(0);
    tk.SetVar(1,cos);
    claspho->SetKinematics(&tk);
  } 
  claspho->ResetBranchAddresses();

  //________________________________________________________
  // Write datasets to data_iso5.root

  TFile *file = new TFile("data_iso5.root","RECREATE");
  claspho->Write();
  file->Close();
  delete claspho;
  delete file;
}
