{
  TDataset::SetDataFolder("./");
  TKinematics tk("tk","iso.6 kinematics",6,"costhkcm:wlab:qsquared",1,1000,0);
  double cos,wlab;

  //________________________________________________________
  // LEPS 2006  -  Kohri et al., PRL 97, 082003
  // Differential cross section
  TDataset *dcsdata = new TDataset("LEPS_2006__Kohri__PRL_97_082003__dcs",
				   "LEPS-PRL97(2006)082003");

  dcsdata->ImportDataset(6,"101");
  dcsdata->MakeSelections("costhkcm"); // 4 selections
  dcsdata->MakeSelections("wlab");     // 18 selections

  // cos selections
  tk.SetVarRange(2,900,2500,100);
  dcsdata->SetBranchAddress("costhkcm",&cos);
  for(int i=1; i<=4; ++i) {
    dcsdata->SetSelection(i);
    dcsdata->GoTo(0);
    tk.SetVar(1,cos);
    dcsdata->SetKinematics(&tk);
  }
  
  // wlab selections
  tk.FixVariable(2);
  tk.SetVarRange(1,-1,1,50);
  dcsdata->ResetBranchAddresses();
  dcsdata->SetBranchAddress("wlab",&wlab);
  for(int i=5; i<=22; ++i) {
    dcsdata->SetSelection(i);
    dcsdata->GoTo(0);
    tk.SetVar(2,wlab);
    dcsdata->SetKinematics(&tk);
  }

  //________________________________________________________
  // LEPS 2006  -  Kohri et al., PRL 97, 082003
  // Beam asymmetry
  TDataset *phodata = new TDataset("LEPS_2006__Kohri__PRL_97_082003__pho",
				   "LEPS-PRL97(2006)082003");

  phodata->ImportDataset(6,"102");
  phodata->MakeSelections("costhkcm"); // 4 selections
  phodata->MakeSelections("wlab");     // 9 selections

  // cos selections
  tk.FixVariables();
  tk.SetVarRange(2,900,2500,100);
  phodata->SetBranchAddress("costhkcm",&cos);
  for(int i=1; i<=4; ++i) {
    phodata->SetSelection(i);
    phodata->GoTo(0);
    tk.SetVar(1,cos);
    phodata->SetKinematics(&tk);
  }
  
  // wlab selections
  tk.FixVariable(2);
  tk.SetVarRange(1,-1,1,50);
  phodata->ResetBranchAddresses();
  phodata->SetBranchAddress("wlab",&wlab);
  for(int i=5; i<=13; ++i) {
    phodata->SetSelection(i);
    phodata->GoTo(0);
    tk.SetVar(2,wlab);
    phodata->SetKinematics(&tk);
  }

  //________________________________________________________
  // CLAS 2010  -  Anefalos Pereira et al., arXiv:0912.4833[nucl-ex]
  // Differential cross section
  TDataset *clasdcs = new TDataset("CLAS_2010__Anefalos__arXiv-0912-4833__dcs",
				   "CLAS-PLB688(2010)289");

  clasdcs->ImportDataset(6,"103");
  clasdcs->MakeSelections("wlab");     // 25 selections
  clasdcs->MakeSelections("costhkcm"); // 18 selections

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(1,-1.,1.,50);
  clasdcs->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=25; ++i) {
    clasdcs->SetSelection(i);
    clasdcs->GoTo(0);
    tk.SetVar(2,wlab);
    clasdcs->SetKinematics(&tk);
  }

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(2,1040.,3550.,200);
  clasdcs->SetBranchAddress("costhkcm",&cos);
  for(int i=26; i<=43; ++i) {
    clasdcs->SetSelection(i);
    clasdcs->GoTo(0);
    tk.SetVar(1,cos);
    clasdcs->SetKinematics(&tk);
  }
  clasdcs->ResetBranchAddresses();

  //________________________________________________________
  TFile *file = new TFile("data_iso6.root","RECREATE");
  dcsdata->Write();
  phodata->Write();
  clasdcs->Write();
  file->Close();
  delete dcsdata;
  delete phodata;
  delete clasdcs;
  delete file;
}
