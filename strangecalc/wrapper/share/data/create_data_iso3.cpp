{
  TDataset::SetDataFolder("./");
 
  //________________________________________________________
  // SAPHIR 2005  -  Lawall et al., Eur.Phys.J. A24(2005)275
  // Differential cross section

  TDataset *dcs = new TDataset("SAPHIR_2005__Lawall__EPJ_A24_275__dcs",
			       "SAPHIR-EPJA24(2005)275");

  dcs->ImportDataset(3,"104");
  dcs->MakeSelections("wlab"); // 12 selections
  dcs->MakeSelections("costhkcm"); // 10 selections

  TKinematics tk("tk","iso.3 kinematics",3,"costhkcm:wlab",1,1000);
  double cos,wlab;

  // wlab selections
  tk.SetVarRange(1,-1,1,50);
  dcs->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=12; ++i) {
    dcs->SetSelection(i);
    dcs->GoTo(0);
    tk.SetVar(2,wlab);
    dcs->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,2600.,100);
  dcs->ResetBranchAddresses();
  dcs->SetBranchAddress("costhkcm",&cos);
  for(int i=13; i<=22; ++i) {
    dcs->SetSelection(i);
    dcs->GoTo(0);
    tk.SetVar(1,cos);
    dcs->SetKinematics(&tk);
  }

  //________________________________________________________
  // SAPHIR 2005  -  Lawall et al., Eur.Phys.J. A24(2005)275
  // Total cross section

  TDataset *tcs = new TDataset("SAPHIR_2005__Lawall__EPJ_A24_275__tcs",
			       "SAPHIR-EPJA24(2005)275");
  
  tcs->ImportDataset(3,"108");  
  tcs->SetKinematics(&tk);

  //________________________________________________________
  // SAPHIR 2005  -  Lawall et al., Eur.Phys.J. A24(2005)275
  // Recoil asymmetry
  
  TDataset *rec = new TDataset("SAPHIR_2005__Lawall__EPJ_A24_275__rec",
			       "SAPHIR-EPJA24(2005)275");
  
  rec->ImportDataset(3,"106");  
  rec->MakeSelections("wlab"); // 2 selections
  rec->MakeSelections("costhkcm"); // 5 selections

  tk.FixVariables();
  double cos,wlab;

  // wlab selections
  tk.SetVarRange(1,-1,1,50);
  rec->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=2; ++i) {
    rec->SetSelection(i);
    rec->GoTo(0);
    tk.SetVar(2,wlab);
    rec->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,2600.,100);
  rec->ResetBranchAddresses();
  rec->SetBranchAddress("costhkcm",&cos);
  for(int i=3; i<=7; ++i) {
    rec->SetSelection(i);
    rec->GoTo(0);
    tk.SetVar(1,cos);
    rec->SetKinematics(&tk);
  }

  //________________________________________________________
  // CB-ELSA/TAPS 2007  -  Castelijns et al., Eur.Phys.J. A31(2007)61
  // Differential cross section

  TDataset *elsadcs = new TDataset("CBELSA-TAPS_2007__Castelijns__EPJ_A31_61__dcs",
				   "CB/ELSA-TAPS-EPJA31(2007)61");

  elsadcs->ImportDataset(3,"110");
  elsadcs->MakeSelections("wlab"); // 12 selections
  elsadcs->MakeSelections("costhkcm"); // 6 selections

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(1,-1,1,50);
  elsadcs->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=12; ++i) {
    elsadcs->SetSelection(i);
    elsadcs->GoTo(0);
    tk.SetVar(2,wlab);
    elsadcs->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,2600.,100);
  elsadcs->ResetBranchAddresses();
  elsadcs->SetBranchAddress("costhkcm",&cos);
  for(int i=13; i<=18; ++i) {
    elsadcs->SetSelection(i);
    elsadcs->GoTo(0);
    tk.SetVar(1,cos);
    elsadcs->SetKinematics(&tk);
  }  

  //________________________________________________________
  // CB-ELSA/TAPS 2007  -  Castelijns et al., Eur.Phys.J. A31(2007)61
  // Total cross section

  TDataset *elsatcs = new TDataset("CBELSA-TAPS_2007__Castelijns__EPJ_A31_61__tcs",
				   "CB/ELSA-TAPS-EPJA31(2007)61");
  
  elsatcs->ImportDataset(3,"111");
  elsatcs->SetKinematics(&tk);

  //________________________________________________________
  // CB-ELSA/TAPS 2007  -  Castelijns et al., Eur.Phys.J. A31(2007)61
  // Recoil asymmetry
  
  TDataset *elsarec = new TDataset("CBELSA-TAPS_2007__Castelijns__EPJ_A31_61__rec",
				   "CB/ELSA-TAPS-EPJA31(2007)61");
  
  elsarec->ImportDataset(3,"112");  
  elsarec->MakeSelections("wlab"); // 12 selections
  elsarec->MakeSelections("costhkcm"); // 6 selections

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(1,-1,1,50);
  elsarec->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=12; ++i) {
    elsarec->SetSelection(i);
    elsarec->GoTo(0);
    tk.SetVar(2,wlab);
    elsarec->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,2600.,100);
  elsarec->ResetBranchAddresses();
  elsarec->SetBranchAddress("costhkcm",&cos);
  for(int i=13; i<=18; ++i) {
    elsarec->SetSelection(i);
    elsarec->GoTo(0);
    tk.SetVar(1,cos);
    elsarec->SetKinematics(&tk);
  }

  //________________________________________________________
  // CLAS 2003  -  Carnahan, Ph.D. CUA
  // Differential cross section

  TDataset *clasdcs = new TDataset("CLAS_2003__Carnahan__PhD__dcs",
				   "CLAS-Ph.D. Carnahan (2003)");

  clasdcs->ImportDataset(3,"109");
  clasdcs->MakeSelections("wlab"); // 6 selections
  clasdcs->MakeSelections("costhkcm"); // 8 selections

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(1,-1,1,50);
  clasdcs->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=6; ++i) {
    clasdcs->SetSelection(i);
    clasdcs->GoTo(0);
    tk.SetVar(2,wlab);
    clasdcs->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,2600.,100);
  clasdcs->ResetBranchAddresses();
  clasdcs->SetBranchAddress("costhkcm",&cos);
  for(int i=7; i<=14; ++i) {
    clasdcs->SetSelection(i);
    clasdcs->GoTo(0);
    tk.SetVar(1,cos);
    clasdcs->SetKinematics(&tk);
  }  

  //________________________________________________________
  // MAMI-C 2012  -  Sergey Prakhov (unpublished)
  // Differential cross section + recoil asymmetries

  TDataset *mamidcs = new TDataset("MAMI_2012__Prakhov__unpublished__dcs",
				   "MAMI-C unpublished (2012)");
  TDataset *mamirec = new TDataset("MAMI_2012__Prakhov__unpublished__rec",
				   "MAMI-C unpublished (2012)");

  mamidcs->ImportDataset(3,"113");
  mamidcs->MakeSelections("wlab"); // 14 selections
  mamidcs->MakeSelections("costhkcm"); // 14 selections
  mamirec->ImportDataset(3,"114");
  mamirec->MakeSelections("wlab"); // 14 selections
  mamirec->MakeSelections("costhkcm"); // 14 selections

  // wlab selections
  tk.FixVariables();
  tk.SetVarRange(1,-1,1,50);
  mamidcs->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=14; ++i) {
    mamidcs->SetSelection(i);
    mamirec->SetSelection(i);
    mamidcs->GoTo(0);
    tk.SetVar(2,wlab);
    mamidcs->SetKinematics(&tk);
    mamirec->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(1);
  tk.SetVarRange(2,1047.,1600.,50);
  mamidcs->ResetBranchAddresses();
  mamidcs->SetBranchAddress("costhkcm",&cos);
  for(int i=15; i<=28; ++i) {
    mamidcs->SetSelection(i);
    mamirec->SetSelection(i);
    mamidcs->GoTo(0);
    tk.SetVar(1,cos);
    mamidcs->SetKinematics(&tk);
    mamirec->SetKinematics(&tk);
  }  

  //________________________________________________________
  // Write datasets to data_iso3.root

  TFile *file = new TFile("data_iso3.root","RECREATE");
  dcs->Write();
  tcs->Write();
  rec->Write();
  elsadcs->Write();
  elsatcs->Write();
  elsarec->Write();
  clasdcs->Write();
  mamidcs->Write();
  mamirec->Write();
  file->Close();
  delete dcs;
  delete tcs;
  delete rec;
  delete elsadcs;
  delete elsatcs;
  delete elsarec;
  delete clasdcs;
  delete mamidcs;
  delete mamirec;
}
