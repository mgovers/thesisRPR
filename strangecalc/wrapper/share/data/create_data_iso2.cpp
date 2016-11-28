// Generate data_iso2.root with
// strangeViewer -q create_data_iso2.cpp

{
  TDataset::SetDataFolder("./");
  TDataset *saphir = new TDataset("SAPHIR_2004__Glander__EPJ_A19_251__dcs",
				  "SAPHIR-EPJA19(2004)251");
  saphir->ImportDataset(2,"104");
  saphir->MakeSelections("wlab");
  saphir->MakeSelections("costhkcm");

  TKinematics tk("tk","iso.2 kinematics",2,"wlab:costhkcm:qsquared",1000,1,0);
  double cos,wlab;

  // wlab selections
  tk.SetVarRange(2,-1.,1.,50);
  saphir->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=33; ++i) {
    saphir->SetSelection(i);
    saphir->GoTo(0);
    tk.SetVar(1,wlab);
    saphir->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1046.,2600.,100);
  saphir->ResetBranchAddresses();
  saphir->SetBranchAddress("costhkcm",&cos);
  for(int i=34; i<=53; ++i) {
    saphir->SetSelection(i);
    saphir->GoTo(0);
    tk.SetVar(2,cos);
    saphir->SetKinematics(&tk);
  }
  saphir->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *clas04 = new TDataset("CLAS_2004__McNabb__PRC69_042201R__dcs",
				  "CLAS-PRC69(2004)042201R");
  clas04->ImportDataset(2,"108");
  clas04->MakeSelections("wlab");
  clas04->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  clas04->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=50; ++i) {
    clas04->SetSelection(i);
    clas04->GoTo(0);
    tk.SetVar(1,wlab);
    clas04->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1046.,2600.,100);
  clas04->ResetBranchAddresses();
  clas04->SetBranchAddress("costhkcm",&cos);
  for(int i=51; i<=67; ++i) {
    clas04->SetSelection(i);
    clas04->GoTo(0);
    tk.SetVar(2,cos);
    clas04->SetKinematics(&tk);
  }
  clas04->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *clas06 = new TDataset("CLAS_2006__Bradford__PRC73_035202__dcs",
				  "CLAS-PRC73(2006)035202");
  clas06->ImportDataset(2,"113");
  clas06->MakeSelections("wlab");
  clas06->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  clas06->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=73; ++i) {
    clas06->SetSelection(i);
    clas06->GoTo(0);
    tk.SetVar(1,wlab);
    clas06->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1046.,3000.,100);
  clas06->ResetBranchAddresses();
  clas06->SetBranchAddress("costhkcm",&cos);
  for(int i=74; i<=91; ++i) {
    clas06->SetSelection(i);
    clas06->GoTo(0);
    tk.SetVar(2,cos);
    clas06->SetKinematics(&tk);
  }
  clas06->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *leps06 = new TDataset("LEPS_2006__Sumihama__PRC73_035214__dcs",
				  "LEPS-PRC73(2006)035214");
  leps06->ImportDataset(2,"118");
  leps06->MakeSelections("wlab");
  leps06->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  leps06->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=18; ++i) {
    leps06->SetSelection(i);
    leps06->GoTo(0);
    tk.SetVar(1,wlab);
    leps06->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1046.,2600.,100);
  leps06->ResetBranchAddresses();
  leps06->SetBranchAddress("costhkcm",&cos);
  for(int i=19; i<=21; ++i) {
    leps06->SetSelection(i);
    leps06->GoTo(0);
    tk.SetVar(2,cos);
    leps06->SetKinematics(&tk);
  }
  leps06->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *lepsPRL = new TDataset("LEPS_2006__Kohri__PRL_082003__dcs",
				   "LEPS-PRL97(2006)082003");
  lepsPRL->ImportDataset(2,"123");
  lepsPRL->MakeSelections("wlab");
  lepsPRL->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  lepsPRL->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=18; ++i) {
    lepsPRL->SetSelection(i);
    lepsPRL->GoTo(0);
    tk.SetVar(1,wlab);
    lepsPRL->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1500.,2400.,100);
  lepsPRL->ResetBranchAddresses();
  lepsPRL->SetBranchAddress("costhkcm",&cos);
  for(int i=19; i<=22; ++i) {
    lepsPRL->SetSelection(i);
    lepsPRL->GoTo(0);
    tk.SetVar(2,cos);
    lepsPRL->SetKinematics(&tk);
  }
  lepsPRL->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *recSaphir = new TDataset("SAPHIR_2004__Glander__EPJ_A19_251__rec",
				     "SAPHIR-EPJA19(2004)251");
  recSaphir->ImportDataset(2,"106");
  recSaphir->MakeSelections("wlab");
  recSaphir->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  recSaphir->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=4; ++i) {
    recSaphir->SetSelection(i);
    recSaphir->GoTo(0);
    tk.SetVar(1,wlab);
    recSaphir->SetKinematics(&tk);
  }

  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,2400.,100);
  recSaphir->ResetBranchAddresses();
  recSaphir->SetBranchAddress("costhkcm",&cos);
  for(int i=5; i<=8; ++i) {
    recSaphir->SetSelection(i);
    recSaphir->GoTo(0);
    tk.SetVar(2,cos);
    recSaphir->SetKinematics(&tk);
  }
  recSaphir->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *recClas04 = new TDataset("CLAS_2004__McNabb__PRC69_042201R__rec",
				     "CLAS-PRC69(2004)042201R");
  recClas04->ImportDataset(2,"110");
  recClas04->MakeSelections("wlab");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  recClas04->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=25; ++i) {
    recClas04->SetSelection(i);
    recClas04->GoTo(0);
    tk.SetVar(1,wlab);
    recClas04->SetKinematics(&tk);
  }
  recClas04->ResetBranchAddresses();

  //________________________________________________________________________
  TDataset *graal = new TDataset("GRAAL_2007__Lleres__EPJA31_79__rec",
				 "GRAAL-EPJA31(2007)79");
  graal->ImportDataset(2,"121");
  graal->MakeSelections("wlab");

  // wlab selections
  graal->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=2; ++i) {
    graal->SetSelection(i);
    graal->GoTo(0);
    tk.SetVar(1,wlab);
    graal->SetKinematics(&tk);
  }
  graal->ResetBranchAddresses();

  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,1500.,100);
  graal->SetSelection(0);
  tk.SetVar(2,0.78);
  graal->AddSelection("costhkcm>.7",&tk);
  tk.SetVar(2,0.21);
  graal->AddSelection("costhkcm>.2 && costhkcm<.23",&tk);
  tk.SetVar(2,-.241);
  graal->AddSelection("costhkcm==-.241",&tk);
  tk.SetVar(2,-.65);
  graal->AddSelection("costhkcm<-.6",&tk);

  //________________________________________________________________________
  TDataset *phoGraal = new TDataset("GRAAL_2007__Lleres__EPJA31_79__pho",
				    "GRAAL-EPJA31(2007)79");
  phoGraal->ImportDataset(2,"119");
  phoGraal->MakeSelections("wlab");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  phoGraal->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=7; ++i) {
    phoGraal->SetSelection(i);
    phoGraal->GoTo(0);
    tk.SetVar(1,wlab);
    phoGraal->SetKinematics(&tk);
  }
  phoGraal->ResetBranchAddresses();

  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,1500.,100);
  phoGraal->SetSelection(0);
  tk.SetVar(2,.95);
  phoGraal->AddSelection("costhkcm<1. && costhkcm>.9",&tk);
  tk.SetVar(2,.77);
  phoGraal->AddSelection("costhkcm<.8 && costhkcm>.7",&tk);
  tk.SetVar(2,-.01);
  phoGraal->AddSelection("costhkcm<0. && costhkcm>-.1",&tk);
  tk.SetVar(2,-.4);
  phoGraal->AddSelection("costhkcm<-.3 && costhkcm>-.5",&tk);

  //________________________________________________________________________
  TDataset *phoLeps03 = new TDataset("LEPS_2003__Zegers__PRL91_092001__pho",
				     "LEPS-PRL91(2003)092001");
  phoLeps03->ImportDataset(2,"117");
  phoLeps03->MakeSelections("wlab");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  phoLeps03->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=9; ++i) {
    phoLeps03->SetSelection(i);
    phoLeps03->GoTo(0);
    tk.SetVar(1,wlab);
    phoLeps03->SetKinematics(&tk);
  }
  phoLeps03->ResetBranchAddresses();

  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,1500.,100);
  phoLeps03->SetSelection(0);
  tk.SetVar(2,.97);
  phoLeps03->AddSelection("costhkcm<1. && costhkcm>.96",&tk);
  tk.SetVar(2,.93);
  phoLeps03->AddSelection("costhkcm<.96 && costhkcm>.9",&tk);
  tk.SetVar(2,.86);
  phoLeps03->AddSelection("costhkcm<.9 && costhkcm>.8",&tk);
  tk.SetVar(2,.75);
  phoLeps03->AddSelection("costhkcm<.8 && costhkcm>.7",&tk);
  tk.SetVar(2,.65);
  phoLeps03->AddSelection("costhkcm<.7 && costhkcm>.5",&tk);

  //________________________________________________________________________
  TDataset *phoLeps06 = new TDataset("LEPS_2006__Kohri__PRL_082003__pho",
				     "LEPS-PRL97(2006)082003");
  phoLeps06->ImportDataset(2,"124");
  phoLeps06->MakeSelections("wlab");
  phoLeps06->MakeSelections("costhkcm");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  phoLeps06->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=9; ++i) {
    phoLeps06->SetSelection(i);
    phoLeps06->GoTo(0);
    tk.SetVar(1,wlab);
    phoLeps06->SetKinematics(&tk);
  }
 
  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,2400.,100);
  phoLeps06->ResetBranchAddresses();
  phoLeps06->SetBranchAddress("costhkcm",&cos);
  for(int i=10; i<=13; ++i) {
    phoLeps06->SetSelection(i);
    phoLeps06->GoTo(0);
    tk.SetVar(2,cos);
    phoLeps06->SetKinematics(&tk);
  }
  phoLeps06->ResetBranchAddresses();
  
  //________________________________________________________________________
  TDataset *cxz = new TDataset("CLAS_2007__Bradford__PRC75_035205__cxz",
			       "CLAS-PRC75(2007)035205");
  cxz->ImportDataset(2,"116");
  cxz->MakeSelections("wlab:observable");
  cxz->MakeSelections("costhkcm:observable");

  // wlab selections
  tk.FixVariable(1);
  tk.SetVarRange(2,-1.,1.,50);
  cxz->SetBranchAddress("wlab",&wlab);
  for(int i=1; i<=32; ++i) {
    cxz->SetSelection(i);
    cxz->GoTo(0);
    tk.SetVar(1,wlab);
    cxz->SetKinematics(&tk);
  }
 
  // cos selections
  tk.FixVariable(2);
  tk.SetVarRange(1,1047.,3000.,100);
  cxz->ResetBranchAddresses();
  cxz->SetBranchAddress("costhkcm",&cos);
  for(int i=33; i<=44; ++i) {
    cxz->SetSelection(i);
    cxz->GoTo(0);
    tk.SetVar(2,cos);
    cxz->SetKinematics(&tk);
  }
  cxz->ResetBranchAddresses();

  //________________________________________________________________________
  diffClas10 = new TDataset("CLAS_2010__Dey__PRC_82_025202__dcs", "CLAS-PRC82(2010)025202");

  diffClas10->ImportDataset(2,"132"); // 2089 datapoints
  
  // Wlab bins
  int nsel= diffClas10->MakeSelections("wlab"); // 112 selections
  diffClas10->SetBranchAddress("wlab",&wlab);
  tk.FixVariables();
  tk.SetVarRange(2,-1.0,1.0,100);
  for(int i=1; i<=nsel; ++i) {
    diffClas10->SetSelection(i);
    diffClas10->GoTo(0);
    tk.SetVar(1,wlab);
    diffClas10->SetKinematics(&tk);
  }
  tk.FixVariable(2);
  
  // Cos bin: 19 bins (centred -0.9 to 0.9 in steps of 0.1)
  tk.SetVarRange(1,938.0,3815.0,200);
  for (int i= 0; i<19 ;i++) 
  {
    diffClas10->SetSelection(0);
    double c = -0.9+(i*0.1);
    char query[128];
    sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
    diffClas10->SetSelection(diffClas10->AddSelection(query)); 
    tk.SetVar(2,c);
    diffClas10->SetKinematics(&tk);
  }
  diffClas10->ViewSelections();

  //________________________________________________________________________
  recClas10 = new TDataset("CLAS_2010__Dey__PRC_82_025202__rec", "CLAS-PRC82(2010)025202");

  recClas10->ImportDataset(2,"134"); // 455 datapoints
  
  // Cos bin: 18 bins (centred -0.8 to 0.9 in steps of 0.1)
  tk.SetVarRange(1,938.0,3815.0,200);
  for (int i= 0; i<19 ;i++) 
  {
    recClas10->SetSelection(0);
    double c = -0.9+(i*0.1);
    char query[128];
    sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
    recClas10->SetSelection(recClas10->AddSelection(query)); 
    tk.SetVar(2,c);
    recClas10->SetKinematics(&tk);
  }
  recClas10->ViewSelections();
  
  //________________________________________________________________________
  TFile *file = new TFile("data_iso2.root","RECREATE");
  saphir->Write();
  clas04->Write();
  clas06->Write();
  leps06->Write();
  lepsPRL->Write();
  recSaphir->Write();
  recClas04->Write();
  graal->Write();
  phoGraal->Write();
  phoLeps03->Write();
  phoLeps06->Write();
  cxz->Write();
  diffClas10->Write();
  recClas10->Write();
  file->Close();
  delete saphir;
  delete clas04;
  delete clas06;
  delete leps06;
  delete lepsPRL;
  delete recSaphir;
  delete recClas04;
  delete graal;
  delete phoGraal;
  delete phoLeps03;
  delete phoLeps06;
  delete cxz;
  delete diffClas10;
  delete recClas10;
  delete file;

}
