//$\gamma p \rightarrow K^{+} \Lambda$ results from McCracken, et al. (CLAS Collaboration) 2009
//for submission to PRC 
//preprint found at http://arxiv.org/abs/0912.4274
//published Phys.Rev.C81:025201,2010 
{
  TFile * file = new TFile("/home/lesley/strangecalc-data/data_iso1.root","UPDATE");
  TDataset::SetDataFolder("/home/lesley/strangecalc-data/");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","iso.1 kinematics",1,"wlab:costhkcm:qsquared",1000,1,0);
  double wlab=0.0;
  int nsel=0;
//*********************************************************************************************
  /*
  NB:No differential cross section results are given for the 1.955, 2.735, or 2.745 GeV bins.
  */
  d = new TDataset("CLAS_2009__McCracken__PRC_81_025201__dcs", "Diff.cross sections by CLAS");

  d->ImportDataset(1,"134"); // 2066 datapoints
  
  nsel= d->MakeSelections("wlab"); // 119 selections
  d->SetBranchAddress("wlab",&wlab);
  tk->SetVarRange(2,-1.0,1.0,100);
  for(int i=1; i<=nsel; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(1,wlab);
    d->SetKinematics(tk);
  }
  tk->FixVariable(2);
  tk->SetVarRange(1,938.0,3815.0,200);
  
  // coarser cos binning: 19 bins (centred -0.9 to 0.9 in steps of 0.1)
  for (int i= 0; i<19 ;i++) 
  {
    d->SetSelection(0);
    double c = -0.9+(i*0.1);
    char query[128];
    sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
    d->SetSelection(d->AddSelection(query)); 
    tk.SetVar(2,c);
    d->SetKinematics(tk);
  }
  d->ViewSelections();
  d->Write();
  delete d;
  
//*********************************************************************************************
  d = new TDataset("CLAS_2009__McCracken__PRC_81_025201__rec", "Recoil polarization by CLAS");

  d->ImportDataset(1,"135"); // 1707 datapoints
  
  nsel= d->MakeSelections("wlab"); // 122 selections
  d->SetBranchAddress("wlab",&wlab);
  tk->FixVariable(1);
  tk->SetVarRange(2,-1.0,1.0,100);
  for(int i=1; i<=nsel; ++i) {
    d->SetSelection(i);
    d->GoTo(0);
    tk.SetVar(1,wlab);
    d->SetKinematics(tk);
  }
  
  tk->FixVariable(2);
  tk->SetVarRange(1,938.0,3815.0,200);
  
  // coarser cos binning: 18 bins (centred -0.8 to 0.9 in steps of 0.1)
  for (int i= 0; i<18 ;i++) 
  {
    d->SetSelection(0);
    double c = -0.8+ (i*0.1);
    char query[128];
    sprintf(query, "costhkcm>%g&&costhkcm<%g",c-0.05,c+0.05);
    d->SetSelection(d->AddSelection(query)); 
    tk.SetVar(2,c);
    d->SetKinematics(tk);
  }
  d->ViewSelections();
  d->Write();
  delete d;
   
//**********************************************************************************************/
  
 file->Close();
 delete tk;
  
}
