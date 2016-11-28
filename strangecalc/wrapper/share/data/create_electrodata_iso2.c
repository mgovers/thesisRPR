//TCalcInfo is actually a wrapper for Data struct (see fitting.h). A datapoint knows its TCalcInfo, but not its TKinematics.
//TKinematics determines the model graph.
//Run this script with strangeViewer and enter the following datasets: 200 201 202 203 204 205 207 208 209.
{
  TFile * storeFile = new TFile("electrodata_iso2.root","RECREATE");
  TDataset::SetDataFolder("./");
  TDataset* d;
  TKinematics* tk = new TKinematics("tk","tk",2,"w:costhkcm:qsquared",2.16,1.0,290000.0);
  int offset = 0;
  double costhkcm,qsquared,w;
//*********************************************************************************************

  d = new TDataset("CLAS_2000__Cha__PhD-unpublished__sigmaL+T","unpolarized cross sec. by CLAS");
  // 11 datapoints
  d->ImportDataset(2,"210");

  tk->FixVariables();
  tk->SetVarRange(2,-1.,1.,50);
  
  tk->SetVar(1,1837.);
  tk->SetVar(3,.477e6);
  d->AddSelection("qsquared>400000",tk);
  
  tk->SetVar(1,1917.);
  tk->SetVar(3,.365e6);
  d->AddSelection("qsquared<400000",tk);

  d->Write();
  delete d;


  //*********************************************************************************************
  delete tk;
  storeFile->Close();
}
