void HWWCut(){
  double xmin = 0.;
  double xmax = 1000.;
  double nbins = 100;
  TString fileInName="yields.root";
  TString outFileName = "cut.root";
  

  TString masses[8]={"300", "400", "500", "600", "700", "800", "900", "1000"};
  

  TFile* fileIn =  new TFile(fileInName);
  TH1F* h_data   = (TH1F*) fileIn->Get("h_data");
  TH1F* h_top    = (TH1F*) fileIn->Get("h_top");
  TH1F* h_ww     = (TH1F*) fileIn->Get("h_ww");
  TH1F* h_dytt   = (TH1F*) fileIn->Get("h_dytt");
  TH1F* h_vv     = (TH1F*) fileIn->Get("h_vv");

  TH1F* h_bkg = h_top->Clone();
  h_bkg->SetNameTitle("background", "background");
  h_bkg->Add(h_ww);
  h_bkg->Add(h_dytt);
  h_bkg->Add(h_vv);

  TFile* fileout = new TFile(outFileName, "recreate");
  fileout->cd();
  double step = (xmax-xmin)/nbins;
  for (unsigned int i = 0; i < 8; ++i){
    TH1F* h_signal = (TH1F*) fileIn->Get("h_signal"+masses[i]);
    TString histoName="significance"+masses[i];
    TH1F* hsignif = new TH1F(histoName, histoName, nbins, xmin, xmax);
    histoName="sig"+masses[i];
    TH1F* hsig    = new TH1F(histoName, histoName, nbins, xmin, xmax);
    histoName="bkg"+masses[i];
    TH1F* hbkg    = new TH1F(histoName, histoName, nbins, xmin, xmax);
    histoName="data"+masses[i];
    TH1F* hdta    = new TH1F(histoName, histoName, nbins, xmin, xmax);
    for (unsigned int j = 1; j < nbins+1; ++j){
      double bkg = h_bkg->Integral(j, nbins+1);
      double sig = h_signal->Integral(j, nbins+1);
      double dta = h_data->Integral(j, nbins+1);
      double signif=sig/sqrt(sig+bkg);
      hsignif->SetBinContent(j, signif);
      hsig->SetBinContent(j, sig);
      hbkg->SetBinContent(j, bkg);
      hdta->SetBinContent(j, dta);
    }
    hsignif->Write();
    hsig->Write();
    hbkg->Write();
    hdta->Write();
    double bestSignif      = hsignif->GetMaximum();
    double bestCut         = xmin + step*hsignif->GetMaximumBin();
    double signalAtCut     = hsig->GetBinContent(hsignif->GetMaximumBin()); 
    double backgroundAtCut = hbkg->GetBinContent(hsignif->GetMaximumBin()); 
    double dataAtCut       = hdta->GetBinContent(hsignif->GetMaximumBin());

    cout << "mass " << masses[i] << "--> best cut at: " << bestCut << 
                                    ", nData: " << dataAtCut << 
                                    ", nSig: " << signalAtCut << 
                                    ", nBkg: " <<  backgroundAtCut << 
                                    ", significance: " << bestSignif << endl;
  } 
}


