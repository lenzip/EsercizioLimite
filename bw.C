void bw(){
  // model to generate events
  RooRealVar mll("mll", "mll", 60, 120);
  RooRealVar m0("m0", "m0", 90);
  RooRealVar Gamma("gamma", "gamma", 2.5);

  RooBreitWigner bw("bw", "bw", mll, m0, Gamma);

  RooDataSet* dataset = bw.generate(mll, 100);

  // now fit that dataset
  RooRealVar m0fit("m0fit", "m0fit", 60, 120);
  RooRealVar Gammafit("gammafir", "gammafit", 0, 10);

  RooBreitWigner bwfit("bwfit", "bwfit", mll, m0fit, Gammafit);
  //Gammafit.setVal(2.5);
  //Gammafit.setConstant(true);

  bwfit.fitTo(*dataset, RooFit::Minos(m0fit));

  // now plot
  RooPlot * frame_mll = mll.frame();
  dataset->plotOn(frame_mll);
  bwfit.plotOn(frame_mll);

  frame_mll->Draw();

  RooAbsReal* nll = bwfit.createNLL(*dataset) ;

  RooMinuit(*nll).migrad() ;
  double minNLL=nll->getVal();

  TCanvas* c = new TCanvas();
  c->cd();

  double m0best = m0fit.getVal();
  double gammabest = Gammafit.getVal();

  RooAbsReal* pll = nll->createProfile(m0fit); 

  //draw the nll function around the minimum and use graphical method to estimate errors
  RooPlot* frame_m0 = m0fit.frame(RooFit::Bins(100), RooFit::Range(89,91), RooFit::Title("LL in m0")) ;
  pll->plotOn(frame_m0, RooFit::ShiftToZero()) ;
  TLine* line = new TLine(89, 0.5, 91, 0.5);
  frame_m0->Draw();
  line->Draw("same");

  
  TH2F* h2d = new TH2F("2d", "2d", 100, 88., 92., 100, 2., 4.);//88., 91., 100, 1.8, 3.9);

  for (unsigned int i = 1; i <= h2d->GetXaxis()->GetNbins(); ++i){
    for (unsigned int j = 1; j <= h2d->GetYaxis()->GetNbins(); ++j){
      double m0here    = h2d->GetXaxis()->GetBinCenter(i);
      double gammaHere = h2d->GetYaxis()->GetBinCenter(j);

      m0fit.setVal(m0here);
      Gammafit.setVal(gammaHere);

      h2d->SetBinContent(i, j, 2*(nll->getVal()-minNLL));

    }
  }

  TCanvas* c2 = new TCanvas();
  c2->cd();

  Double_t contours[2];
  contours[0] = 1.;
  contours[1] = 3.84;
  h2d->Draw("COLZ");
  TH2* h2dclone = (TH2*) h2d->Clone();
  h2dclone->SetContour(2, contours);
  h2dclone->SetLineStyle(2);
  h2dclone->Draw("CONT2 LIST same");
  TMarker* trueVal = new TMarker(90., 2.5, 20);
  trueVal->Draw("sames");
  TMarker* fitVal = new TMarker(m0best, gammabest, 29);
  fitVal->Draw("sames");



}
