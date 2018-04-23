#include "TGraphAsymmErrors.h"

void HWWLimit(){
  TString fileInName="cut.root";
  TString outFileName = "limit.root";

  const unsigned int nMasses=7;

  TString masses[nMasses]={"300", "400", "500", "600", "700", "900", "1000"};
  double massesD[nMasses]={300, 400, 500, 600, 700, 900, 1000};

  TFile filein(fileInName);

  double expected[nMasses];
  double expectedUp1[nMasses];
  double expectedUp2[nMasses];
  double expectedDo1[nMasses];
  double expectedDo2[nMasses];
  double observed[nMasses];

  TRandom3 rand;
  
  for (unsigned int i = 0; i < nMasses; ++i){
    TH1F* hdta = (TH1F*) filein.Get("data"+masses[i]);
    TH1F* hbkg = (TH1F*) filein.Get("bkg"+masses[i]);
    TH1F* hsig = (TH1F*) filein.Get("sig"+masses[i]);
    TH1F* hsignif = (TH1F*) filein.Get("significance"+masses[i]);

    double signalAtCut     = hsig->GetBinContent(hsignif->GetMaximumBin());
    double backgroundAtCut = hbkg->GetBinContent(hsignif->GetMaximumBin());
    double dataAtCut       = hdta->GetBinContent(hsignif->GetMaximumBin());

    // compute expected limit with 10% uncertainty on the background
    // observation is 0 in this case
    // assume 10% error on background
    double sigmaGaus = 0.10*backgroundAtCut;
    double sigmaStat = sqrt(backgroundAtCut);
    double sigmaTot = sqrt(sigmaGaus*sigmaGaus + sigmaStat*sigmaStat);
   
   
    double expectedLimit = 0. + sigmaTot*1.645;
    Double_t xq[5] = {0.025, 0.34, 0.5, 0.84, 0.975};  // position where to compute the quantiles in [0,1]
    Double_t yq[5];
    TH1F* hToy = new TH1F("hToy", "hToy", 2000, -1000, 1000); 
  
    for (unsigned k =0; k < 1000; ++k){
      int dataToy = rand.Poisson(backgroundAtCut);
      double observedLimitToy = (dataToy - backgroundAtCut) + sigmaTot*1.645;
      //cout << "observed limit toy " << observedLimitToy << endl;
      hToy->Fill(observedLimitToy);
    }
    hToy->GetQuantiles(5,yq,xq);
    
    double observedLimit = (dataAtCut - backgroundAtCut) + sigmaTot*1.645;

    cout << "mass " << masses[i] << "-->"
         << "expected limit: " << expectedLimit << " or " << yq[2] 
         << ", 2 sigma down: " << yq[0]
         << ", 1 sigma down: " << yq[1]
         << ", 1 sigma up: " << yq[3]
         << ", 2 sigma up: " << yq[4]
         << ", observed: " << observedLimit << endl;
    delete hToy;
    expected[i] = expectedLimit/signalAtCut;
    observed[i] = observedLimit/signalAtCut;
    expectedUp1[i] = (yq[3]-expectedLimit)/signalAtCut;
    expectedUp2[i] = (yq[4]-expectedLimit)/signalAtCut;
    expectedDo1[i] = (expectedLimit-yq[1])/signalAtCut;
    expectedDo2[i] = (expectedLimit-yq[0])/signalAtCut;

  }
 
  TFile out(outFileName, "recreate");
  out.cd();
  TCanvas * c = new TCanvas();
  c->cd();
  TGraphAsymmErrors * twoSigmaBandG = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, expectedUp2, expectedDo2);
  twoSigmaBandG->SetFillColor(kYellow);
  TGraphAsymmErrors * oneSigmaBandG = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, expectedUp1, expectedDo1);
  oneSigmaBandG->SetFillColor(kGreen);
  TGraphAsymmErrors * expectedG   = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, 0, 0);
  expectedG->SetLineStyle(2);
  TGraphAsymmErrors * observedG   = new TGraphAsymmErrors(nMasses, massesD, observed, 0, 0, 0, 0);

  twoSigmaBandG->Draw("A3");
  oneSigmaBandG->Draw("3 same");
  expectedG->Draw("L same");
  observedG->Draw("L same");
  c->Write();
}
