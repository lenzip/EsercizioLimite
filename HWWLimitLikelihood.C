#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "HWWWorkspace.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"

#include <iostream>

using namespace RooFit ;
using namespace std ;

std::pair<double, RooAbsReal*> getLimit(TH1* hdata, RooWorkspace * w, bool doCleanup=true){ 
  RooRealVar * x = w->var("x");
  RooAbsPdf* model = w->pdf("model");
  RooRealVar * mu = w->var("mu");
  const RooArgSet* constraints = w->set("constraints");
 
  RooDataHist data("data", "data", *x, hdata);  


  RooAbsReal* nll = model->createNLL(data, ExternalConstraints(*constraints) );
  RooMinuit m(*nll);
  m.setPrintLevel(-1000);
  m.migrad();
  m.hesse();
  m.migrad();
  // extract the profile likelihood in the signal strength parameter
  // to be used toset the limit
  RooAbsReal* profileLogLikelihood = nll->createProfile(*mu) ;
  double limit = -1;
  // find the limit as the point at which the NLL crosses 2
  for (unsigned int k=0; k < 1000; ++k){
    mu->setVal(k*0.01);
    double deltaNLL = profileLogLikelihood->getVal();
    if (deltaNLL > 2){
      limit = k*0.01;
      break;
    }
  } 
  if (doCleanup) delete nll;
  return std::make_pair(limit, profileLogLikelihood);
}

void HWWLimitLikelihood2(){
  TString fileCutName="cut.root";
  TString filePlotName="yields.root";
  
  TFile outputfile("limitProfileLikelihood.root", "recreate");
  const unsigned int nMasses=7;

  TString masses[nMasses]={"300", "400", "500", "600", "700", "900", "1000"};
  double massesD[nMasses]={300, 400, 500, 600, 700, 900, 1000};

  TFile filePlotIn(filePlotName);
  TFile fileCutIn(fileCutName);

  // take histograms in
  TH1F* h_data   = (TH1F*) filePlotIn.Get("h_data");
  TH1F* h_top    = (TH1F*) filePlotIn.Get("h_top");
  TH1F* h_ww     = (TH1F*) filePlotIn.Get("h_ww");
  TH1F* h_dytt   = (TH1F*) filePlotIn.Get("h_dytt");
  TH1F* h_vv     = (TH1F*) filePlotIn.Get("h_vv"); 
  
  TRandom3 rand;
  // silence roofit
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
 
  // vectors to store limit values
  double expected[nMasses];
  double expectedUp1[nMasses];
  double expectedUp2[nMasses];
  double expectedDo1[nMasses];
  double expectedDo2[nMasses];
  double observed[nMasses];

  // loop over masses
  for (unsigned int i = 0; i < nMasses; ++i){
    TH1F* hdta = (TH1F*) fileCutIn.Get("data"+masses[i]);
    TH1F* hbkg = (TH1F*) fileCutIn.Get("bkg"+masses[i]);
    TH1F* hsig = (TH1F*) fileCutIn.Get("sig"+masses[i]);
    TH1F* hsignif = (TH1F*) fileCutIn.Get("significance"+masses[i]);


    // get process yields cutting at max significance
    double signalAtCut     = hsig->GetBinContent(hsignif->GetMaximumBin());
    double backgroundAtCut = hbkg->GetBinContent(hsignif->GetMaximumBin());
    double dataAtCut       = hdta->GetBinContent(hsignif->GetMaximumBin());

    double topAtCut  = h_top->Integral(hsignif->GetMaximumBin(), 1001);
    double wwAtCut   = h_ww->Integral(hsignif->GetMaximumBin(), 1001);
    double dyttAtCut = h_dytt->Integral(hsignif->GetMaximumBin(), 1001);
    double vvAtCut   = h_vv->Integral(hsignif->GetMaximumBin(), 1001); 

    cout << "top: " << topAtCut << endl;
    cout << "ww: " << wwAtCut << endl;
    cout << "dytt: " << dyttAtCut << endl;
    cout << "vv: " << vvAtCut << endl;
    cout << "signal: " << signalAtCut << endl;
    cout << "background: " << backgroundAtCut << endl;

    TH1F hData("data", "data", 1,0,1);
    hData.SetBinContent(1,dataAtCut);
    TH1F hBkg("hBkg", "hBkg", 1,0,1);
    hBkg.SetBinContent(1,backgroundAtCut);
    TH1F hww("hww", "hww", 1,0,1);
    hww.SetBinContent(1,wwAtCut);
    TH1F htop("htop", "htop", 1,0,1);
    htop.SetBinContent(1,topAtCut);
    TH1F hdytt("hdytt", "hdytt", 1,0,1);
    hdytt.SetBinContent(1,max(dyttAtCut, 0.001));
    TH1F hvv("hvv", "hvv", 1,0,1);
    hvv.SetBinContent(1,vvAtCut);
    TH1F hs("hsig", "hsig", 1,0,1);
    hs.SetBinContent(1,signalAtCut);
    RooWorkspace * w = getWorkspace(&hww, &htop, &hdytt, &hvv, &hs, masses[i]);
    

    Double_t xq[5] = {0.025, 0.34, 0.5, 0.84, 0.975};  // position where to compute the quantiles in [0,1]
    Double_t yq[5];
    // auxiliary histogram to compute quantiles
    TH1F* hLimitToy = new TH1F("hLimitToy", "hLimitToy", 2000, -1000, 1000);
    

    // run 100 toys to get the 1 sigma and 2 sigma bands
    for (unsigned j =0; j < 100; ++j){
      //cout << "Toy>>>" << j << endl;
      int dataToy = rand.Poisson(backgroundAtCut);
      TH1F * hToy = new TH1F("hToy", "hToy", 1,0,1);
      hToy->SetBinContent(1,dataToy);

      // negative log likelyhood for each toy:
      // fit to the data
      std::pair<double, RooAbsReal*> limit = getLimit(hToy, w);
      hLimitToy->Fill(limit.first);
      delete hToy;
    }  
    // find quantiles on the toys to get the 1 and 2 sigma bands
    hLimitToy->GetQuantiles(5,yq,xq);

    // get the observed limit
    std::pair<double, RooAbsReal*> observedLimit = getLimit(&hData, w, false);

    
    
    // get the profile likelihood for expected background (just for comparison purposes)
    std::pair<double, RooAbsReal*> expectedLimit = getLimit(&hBkg, w, false);

    RooRealVar * mu = w->var("mu");

    RooPlot* frame1 = mu->frame(Bins(100),Range(0., 10.),Title("profileLL in mu")) ;
    expectedLimit.second->plotOn(frame1,LineColor(kBlack)) ;
    observedLimit.second->plotOn(frame1,LineColor(kBlue)) ;
    frame1->SetMinimum(0);
    frame1->SetMaximum(5);
    frame1->Draw();
    frame1->SetNameTitle(masses[i], masses[i]);
    outputfile.cd();
    frame1->Write();

   
    cout << "mass " << masses[i] << "-->"
         << "expected limit: " << yq[2]
         << ", 2 sigma down: " << yq[0]
         << ", 1 sigma down: " << yq[1]
         << ", 1 sigma up: " << yq[3]
         << ", 2 sigma up: " << yq[4]
         << ", observed: " << observedLimit.first << endl;
    expected[i] = yq[2];
    observed[i] = observedLimit.first;
    expectedUp1[i] = (yq[3]-yq[2]);
    expectedUp2[i] = (yq[4]-yq[2]);
    expectedDo1[i] = (yq[2]-yq[1]);
    expectedDo2[i] = (yq[2]-yq[0]);

    delete hLimitToy;

  }

  TCanvas * c = new TCanvas();
  c->cd();
  TGraphAsymmErrors * twoSigmaBandG = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, expectedDo2, expectedUp2);
  twoSigmaBandG->SetFillColor(kYellow);
  TGraphAsymmErrors * oneSigmaBandG = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, expectedDo1, expectedUp1);
  oneSigmaBandG->SetFillColor(kGreen);
  TGraphAsymmErrors * expectedG   = new TGraphAsymmErrors(nMasses, massesD, expected, 0, 0, 0, 0);
  expectedG->SetLineStyle(2);
  TGraphAsymmErrors * observedG   = new TGraphAsymmErrors(nMasses, massesD, observed, 0, 0, 0, 0);

  twoSigmaBandG->Draw("A3");
  oneSigmaBandG->Draw("3 same");
  expectedG->Draw("L same");
  observedG->Draw("L same");
  outputfile.cd();
  c->Write();
}
