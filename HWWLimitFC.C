#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/ProbFunc.h"
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
#include "RooStats/ModelConfig.h"

#include <iostream>

using namespace RooFit ;
using namespace RooStats ;
using namespace std ;


RooDataHist* th1toDataHist(RooWorkspace * w, const TH1* h){
  RooRealVar * x = w->var("x");
  RooDataHist* data = new RooDataHist("data", "data", *x, h);
  return data;
}

RooDataHist*  generateOneToy(RooWorkspace * w, const TH1* hbkg){
  w->loadSnapshot("bonly");
  RooAbsPdf* model = w->pdf("model");
  RooRealVar * x = w->var("x");
  RooDataHist* dataset = model->generateBinned(RooArgSet(*x), hbkg->Integral(), Extended());
  
  return dataset;
  
}


double getLimit(RooDataHist* data, RooWorkspace * w, bool doCleanup=true){ 
  double limit = -1;
  RooRealVar * x = w->var("x");
  RooRealVar * mu = w->var("mu");
  
   
  RooAbsPdf* model = w->pdf("model");
  const RooArgSet* constraints = w->set("constraints");
  const RooArgSet* nuisances = w->set("nuisances");
 

  RooAbsReal* nll = model->createNLL(*data, ExternalConstraints(*constraints) );
  RooMinuit m(*nll);
  m.setPrintLevel(-1000); //silence minuit
  m.migrad(); //minimize
  if (mu->getVal() < 0){ //if best fit is negative
    //redo the fit freezing mu to 0
    mu->setVal(0);
    mu->setConstant(true);
    m.migrad();
  }

  double minNLL = nll->getVal(); //value of the negative log likelihood at global minimum
    
  for (unsigned int k = 0; k < 5000; ++k){
    mu->setVal(k*0.01);  //set mu to a given value
    mu->setConstant(true); // make it constant, so the fit does not touch it
    m.migrad(); // minimize
    double minNLLmu = nll->getVal(); // net the minimum of NLL for that choice of mu
    double q_mu = 2*(minNLLmu - minNLL); // compute q_mu
    double CLsb = 1.-ROOT::Math::chisquared_cdf(q_mu,1.);  /// compute its p-value 
    //cout << "q_mu: " << q_mu << " CLsb " << CLsb << endl;
    if (CLsb < 0.05){
       limit = k*0.01;
       break;
    }  
  }
  if (doCleanup) delete nll;
  return limit;
  
}

void HWWLimitLikelihood(){
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
 
  TH1F* h_bkg = (TH1F*) h_top->Clone();
  h_bkg->Add(h_ww);
  h_bkg->Add(h_dytt);
  h_bkg->Add(h_vv);

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
    TH1F* h_sig = (TH1F*) filePlotIn.Get("h_signal"+masses[i]);

    RooWorkspace * w = getWorkspace(h_ww, h_top, h_dytt, h_vv, h_sig, masses[i]);
    

    Double_t xq[5] = {0.025, 0.34, 0.5, 0.84, 0.975};  // position where to compute the quantiles in [0,1]
    Double_t yq[5];
    
    // get the expected limit 
    double expectedLimit = getLimit(th1toDataHist(w, h_bkg), w, false);

    // auxiliary histogram to compute quantiles
    TH1F* hLimitToy = new TH1F("hLimitToy", "hLimitToy", 2000, -1000, 1000);
    // run 100 toys to get the 1 sigma and 2 sigma bands
    for (unsigned j =0; j < 100; ++j){
      
      // negative log likelyhood for each toy:
      // fit to the data
      RooDataHist* hToy = generateOneToy(w, h_bkg);
      double limit = getLimit(hToy, w);
      hLimitToy->Fill(limit);
      delete hToy;
    }  
    // find quantiles on the toys to get the 1 and 2 sigma bands
    hLimitToy->GetQuantiles(5,yq,xq);

    // get the observed limit
    double observedLimit = getLimit(th1toDataHist(w, h_data), w, false);

    
    

    cout << "mass " << masses[i] << "-->"
         << "expected limit: " << expectedLimit 
         << ", 2 sigma down: " << yq[0]
         << ", 1 sigma down: " << yq[1]
         << ", 1 sigma up: " << yq[3]
         << ", 2 sigma up: " << yq[4]
         << ", observed: " << observedLimit << endl;
    expected[i] = yq[2];
    observed[i] = observedLimit;
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
