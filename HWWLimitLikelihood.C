#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooConstVar.h"

#include <iostream>

using namespace RooFit ;
using namespace std ;

std::pair<double, RooAbsReal*> getLimit(RooDataHist& hist, 
                RooAbsPdf& model, 
                RooArgSet& constraints, 
                RooRealVar& variable,
                bool doCleanup=true){ 
  RooAbsReal* nll = model.createNLL(hist, ExternalConstraints(constraints) );
  RooMinuit m(*nll);
  m.setPrintLevel(-1000);
  m.migrad();
  m.hesse();
  m.migrad();
  // extract the profile likelihood in the signal strength parameter
  // to be used toset the limit
  RooAbsReal* profileLogLikelihood = nll->createProfile(variable) ;
  double limit = -1;
  // find the limit as the point at which the NLL crosses 2
  for (unsigned int k=0; k < 1000; ++k){
    variable.setVal(k*0.01);
    double deltaNLL = profileLogLikelihood->getVal();
    if (deltaNLL > 2){
      limit = k*0.01;
      break;
    }
  } 
  if (doCleanup) delete nll;
  return std::make_pair(limit, profileLogLikelihood);
}

void HWWLimitLikelihood(){
  TString fileCutName="cut.root";
  TString filePlotName="yields.root";
  
  TFile outputfile("limitProfileLikelihood.root", "recreate");
  const unsigned int nMasses=8;

  TString masses[nMasses]={"300", "400", "500", "600", "700", "800", "900", "1000"};
  double massesD[nMasses]={300, 400, 500, 600, 700, 800, 900, 1000};

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

    // convert single-binned histograms (representing event count) 
    // to RooDataHist to be used with roofit
    TH1F hData("data", "data", 1,0,1);
    hData.SetBinContent(1,dataAtCut);
    TH1F hBkg("hBkg", "hBkg", 1,0,1);
    hBkg.SetBinContent(1,backgroundAtCut);

    // this is a "dummy" variable used to just count events
    RooRealVar events("events", "events", 0, 1);

    RooDataHist bkgData("bkgData", "bkgData", RooArgList(events), &hBkg);
    RooDataHist dhData("dhData", "dhData", RooArgList(events), &hData);


    // these represent single-binned PDFs (a workaround to make roofit work with event counts rather than PDFs)
    RooHistPdf count_ww("count_ww", "count_ww", events, bkgData);
    RooHistPdf count_top("count_top", "count_top",events, bkgData);
    RooHistPdf count_dytt("count_dytt", "count_dytt",events, bkgData);
    RooHistPdf count_vv("count_vv", "count_vv", events, bkgData);
    RooHistPdf count_signal("count_signal", "count_signal", events, bkgData);
    
    // the expected normalizations of each process
    RooRealVar top("top", "top", topAtCut);
    RooRealVar ww("ww", "ww", wwAtCut);
    RooRealVar dytt("dytt", "dytt", dyttAtCut);
    RooRealVar vv("vv", "vv", vvAtCut);
    RooRealVar signal("signal", "signal", signalAtCut);

    // a multiplicative factor for each process
    RooRealVar r_top("r_top", "r_top", 1, 0.1, 2);
    RooRealVar r_ww("r_ww", "r_ww", 1, 0.1, 2);
    RooRealVar r_dytt("r_dytt", "r_dytt", 1, 0.1, 2);
    RooRealVar r_vv("r_vv", "r_vv", 1, 0.1, 2);
    RooRealVar r_signal("r_signal", "r_signal", 1, -20, 20);

    // the actual normalization for each process
    RooFormulaVar top_norm("top_norm", "top*r_top", RooArgSet(top, r_top));
    RooFormulaVar ww_norm("ww_norm", "ww*r_ww", RooArgSet(ww, r_ww));
    RooFormulaVar dytt_norm("dytt_norm", "dytt*r_dytt", RooArgSet(dytt, r_dytt));
    RooFormulaVar vv_norm("vv_norm", "vv*r_vv", RooArgSet(vv, r_vv));
    RooFormulaVar signal_norm("signal_norm", "signal*r_signal", RooArgSet(signal, r_signal));

    // gaussian constraints on the normalizations of each background process
    // These are nuisance parameters
    RooGaussian ww_constraint("ww_constraint", "ww_constraint", r_ww, RooConst(1.), RooConst(0.30));
    RooGaussian top_constraint("top_constraint", "top_constraint", r_top, RooConst(1.), RooConst(0.30));
    RooGaussian dytt_constraint("dytt_constraint", "dytt_constraint", r_dytt, RooConst(1.), RooConst(0.30));
    RooGaussian vv_constraint("vv_constraint", "vv_constraint", r_vv, RooConst(1.), RooConst(0.30)); 
  
    // model representing the data as signal + background
    RooAddPdf model("model", "model", RooArgList(count_ww, count_top, count_dytt, count_vv, count_signal), RooArgList(ww_norm, top_norm, dytt_norm, vv_norm, signal_norm));

    RooArgSet constraints(ww_constraint, top_constraint, dytt_constraint, vv_constraint);

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
      RooDataHist dhToy("dhToy", "dhToy", RooArgList(events), hToy);

      // negative log likelyhood for each toy:
      // fit to the data
      std::pair<double, RooAbsReal*> limit = getLimit(dhToy, model, constraints, r_signal);
      hLimitToy->Fill(limit.first);
      delete hToy;
    }  
    // find quantiles on the toys to get the 1 and 2 sigma bands
    hLimitToy->GetQuantiles(5,yq,xq);

    // get the observed limit
    std::pair<double, RooAbsReal*> observedLimit = getLimit(dhData, model, constraints, r_signal, false);

    
    
    // get the profile likelihood for expected background (just for comparison purposes)
    std::pair<double, RooAbsReal*> expectedLimit = getLimit(bkgData, model, constraints, r_signal, false);

    RooPlot* frame1 = r_signal.frame(Bins(100),Range(0., 10.),Title("profileLL in r_signal")) ;
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

  }


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
  outputfile.cd();
  c->Write();

}
