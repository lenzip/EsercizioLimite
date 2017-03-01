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
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"

#include <iostream>

using namespace RooFit ;
using namespace RooStats ;
using namespace std ;


std::pair<double, RooAbsReal*> getSignif(TH1* hdata, RooWorkspace * w, bool doCleanup=true){ 
  double q0 = 0;
  RooRealVar * x = w->var("x");
  RooDataHist data("data", "data", *x, hdata);  
  RooRealVar * mu = w->var("mu");
  
  RooAbsPdf* model = w->pdf("model");
  const RooArgSet* constraints = w->set("constraints");
 


  RooAbsReal* nll = model->createNLL(data, ExternalConstraints(*constraints) );
  RooMinuit m(*nll);
  m.setPrintLevel(-1000);
  m.migrad();
  m.hesse();
  m.migrad();
  double muhat = mu->getVal();
  RooAbsReal* profileLogLikelihood = nll->createProfile(*mu) ;
  if (muhat>0){
    mu->setVal(0);
     q0 = 2*profileLogLikelihood->getVal();
  }
  
  
  if (doCleanup) delete nll;
  return std::make_pair(sqrt(q0), profileLogLikelihood);
}

void HWWSignificance(){
  TString fileCutName="cut.root";
  TString filePlotName="yields.root";
  
  TFile outputfile("significanceProfileLikelihood.root", "recreate");
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
    
    // get the observed signif
    std::pair<double, RooAbsReal*> observedSignif = getSignif(&hData, w, false);

    
    
    cout << "mass " << masses[i] << "-->" << observedSignif.first << " sigmas" << endl;
    observed[i] = observedSignif.first;  
  }

  TCanvas * c = new TCanvas();
  c->cd();
  TGraphAsymmErrors * observedG   = new TGraphAsymmErrors(nMasses, massesD, observed, 0, 0, 0, 0);

  observedG->Draw("AL");
  outputfile.cd();
  c->Write();



}
