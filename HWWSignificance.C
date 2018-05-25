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
 
 
  // silence roofit
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
 
  // vectors to store limit values
  double observed[nMasses];

  // loop over masses
  for (unsigned int i = 0; i < nMasses; ++i){
    TH1F* h_sig = (TH1F*) filePlotIn.Get("h_signal"+masses[i]);
    RooWorkspace * w = getWorkspace(h_ww, h_top, h_dytt, h_vv, h_sig, masses[i]);
    
    // get the observed signif
    std::pair<double, RooAbsReal*> observedSignif = getSignif(h_data, w, false);

    
    
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
