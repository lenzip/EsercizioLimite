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
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestInverter.h"

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
  
/*   
  RooAbsPdf* model = w->pdf("model");
  const RooArgSet* constraints = w->set("constraints");
  const RooArgSet* nuisances = w->set("nuisances");
 


  RooAbsReal* nll = model->createNLL(*data, ExternalConstraints(*constraints) );
  RooMinuit m(*nll);
  m.setPrintLevel(-1000);
  m.migrad();
  //m.hesse();
  //m.migrad();
  double muhat       = mu->getVal();
  double delta_muhat = mu->getError();
  cout << "muhat " << muhat << "+/-" << delta_muhat << endl;
  // extract the profile likelihood in the signal strength parameter
  // to be used toset the limit
  RooAbsReal* profileLogLikelihood = nll->createProfile(*mu) ;
  // find the limit as the point at which the NLL crosses 2
  double cumulativeIntegral = 0.;
  int k = 0;
  for (unsigned int k = 0; k < 5000; ++k){
    mu->setVal(k*0.01);
    double t_mu_tilde = profileLogLikelihood->getVal() ;
    //cout << "mutilde " << t_mu_tilde << endl;
    // the asumptotic cumulative distribution for t_mu_tilde is:
    // https://arxiv.org/pdf/1007.1727.pdf eq (46)
    if (t_mu_tilde < muhat*muhat/(delta_muhat*delta_muhat))
      cumulativeIntegral = 2*ROOT::Math::normal_cdf(sqrt(t_mu_tilde)) -1.;
    else
      cumulativeIntegral = ROOT::Math::normal_cdf(sqrt(t_mu_tilde))+ROOT::Math::normal_cdf((t_mu_tilde + muhat*muhat/(delta_muhat*delta_muhat))/2/muhat/delta_muhat) - 1.;
    //cout << k*0.01 << " " << cumulativeIntegral << endl; 
    if (cumulativeIntegral > 0.95){
       limit = k*0.01;
       break;
    }  
  }
  if (doCleanup) delete nll;
  return limit;
*/  
  
    
  w->var("mu")->setVal(1.);
  ModelConfig config("config", w);
  config.SetParametersOfInterest("mu");
  config.SetPdf("constrained_model");
  config.SetNuisanceParameters("mu_ww,mu_top,mu_dytt,mu_vv");
  config.SetObservables("x");
  config.SetSnapshot(*(w->var("mu")));

  ModelConfig* config_bonly = (ModelConfig*) config.Clone();
  config_bonly->SetName("config_bonly");
  RooRealVar* poi = (RooRealVar*) config_bonly->GetParametersOfInterest()->first();
  poi->setVal(0);
  config_bonly->SetSnapshot( *poi  );
  poi->setVal(1);
  
  AsymptoticCalculator ac(*data, *config_bonly, config);
  ac.SetOneSided(true);
  // profile likelihood test statistics 
  ProfileLikelihoodTestStat profll(*config.GetPdf());
  // use one-sided profile likelihood
  profll.SetOneSided(true);
  AsymptoticCalculator::SetPrintLevel(-1);
  HypoTestInverter calc(ac);
  calc.SetTestStatistic(profll);
  // set confidence level (e.g. 95% upper limits)
  calc.SetConfidenceLevel(0.95); 
  int npoints = 100;  // number of points to scan
  double poimin = poi->getMin();
  double poimax = poi->getMax();

  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
  calc.SetFixedScan(npoints,poimin,poimax);	

  HypoTestInverterResult * r = calc.GetInterval();
  
  

  double limit2 = r->UpperLimit();	 
  delete config_bonly; 

  return limit2;
  
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
    // auxiliary histogram to compute quantiles
    
    // get the profile likelihood for expected background (just for comparison purposes)
    double expectedLimit = getLimit(th1toDataHist(w, h_bkg), w, false);

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
