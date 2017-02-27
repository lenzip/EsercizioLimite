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

RooWorkspace* getWorkspace(TH1* ww,
                           TH1* top,
                           TH1* dytt,
                           TH1* vv, 
                           TH1* sig){

  // get x axis histogram from data (assume it is the same for all)
  // x is a variable running over the bins that have been measured
  RooRealVar x("x", "x", ww->GetXaxis()->GetXmin(),  ww->GetXaxis()->GetXmax());

  // transform the normalizations of each background in a RooRealVar, fixed
  RooRealVar nu_ww  ("nu_ww",   "nu_ww",   ww->Integral());
  RooRealVar nu_top ("nu_top",  "nu_top",  top->Integral());
  RooRealVar nu_dytt("nu_dytt", "nu_dytt", dytt->Integral());
  RooRealVar nu_vv  ("nu_vv",   "nu_vv",   vv->Integral());
  // same for the signal
  RooRealVar nu_s   ("nu_s",    "nu_s",   sig->Integral());

  // Define a multiplier for each normalization, including the signal
  RooRealVar mu     ("mu",      "mu",      -10, 10); // this is the POI
  // the following are nuisances
  RooRealVar mu_ww  ("mu_ww",   "mu_ww",   0.1, 10); 
  RooRealVar mu_top ("mu_top",  "mu_top",  0.1, 10); 
  RooRealVar mu_dytt("mu_dytt", "mu_dytt", 0.1, 10); 
  RooRealVar mu_vv  ("mu_vv",   "mu_vv",   0.1, 10); 

  // we now write the normalization of each background and the signal as
  // the relevant mu_ times teh relevant nu_
  RooFormulaVar norm_s   ("norm_s",    "mu*nu_s",         RooArgSet(mu, nu_s));
  RooFormulaVar norm_ww  ("norm_ww",   "mu_ww*nu_ww",     RooArgSet(mu_ww, nu_ww));
  RooFormulaVar norm_top ("norm_top",  "mu_top*nu_top",   RooArgSet(mu_top, nu_top));
  RooFormulaVar norm_dytt("norm_dytt", "mu_dytt*nu_dytt", RooArgSet(mu_dytt, nu_dytt));
  RooFormulaVar norm_vv  ("norm_vv",   "mu_vv*nu_vv",     RooArgSet(mu_vv, nu_vv));

  // Define PDFs in x for each background and for the signal:
  // to do so we need to transform each input histogram
  // to a RooDataHist and then to a RooHistPdf
  RooDataHist h_s    ("h_s",    "h_s",    RooArgSet(x), sig);
  RooDataHist h_ww   ("h_ww",   "h_ww",   RooArgSet(x), ww);
  RooDataHist h_top  ("h_top",  "h_top",  RooArgSet(x), top);
  RooDataHist h_dytt ("h_dytt", "h_dytt", RooArgSet(x), dytt);
  RooDataHist h_vv   ("h_vv",   "h_vv",   RooArgSet(x), vv);

  RooHistPdf pdf_s    ("pdf_s",    "pdf_s",    RooArgSet(x), h_s);
  RooHistPdf pdf_ww   ("pdf_ww",   "pdf_ww",   RooArgSet(x), h_ww);
  RooHistPdf pdf_top  ("pdf_top",  "pdf_top",  RooArgSet(x), h_top);
  RooHistPdf pdf_dytt ("pdf_dytt", "pdf_dytt", RooArgSet(x), h_dytt);
  RooHistPdf pdf_vv   ("pdf_vv",   "pdf_vv",   RooArgSet(x), h_vv);

  // we now build a model.
  // data are represented by the sum of each of the above PDFs, 
  // each one multiplied by its integral 
  RooAddPdf model("model", "model", RooArgSet(pdf_s,  pdf_ww,  pdf_top,  pdf_dytt,  pdf_vv),
                                    RooArgSet(norm_s, norm_ww, norm_top, norm_dytt, norm_vv));

  // we now impose Gaussian constraints on the nuisances
  RooGaussian constraint_ww   ("constraint_ww",   "constraint_ww",   mu_ww,   RooConst(1.), RooConst(0.1));
  RooGaussian constraint_top  ("constraint_top",  "constraint_top",  mu_top,  RooConst(1.), RooConst(0.1));
  RooGaussian constraint_dytt ("constraint_dytt", "constraint_dytt", mu_dytt, RooConst(1.), RooConst(0.1));
  RooGaussian constraint_vv   ("constraint_vv",   "constraint_vv",   mu_vv,   RooConst(1.), RooConst(0.1));

  RooArgSet constraints(constraint_ww, constraint_top, constraint_dytt, constraint_vv);

  RooWorkspace* w = new RooWorkspace("w", "w");
  w->import(model);
  w->import(constraints);

  return w;

}
