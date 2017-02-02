#include "TGraphAsymmErrors.h"

using namespace RooFit ;

void HWWLimitLikelihood(){
  TString fileCutName="cut.root";
  TString filePlotName="output.root";
  
  const unsigned int nMasses=1;

  TString masses[nMasses]={"300"};//, "400", "500", "600", "700", "800", "900", "1000"};
  double massesD[nMasses]={300};//, 400, 500, 600, 700, 800, 900, 1000};

  TFile filePlotIn(filePlotName);

  TH1F* h_data   = (TH1F*) filePlotIn.Get("h_data");
  TH1F* h_top    = (TH1F*) filePlotIn.Get("h_top");
  TH1F* h_ww     = (TH1F*) filePlotIn.Get("h_ww");
  TH1F* h_dytt   = (TH1F*) filePlotIn.Get("h_dytt");
  TH1F* h_vv     = (TH1F*) filePlotIn.Get("h_vv"); 
  
  TFile fileCutIn(fileCutName);

  TRandom3 rand;
  
  for (unsigned int i = 0; i < nMasses; ++i){
    TH1F* hdta = (TH1F*) fileCutIn.Get("data"+masses[i]);
    TH1F* hbkg = (TH1F*) fileCutIn.Get("bkg"+masses[i]);
    TH1F* hsig = (TH1F*) fileCutIn.Get("sig"+masses[i]);
    TH1F* hsignif = (TH1F*) fileCutIn.Get("significance"+masses[i]);


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
    cout << "background" << backgroundAtCut << endl;


    TH1F data("data", "data", 1,0,1);
    data.SetBinContent(1,dataAtCut);
    TH1F bkg("bkg", "bkg", 1,0,1);
    bkg.SetBinContent(1,backgroundAtCut);

    RooRealVar events("events", "events", 0, 1);
    RooDataHist bkgData("bkgData", "bkgData", RooArgList(events), &bkg);


    RooHistPdf count_ww("count_ww", "count_ww", events, bkgData);
    RooHistPdf count_top("count_top", "count_top",events, bkgData);
    RooHistPdf count_dytt("count_dytt", "count_dytt",events, bkgData);
    RooHistPdf count_vv("count_vv", "count_vv", events, bkgData);
    RooHistPdf count_signal("count_signal", "count_signal", events, bkgData);
   /* 
    RooRealVar top("top", "top", topAtCut, topAtCut/100., topAtCut*100);
    RooRealVar ww("ww", "ww", wwAtCut, wwAtCut/100., wwAtCut*100);
    RooRealVar dytt("dytt", "dytt", dyttAtCut, dyttAtCut/100., dyttAtCut*100);
    RooRealVar vv("vv", "vv", vvAtCut, vvAtCut/100., vvAtCut*100);
    RooRealVar signal("signal", "signal", signalAtCut/100., signalAtCut*100);

    RooGaussian ww_constraint("ww_constraint", "ww_constraint", ww, RooConst(wwAtCut), RooConst(wwAtCut*0.01));
    RooGaussian top_constraint("top_constraint", "top_constraint", top, RooConst(topAtCut), RooConst(topAtCut*0.01));
    RooGaussian dytt_constraint("dytt_constraint", "dytt_constraint", dytt, RooConst(dyttAtCut), RooConst(dyttAtCut*0.01));
    RooGaussian vv_constraint("vv_constraint", "vv_constraint", vv, RooConst(vvAtCut), RooConst(vvAtCut*0.01));

    RooAddPdf model("model", "model", RooArgList(count_ww, count_top, count_dytt, count_vv, count_signal), RooArgList(ww, top, dytt, vv, signal));
    //RooProdPdf modelc("modelc", "modelc", RooArgSet(model, ww_constraint, top_constraint, dytt_constraint, vv_constraint));

    model.fitTo(bkgData,ExternalConstraints(RooArgSet(ww_constraint, top_constraint, dytt_constraint, vv_constraint) ), Minos(signal));
  */
    RooRealVar top("top", "top", topAtCut);
    RooRealVar ww("ww", "ww", wwAtCut);
    RooRealVar dytt("dytt", "dytt", dyttAtCut);
    RooRealVar vv("vv", "vv", vvAtCut);
    RooRealVar signal("signal", "signal", signalAtCut);

    RooRealVar r_top("r_top", "r_top", 1, 0.1, 2);
    RooRealVar r_ww("r_ww", "r_ww", 1, 0.1, 2);
    RooRealVar r_dytt("r_dytt", "r_dytt", 1, 0.1, 2);
    RooRealVar r_vv("r_vv", "r_vv", 1, 0.1, 2);
    RooRealVar r_signal("r_signal", "r_signal", 1, 0.1, 2);

    RooFormulaVar top_norm("top_norm", "top*r_top", RooArgSet(top, r_top));
    RooFormulaVar ww_norm("ww_norm", "ww*r_ww", RooArgSet(ww, r_ww));
    RooFormulaVar dytt_norm("dytt_norm", "dytt*r_dytt", RooArgSet(dytt, r_dytt));
    RooFormulaVar vv_norm("vv_norm", "vv*r_vv", RooArgSet(vv, r_vv));
    RooFormulaVar signal_norm("signal_norm", "signal*r_signal", RooArgSet(signal, r_signal));

    RooLognormal ww_constraint("ww_constraint", "ww_constraint", r_ww, RooConst(1.), RooConst(1.10));
    RooLognormal top_constraint("top_constraint", "top_constraint", r_top, RooConst(1.), RooConst(1.10));
    RooLognormal dytt_constraint("dytt_constraint", "dytt_constraint", r_dytt, RooConst(1.), RooConst(1.10));
    RooLognormal vv_constraint("vv_constraint", "vv_constraint", r_vv, RooConst(1.), RooConst(1.10)); 
  
    RooAddPdf model("model", "model", RooArgList(count_ww, count_top, count_dytt, count_vv, count_signal), RooArgList(ww_norm, top_norm, dytt_norm, vv_norm, signal_norm));
    
    model.fitTo(bkgData,ExternalConstraints(RooArgSet(ww_constraint, top_constraint, dytt_constraint, vv_constraint) ), Minos(r_signal));


  }

}
