void HWWYields(){
  double xmin = 0.;
  double xmax = 1000.;
  double nbins = 100;
  TString var = "mTi*(mTi<1000)+999.9*(mTi>=1000)";  // this is just mTi, bounded to 1000
  TString DirMC   = "/tmp/mc/";
  TString DirData = "/tmp/data/";
  TString outFileName="yields.root";
  
  TString Lumi = "*2.3" ;


  TString cut="";
  cut += "(   mll>50 && mth>60";                                                                 // cut SM Higgs away
  cut += " && std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>20";               // trigger safe lepton cuts
  cut += " && std_vector_lepton_pt[2]<10 ";                                            // third lepton veto against WZ
  cut += " && ptll > 30 ";                                                             // pt of the lepton pair, against DY->TauTau
  cut += " && (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -11*13)";// different flavor, most sensitive category
  cut += " && ( std_vector_jet_pt[0] < 20 || std_vector_jet_cmvav2[0] < -0.715 ) \
           && ( std_vector_jet_pt[1] < 20 || std_vector_jet_cmvav2[1] < -0.715 ) \
           && ( std_vector_jet_pt[2] < 20 || std_vector_jet_cmvav2[2] < -0.715 ) \
           && ( std_vector_jet_pt[3] < 20 || std_vector_jet_cmvav2[3] < -0.715 ) \
           && ( std_vector_jet_pt[4] < 20 || std_vector_jet_cmvav2[4] < -0.715 ) \
           && ( std_vector_jet_pt[5] < 20 || std_vector_jet_cmvav2[5] < -0.715 ) \
           && ( std_vector_jet_pt[6] < 20 || std_vector_jet_cmvav2[6] < -0.715 ) \
           && ( std_vector_jet_pt[7] < 20 || std_vector_jet_cmvav2[7] < -0.715 ) \
           && ( std_vector_jet_pt[8] < 20 || std_vector_jet_cmvav2[8] < -0.715 ) \
           && ( std_vector_jet_pt[9] < 20 || std_vector_jet_cmvav2[9] < -0.715 ) )";   // veto events with any b-tagged jet above 20 GeV

  TString weightMC="";
  weightMC += "*baseW";                               // normalization to 1/fb
  weightMC += "*metFilter";                           // filter to reject instrumental met background 
  weightMC += "*puW";                                 // weight to reproduce the PU distribution observed in data
  weightMC += "*bPogSF";                              // correction for different data vs MC btagging efficiency
  weightMC += "*effTrigW";                            // correction for different data vs MC trigger efficiency
  weightMC += "*std_vector_lepton_idisoW[0]";         // correction for different data vs MC lepton id and isolation efficiencyes (leading lepton)
  weightMC += "*std_vector_lepton_idisoW[1]";         // same ad above (subleading lepton)
  weightMC += "*std_vector_lepton_genmatched[0]";     // auxiliary information for data/mc scale factors (is the leading lepton a true lepton)
  weightMC += "*std_vector_lepton_genmatched[1]";     // auxiliary information for data/mc scale factors (is the subleading lepton a true lepton)

  TString weightData = "*trigger*metFilter";


  TFile output(outFileName, "recreate");
  output.cd();


  // DATA

  TChain * Data = new TChain("Data");
  Data->Add(DirData+"latino_Run2015C_16Dec2015_DoubleEG.root/latino") ;
  Data->Add(DirData+"latino_Run2015C_16Dec2015_DoubleMuon.root/latino") ;
  Data->Add(DirData+"latino_Run2015C_16Dec2015_MuonEG.root/latino") ;
  Data->Add(DirData+"latino_Run2015C_16Dec2015_SingleMuon.root/latino") ;
  Data->Add(DirData+"latino_Run2015D_16Dec2015_DoubleEG.root/latino") ;
  Data->Add(DirData+"latino_Run2015D_16Dec2015_DoubleMuon.root/latino") ;
  Data->Add(DirData+"latino_Run2015D_16Dec2015_MuonEG.root/latino") ;
  Data->Add(DirData+"latino_Run2015D_16Dec2015_SingleMuon.root/latino") ;
  TH1F* hdata = new TH1F("h_data","Data",nbins, xmin, xmax);
  h_data->Sumw2();
  Data->Draw(var+">> h_data",cut+weightData);
  h_data->Write();
  cout << "DATA: " << h_data->Integral() << endl;

  // WW
  TChain * WW =  new TChain("WW");
  WW->Add(DirMC+"/latino_WWTo2L2Nu.root/latino") ;
  WW->Add(DirMC+"/latino_GluGluWWTo2L2Nu_MCFM.root/latino");  
  TH1F* h_ww = new TH1F("h_ww","ww",nbins, xmin, xmax);
  h_ww->Sumw2();
  WW->Draw(var+">> h_ww",cut+weightMC+Lumi);
  h_ww->Write();
  cout << "WW: " << h_ww->Integral() << endl;

  // Top
  TChain * Top =  new TChain("Top");
  Top->Add(DirMC+"/latino_TTTo2L2Nu.root/latino");
  Top->Add(DirMC+"/latino_ST_tW_antitop.root/latino");
  Top->Add(DirMC+"/latino_ST_tW_top.root/latino");
  TH1F* h_top = new TH1F("h_top","Top",nbins, xmin, xmax);
  h_top->Sumw2();
  Top->Draw(var+">> h_top",cut+weightMC+Lumi);
  h_top->Write();
  cout << "Top: " << h_top->Integral() << endl;

  // DYtt
  TChain * DYtt =  new TChain("DYtt");
  DYtt->Add(DirMC+"/latino_DYJetsToLL_M-50-LO__part0.root/latino");
  DYtt->Add(DirMC+"/latino_DYJetsToLL_M-50-LO__part1.root/latino");
  TH1F* h_dytt = new TH1F("h_dytt","DYTT",nbins, xmin, xmax);
  h_dytt->Sumw2();
  DYtt->Draw(var+">> h_dytt",cut+weightMC+Lumi);
  h_dytt->Write();
  cout << "DY->tautau: "  << h_dytt->Integral() << endl;

  // VV = WZ ZZ  WGsToElNuM WGsToMuNu WGsToTauNu WG WZ2L2Q ZZ2L2Q ZGToLLuG
  TChain * VV =  new TChain("VV");
  VV->Add(DirMC+"/latino_WZTo3LNu.root/latino");
  VV->Add(DirMC+"/latino_ZZTo2L2Nu.root/latino");
  VV->Add(DirMC+"/latino_WZZ.root/latino");
  VV->Add(DirMC+"/latino_Wg_MADGRAPHMLM.root/latino");
  VV->Add(DirMC+"/latino_WgStarLNuEE.root/latino");
  VV->Add(DirMC+"/latino_WgStarLNuMuMu.root/latino");
  TH1F* hvv = new TH1F("h_vv","WZ/ZZ/W#gamma/Z#gamma",nbins, xmin, xmax);
  h_vv->Sumw2();
  VV->Draw(var+">> h_vv",cut+weightMC+Lumi);
  h_vv->Write();
  cout << "MultiBosons: " << h_vv->Integral() << endl;

  TString masses[7]={"300", "400", "500", "600", "700", "900", "1000"};
  for (unsigned int i = 0; i < 7; ++i){
    TChain * Signal = new TChain("Signal");
    Signal->Add(DirMC+"/latino_GluGluHToWWTo2L2Nu_M"+masses[i]+".root/latino");  
    Signal->Add(DirMC+"/latino_VBFHToWWTo2L2Nu_M"+masses[i]+".root/latino");  
    TH1F* hsig = new TH1F("h_signal"+masses[i],"Signal "+masses[i]+" GeV",nbins, xmin, xmax);
    hsig->Sumw2();
    Signal->Draw(var+">> h_signal"+masses[i], cut+weightMC+Lumi);
    cout << "Signal M"+masses[i]+": " << hsig->Integral() << endl;
    hsig->Write();
    delete Signal;
    delete hsig; 
  }

/*
  // VVV = WWGJets WZZJets ZZZJets WWZJets WWWJets TTWJets TTZJets TTWWJets TTGJets

  TLegend * leg = new TLegend(0.2, 0.6, 0.8, 0.9);
  leg->SetFillColor(0);
  leg->SetNColumns(3);
  //leg->AddEntry(p1, "phantom", "l");
  leg->AddEntry(hvv, "WZ,ZZ,VVV", "f");
  leg->AddEntry(hwjet, "W+jets", "f");
  leg->AddEntry(hvg, "V+#gamma,V+#gamma^{*}", "f");
  leg->AddEntry(htop, "Top", "f");
  leg->AddEntry(hdytt, "DY, DYTT", "f");
  leg->AddEntry(hww, "WW", "f");        
  leg->AddEntry(hh, "On shell H", "f");


   
  THStack* stack =  new THStack();
  stack->SetNameTitle("bkg", "bkg");
  TH1F* bkgByHand = (TH1F*) hwjet->Clone();
  stack->Add(hvv);
  stack->Add(hwjet);
  stack->Add(hvg);
  stack->Add(htop);
  //stack->Add(hdy);
  stack->Add(hdytt);
  stack->Add(hww);
  stack->Add(hh);
  
*/  

}
