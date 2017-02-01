void HWWStack(){
  TString fileInName="output.root";
  TString mass="1000";
  double scaleSignal=100.;
  TFile* fileIn =  new TFile(fileInName);


  TH1F* h_data   = (TH1F*) fileIn->Get("h_data");
  h_data->SetLineColor(1);
  h_data->SetLineWidth(2);
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(20);
  TH1F* h_top    = (TH1F*) fileIn->Get("h_top");
  h_top->SetFillColor(kAzure);
  h_top->SetLineColor(kAzure+3);
  TH1F* h_ww     = (TH1F*) fileIn->Get("h_ww");
  h_ww->SetFillColor(kYellow);
  h_ww->SetLineColor(kYellow+3);
  TH1F* h_dytt   = (TH1F*) fileIn->Get("h_dytt");
  h_dytt->SetFillColor(kGreen);
  h_dytt->SetLineColor(kGreen+3);
  TH1F* h_vv     = (TH1F*) fileIn->Get("h_vv");
  h_vv->SetFillColor(kViolet);
  h_vv->SetLineColor(kViolet+3);
  TH1F* h_signal = (TH1F*) fileIn->Get("h_signal"+mass);
  h_signal->Scale(scaleSignal);
  h_signal->SetLineColor(1);

  THStack* stack =  new THStack();
  stack->SetNameTitle("stack", "stack");
  stack->Add(h_ww, "hist");
  stack->Add(h_top, "hist");
  stack->Add(h_dytt, "hist");
  stack->Add(h_vv, "hist");
  //stack->Add(h_signal, "hist");

  TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetNColumns(2);
  leg->AddEntry(h_vv, "WZ,ZZ,VVV", "f");
  leg->AddEntry(h_top, "Top", "f");
  leg->AddEntry(h_dytt, "DY, DYTT", "f");
  leg->AddEntry(h_ww, "WW", "f");        
  leg->AddEntry(h_signal, "signal "+mass+" GeV", "l");
  leg->AddEntry(h_data, "Data", "lp"); 


  stack->Draw();
  h_signal->Draw("HISTsame");
  h_data->Draw("same");
  leg->Draw("same");
}
