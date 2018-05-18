#include "TMath.h"
#include "TH1F.h"
#include "TRandom3.h"

void limit(double b = 5, int nobs=5, double signalToTest=5){
  double shat_obs = TMath::Max(nobs-b, 0.);

  double L_obs    = TMath::Poisson(nobs, b+signalToTest);
  double Lhat_obs = TMath::Poisson(nobs, b+shat_obs);

  double l_obs = -2*log(L_obs/Lhat_obs);

  TH1F * hsb = new TH1F("hsb", "hsb", 1000000, 0, 100.);

  TRandom3 rand;
  for (unsigned int i =0; i < 1000000; ++i){
    int n = rand.Poisson(b+signalToTest);
    double shat = TMath::Max(n-b, 0.);
    double L    = TMath::Poisson(n, b+signalToTest);
    double Lhat = TMath::Poisson(n, b+shat);
    double l = -2*log(L/Lhat);
    hsb->Fill(l);
  }

  TH1F * hb = new TH1F("hb", "hb", 1000000, 0, 100.);	

  for (unsigned int i =0; i < 1000000; ++i){
    int n = rand.Poisson(b);
    double shat = TMath::Max(n-b, 0.);
    double L    = TMath::Poisson(n, b+signalToTest);
    double Lhat = TMath::Poisson(n, b+shat);
    double l = -2*log(L/Lhat);
    hb->Fill(l);
  }	

  TLine * line = new TLine(l_obs, 0, l_obs, 100000);
  line->SetLineWidth(3);
  hsb->SetLineColor(kRed);
  hsb->Draw();
  hb->SetLineColor(kBlue);
  hb->Draw("sames");
  line->Draw("same");
  TH1* hsbclone = (TH1*) hsb->Clone();
  hsbclone->GetXaxis()->SetRangeUser(l_obs, 100.);
  TH1* hbclone = (TH1*) hb->Clone();
  hbclone->GetXaxis()->SetRangeUser(l_obs, 100.);

  double CLsb=hsbclone->Integral()/hsb->Integral();
  double CLb=hbclone->Integral()/hb->Integral();
  double CLs = CLsb/CLb;

  cout << "lObs " << l_obs << " pvalue s " << CLsb << " pvalue b " << CLb  << " CLs: " << CLs << endl;
  
}
