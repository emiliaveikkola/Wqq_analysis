// Purpose: Draw Wqq variations for ISR, FSR and JER from Emilia with data
//          Perform TF1 fit to Data/MC ratio given these variations
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include "../minitools/tdrstyle_mod22.C"

TH1D *hisr_up(0), *hisr_dw(0);
TH1D *hfsr_up(0), *hfsr_dw(0);
TH1D *hjer_up(0), *hjer_dw(0);
TH1D *hist_up(0), *hist_dw(0);

TH1D *hisr_up_scaled(0), *hisr_dw_scaled(0);
TH1D *hfsr_up_scaled(0), *hfsr_dw_scaled(0);
TH1D *hjer_up_scaled(0), *hjer_dw_scaled(0);

Double_t hist_func(Double_t *x, Double_t *p) {

  assert(hist_up);
  assert(hist_dw);
  double mass = x[0];
  int i = hisr_up->FindBin(mass);
  double up = hist_up->GetBinContent(i);
  double dw = hist_dw->GetBinContent(i);
  // Line going through zero to up and dw
  double slope = (up-dw)/2.;

  return p[0]*slope;
} // hist_func

// Sum of ISR, FSR and JER variations
Double_t fit_func(Double_t *x, Double_t *p) {

  hist_up = hisr_up; hist_dw = hisr_dw;
  double isr = hist_func(x,&p[0]);
  
  hist_up = hfsr_up; hist_dw = hfsr_dw;
  double fsr = hist_func(x,&p[1]);

  hist_up = hjer_up; hist_dw = hjer_dw;
  double jer = hist_func(x,&p[2]);

  // Use co-opted data bins 98,99,100 to impose parameter chi2 penalties
  double mass = x[0];
  int i = hist_up->FindBin(mass);
  if (i==98)  { isr=fsr=jer=0; isr=p[0]; }
  if (i==99)  { isr=fsr=jer=0; fsr=p[1]; }
  if (i==100) { isr=fsr=jer=0; jer=p[2]; }
  
  return (1 + isr + fsr + jer);
}

Double_t fit_func_scaled(Double_t *x, Double_t *p) {
  // For ISR:
  hist_up = hisr_up_scaled;
  hist_dw = hisr_dw_scaled;
  double isr = hist_func(x, &p[0]);

  // For FSR:
  hist_up = hfsr_up_scaled;
  hist_dw = hfsr_dw_scaled;
  double fsr = hist_func(x, &p[1]);

  // For JER:
  hist_up = hjer_up_scaled;
  hist_dw = hjer_dw_scaled;
  double jer = hist_func(x, &p[2]);

  // Use co-opted data bins 98,99,100 to impose parameter chi2 penalties
  double mass = x[0];
  int i = hist_up->FindBin(mass);
  if (i==98)  { isr=fsr=jer=0; isr=p[0]; }
  if (i==99)  { isr=fsr=jer=0; fsr=p[1]; }
  if (i==100) { isr=fsr=jer=0; jer=p[2]; }

  // same final expression
  return (1 + isr + fsr + jer);
}

// Generic error calculation using differentials
void drawErrBand(TF1 *f1, TFitResultPtr &fp, double xmin, double xmax,
			double scale) {

  TMatrixDSym emat = fp->GetCovarianceMatrix();

  //TH1D *he = new TH1D(Form("he_%d",(int)fp),"",1,xmin,xmax);
  TGraphErrors *ge = new TGraphErrors(int(xmax-xmin)+1);
  TGraph *ge_up = new TGraph(int(xmax-xmin)+1);
  TGraph *ge_dw = new TGraph(int(xmax-xmin)+1);
  vector<double> dfdp(f1->GetNpar());
  //vector<double> dp(f1->GetNpar());
  for (int i = 0; i != ge->GetN(); ++i) {

    double x = xmin + (xmax-xmin)*i/(ge->GetN()-1);
    double y = f1->Eval(x);
    for (int j = 0; j != f1->GetNpar(); ++j) {

      double p = f1->GetParameter(j);
      double ep = f1->GetParError(j);
      double dp = 0.1*ep;
      f1->SetParameter(j, p+0.1*ep);
      double yup = f1->Eval(x);
      f1->SetParameter(j, p-0.1*ep);
      double ydw = f1->Eval(x);
      f1->SetParameter(j, p);
      dfdp[j] = (dp!=0 ? 0.5*(yup-ydw)/dp : 0);
      //dp[j] = ep;
    } // for j
    
    double e2(0);
    for (int j = 0; j != f1->GetNpar(); ++j) {
      for (int k = 0; k != f1->GetNpar(); ++k) {
	e2 += dfdp[j]*dfdp[k]*emat[j][k];
      } // for k
    } //  for j
    
    double err = sqrt(e2);
    if (scale <0) err *= sqrt(f1->GetChisquare()/f1->GetNDF());
    if (scale!=0) err *= fabs(scale);
    ge->SetPoint(i, x, y);
    ge->SetPointError(i, 0, err);
    ge_up->SetPoint(i, x, y+err);
    ge_dw->SetPoint(i, x, y-err);
  } // for i

  ge->SetFillColorAlpha(kBlue-9,0.5);
  ge->SetFillStyle(1001);
  ge->Draw("SAMEE3");
  ge_up->SetLineStyle(kDotted);
  ge_up->Draw("SAMEL");
  ge_dw->SetLineStyle(kDotted);
  ge_dw->Draw("SAMEL");

} // drawErr

void drawWqqFSR() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // --- Open original files as before ---
  TFile *fisr = new TFile("histograms_ISR.root","READ");
  assert(fisr && !fisr->IsZombie());

  TFile *ffsr = new TFile("histograms_FSR.root","READ");
  assert(ffsr && !ffsr->IsZombie());

  TFile *fjer = new TFile("histograms_JER.root","READ");
  assert(fjer && !fjer->IsZombie());

  TFile *ffsr_fsr = new TFile("histograms_FSR.root","READ");
  assert(ffsr_fsr && !ffsr_fsr->IsZombie());

  curdir->cd();
  
  // 40,160 or 70,130
  hisr_dw = (TH1D*)fisr->Get("RatioHist_ISR070"); assert(hisr_dw);
  hisr_up = (TH1D*)fisr->Get("RatioHist_ISR130"); assert(hisr_up);

  // 0990,1010 or 0995,1005 
  hfsr_dw = (TH1D*)ffsr->Get("RatioHist_FSR0995"); assert(hfsr_dw);
  hfsr_up = (TH1D*)ffsr->Get("RatioHist_FSR1005"); assert(hfsr_up);

  // 80,120 or 90,110
  hjer_dw = (TH1D*)fjer->Get("RatioHist_JER090"); assert(hjer_dw);
  hjer_up = (TH1D*)fjer->Get("RatioHist_JER110"); assert(hjer_up);

  TH1D *hdata = (TH1D*)ffsr_fsr->Get("RatioHistogram_down2"); assert(hdata);


  hisr_up_scaled  = (TH1D*)hisr_up->Clone("hisr_up_scaled");
  hisr_dw_scaled  = (TH1D*)hisr_dw->Clone("hisr_dw_scaled");

  hfsr_up_scaled  = (TH1D*)hfsr_up->Clone("hfsr_up_scaled");
  hfsr_dw_scaled  = (TH1D*)hfsr_dw->Clone("hfsr_dw_scaled");

  hjer_up_scaled  = (TH1D*)hjer_up->Clone("hjer_up_scaled");
  hjer_dw_scaled  = (TH1D*)hjer_dw->Clone("hjer_dw_scaled");
  

  // Co-opt bins 98,99,100 at [97,100] to hold parameter constraints
  // (bit of a hack, but I'm too lazy to write new chi2 function)
  hdata->SetBinContent(98,1);  hdata->SetBinError(98,1);
  hdata->SetBinContent(99,1);  hdata->SetBinError(99,1);
  hdata->SetBinContent(100,1); hdata->SetBinError(100,1);

  // TH1D *h = tdrHist("h","Ratio to MC",0.68,1.36,"Mass (GeV)",65.,100.);
  // [65,100] => [66,97] to avoid edge effects in JER
  TH1D *h1 = tdrHist("h1","Ratio to MC",0.86,1.21,"Mass (GeV)",66.,97.);
  //lumi_13TeV = "Run2";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,4,8,kSquare);

  TH1D *h2 = tdrHist("h2","Data/fit -1 (%)",-2,2,"Mass (GeV)",66.,97.);
  TCanvas *c2 = tdrCanvas("c2",h2,4,8,kSquare);

  TCanvas *c = tdrDiCanvas("c",h1,h2,4,11);

  c->cd(1);
  

  //TF1 *fit = new TF1("fit",fit_func,66,97,3);
  TF1 *fit1 = new TF1("fit1",fit_func,70,100,3); // use co-opted range 97-100 GeV
  fit1->SetParName(0,"isr_up (+30%)");
  fit1->SetParName(1,"fsr_up (+0.5%)");
  fit1->SetParName(2,"jer_up (+10%)");
  fit1->FixParameter(0,0);
  //fit1->FixParameter(1,0);
  //fit1->FixParameter(2,0);
  //fit1->SetParameters(1,-1,0);
  fit1->SetParLimits(2,-1,0);
  hdata->Fit(fit1,"RN");
  fit1->SetRange(66,100);

  double chi2_1 = fit1->GetChisquare();
  int    ndf_1  = fit1->GetNDF();
  double scaleFactor = (ndf_1 > 0) ? TMath::Sqrt(chi2_1 / ndf_1) : 1.0;

  auto scaleBinErrors = [&](TH1D* h, double sf)
  {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      double oldErr = h->GetBinError(i);
      double newErr = oldErr * sf;
      h->SetBinError(i, newErr);
    }
  };

  scaleBinErrors(hisr_up_scaled, scaleFactor);
  scaleBinErrors(hisr_dw_scaled, scaleFactor);
  scaleBinErrors(hfsr_up_scaled, scaleFactor);
  scaleBinErrors(hfsr_dw_scaled, scaleFactor);
  scaleBinErrors(hjer_up_scaled, scaleFactor);
  scaleBinErrors(hjer_dw_scaled, scaleFactor);

  // Scale factor for errors
  double errorScale = 1.0;
  if (ndf_1 > 0) {
    errorScale = TMath::Sqrt(chi2_1 / (double)ndf_1);
  }

  // ------------------------------------------------------------------
  // 2) Create a new histogram with scaled errors
  // ------------------------------------------------------------------
  TH1D *hdata_scaled = (TH1D*)hdata->Clone("hdata_scaled");
  for (int i = 1; i <= hdata_scaled->GetNbinsX(); ++i)
  {
    double oldErr = hdata_scaled->GetBinError(i);
    double newErr = oldErr * errorScale;
    hdata_scaled->SetBinError(i, newErr);
  }

    TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(66,1,97,1);
  
  tdrDraw(hisr_up_scaled,"Pz",kOpenSquare,kBlue-9);
  tdrDraw(hisr_dw_scaled,"Pz",kFullSquare,kBlue);

  tdrDraw(hfsr_up_scaled,"Pz",kOpenCircle,kRed-9);
  tdrDraw(hfsr_dw_scaled,"Pz",kFullCircle,kRed);

  tdrDraw(hjer_up_scaled,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(hjer_dw_scaled,"Pz",kOpenDiamond,kGreen-9);

  tdrDraw(hdata_scaled,"Pz",kFullCircle,kBlack);

  // ------------------------------------------------------------------
  // 3) Second fit with new histogram (bigger errors)
  // ------------------------------------------------------------------
  TF1 *fit2 = new TF1("fit2", fit_func_scaled, 70, 100, 3);
  fit2->SetParName(0, "isr_up (+30%)");
  fit2->SetParName(1, "fsr_up (+0.5%)");
  fit2->SetParName(2, "jer_up (+10%)");
  fit2->FixParameter(0,0);
  //fit2->FixParameter(1,0);
  //fit2->FixParameter(2,0);
  //fit2->SetParameters(1, -1, 0);
  //fit2->SetParLimits(2, 0, 100);
  fit2->SetRange(66,100);
  fit2->SetParLimits(2,-1,0);
  fit2->SetLineColor(kBlack);
  fit2->SetLineWidth(3);
  fit2->Draw("SAME");

  // Important: We now fit hdata_scaled, which has bigger errors
  hdata_scaled->Fit(fit2, "RN");

    // Get fit results
  double chi2 = fit2->GetChisquare();
  int ndf = fit2->GetNDF();
  double par0 = fit2->GetParameter(0);
  double par1 = fit2->GetParameter(1);
  double par2 = fit2->GetParameter(2);
  double err0 = fit2->GetParError(0);
  double err1 = fit2->GetParError(1);
  double err2 = fit2->GetParError(2);

  // Add chi2 and parameter values to the canvas
  TLatex latex;
  latex.SetTextSize(0.045);
  latex.SetNDC(); // Use normalized device coordinates for positioning
  latex.DrawLatex(0.32, 0.85, Form("#chi^{2}/NDF = %.2f / %d", chi2, ndf));
  latex.DrawLatex(0.32, 0.80, Form("ISR: %.6f #pm %.6f", par0, err0));
  latex.DrawLatex(0.32, 0.75, Form("FSR: %.6f #pm %.6f", par1, err1));
  latex.DrawLatex(0.32, 0.70, Form("JER: %.6f #pm %.6f", par2, err2));


  // Add the legend after all drawings
  TLegend *legend = tdrLeg(0.65, 0.88-0.04*8, 0.95, 0.88); 

  legend->AddEntry(hisr_up_scaled, "ISR Up", "P");
  legend->AddEntry(hisr_dw_scaled, "ISR Down", "P");
  legend->AddEntry(hfsr_up_scaled, "FSR Up", "P");
  legend->AddEntry(hfsr_dw_scaled, "FSR Down", "P");
  legend->AddEntry(hjer_up_scaled, "JER Up", "P");
  legend->AddEntry(hjer_dw_scaled, "JER Down", "P");
  legend->AddEntry(hdata_scaled, "FSR up x2", "P");
  legend->AddEntry(fit2, "Fit", "L");
  //legend->AddEntry(fit_mod, "Fit Mod", "L");

  legend->Draw("SAME");

  c->cd(2);

  // Create a histogram for the fit function values
  TH1D *h_fit = (TH1D*)hdata_scaled->Clone("h_fit");
  h_fit->Reset(); // Clear any existing content

  int nBins = h_fit->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
      double x = h_fit->GetBinCenter(i);
      double fit_value = fit2->Eval(x);
      h_fit->SetBinContent(i, fit_value);
      h_fit->SetBinError(i, 0); // Assuming negligible error on fit values for now
  }

  // Clone the data histogram to create the ratio histogram
  TH1D *h_ratio = (TH1D*)hdata_scaled->Clone("h_ratio");

  // ---- Get the fit results ----
  Double_t chi2_fit = fit2->GetChisquare(); // total chi-squared
  Int_t    ndf_fit  = fit2->GetNDF();       // number of degrees of freedom

  // Divide data histogram by fit histogram
  h_ratio->Divide(h_fit);

  // Transform the ratio histogram to (data/fit - 1) * 100
  for (int i = 1; i <= nBins; ++i) {
      double ratio = h_ratio->GetBinContent(i);
      double ratio_error = h_ratio->GetBinError(i)*TMath::Sqrt(chi2_fit/ndf_fit);//h_ratio->GetBinError(i);

      // Compute percentage difference
      double percent_diff = (ratio - 1.0) * 100.0;
      double percent_error = ratio_error * 100.0;

      h_ratio->SetBinContent(i, percent_diff);
      h_ratio->SetBinError(i, percent_error);
  }

    // Add a horizontal line at zero
  TLine *line = new TLine(66., 0.0, 97., 0.0);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2); // Dashed line
  line->Draw("SAME");

  h_ratio->Draw("same");
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(kFullCircle);
  h_ratio->SetLineWidth(2);
  h_ratio->SetLineColor(kBlack);



   // Save the canvas
  c->SaveAs("pdf/drawWqqFSR_down2.pdf");   // Save as PDF
  
} // drawWqqFromEmilia
