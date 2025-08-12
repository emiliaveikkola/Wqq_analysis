#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <TLatex.h>
#include "../minitools/tdrstyle_mod22.C"
#include "TBox.h"

void SignalvsBackground() {
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TFile *file1 = new TFile("rootfiles/Winter24_TTtoLNu2Q.root", "READ");

  TH2D* h2_control_raw    = (TH2D*)file1->Get("h2_Wmass_ptpair_qg");
  TH2D* h2_control2_raw    = (TH2D*)file1->Get("h2_Wmass_ptpair_qq_others");
  TH2D* h2_background_raw = (TH2D*)file1->Get("h2_Wmass_ptpair_gg");
  TH2D* h2_signal_raw     = (TH2D*)file1->Get("h2_Wmass_ptpair_qq");

  h2_control_raw->Add(h2_control2_raw,1);

  TH2D* h2_control    = (TH2D*)h2_control_raw->Clone("h2_control_raw");
  TH2D* h2_background = (TH2D*)h2_background_raw->Clone("h2_background_raw");
  TH2D* h2_signal     = (TH2D*)h2_signal_raw->Clone("h2_signal_raw");
  
  TH1D *h1 = tdrHist("h1","Mass (GeV)",10,250.,"ptpair",15,229);
  //lumi_13TeV = "Run2";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,8,0,kSquare);

  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetLabelSize(0.045);
  h1->GetYaxis()->SetLabelSize(0.045);


// ——— Custom axis ranges ———
// User can modify these to zoom in on a region of interest
double xMin = 15.0, xMax = 229.0;
double yMin = 15.0, yMax = 249.0;
// apply to the 2D “signal” axis, which is used to draw the frame
h2_control->GetXaxis()->SetRangeUser(xMin, xMax);
h2_control->GetYaxis()->SetRangeUser(yMin, yMax);
// ————————————————

  // Get bin dimensions
  int nx = h2_control->GetNbinsX();
  int ny = h2_control->GetNbinsY();

  // Normalize each bin so its content is fraction of total events in that bin
  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double c = h2_control->GetBinContent(ix, iy);
      double b = h2_background->GetBinContent(ix, iy);
      double s = h2_signal->GetBinContent(ix, iy);
      double tot = c + b + s;
      if (tot > 0) {
        h2_control   ->SetBinContent(ix, iy, c/tot);
        h2_background->SetBinContent(ix, iy, b/tot);
        h2_signal    ->SetBinContent(ix, iy, s/tot);
      }
    }
  }

  // Find maxima to normalize alpha
  double maxSig = h2_signal->GetMaximum();
  double maxBkg = h2_background->GetMaximum();
  double maxCR  = h2_control->GetMaximum();

  auto paint = [&](TH2D* h, Int_t color, double max, double globalAlpha){
    for(int ix=1; ix<=nx; ++ix){
      double x1 = h->GetXaxis()->GetBinLowEdge(ix);
      double x2 = h->GetXaxis()->GetBinUpEdge(ix);
      if (x2 < xMin || x1 > xMax) continue;
      for(int iy=1; iy<=ny; ++iy){
        double y1 = h->GetYaxis()->GetBinLowEdge(iy);
        double y2 = h->GetYaxis()->GetBinUpEdge(iy);
        if (y2 < yMin || y1 > yMax) continue;
        double val = h->GetBinContent(ix,iy);
        if(val<=0) continue;
        double alpha = std::min(val/max*globalAlpha, 1.0);
        TBox *box = new TBox(x1,y1,x2,y2);
        box->SetFillColorAlpha(color, alpha);
        box->SetLineColor(color);
        box->Draw("same");
      }
    }
  };

  c1->SetLogz();

  // Paint control (green), background (blue), and signal (red)
  paint(h2_control,  kGreen+2, maxCR,  1);
  paint(h2_signal,    kRed+1,   maxSig, 1);
  paint(h2_background, kBlue-4,  maxBkg, 1);

  TH1D* dummyControl = new TH1D("dummyControl", "", 1, 0, 1);
  dummyControl->SetFillColorAlpha(kGreen+2, 1);
  dummyControl->SetLineColor(kGreen+2);
  
  TLegend *legend = tdrLeg(0.72, 0.75-0.015*3, 0.9, 0.9);
  TH1D* dummySignal = new TH1D("dummySignal", "", 1, 0, 1);
  dummySignal->SetFillColorAlpha(kRed+1, 0.7);
  dummySignal->SetLineColor(kRed+1);

  TH1D* dummyBackground = new TH1D("dummyBackground", "", 1, 0, 1);
  dummyBackground->SetFillColorAlpha(kBlue-4, 0.7);
  dummyBackground->SetLineColor(kBlue-4);



  legend->AddEntry(dummySignal,    "Signal",    "f");
  legend->AddEntry(dummyBackground,"Background","f");
  legend->AddEntry(dummyControl,   "Control",   "f");
  legend->SetTextSize(0.03);

  c1->RedrawAxis();

  TFile* fout = new TFile("Combined.root", "RECREATE");
  fout->cd();
  // write the canvas so that the overlaid colors are preserved
  c1->Write("c_SignalvsBackground");
  fout->Close();
  c1->Update();
  c1->SaveAs("pdf/SignalvsBackground.pdf");   // Save as PDF

  // === 1D Projections onto ptpair (X-axis) ===
  TH1D *h2 = tdrHist("h2","N frac",0,1,"ptpair",0,230);
  TCanvas *c2 = tdrCanvas("c2", h2, 8, 11, kSquare);

  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.045);
  h2->GetYaxis()->SetLabelSize(0.045);

  TH1D *projControlX   = (TH1D*)h2_control_raw->ProjectionX("projControlX");
  TH1D *projBackgroundX= (TH1D*)h2_background_raw->ProjectionX("projBackgroundX");
  TH1D *projSignalX    = (TH1D*)h2_signal_raw->ProjectionX("projSignalX");

  // Bin-by-bin normalization so control+background+signal = 1 per bin
  int nbinsX = projControlX->GetNbinsX();
  TH1D *hFracControl   = (TH1D*)projControlX   ->Clone("hFracControl");
  TH1D *hFracBackground= (TH1D*)projBackgroundX->Clone("hFracBackground");
  TH1D *hFracSignal    = (TH1D*)projSignalX    ->Clone("hFracSignal");
  for (int ib = 1; ib <= nbinsX; ++ib) {
    double c = projControlX   ->GetBinContent(ib);
    double b = projBackgroundX->GetBinContent(ib);
    double s = projSignalX    ->GetBinContent(ib);
    double tot = c + b + s;
    double fc = (tot>0 ? c/tot : 0);
    double fb = (tot>0 ? b/tot : 0);
    double fs = (tot>0 ? s/tot : 0);
    hFracControl   ->SetBinContent(ib, fc);
    hFracBackground->SetBinContent(ib, fb);
    hFracSignal    ->SetBinContent(ib, fs);
    // binomial errors
    double ec = (tot>0 ? sqrt(fc*(1-fc)/tot) : 0);
    double eb = (tot>0 ? sqrt(fb*(1-fb)/tot) : 0);
    double es = (tot>0 ? sqrt(fs*(1-fs)/tot) : 0);
    hFracControl   ->SetBinError(ib, ec);
    hFracBackground->SetBinError(ib, eb);
    hFracSignal    ->SetBinError(ib, es);
  }

  hFracControl   ->SetLineColor(kGreen+2);
  hFracBackground->SetLineColor(kBlue-4);
  hFracSignal    ->SetLineColor(kRed+1);
  hFracControl   ->SetLineWidth(2);
  hFracBackground->SetLineWidth(2);
  hFracSignal    ->SetLineWidth(2);

  hFracControl   ->Draw("E same");
  hFracBackground->Draw("E SAME");
  hFracSignal    ->Draw("E SAME");
  TLegend *leg2 = tdrLeg(0.62, 0.75-0.015*3, 0.9, 0.9);
  leg2->AddEntry(hFracSignal,    "Signal (qq)",    "l");
  leg2->AddEntry(hFracBackground,"Background (gg)","l");
  leg2->AddEntry(hFracControl,   "Control (qq* + qg)",   "l");
  leg2->SetTextSize(0.03);
  c2->RedrawAxis();
  c2->SaveAs("pdf/SignalvsBackground_projX.pdf");

  // === 1D Projections onto W mass (Y-axis) ===
  TH1D *h3 = tdrHist("h3","N",0,1, "Mass (GeV)",0,200);
  TCanvas *c3 = tdrCanvas("c3", h3, 8, 11, kSquare);

  h3->GetXaxis()->SetTitleSize(0.05);
  h3->GetYaxis()->SetTitleSize(0.05);
  h3->GetXaxis()->SetLabelSize(0.045);
  h3->GetYaxis()->SetLabelSize(0.045);

  TH1D *projControlY   = (TH1D*)h2_control_raw->ProjectionY("projControlY");
  TH1D *projBackgroundY= (TH1D*)h2_background_raw->ProjectionY("projBackgroundY");
  TH1D *projSignalY    = (TH1D*)h2_signal_raw->ProjectionY("projSignalY");

  // Bin-by-bin normalization so control+background+signal = 1 per W-mass bin
  int nbinsY = projControlY->GetNbinsX();
  TH1D *hFracControlY   = (TH1D*)projControlY   ->Clone("hFracControlY");
  TH1D *hFracBackgroundY= (TH1D*)projBackgroundY->Clone("hFracBackgroundY");
  TH1D *hFracSignalY    = (TH1D*)projSignalY    ->Clone("hFracSignalY");
  for (int ib = 1; ib <= nbinsY; ++ib) {
    double c = projControlY   ->GetBinContent(ib);
    double b = projBackgroundY->GetBinContent(ib);
    double s = projSignalY    ->GetBinContent(ib);
    double tot = c + b + s;
    double fc = (tot>0 ? c/tot : 0);
    double fb = (tot>0 ? b/tot : 0);
    double fs = (tot>0 ? s/tot : 0);
    hFracControlY   ->SetBinContent(ib, fc);
    hFracBackgroundY->SetBinContent(ib, fb);
    hFracSignalY    ->SetBinContent(ib, fs);
    // binomial errors
    double ec = (tot>0 ? sqrt(fc*(1-fc)/tot) : 0);
    double eb = (tot>0 ? sqrt(fb*(1-fb)/tot) : 0);
    double es = (tot>0 ? sqrt(fs*(1-fs)/tot) : 0);
    hFracControlY   ->SetBinError(ib, ec);
    hFracBackgroundY->SetBinError(ib, eb);
    hFracSignalY    ->SetBinError(ib, es);
  }

  // Draw normalized W-mass fractions
  hFracControlY   ->SetLineColor(kGreen+2);
  hFracBackgroundY->SetLineColor(kBlue-4);
  hFracSignalY    ->SetLineColor(kRed+1);
  hFracControlY   ->SetLineWidth(2);
  hFracBackgroundY->SetLineWidth(2);
  hFracSignalY    ->SetLineWidth(2);
  hFracControlY   ->Draw("E same");
  hFracBackgroundY->Draw("E SAME");
  hFracSignalY    ->Draw("E SAME");
  TLegend *leg3 = tdrLeg(0.62, 0.45-0.015*3, 0.9, 0.65);
  leg3->AddEntry(hFracSignalY,    "Signal (qq)",    "l");
  leg3->AddEntry(hFracBackgroundY,"Background (gg)","l");
  leg3->AddEntry(hFracControlY,   "Control (qq* + qg)",   "l");
  leg3->SetTextSize(0.03);
  c3->RedrawAxis();
  c3->SaveAs("pdf/SignalvsBackground_projY.pdf");
  // === Purity plots ===
  // Clone and compute bin-by-bin purity = category/(control+background+signal)
  TH2D *h2_puritySignal = (TH2D*)h2_signal_raw->Clone("h2_puritySignal");
  TH2D *h2_purityBkg    = (TH2D*)h2_background_raw->Clone("h2_purityBkg");

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double c = h2_control_raw   ->GetBinContent(ix, iy);
      double b = h2_background_raw->GetBinContent(ix, iy);
      double s = h2_signal_raw    ->GetBinContent(ix, iy);
      double tot = c + b + s;
      double ps = (tot>0 ? s/tot : 0);
      double pb = (tot>0 ? b/tot : 0);
      h2_puritySignal->SetBinContent(ix, iy, ps);
      h2_purityBkg   ->SetBinContent(ix, iy, pb);
    }
  }

  // Draw and save signal purity heatmap
  TH1D *h4 = tdrHist("h1","Mass (GeV)",10,250.,"ptpair",15,229);
  TCanvas *c4 = tdrCanvas("c4", h4, 8, 0, kSquare);
  h4->GetXaxis()->SetTitleSize(0.05);
  h4->GetYaxis()->SetTitleSize(0.05);
  h4->GetXaxis()->SetLabelSize(0.045);
  h4->GetYaxis()->SetLabelSize(0.045);

  gPad->SetRightMargin(1.01);

  h2_puritySignal->Draw("COLZ same");
  c4->RedrawAxis();
  c4->SaveAs("pdf/SignalPurity.pdf");

  // Draw and save background purity heatmap
  TH1D *h5 = tdrHist("h5","Mass (GeV)",10,250.,"ptpair",15,229);
  TCanvas *c5 = tdrCanvas("c5", h1, 8, 0, kSquare);

  h5->GetXaxis()->SetTitleSize(0.05);
  h5->GetYaxis()->SetTitleSize(0.05);
  h5->GetXaxis()->SetLabelSize(0.045);
  h5->GetYaxis()->SetLabelSize(0.045);

  gPad->SetRightMargin(1.01);

  h2_purityBkg->Draw("COLZ same");
  c5->RedrawAxis();
  c5->SaveAs("pdf/BackgroundPurity.pdf");
}