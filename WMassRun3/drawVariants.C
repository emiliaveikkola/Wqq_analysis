// Purpose: Draw Wqq variations for ISR, FSR and JER from Emilia with data
//          Perform TF1 fit to Data/MC ratio given these variations
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <TLatex.h>
#include "../minitools/tdrstyle_mod22.C"


void drawVariants() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open the output files
  TFile *file1 = new TFile("rootfiles/Winter24_TTtoLNu2Q_Wfinal.root", "READ");
  TFile *file2 = new TFile("rootfiles/Muon_Run2024CDE_Reprocessing_Wfinal.root", "READ");
  TFile *file3 = new TFile("rootfiles/Muon_Run2024FGHI_Prompt.root", "READ");


  if (!file1 || file1->IsZombie()) {
    std::cerr << "Error: Cannot open file 'Winter24_TTtoLNu2Q.root'" << std::endl;
    return;
}
if (!file2 || file2->IsZombie()) {
    std::cerr << "Error: Cannot open file 'Muon_Run2024CDE_Reprocessing.root'" << std::endl;
    return;
}
if (!file3 || file3->IsZombie()) {
    std::cerr << "Error: Cannot open file 'Muon_Run2024FGHI_Prompt.root'" << std::endl;
    return;
}

  curdir->cd();
  
  // 40,160 or 70,130
  TProfile* prof_W_ptpair_MC = (TProfile*)file1->Get("prof_W_ptpair");
  TProfile* prof_W_ptpair_DATA = (TProfile*)file2->Get("prof_W_ptpair");

  TProfile* prof_W_ptpair_inWindow_MC = (TProfile*)file1->Get("prof_W_inWindow_ptpair");
  TProfile* prof_W_ptpair_inWindow_DATA = (TProfile*)file2->Get("prof_W_inWindow_ptpair");

  TProfile* prof_W_ptpair_improved_MC = (TProfile*)file1->Get("prof_W_ptpair_improved");
  TProfile* prof_W_ptpair_improved_DATA = (TProfile*)file2->Get("prof_W_ptpair_improved");

  TProfile* prof_W_ptpair_inWindow_improved_MC = (TProfile*)file1->Get("prof_W_inWindow_ptpair_improved");
  TProfile* prof_W_ptpair_inWindow_improved_DATA = (TProfile*)file2->Get("prof_W_inWindow_ptpair_improved");

if (!prof_W_ptpair_inWindow_MC) {
    std::cerr << "Error: Histogram 'prof_W_ptpair_inWindow_MC' not found in file1" << std::endl;
    return;
}
if (!prof_W_ptpair_DATA) {
    std::cerr << "Error: Histogram 'prof_W_ptpair_DATA' not found in file2" << std::endl;
    return;
}

  // Convert profiles to TH1D via ProjectionX
  TH1D *h_W_MC              = prof_W_ptpair_MC->ProjectionX("h_W_MC");
  TH1D *h_W_DATA            = prof_W_ptpair_DATA->ProjectionX("h_W_DATA");
  TH1D *h_W_MC_improved     = prof_W_ptpair_improved_MC->ProjectionX("h_W_MC_improved");
  TH1D *h_W_DATA_improved   = prof_W_ptpair_improved_DATA->ProjectionX("h_W_DATA_improved");
  TH1D *h_W_MC_inWindow     = prof_W_ptpair_inWindow_MC->ProjectionX("h_W_MC_inWindow");
  TH1D *h_W_DATA_inWindow   = prof_W_ptpair_inWindow_DATA->ProjectionX("h_W_DATA_inWindow");
  TH1D *h_W_MC_inWindow_improved   = prof_W_ptpair_inWindow_improved_MC->ProjectionX("h_W_MC_inWindow_improved");
  TH1D *h_W_DATA_inWindow_improved = prof_W_ptpair_inWindow_improved_DATA->ProjectionX("h_W_DATA_inWindow_improved");

  TH1D *h1 = tdrHist("h1","Mass (GeV)",18.,180.,"#sqrt{p_{T1}*p_{T2}}",20.,130.);
  //lumi_13TeV = "Run2";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,4,8,kSquare);

  TH1D *h2 = tdrHist("h2","MC/DATA",0.7,1.1,"#sqrt{p_{T1}*p_{T2}}",20.,130);
  TCanvas *c2 = tdrCanvas("c2",h2,4,8,kSquare);

  TCanvas *c = tdrDiCanvas("c",h1,h2,8,11);

  c->cd(1);

  // Draw original MC (projection)
  tdrDraw(h_W_MC,            "hist", kNone,kBlue,kSolid,-1,1001,kNone);
  // Draw original Data
  tdrDraw(h_W_DATA,          "Pz",    kFullCircle, kAzure+7);
  // Draw improved MC (projection)
  tdrDraw(h_W_MC_improved,   "hist", kNone,kGreen+3,kSolid,-1,1001,kNone);
  // Draw improved Data
  tdrDraw(h_W_DATA_improved, "Pz",    kFullSquare, kSpring-5);
  // Draw inWindow MC (projection)
  tdrDraw(h_W_MC_inWindow,   "hist", kNone, kRed-4, kSolid, -1, 1001, kNone);
  // Draw inWindow Data
  tdrDraw(h_W_DATA_inWindow, "Pz",    kFullTriangleUp, kOrange-2);
  // Draw inWindow improved MC (projection)
  tdrDraw(h_W_MC_inWindow_improved,   "hist", kNone, kViolet-5, kSolid, -1, 1001, kNone);
  // Draw inWindow improved Data
  tdrDraw(h_W_DATA_inWindow_improved, "Pz",    kFullSquare,    kMagenta-4);
  // Add W mass PDG line at 80.4 GeV (unchanged)
  TLine *lineW = new TLine(20., 80.4, 130., 80.4);
  lineW->SetLineColor(kGray+1);
  lineW->SetLineStyle(2);
  lineW->Draw("SAME");
  // Legend
  TLegend *legend = tdrLeg(0.32, 0.75-0.015*8, 0.6, 0.9);
  legend->AddEntry(h_W_MC,            "MC",            "l");
  legend->AddEntry(h_W_DATA,          "Data (CDE)",          "ple");
  legend->AddEntry(h_W_MC_improved,   "MC improved",   "l");
  legend->AddEntry(h_W_DATA_improved, "Data improved", "ple");
  legend->AddEntry(h_W_MC_inWindow,   "MC inWindow",   "l");
  legend->AddEntry(h_W_DATA_inWindow, "Data inWindow", "ple");
  legend->AddEntry(h_W_MC_inWindow_improved,   "MC inWindow improved",   "l");
  legend->AddEntry(h_W_DATA_inWindow_improved, "Data inWindow improved", "ple");
  legend->SetTextSize(0.03);
  legend->Draw();

  c->cd(2);

  // Ratio Data/MC using TProfile projection
  TH1D *p1 = prof_W_ptpair_DATA->ProjectionX("p1");
  TH1D *p2 = prof_W_ptpair_MC->ProjectionX("p2");
  p1->Divide(p2);
  tdrDraw(p1, "Pz", kFullCircle, kAzure+7);

  // Ratio improved Data/MC via Projection
  TH1D *p1_impr = prof_W_ptpair_improved_DATA->ProjectionX("p1_impr");
  TH1D *p2_impr = prof_W_ptpair_improved_MC->ProjectionX("p2_impr");
  p1_impr->Divide(p2_impr);
  tdrDraw(p1_impr, "Pz", kFullSquare, kSpring-5);
  // Ratio inWindow Data/MC via Projection
  TH1D *p1_inWindow = prof_W_ptpair_inWindow_DATA->ProjectionX("p1_inWindow");
  TH1D *p2_inWindow = prof_W_ptpair_inWindow_MC->ProjectionX("p2_inWindow");
  p1_inWindow->Divide(p2_inWindow);
  tdrDraw(p1_inWindow, "Pz", kFullTriangleUp, kOrange-2);
  // Ratio inWindow improved Data/MC via Projection
  TH1D *p1_inWindow_impr = prof_W_ptpair_inWindow_improved_DATA->ProjectionX("p1_inWindow_impr");
  TH1D *p2_inWindow_impr = prof_W_ptpair_inWindow_improved_MC->ProjectionX("p2_inWindow_impr");
  p1_inWindow_impr->Divide(p2_inWindow_impr);
  tdrDraw(p1_inWindow_impr, "Pz", kFullTriangleUp, kMagenta-4);

  // Add a horizontal line at zero
  TLine *line = new TLine(20., 1, 130., 1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2); // Dashed line
  line->Draw("SAME");


  c->RedrawAxis();
  c->Update();
   // Save the canvas
  c->SaveAs("pdf/drawVariants_Wfinal.pdf");   // Save as PDF
  
} // drawWqqFromEmilia
