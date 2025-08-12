#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <algorithm>
#include "../minitools/tdrstyle_mod22.C"

void MCvsDATA_pt_Run2() {
    TFile *file = new TFile("../output/output_MCSCALEDRun2_tagprobe_chi2_FSR_before.root", "READ");
    TH3D* h3PtFlavorPairs_DATAMC_MC = (TH3D*)file->Get("h3PtFlavorPairs_DATAMC_CHI_1.000_1.000");
      // Project Z for different indices and draw
      TH1D* h3all = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3all", 1, 3, 1, 3);

      TH1D* h3gencsall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3gencsall", 1, 1, 1, 3);
      TH1D* h3genudall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genudall", 2, 2, 1, 3);
      TH1D* h3genxall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genxall", 3, 3, 1, 3);

      TH1D* h3tagcsall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagcsall", 1, 3, 1, 1);
      TH1D* h3tagudall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagudall", 1, 3, 2, 2);
      TH1D* h3tagxall = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagxall", 1, 3, 3, 3);

      TH1D* h3gencstagcs = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagcs", 1, 1, 1, 1);
      TH1D* h3gencstagud = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagud", 1, 1, 2, 2);
      TH1D* h3gencstagx = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagx", 1, 1, 3, 3);
      
      TH1D* h3genudtagcs = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagcs", 2, 2, 1, 1);
      TH1D* h3genudtagud = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagud", 2, 2, 2, 2);
      TH1D* h3genudtagx = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagx", 2, 2, 3, 3);
      
      TH1D* h3genxtagcs = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagcs", 3, 3, 1, 1);
      TH1D* h3genxtagud = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagud", 3, 3, 2, 2);
      TH1D* h3genxtagx = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagx", 3, 3, 3, 3);

cout << "here" << flush;

      TFile *file3 = new TFile("../output/output_MCSCALEDRun2_tagprobe_chi2_FSR_after.root", "READ");
      TH3D* h3PtFlavorPairs_DATAMC_SCALEDMASS = (TH3D*)file3->Get("h3PtFlavorPairs_DATAMC_CHI_1.000_1.000");

      // Project Z for different indices and draw
      TH1D* h3all_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3all_scaled", 1, 3, 1, 3);

      TH1D* h3gencsall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3gencsall_scaled", 1, 1, 1, 3);
      TH1D* h3genudall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genudall_scaled", 2, 2, 1, 3);
      TH1D* h3genxall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genxall_scaled", 3, 3, 1, 3);

      TH1D* h3tagcsall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3tagcsall_scaled", 1, 3, 1, 1);
      TH1D* h3tagudall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3tagudall_scaled", 1, 3, 2, 2);
      TH1D* h3tagxall_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3tagxall_scaled", 1, 3, 3, 3);

      TH1D* h3gencstagcs_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3gencstagcs_scaled", 1, 1, 1, 1);
      TH1D* h3gencstagud_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3gencstagud_scaled", 1, 1, 2, 2);
      TH1D* h3gencstagx_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3gencstagx_scaled", 1, 1, 3, 3);
      
      TH1D* h3genudtagcs_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genudtagcs_scaled", 2, 2, 1, 1);
      TH1D* h3genudtagud_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genudtagud_scaled", 2, 2, 2, 2);
      TH1D* h3genudtagx_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genudtagx_scaled", 2, 2, 3, 3);
      
      TH1D* h3genxtagcs_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genxtagcs_scaled", 3, 3, 1, 1);
      TH1D* h3genxtagud_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genxtagud_scaled", 3, 3, 2, 2);
      TH1D* h3genxtagx_scaled = h3PtFlavorPairs_DATAMC_SCALEDMASS->ProjectionZ("h3genxtagx_scaled", 3, 3, 3, 3);


      TFile *file2 = new TFile("../output/output_DATARun2_tagprobe_before.root", "READ");
      TH3D* h3PtFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3PtFlavorPairs_DATAMC");

      TH1D* h3all_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

      TH1D* h3tagcsall_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
      TH1D* h3tagudall_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
      TH1D* h3tagxall_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);

      TH1D* h3corr = new TH1D("h3corr",";Mass;N",200,0,200);

      TH1D* h3corrdata = new TH1D("h3corrdata",";Mass;N",200,0,200);


      // Cloning and renaming histograms
TH1D* h3all_mass = (TH1D*)h3all->Clone("h3all_mass");
TH1D* h3all_mass2 = (TH1D*)h3all->Clone("h3all_mass2");

TH1D* h3gencsall_mass = (TH1D*)h3gencsall->Clone("h3gencsall_mass");
TH1D* h3genudall_mass = (TH1D*)h3genudall->Clone("h3genudall_mass");
TH1D* h3genxall_mass = (TH1D*)h3genxall->Clone("h3genxall_mass");

TH1D* h3tagcsall_mass = (TH1D*)h3tagcsall->Clone("h3tagcsall_mass");
TH1D* h3tagudall_mass = (TH1D*)h3tagudall->Clone("h3tagudall_mass");
TH1D* h3tagxall_mass = (TH1D*)h3tagxall->Clone("h3tagxall_mass");

TH1D* h3gencstagcs_mass = (TH1D*)h3gencstagcs->Clone("h3gencstagcs_mass");
TH1D* h3gencstagud_mass = (TH1D*)h3gencstagud->Clone("h3gencstagud_mass");
TH1D* h3gencstagx_mass = (TH1D*)h3gencstagx->Clone("h3gencstagx_mass");

TH1D* h3genudtagcs_mass = (TH1D*)h3genudtagcs->Clone("h3genudtagcs_mass");
TH1D* h3genudtagud_mass = (TH1D*)h3genudtagud->Clone("h3genudtagud_mass");
TH1D* h3genudtagx_mass = (TH1D*)h3genudtagx->Clone("h3genudtagx_mass");

TH1D* h3genxtagcs_mass = (TH1D*)h3genxtagcs->Clone("h3genxtagcs_mass");
TH1D* h3genxtagud_mass = (TH1D*)h3genxtagud->Clone("h3genxtagud_mass");
TH1D* h3genxtagx_mass = (TH1D*)h3genxtagx->Clone("h3genxtagx_mass");


    
    h3gencstagcs->Scale(1./h3all->Integral());
    h3gencstagud->Scale(1./h3all->Integral());
    h3gencstagx->Scale(1./h3all->Integral());
    
    h3genudtagcs->Scale(1./h3all->Integral());
    h3genudtagud->Scale(1./h3all->Integral());
    h3genudtagx->Scale(1./h3all->Integral());

    h3genxtagcs->Scale(1./h3all->Integral());
    h3genxtagud->Scale(1./h3all->Integral());
    h3genxtagx->Scale(1./h3all->Integral());

    h3gencsall->Scale(1./h3all->Integral());
    h3genudall->Scale(1./h3all->Integral());
    h3genxall->Scale(1./h3all->Integral());

    h3tagcsall->Scale(1./h3all->Integral());
    h3tagudall->Scale(1./h3all->Integral());
    h3tagxall->Scale(1./h3all->Integral());

    h3all->Scale(1./h3all->Integral());

    h3gencstagcs_scaled->Scale(1./h3all_scaled->Integral());
    h3gencstagud_scaled->Scale(1./h3all_scaled->Integral());
    h3gencstagx_scaled->Scale(1./h3all_scaled->Integral());
    
    h3genudtagcs_scaled->Scale(1./h3all_scaled->Integral());
    h3genudtagud_scaled->Scale(1./h3all_scaled->Integral());
    h3genudtagx_scaled->Scale(1./h3all_scaled->Integral());

    h3genxtagcs_scaled->Scale(1./h3all_scaled->Integral());
    h3genxtagud_scaled->Scale(1./h3all_scaled->Integral());
    h3genxtagx_scaled->Scale(1./h3all_scaled->Integral());

    h3gencsall_scaled->Scale(1./h3all_scaled->Integral());
    h3genudall_scaled->Scale(1./h3all_scaled->Integral());
    h3genxall_scaled->Scale(1./h3all_scaled->Integral());

    h3tagcsall_scaled->Scale(1./h3all_scaled->Integral());
    h3tagudall_scaled->Scale(1./h3all_scaled->Integral());
    h3tagxall_scaled->Scale(1./h3all_scaled->Integral());

    h3all_scaled->Scale(1./h3all_scaled->Integral());

    h3tagcsall_data->Scale(1./h3all_data->Integral());
    h3tagudall_data->Scale(1./h3all_data->Integral());
    h3tagxall_data->Scale(1./h3all_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

    // Scaling histograms
/*
h3gencstagcs_scaled->Scale(0.789996);
h3gencstagx_scaled->Scale(1.09397);
h3gencstagud_scaled->Scale(1.);
h3genudtagcs_scaled->Scale(1.);
h3genudtagx_scaled->Scale(1.07117);
h3genudtagud_scaled->Scale(0.782763);
h3genxtagcs_scaled->Scale(0.817742);
h3genxtagx_scaled->Scale(1.07278);
h3genxtagud_scaled->Scale(0.812845);

// Reset h3all_scaled to zero
h3all_scaled->Reset();
// Summing the scaled histograms
h3all_scaled->Add(h3gencstagcs_scaled);
h3all_scaled->Add(h3gencstagud_scaled);
h3all_scaled->Add(h3gencstagx_scaled);
h3all_scaled->Add(h3genudtagcs_scaled);
h3all_scaled->Add(h3genudtagud_scaled);
h3all_scaled->Add(h3genudtagx_scaled);
h3all_scaled->Add(h3genxtagcs_scaled);
h3all_scaled->Add(h3genxtagud_scaled);
h3all_scaled->Add(h3genxtagx_scaled);

h3tagcsall_scaled->Reset();
h3tagcsall_scaled->Add(h3gencstagcs_scaled);
h3tagcsall_scaled->Add(h3genudtagcs_scaled);
h3tagcsall_scaled->Add(h3genxtagcs_scaled);

h3tagudall_scaled->Reset();
h3tagudall_scaled->Add(h3gencstagud_scaled);
h3tagudall_scaled->Add(h3genudtagud_scaled);
h3tagudall_scaled->Add(h3genxtagud_scaled);

h3tagxall_scaled->Reset();
h3tagxall_scaled->Add(h3gencstagx_scaled);
h3tagxall_scaled->Add(h3genudtagx_scaled);
h3tagxall_scaled->Add(h3genxtagx_scaled);

h3gencsall_scaled->Reset();
h3gencsall_scaled->Add(h3gencstagcs_scaled);
h3gencsall_scaled->Add(h3gencstagud_scaled);
h3gencsall_scaled->Add(h3gencstagx_scaled);

h3genudall_scaled->Reset();
h3genudall_scaled->Add(h3genudtagcs_scaled);
h3genudall_scaled->Add(h3genudtagud_scaled);
h3genudall_scaled->Add(h3genudtagx_scaled);

h3genxall_scaled->Reset();
h3genxall_scaled->Add(h3genxtagcs_scaled);
h3genxall_scaled->Add(h3genxtagud_scaled);
h3genxall_scaled->Add(h3genxtagx_scaled);

// Final scaling

double scale_factor2 = 1.0 / h3all_scaled->Integral();

h3gencstagx_scaled->Scale(scale_factor2);
h3genudtagx_scaled->Scale(scale_factor2);
h3genxtagx_scaled->Scale(scale_factor2);

h3gencstagud_scaled->Scale(scale_factor2);
h3genudtagud_scaled->Scale(scale_factor2);
h3genxtagud_scaled->Scale(scale_factor2);

h3gencstagcs_scaled->Scale(scale_factor2);
h3genudtagcs_scaled->Scale(scale_factor2);
h3genxtagcs_scaled->Scale(scale_factor2);

h3gencsall_scaled->Scale(scale_factor2);
h3genudall_scaled->Scale(scale_factor2);
h3genxall_scaled->Scale(scale_factor2);

h3tagcsall_scaled->Scale(scale_factor2);
h3tagudall_scaled->Scale(scale_factor2);
h3tagxall_scaled->Scale(scale_factor2);

h3all_scaled->Scale(scale_factor2);






h3gencstagx_mass->Scale(1);
h3genudtagx_mass->Scale(1);
h3genxtagx_mass->Scale(1);

h3gencstagud_mass->Scale(0.9);
h3genudtagud_mass->Scale(0.8);
h3genxtagud_mass->Scale(0.8);

h3gencstagcs_mass->Scale(0.9);
h3genudtagcs_mass->Scale(0.8);
h3genxtagcs_mass->Scale(0.8);


// Reset h3all_mass to zero
h3all_mass->Reset();
// Summing the scaled histograms
h3all_mass->Add(h3gencstagcs_mass);
h3all_mass->Add(h3gencstagud_mass);
h3all_mass->Add(h3gencstagx_mass);
h3all_mass->Add(h3genudtagcs_mass);
h3all_mass->Add(h3genudtagud_mass);
h3all_mass->Add(h3genudtagx_mass);
h3all_mass->Add(h3genxtagcs_mass);
h3all_mass->Add(h3genxtagud_mass);
h3all_mass->Add(h3genxtagx_mass);

h3tagcsall_mass->Reset();
h3tagcsall_mass->Add(h3gencstagcs_mass);
h3tagcsall_mass->Add(h3genudtagcs_mass);
h3tagcsall_mass->Add(h3genxtagcs_mass);

h3tagudall_mass->Reset();
h3tagudall_mass->Add(h3gencstagud_mass);
h3tagudall_mass->Add(h3genudtagud_mass);
h3tagudall_mass->Add(h3genxtagud_mass);

h3tagxall_mass->Reset();
h3tagxall_mass->Add(h3gencstagx_mass);
h3tagxall_mass->Add(h3genudtagx_mass);
h3tagxall_mass->Add(h3genxtagx_mass);

h3gencsall_mass->Reset();
h3gencsall_mass->Add(h3gencstagcs_mass);
h3gencsall_mass->Add(h3gencstagud_mass);
h3gencsall_mass->Add(h3gencstagx_mass);

h3genudall_mass->Reset();
h3genudall_mass->Add(h3genudtagcs_mass);
h3genudall_mass->Add(h3genudtagud_mass);
h3genudall_mass->Add(h3genudtagx_mass);

h3genxall_mass->Reset();
h3genxall_mass->Add(h3genxtagcs_mass);
h3genxall_mass->Add(h3genxtagud_mass);
h3genxall_mass->Add(h3genxtagx_mass);

// Final scaling

double scale_factor2_mass = 1.0 / h3all_mass->Integral();

h3gencstagx_mass->Scale(scale_factor2);
h3genudtagx_mass->Scale(scale_factor2);
h3genxtagx_mass->Scale(scale_factor2);

h3gencstagud_mass->Scale(scale_factor2);
h3genudtagud_mass->Scale(scale_factor2);
h3genxtagud_mass->Scale(scale_factor2);

h3gencstagcs_mass->Scale(scale_factor2);
h3genudtagcs_mass->Scale(scale_factor2);
h3genxtagcs_mass->Scale(scale_factor2);

h3gencsall_mass->Scale(scale_factor2);
h3genudall_mass->Scale(scale_factor2);
h3genxall_mass->Scale(scale_factor2);

h3tagcsall_mass->Scale(scale_factor2);
h3tagudall_mass->Scale(scale_factor2);
h3tagxall_mass->Scale(scale_factor2);

h3all_mass->Scale(scale_factor2);

*/
    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

        // Create a canvas
    TCanvas *c3 = new TCanvas("c3", "Canvas with 3x3 Grid", 2560, 1440);

    // Divide the canvas into a 3x3 grid
    c3->Divide(4, 4);
    gStyle->SetOptStat(0);

    for (int i = 1; i <= 16; ++i) {
    c3->cd(i);
    
    // Set minimal margins inside each pad
    gPad->SetLeftMargin(0.12);  // Leave enough space for the y-axis title
    gPad->SetRightMargin(0.02);  // Minimal right margin
    gPad->SetTopMargin(0.02);  // Minimal top margin
    gPad->SetBottomMargin(0.12);
    gPad->Update();
    }
    c3->cd(1);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagx_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagx_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagx->GetYaxis()->SetTitle("N frac");

    // Set the x-axis range
    //h3gencstagx->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3gencstagx->GetYaxis()->SetRangeUser(0, 0.0101);
    h3gencstagx->SetTitle("");

    TLegend *leg1= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg1->AddEntry(h3gencstagx, "FSR, no pgof", "FPLE");
    leg1->AddEntry(h3gencstagx_scaled, "FSR, pgof", "FPLE");
    //leg1->AddEntry(h3gencstagx_mass, "MC (init.+JES)", "FPLE");

     // Update the canvas to reflect the changes
    TLatex *tex1 = new TLatex();
    tex1->SetNDC(); tex1->SetTextSize(0.045);
    tex1->DrawLatex(0.17,0.85,"gencs,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(2);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);


        // Draw the histogram
    tdrDraw(h3genudtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagx_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagx_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagx->GetYaxis()->SetTitle("N frac");

    // Set the x-axis range
    h3genudtagx->GetYaxis()->SetRangeUser(0, 0.0135); // Adjust the range as needed
    h3genudtagx->SetTitle("");

    TLegend *leg2= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg2->AddEntry(h3genudtagx, "FSR, no pgof", "FPLE");
    leg2->AddEntry(h3genudtagx_scaled, "FSR, pgof", "FPLE");
    //leg2->AddEntry(h3genudtagx_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex2 = new TLatex();
    tex2->SetNDC(); tex2->SetTextSize(0.045);
    tex2->DrawLatex(0.17,0.85,"genud,tagx");
    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(3);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagx_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagx_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagx->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genxtagx->GetYaxis()->SetRangeUser(0, 0.014); // Adjust the range as needed
    h3genxtagx->SetTitle("");

    TLegend *leg3= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg3->AddEntry(h3genxtagx, "FSR, no pgof", "FPLE");
    leg3->AddEntry(h3genxtagx_scaled, "FSR, pgof", "FPLE");
    //leg3->AddEntry(h3genxtagx_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex3 = new TLatex();
    tex3->SetNDC(); tex3->SetTextSize(0.045);
    tex3->DrawLatex(0.17,0.85,"genx,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(4);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagxall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagxall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagxall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagxall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    //h3tagxall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.041); // Adjust the range as needed
    h3tagxall->SetTitle("");

    TLegend *leg4 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg4->AddEntry(h3tagxall, "FSR, no pgof", "FPLE");
    leg4->AddEntry(h3tagxall_scaled, "FSR, pgof", "FPLE");
    //leg4->AddEntry(h3tagxall_mass, "MC (init.+JES)", "FPLE");
    leg4->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex4 = new TLatex();
    tex4->SetNDC(); tex4->SetTextSize(0.045);
    tex4->DrawLatex(0.17,0.85,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();



    c3->cd(5);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagud_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagud_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagud->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3gencstagud->GetYaxis()->SetRangeUser(0, 0.00125); // Adjust the range as needed
    h3gencstagud->SetTitle("");

    TLegend *leg5= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg5->AddEntry(h3gencstagud, "FSR, no pgof", "FPLE");
    leg5->AddEntry(h3gencstagud_scaled, "FSR, pgof", "FPLE");
    //leg5->AddEntry(h3gencstagud_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex5 = new TLatex();
    tex5->SetNDC(); tex5->SetTextSize(0.045);
    tex5->DrawLatex(0.17,0.85,"gencs,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(6);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagud_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagud_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagud->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genudtagud->GetYaxis()->SetRangeUser(0, 0.0045); // Adjust the range as needed
    h3genudtagud->SetTitle("");

    TLatex *tex6 = new TLatex();
    tex6->SetNDC(); tex6->SetTextSize(0.045);
    tex6->DrawLatex(0.17,0.85,"genud,tagud");

    TLegend *leg6= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg6->AddEntry(h3genudtagud, "FSR, no pgof", "FPLE");
    leg6->AddEntry(h3genudtagud_scaled, "FSR, pgof", "FPLE");
    //leg6->AddEntry(h3genudtagud_mass, "MC (init.+JES)", "FPLE");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(7);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagud_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagud_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagud->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genxtagud->GetYaxis()->SetRangeUser(0, 0.004); // Adjust the range as needed
    h3genxtagud->SetTitle("");

    TLegend *leg7 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg7->AddEntry(h3genxtagud, "FSR, no pgof", "FPLE");
    leg7->AddEntry(h3genxtagud_scaled, "FSR, pgof", "FPLE");
    //leg7->AddEntry(h3genxtagud_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex7 = new TLatex();
    tex7->SetNDC(); tex7->SetTextSize(0.045);
    tex7->DrawLatex(0.17,0.85,"genx,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(8);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagudall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagudall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagudall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagudall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagudall->GetYaxis()->SetRangeUser(0, 0.0095); // Adjust the range as needed
    h3tagudall->SetTitle("");

    TLegend *leg8 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg8->AddEntry(h3tagudall, "FSR, no pgof", "FPLE");
    leg8->AddEntry(h3tagudall_scaled, "FSR, pgof", "FPLE");
    //leg8->AddEntry(h3tagudall_mass, "MC (init.+JES)", "FPLE");
    leg8->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex8 = new TLatex();
    tex8->SetNDC(); tex8->SetTextSize(0.045);
    tex8->DrawLatex(0.17,0.85,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(9);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagcs_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagcs_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagcs->GetYaxis()->SetTitle("N frac");

    h3gencstagcs->GetYaxis()->SetRangeUser(0, 0.004); // Adjust the range as needed
    h3gencstagcs->SetTitle("");

    TLegend *leg9 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg9->AddEntry(h3gencstagcs, "FSR, no pgof", "FPLE");
    leg9->AddEntry(h3gencstagcs_scaled, "FSR, pgof", "FPLE");
    //leg9->AddEntry(h3gencstagcs_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex9 = new TLatex();
    tex9->SetNDC(); tex9->SetTextSize(0.045);
    tex9->DrawLatex(0.17,0.85,"gencs,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(10);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagcs_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagcs_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagcs->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    //h3genudtagcs->SetMinimum(0);
    h3genudtagcs->GetYaxis()->SetRangeUser(0, 0.00028); // Adjust the range as needed
    h3genudtagcs->SetTitle("");
    h3genudtagcs->GetYaxis()->SetNoExponent();


    TLegend *leg10 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg10->AddEntry(h3genudtagcs, "FSR, no pgof", "FPLE");
    leg10->AddEntry(h3genudtagcs_scaled, "FSR, pgof", "FPLE");
    //leg10->AddEntry(h3genudtagcs_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex10 = new TLatex();
    tex10->SetNDC(); tex10->SetTextSize(0.045);
    tex10->DrawLatex(0.17,0.85,"genud,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(11);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagcs_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagcs_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagcs->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genxtagcs->GetYaxis()->SetRangeUser(0, 0.0016); // Adjust the range as needed
    h3genxtagcs->SetTitle("");

    TLegend *leg11 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg11->AddEntry(h3genxtagcs, "FSR, no pgof", "FPLE");
    leg11->AddEntry(h3genxtagcs_scaled, "FSR, pgof", "FPLE");
    //leg11->AddEntry(h3genxtagcs_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex11 = new TLatex();
    tex11->SetNDC(); tex11->SetTextSize(0.045);
    tex11->DrawLatex(0.17,0.85,"genx,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(12);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagcsall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagcsall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagcsall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.0055); // Adjust the range as needed
    h3tagcsall->SetTitle("");

    TLegend *leg12 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg12->AddEntry(h3tagcsall, "FSR, no pgof", "FPLE");
    leg12->AddEntry(h3tagcsall_scaled, "FSR, pgof", "FPLE");
    //leg12->AddEntry(h3tagcsall_mass, "MC (init.+JES)", "FPLE");
    leg12->AddEntry(h3tagcsall_data, "DATA", "PLE");

    TLatex *tex12 = new TLatex();
    tex12->SetNDC(); tex12->SetTextSize(0.045);
    tex12->DrawLatex(0.17,0.85,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(13);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencsall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencsall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencsall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencsall_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3gencsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3gencsall->GetYaxis()->SetRangeUser(0, 0.0135); // Adjust the range as needed
    h3gencsall->SetTitle("");

    TLegend *leg13= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg13->AddEntry(h3gencsall, "FSR, no pgof", "FPLE");
    leg13->AddEntry(h3gencsall_scaled, "FSR, pgof", "FPLE");
    //leg13->AddEntry(h3gencsall_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex13 = new TLatex();
    tex13->SetNDC(); tex13->SetTextSize(0.045);
    tex13->DrawLatex(0.17,0.85,"gencs,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(14);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudall_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genudall->GetYaxis()->SetRangeUser(0, 0.017); // Adjust the range as needed
    h3genudall->SetTitle("");

        TLegend *leg14= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg14->AddEntry(h3genudall, "FSR, no pgof", "FPLE");
    leg14->AddEntry(h3genudall_scaled, "FSR, pgof", "FPLE");
    //leg14->AddEntry(h3genudall_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex14 = new TLatex();
    tex14->SetNDC(); tex14->SetTextSize(0.045);
    tex14->DrawLatex(0.17,0.85,"genud,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(15);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxall_mass->SetFillColorAlpha(kSpring-5,0.25);

    h3genxall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3genxall->GetYaxis()->SetRangeUser(0, 0.017); // Adjust the range as needed
    h3genxall->SetTitle("");

    TLegend *leg15= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg15->AddEntry(h3genxall, "FSR, no pgof", "FPLE");
    leg15->AddEntry(h3genxall_scaled, "FSR, pgof", "FPLE");
    //leg15->AddEntry(h3genxall_mass, "MC (init.+JES)", "FPLE");

    TLatex *tex15 = new TLatex();
    tex15->SetNDC(); tex15->SetTextSize(0.045);
    tex15->DrawLatex(0.17,0.85,"genx,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(16);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3all_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3all_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3all_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3all->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3all->GetYaxis()->SetRangeUser(0, 0.048); // Adjust the range as needed
    h3all->SetTitle("");

    TLegend *leg16 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg16->AddEntry(h3all, "FSR, no pgof", "FPLE");
    leg16->AddEntry(h3all_scaled, "FSR, pgof", "FPLE");
    //leg16->AddEntry(h3all_mass, "MC (init.+JES)", "FPLE");
    leg16->AddEntry(h3all_data, "DATA", "PLE");
    

    TLatex *tex16 = new TLatex();
    tex16->SetNDC(); tex16->SetTextSize(0.045);
    tex16->DrawLatex(0.17,0.85,"all");

h3gencstagx->GetXaxis()->SetTitleSize(0.06);  
h3gencstagx->GetXaxis()->SetLabelSize(0.05);  
h3gencstagx->GetYaxis()->SetTitleSize(0.06);  
h3gencstagx->GetYaxis()->SetLabelSize(0.05);  

h3genudtagx->GetXaxis()->SetTitleSize(0.06);  
h3genudtagx->GetXaxis()->SetLabelSize(0.05);  
h3genudtagx->GetYaxis()->SetTitleSize(0.06);  
h3genudtagx->GetYaxis()->SetLabelSize(0.05);  

h3genxtagx->GetXaxis()->SetTitleSize(0.06);  
h3genxtagx->GetXaxis()->SetLabelSize(0.05);  
h3genxtagx->GetYaxis()->SetTitleSize(0.06);  
h3genxtagx->GetYaxis()->SetLabelSize(0.05);  

h3tagxall->GetXaxis()->SetTitleSize(0.06);  
h3tagxall->GetXaxis()->SetLabelSize(0.05);  
h3tagxall->GetYaxis()->SetTitleSize(0.06);  
h3tagxall->GetYaxis()->SetLabelSize(0.05);  

h3gencstagud->GetXaxis()->SetTitleSize(0.06);  
h3gencstagud->GetXaxis()->SetLabelSize(0.05);  
h3gencstagud->GetYaxis()->SetTitleSize(0.06);  
h3gencstagud->GetYaxis()->SetLabelSize(0.05);  

h3genudtagud->GetXaxis()->SetTitleSize(0.06);  
h3genudtagud->GetXaxis()->SetLabelSize(0.05);  
h3genudtagud->GetYaxis()->SetTitleSize(0.06);  
h3genudtagud->GetYaxis()->SetLabelSize(0.05);  

h3genxtagud->GetXaxis()->SetTitleSize(0.06);  
h3genxtagud->GetXaxis()->SetLabelSize(0.05);  
h3genxtagud->GetYaxis()->SetTitleSize(0.06);  
h3genxtagud->GetYaxis()->SetLabelSize(0.05);  

h3tagudall->GetXaxis()->SetTitleSize(0.06);  
h3tagudall->GetXaxis()->SetLabelSize(0.05);  
h3tagudall->GetYaxis()->SetTitleSize(0.06);  
h3tagudall->GetYaxis()->SetLabelSize(0.05);  

h3gencstagcs->GetXaxis()->SetTitleSize(0.06);  
h3gencstagcs->GetXaxis()->SetLabelSize(0.05);  
h3gencstagcs->GetYaxis()->SetTitleSize(0.06);  
h3gencstagcs->GetYaxis()->SetLabelSize(0.05);  

h3genudtagcs->GetXaxis()->SetTitleSize(0.06);  
h3genudtagcs->GetXaxis()->SetLabelSize(0.05);  
h3genudtagcs->GetYaxis()->SetTitleSize(0.06);  
h3genudtagcs->GetYaxis()->SetLabelSize(0.05);  

h3genxtagcs->GetXaxis()->SetTitleSize(0.06);  
h3genxtagcs->GetXaxis()->SetLabelSize(0.05);  
h3genxtagcs->GetYaxis()->SetTitleSize(0.06);  
h3genxtagcs->GetYaxis()->SetLabelSize(0.05);  

h3tagcsall->GetXaxis()->SetTitleSize(0.06);  
h3tagcsall->GetXaxis()->SetLabelSize(0.05);  
h3tagcsall->GetYaxis()->SetTitleSize(0.06);  
h3tagcsall->GetYaxis()->SetLabelSize(0.05);  

h3gencsall->GetXaxis()->SetTitleSize(0.06);  
h3gencsall->GetXaxis()->SetLabelSize(0.05);  
h3gencsall->GetYaxis()->SetTitleSize(0.06);  
h3gencsall->GetYaxis()->SetLabelSize(0.05);  

h3genudall->GetXaxis()->SetTitleSize(0.06);  
h3genudall->GetXaxis()->SetLabelSize(0.05);  
h3genudall->GetYaxis()->SetTitleSize(0.06);  
h3genudall->GetYaxis()->SetLabelSize(0.05);  

h3genxall->GetXaxis()->SetTitleSize(0.06);  
h3genxall->GetXaxis()->SetLabelSize(0.05);  
h3genxall->GetYaxis()->SetTitleSize(0.06);  
h3genxall->GetYaxis()->SetLabelSize(0.05);  

h3all->GetXaxis()->SetTitleSize(0.06);  
h3all->GetXaxis()->SetLabelSize(0.05);  
h3all->GetYaxis()->SetTitleSize(0.06);  
h3all->GetYaxis()->SetLabelSize(0.05);


    
    h3gencstagcs->GetYaxis()->SetTitleOffset(1.1);
    h3gencstagud->GetYaxis()->SetTitleOffset(1.1);
    h3gencstagx->GetYaxis()->SetTitleOffset(1.1);
    h3all->GetYaxis()->SetTitleOffset(1.1);
    h3genudtagcs->GetYaxis()->SetTitleOffset(1.15);
    h3genudtagud->GetYaxis()->SetTitleOffset(1.1);
    h3genudtagx->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagcs->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagud->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagx->GetYaxis()->SetTitleOffset(1.1);

    h3gencsall->GetYaxis()->SetTitleOffset(1.1);
    h3genudall->GetYaxis()->SetTitleOffset(1.1);
    h3genxall->GetYaxis()->SetTitleOffset(1.1);
    h3tagcsall->GetYaxis()->SetTitleOffset(1.1);
    h3tagudall->GetYaxis()->SetTitleOffset(1.1);
    h3tagxall->GetYaxis()->SetTitleOffset(1.1);

    // Update the canvas to reflect the changes
    gPad->Update();


    // Save the canvas as a .pdf file
    c3->SaveAs("pdf/MCvsDATA_Pt_Run2_FSR.pdf");










setTDRStyle();
        // Create a canvas
    TCanvas *c100 = new TCanvas("c100", "Canvas with 3x3 Grid", 2560, 1140);
    

    // Divide the canvas into a 3x3 grid
    c100->Divide(4, 2);

    c100->cd(1);
    
    //TH1D *h101 = tdrHist("h101","Average mass (GeV)",0,0.1,"Jet pair, reco",55,115);
    //TCanvas *c101 = tdrCanvas("c101",h101,8,11,kSquare);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagcsall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagcsall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagcsall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.0059);
    h3tagcsall->GetYaxis()->SetTitleSize(0.045);
    h3tagcsall->GetXaxis()->SetTitleSize(0.045);
    h3tagcsall->GetYaxis()->SetLabelSize(0.045);
    h3tagcsall->GetXaxis()->SetLabelSize(0.045);
    h3tagcsall->GetYaxis()->SetNdivisions(505);
    // Adjust the y-axis title location
    h3tagcsall->GetYaxis()->SetTitleOffset(1.5);

    TLegend *leg103 = tdrLeg(0.68,0.9-0.05*3,0.85,0.9);
    leg103->AddEntry(h3tagcsall, "FSR, no pgof", "FPLE");
    leg103->AddEntry(h3tagcsall_scaled, "FSR, pgof", "FPLE");
    //leg103->AddEntry(h3tagcsall_mass, "MC (init.+JES)", "FPLE");
    leg103->AddEntry(h3tagcsall_data, "DATA", "PLE");
    //leg103->SetTextSize(0.037);

    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.045);
    tex103->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(2);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagudall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagudall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagudall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagudall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagudall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetTitleSize(0.045);
    h3tagudall->GetXaxis()->SetTitleSize(0.045);
    h3tagudall->GetYaxis()->SetLabelSize(0.045);
    h3tagudall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagudall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg102 = tdrLeg(0.68,0.9-0.05*3,0.85,0.9);
    leg102->AddEntry(h3tagudall, "FSR, no pgof", "FPLE");
    leg102->AddEntry(h3tagudall_scaled, "FSR, pgof", "FPLE");
    //leg102->AddEntry(h3tagudall_mass, "MC (init.+JES)", "FPLE");
    leg102->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex102 = new TLatex();
    tex102->SetNDC(); tex102->SetTextSize(0.045);
    tex102->DrawLatex(0.17,0.92,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(3);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagxall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagxall_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagxall_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagxall->GetYaxis()->SetTitle("N frac");
    h3tagxall->GetYaxis()->SetTitleSize(0.045);
    h3tagxall->GetXaxis()->SetTitleSize(0.045);
    h3tagxall->GetYaxis()->SetLabelSize(0.045);
    h3tagxall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagxall->GetYaxis()->SetTitleOffset(1.65);
    

    // Set the x-axis range
    h3tagxall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.039); // Adjust the range as needed

    TLegend *leg101 = tdrLeg(0.68,0.9-0.05*3,0.85,0.9);
    leg101->AddEntry(h3tagxall, "FSR, no pgof", "FPLE");
    leg101->AddEntry(h3tagxall_scaled, "FSR, pgof", "FPLE");
    //leg101->AddEntry(h3tagxall_mass, "MC (init.+JES)", "FPLE");
    leg101->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex101 = new TLatex();
    tex101->SetNDC(); tex101->SetTextSize(0.045);
    tex101->DrawLatex(0.17,0.92,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(4);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3all_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3all_mass,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3all_mass->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3all->GetYaxis()->SetTitle("N frac");
    h3all->GetYaxis()->SetTitleSize(0.045);
    h3all->GetXaxis()->SetTitleSize(0.045);
    h3all->GetYaxis()->SetLabelSize(0.045);
    h3all->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3all->GetYaxis()->SetTitleOffset(1.65);

    // Set the x-axis range
    h3all->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3all->GetYaxis()->SetRangeUser(0, 0.049); // Adjust the range as needed

    TLegend *leg104 = tdrLeg(0.68,0.9-0.05*3,0.85,0.9);
    leg104->AddEntry(h3all, "FSR, no pgof", "FPLE");
    leg104->AddEntry(h3all_scaled, "FSR, pgof", "FPLE");
    //leg104->AddEntry(h3all_mass, "MC (init.+JES)", "FPLE");
    leg104->AddEntry(h3all_data, "DATA", "PLE");
    

    TLatex *tex104 = new TLatex();
    tex104->SetNDC(); tex104->SetTextSize(0.045);
    tex104->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(5);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagcsvsdata = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata");
    h_tagcsvsdata->Divide(h3tagcsall_data);

    TH1D *h_tagcsvsdata_mass = (TH1D*)h3tagcsall_mass->Clone("h_tagcsvsdata_mass");
    h_tagcsvsdata_mass->Divide(h3tagcsall_data);

    TH1D *h_tagcsvsdata_scaled = (TH1D*)h3tagcsall_scaled->Clone("h_tagcsvsdata_scaled");
    h_tagcsvsdata_scaled->Divide(h3tagcsall_data);

    tdrDraw(h_tagcsvsdata,"Pz",kFullCircle,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1);

    tdrDraw(h_tagcsvsdata_mass,"Pz",kFullCircle,kSpring-5);
    h_tagcsvsdata_mass->SetMarkerSize(1);

    //tdrDraw(h_tagcsvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    //h_tagcsvsdata_scaled->SetMarkerSize(1);

    h_tagcsvsdata->GetYaxis()->SetRangeUser(0, 3.5); // Adjust the range as needed
    h_tagcsvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagcsvsdata->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h_tagcsvsdata->GetYaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetYaxis()->SetLabelSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h_tagcsvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line107 = new TLine(-1, 1, 1, 1);
    line107->SetLineColor(kGray);
    line107->SetLineStyle(kDashed);
    line107->Draw();

    tdrDraw(h_tagcsvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagcsvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagcsvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_tagcsvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_tagcsvsdata_mass->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f7 = new TF1("f7", "[0]", -0.7, 0.7);
    f7->FixParameter(0, 1);
    TF1* f7m = new TF1("f7m", "[0]", -0.7, 0.7);
    f7m->FixParameter(0, 1);
    TF1* f7f = new TF1("f7f", "[0]", -0.7, 0.7);
    f7f->FixParameter(0, 1);

    h_tagcsvsdata->Fit(f7, "RN");
    double chi27 = f7->GetChisquare();
    int ndf7 = f7->GetNDF();
    h_tagcsvsdata_scaled->Fit(f7m, "RN");
    double chi27m = f7m->GetChisquare();
    int ndf7m = f7m->GetNDF();
    h_tagcsvsdata_mass->Fit(f7f, "RN");
    double chi27f = f7f->GetChisquare();
    int ndf7f = f7f->GetNDF();

    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi27, ndf7));
    tex103->DrawLatex(0.17, 0.78, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi27m, ndf7m));
    //tex103->DrawLatex(0.17, 0.71, Form("#chi_{MC (init.+JES)}^{2} / NDF = %1.1f / %d", chi27f, ndf7f));

    TLegend *leg107 = tdrLeg(0.62,0.9-0.05*2,0.8,0.9);
    leg107->AddEntry(h_tagcsvsdata, "FSR, no pgof", "PLE");
    leg107->AddEntry(h_tagcsvsdata_scaled, "FSR, pgof", "PLE");
    //leg107->AddEntry(h_tagcsvsdata_mass, "MC (init.+JES)", "PLE");

    TLatex *tex107 = new TLatex();
    tex107->SetNDC(); tex107->SetTextSize(0.045);
    tex107->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();
    gPad->RedrawAxis();

    c100->cd(6);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagudvsdata = (TH1D*)h3tagudall->Clone("h_tagudvsdata");
    h_tagudvsdata->Divide(h3tagudall_data);

    TH1D *h_tagudvsdata_mass = (TH1D*)h3tagudall_mass->Clone("h_tagudvsdata_mass");
    h_tagudvsdata_mass->Divide(h3tagudall_data);

    TH1D *h_tagudvsdata_scaled = (TH1D*)h3tagudall_scaled->Clone("h_tagudvsdata_scaled");
    h_tagudvsdata_scaled->Divide(h3tagudall_data);

    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagudvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagudvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_tagudvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_tagudvsdata_mass->SetMarkerSize(1.3);
    

    h_tagudvsdata->GetYaxis()->SetRangeUser(0, 4.1); // Adjust the range as needed
    h_tagudvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagudvsdata->GetXaxis()->SetTitle("(s-c)/(s+c)");
    // Adjust the y-axis title location
    h_tagudvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line106 = new TLine(-1, 1, 1, 1);
    line106->SetLineColor(kGray);
    line106->SetLineStyle(kDashed);
    line106->Draw();

    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagudvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagudvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_tagudvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_tagudvsdata_mass->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f6 = new TF1("f6", "[0]", -0.7, 0.7);
    f6->FixParameter(0, 1);
    TF1* f6m = new TF1("f6m", "[0]", -0.7, 0.7);
    f6m->FixParameter(0, 1);
    TF1* f6f = new TF1("f6f", "[0]", -0.7, 0.7);
    f6f->FixParameter(0, 1);

    h_tagudvsdata->Fit(f6, "RN");
    double chi26 = f6->GetChisquare();
    int ndf6 = f6->GetNDF();
    h_tagudvsdata_scaled->Fit(f6m, "RN");
    double chi26m = f6m->GetChisquare();
    int ndf6m = f6m->GetNDF();
    h_tagudvsdata_mass->Fit(f6f, "RN");
    double chi26f = f6f->GetChisquare();
    int ndf6f = f6f->GetNDF();


    TLegend *leg106 = tdrLeg(0.62,0.9-0.05*2,0.8,0.9);
    leg106->AddEntry(h_tagudvsdata, "FSR, no pgof", "PLE");
    leg106->AddEntry(h_tagudvsdata_scaled, "FSR, pgof", "PLE");
    //leg106->AddEntry(h_tagudvsdata_mass, "MC (init.+JES)", "PLE");
    

    TLatex *tex106 = new TLatex();
    tex106->SetNDC(); tex106->SetTextSize(0.045);
    tex106->DrawLatex(0.17,0.92,"genall,tagud");

    // Draw the chi2/ndf value on the plot
    tex106->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi26, ndf6));
    tex106->DrawLatex(0.17, 0.78, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi26m, ndf6m));
    //tex106->DrawLatex(0.17, 0.71, Form("#chi_{MC (init.+JES)}^{2} / NDF = %1.1f / %d", chi26f, ndf6f));

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(7);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagxvsdata = (TH1D*)h3tagxall->Clone("h_tagxvsdata");
    h_tagxvsdata->Divide(h3tagxall_data);

    TH1D *h_tagxvsdata_mass = (TH1D*)h3tagxall_mass->Clone("h_tagxvsdata_mass");
    h_tagxvsdata_mass->Divide(h3tagxall_data);

    TH1D *h_tagxvsdata_scaled = (TH1D*)h3tagxall_scaled->Clone("h_tagxvsdata_scaled");
    h_tagxvsdata_scaled->Divide(h3tagxall_data);

    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagxvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagxvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_tagxvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_tagxvsdata_mass->SetMarkerSize(1.3);

    h_tagxvsdata->GetYaxis()->SetRangeUser(0, 2.55); // Adjust the range as needed
    h_tagxvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagxvsdata->GetXaxis()->SetTitle("(s-c)/(s+c)");
    // Adjust the y-axis title location
    h_tagxvsdata->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line105 = new TLine(-1, 1, 1, 1);
    line105->SetLineColor(kGray);
    line105->SetLineStyle(kDashed);
    line105->Draw();

    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagxvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagxvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_tagxvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_tagxvsdata_mass->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f5 = new TF1("f5", "[0]", -0.7, 0.7);
    f5->FixParameter(0, 1);
    TF1* f5m = new TF1("f5m", "[0]", -0.7, 0.7);
    f5m->FixParameter(0, 1);
    TF1* f5f = new TF1("f5f", "[0]", -0.7, 0.7);
    f5f->FixParameter(0, 1);

    h_tagxvsdata->Fit(f5, "RN");
    double chi25 = f5->GetChisquare();
    int ndf5 = f5->GetNDF();
    h_tagxvsdata_scaled->Fit(f5m, "RN");
    double chi25m = f5m->GetChisquare();
    int ndf5m = f5m->GetNDF();
    h_tagxvsdata_mass->Fit(f5f, "RN");
    double chi25f = f5f->GetChisquare();
    int ndf5f = f5f->GetNDF();

    TLegend *leg105 = tdrLeg(0.62,0.9-0.05*2,0.8,0.9);
    leg105->AddEntry(h_tagxvsdata, "FSR, no pgof", "PLE");
    leg105->AddEntry(h_tagxvsdata_scaled, "FSR, pgof", "PLE");
    //leg105->AddEntry(h_tagxvsdata_mass, "MC (init.+JES)", "PLE");
    

    TLatex *tex105 = new TLatex();
    tex105->SetNDC(); tex105->SetTextSize(0.045);
    tex105->DrawLatex(0.17,0.92,"genall,tagx");

    // Draw the chi2/ndf value on the plot
    tex105->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi25, ndf5));
    tex105->DrawLatex(0.17, 0.78, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi25m, ndf5m));
    //tex105->DrawLatex(0.17, 0.71, Form("#chi_{MC (init.+JES)}^{2} / NDF = %1.1f / %d", chi25f, ndf5f));

    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(8);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_allvsdata = (TH1D*)h3all->Clone("h_allvsdata");
    h_allvsdata->Divide(h3all_data);

    TH1D *h_allvsdata_mass = (TH1D*)h3all_mass->Clone("h_allvsdata_mass");
    h_allvsdata_mass->Divide(h3all_data);

    TH1D *h_allvsdata_scaled = (TH1D*)h3all_scaled->Clone("h_allvsdata_scaled");
    h_allvsdata_scaled->Divide(h3all_data);

    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);

    tdrDraw(h_allvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_allvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_allvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_allvsdata_mass->SetMarkerSize(1.3);

    h_allvsdata->GetYaxis()->SetRangeUser(0, 2.6); // Adjust the range as needed
    h_allvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_allvsdata->GetXaxis()->SetTitle("(s-c)/(s+c)");
    // Adjust the y-axis title location
    h_allvsdata->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line108 = new TLine(-1, 1, 1, 1);
    line108->SetLineColor(kGray);
    line108->SetLineStyle(kDashed);
    line108->Draw();

    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);

    tdrDraw(h_allvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_allvsdata_scaled->SetMarkerSize(1.2);

    //tdrDraw(h_allvsdata_mass,"Pz",kFullDiamond,kSpring-5);
    //h_allvsdata_mass->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f8 = new TF1("f8", "[0]", -0.7, 0.7);
    f8->FixParameter(0, 1);
    TF1* f8m = new TF1("f8m", "[0]", -0.7, 0.7);
    f8m->FixParameter(0, 1);
    TF1* f8f = new TF1("f8f", "[0]", -0.7, 0.7);
    f8f->FixParameter(0, 1);

    h_allvsdata->Fit(f8, "RN");
    h_allvsdata_scaled->Fit(f8m, "RN");
    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();
    double chi28m = f8m->GetChisquare();
    int ndf8m = f8m->GetNDF();
    h_allvsdata_mass->Fit(f8f, "RN");
    double chi28f = f8f->GetChisquare();
    int ndf8f = f8f->GetNDF();

    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi28, ndf8));
    tex103->DrawLatex(0.17, 0.78, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi28m, ndf8m));
    //tex103->DrawLatex(0.17, 0.71, Form("#chi_{MC (init.+JES)}^{2} / NDF = %1.1f / %d", chi28f, ndf8f));

    TLegend *leg108 = tdrLeg(0.62,0.9-0.05*2,0.8,0.9);
    leg108->AddEntry(h_allvsdata, "FSR, no pgof", "PLE");
    leg108->AddEntry(h_allvsdata_scaled, "FSR, pgof", "PLE");
    //leg108->AddEntry(h_allvsdata_mass, "MC (init.+JES)", "PLE");
    

    TLatex *tex108 = new TLatex();
    tex108->SetNDC(); tex108->SetTextSize(0.045);
    tex108->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();



    // Save the canvas as a .pdf file
    c100->SaveAs("pdf/MCvsDATA_Pt2_Run2_FSR.pdf");

}