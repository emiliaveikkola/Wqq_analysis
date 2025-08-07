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
#include "tdrstyle_mod22.C"

void MCvsDATA_mass_Run2() {
    TFile *file = new TFile("output_MCRun2_tagprobe_test4.root", "READ");
    TH3D* h3MassFlavorPairs_DATAMC_MC = (TH3D*)file->Get("h3MassFlavorPairs_DATAMC");
    
      // Project Z for different indices and draw
      TH1D* h3all = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3all", 1, 3, 1, 3);

      TH1D* h3gencsall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencsall", 1, 1, 1, 3);
      TH1D* h3genudall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudall", 2, 2, 1, 3);
      TH1D* h3genxall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxall", 3, 3, 1, 3);

      TH1D* h3tagcsall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagcsall", 1, 3, 1, 1);
      TH1D* h3tagudall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagudall", 1, 3, 2, 2);
      TH1D* h3tagxall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagxall", 1, 3, 3, 3);

      TH1D* h3gencstagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagcs", 1, 1, 1, 1);
      TH1D* h3gencstagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagud", 1, 1, 2, 2);
      TH1D* h3gencstagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagx", 1, 1, 3, 3);
      
      TH1D* h3genudtagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagcs", 2, 2, 1, 1);
      TH1D* h3genudtagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagud", 2, 2, 2, 2);
      TH1D* h3genudtagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagx", 2, 2, 3, 3);
      
      TH1D* h3genxtagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagcs", 3, 3, 1, 1);
      TH1D* h3genxtagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagud", 3, 3, 2, 2);
      TH1D* h3genxtagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagx", 3, 3, 3, 3);

      TFile *file3 = new TFile("output_MCSCALEDRun2_tagprobe_test5.root", "READ");
      TH3D* h3MassFlavorPairs_DATAMC_RESP = (TH3D*)file3->Get("h3MassFlavorPairs_DATAMC_RESP");

      // Project Z for different indices and draw
      TH1D* h3all_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3all_resp", 1, 3, 1, 3);

      TH1D* h3gencsall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3gencsall_resp", 1, 1, 1, 3);
      TH1D* h3genudall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genudall_resp", 2, 2, 1, 3);
      TH1D* h3genxall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genxall_resp", 3, 3, 1, 3);

      TH1D* h3tagcsall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3tagcsall_resp", 1, 3, 1, 1);
      TH1D* h3tagudall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3tagudall_resp", 1, 3, 2, 2);
      TH1D* h3tagxall_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3tagxall_resp", 1, 3, 3, 3);

      TH1D* h3gencstagcs_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3gencstagcs_resp", 1, 1, 1, 1);
      TH1D* h3gencstagud_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3gencstagud_resp", 1, 1, 2, 2);
      TH1D* h3gencstagx_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3gencstagx_resp", 1, 1, 3, 3);
      
      TH1D* h3genudtagcs_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genudtagcs_resp", 2, 2, 1, 1);
      TH1D* h3genudtagud_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genudtagud_resp", 2, 2, 2, 2);
      TH1D* h3genudtagx_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genudtagx_resp", 2, 2, 3, 3);
      
      TH1D* h3genxtagcs_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genxtagcs_resp", 3, 3, 1, 1);
      TH1D* h3genxtagud_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genxtagud_resp", 3, 3, 2, 2);
      TH1D* h3genxtagx_resp = h3MassFlavorPairs_DATAMC_RESP->ProjectionZ("h3genxtagx_resp", 3, 3, 3, 3);

      TFile *file2 = new TFile("output_DATARun2_tagprobe_test4.root", "READ");
      TH3D* h3MassFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATAMC");

      TH1D* h3all_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

      TH1D* h3tagcsall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
      TH1D* h3tagudall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
      TH1D* h3tagxall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);

      TH1D* h3corr = new TH1D("h3corr",";Mass;N",200,0,200);

      TH1D* h3corrdata = new TH1D("h3corrdata",";Mass;N",200,0,200);


      // Cloning and renaming histograms
TH1D* h3all_scaled = (TH1D*)h3all->Clone("h3all_scaled");
TH1D* h3all_scaled2 = (TH1D*)h3all->Clone("h3all_scaled2");

TH1D* h3gencsall_scaled = (TH1D*)h3gencsall->Clone("h3gencsall_scaled");
TH1D* h3genudall_scaled = (TH1D*)h3genudall->Clone("h3genudall_scaled");
TH1D* h3genxall_scaled = (TH1D*)h3genxall->Clone("h3genxall_scaled");

TH1D* h3tagcsall_scaled = (TH1D*)h3tagcsall->Clone("h3tagcsall_scaled");
TH1D* h3tagudall_scaled = (TH1D*)h3tagudall->Clone("h3tagudall_scaled");
TH1D* h3tagxall_scaled = (TH1D*)h3tagxall->Clone("h3tagxall_scaled");

TH1D* h3gencstagcs_scaled = (TH1D*)h3gencstagcs->Clone("h3gencstagcs_scaled");
TH1D* h3gencstagud_scaled = (TH1D*)h3gencstagud->Clone("h3gencstagud_scaled");
TH1D* h3gencstagx_scaled = (TH1D*)h3gencstagx->Clone("h3gencstagx_scaled");

TH1D* h3genudtagcs_scaled = (TH1D*)h3genudtagcs->Clone("h3genudtagcs_scaled");
TH1D* h3genudtagud_scaled = (TH1D*)h3genudtagud->Clone("h3genudtagud_scaled");
TH1D* h3genudtagx_scaled = (TH1D*)h3genudtagx->Clone("h3genudtagx_scaled");

TH1D* h3genxtagcs_scaled = (TH1D*)h3genxtagcs->Clone("h3genxtagcs_scaled");
TH1D* h3genxtagud_scaled = (TH1D*)h3genxtagud->Clone("h3genxtagud_scaled");
TH1D* h3genxtagx_scaled = (TH1D*)h3genxtagx->Clone("h3genxtagx_scaled");


    
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

    h3gencstagcs_resp->Scale(1./h3all_resp->Integral());
    h3gencstagud_resp->Scale(1./h3all_resp->Integral());
    h3gencstagx_resp->Scale(1./h3all_resp->Integral());
    
    h3genudtagcs_resp->Scale(1./h3all_resp->Integral());
    h3genudtagud_resp->Scale(1./h3all_resp->Integral());
    h3genudtagx_resp->Scale(1./h3all_resp->Integral());

    h3genxtagcs_resp->Scale(1./h3all_resp->Integral());
    h3genxtagud_resp->Scale(1./h3all_resp->Integral());
    h3genxtagx_resp->Scale(1./h3all_resp->Integral());

    h3gencsall_resp->Scale(1./h3all_resp->Integral());
    h3genudall_resp->Scale(1./h3all_resp->Integral());
    h3genxall_resp->Scale(1./h3all_resp->Integral());

    h3tagcsall_resp->Scale(1./h3all_resp->Integral());
    h3tagudall_resp->Scale(1./h3all_resp->Integral());
    h3tagxall_resp->Scale(1./h3all_resp->Integral());

    h3all_resp->Scale(1./h3all_resp->Integral());

    h3tagcsall_data->Scale(1./h3all_data->Integral());
    h3tagudall_data->Scale(1./h3all_data->Integral());
    h3tagxall_data->Scale(1./h3all_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

// Total counts in h3all
double entries_h3all = h3all->Integral();
std::cout << "Total counts in h3all: " << entries_h3all << std::endl;

// Total counts in h3gencsall
double entries_h3gencsall = h3gencsall->Integral();
std::cout << "Total counts in h3gencsall: " << entries_h3gencsall << std::endl;

// Total counts in h3genudall
double entries_h3genudall = h3genudall->Integral();
std::cout << "Total counts in h3genudall: " << entries_h3genudall << std::endl;

// Total counts in h3genxall
double entries_h3genxall = h3genxall->Integral();
std::cout << "Total counts in h3genxall: " << entries_h3genxall << std::endl;

// Total counts in h3tagcsall
double entries_h3tagcsall = h3tagcsall->Integral();
std::cout << "Total counts in h3tagcsall: " << entries_h3tagcsall << std::endl;

// Total counts in h3tagudall
double entries_h3tagudall = h3tagudall->Integral();
std::cout << "Total counts in h3tagudall: " << entries_h3tagudall << std::endl;

// Total counts in h3tagxall
double entries_h3tagxall = h3tagxall->Integral();
std::cout << "Total counts in h3tagxall: " << entries_h3tagxall << std::endl;

// Total counts in h3gencstagcs
double entries_h3gencstagcs = h3gencstagcs->Integral();
std::cout << "Total counts in h3gencstagcs: " << entries_h3gencstagcs << std::endl;

// Continue this pattern for the remaining histograms
double entries_h3gencstagud = h3gencstagud->Integral();
std::cout << "Total counts in h3gencstagud: " << entries_h3gencstagud << std::endl;

double entries_h3gencstagx = h3gencstagx->Integral();
std::cout << "Total counts in h3gencstagx: " << entries_h3gencstagx << std::endl;

double entries_h3genudtagcs = h3genudtagcs->Integral();
std::cout << "Total counts in h3genudtagcs: " << entries_h3genudtagcs << std::endl;

double entries_h3genudtagud = h3genudtagud->Integral();
std::cout << "Total counts in h3genudtagud: " << entries_h3genudtagud << std::endl;

double entries_h3genudtagx = h3genudtagx->Integral();
std::cout << "Total counts in h3genudtagx: " << entries_h3genudtagx << std::endl;

double entries_h3genxtagcs = h3genxtagcs->Integral();
std::cout << "Total counts in h3genxtagcs: " << entries_h3genxtagcs << std::endl;

double entries_h3genxtagud = h3genxtagud->Integral();
std::cout << "Total counts in h3genxtagud: " << entries_h3genxtagud << std::endl;

double entries_h3genxtagx = h3genxtagx->Integral();
std::cout << "Total counts in h3genxtagx: " << entries_h3genxtagx << std::endl;

// For the data histograms
double entries_h3all_data = h3all_data->Integral();
std::cout << "Total counts in h3all_data: " << entries_h3all_data << std::endl;

double entries_h3tagcsall_data = h3tagcsall_data->Integral();
std::cout << "Total counts in h3tagcsall_data: " << entries_h3tagcsall_data << std::endl;

double entries_h3tagudall_data = h3tagudall_data->Integral();
std::cout << "Total counts in h3tagudall_data: " << entries_h3tagudall_data << std::endl;

double entries_h3tagxall_data = h3tagxall_data->Integral();
std::cout << "Total counts in h3tagxall_data: " << entries_h3tagxall_data << std::endl;

// Scaling histograms
/*
h3gencstagcs_scaled->Scale(0.802232);
h3gencstagx_scaled->Scale(1.06251);
h3gencstagud_scaled->Scale(1.2);
h3genudtagcs_scaled->Scale(0.8);
h3genudtagx_scaled->Scale(1.08923);
h3genudtagud_scaled->Scale(0.740574);
h3genxtagcs_scaled->Scale(0.80194);
h3genxtagx_scaled->Scale(1.07729);
h3genxtagud_scaled->Scale(0.80389);
*/
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

// Total counts in h3all
double entries_h3all_scaled = h3all_scaled->Integral();
std::cout << "Total counts in h3all_scaled: " << entries_h3all_scaled << std::endl;

    // Scaling histograms
/*
h3gencstagx_resp->Scale(1);
h3genudtagx_resp->Scale(1);
h3genxtagx_resp->Scale(1);

h3gencstagud_resp->Scale(0.9);
h3genudtagud_resp->Scale(0.8);
h3genxtagud_resp->Scale(0.8);

h3gencstagcs_resp->Scale(0.9);
h3genudtagcs_resp->Scale(0.8);
h3genxtagcs_resp->Scale(0.8);

*/

// Reset h3all_resp to zero


h3all_resp->Reset();
// Summing the resp histograms
h3all_resp->Add(h3gencstagcs_resp);
h3all_resp->Add(h3gencstagud_resp);
h3all_resp->Add(h3gencstagx_resp);
h3all_resp->Add(h3genudtagcs_resp);
h3all_resp->Add(h3genudtagud_resp);
h3all_resp->Add(h3genudtagx_resp);
h3all_resp->Add(h3genxtagcs_resp);
h3all_resp->Add(h3genxtagud_resp);
h3all_resp->Add(h3genxtagx_resp);

h3tagcsall_resp->Reset();
h3tagcsall_resp->Add(h3gencstagcs_resp);
h3tagcsall_resp->Add(h3genudtagcs_resp);
h3tagcsall_resp->Add(h3genxtagcs_resp);

h3tagudall_resp->Reset();
h3tagudall_resp->Add(h3gencstagud_resp);
h3tagudall_resp->Add(h3genudtagud_resp);
h3tagudall_resp->Add(h3genxtagud_resp);

h3tagxall_resp->Reset();
h3tagxall_resp->Add(h3gencstagx_resp);
h3tagxall_resp->Add(h3genudtagx_resp);
h3tagxall_resp->Add(h3genxtagx_resp);

h3gencsall_resp->Reset();
h3gencsall_resp->Add(h3gencstagcs_resp);
h3gencsall_resp->Add(h3gencstagud_resp);
h3gencsall_resp->Add(h3gencstagx_resp);

h3genudall_resp->Reset();
h3genudall_resp->Add(h3genudtagcs_resp);
h3genudall_resp->Add(h3genudtagud_resp);
h3genudall_resp->Add(h3genudtagx_resp);

h3genxall_resp->Reset();
h3genxall_resp->Add(h3genxtagcs_resp);
h3genxall_resp->Add(h3genxtagud_resp);
h3genxall_resp->Add(h3genxtagx_resp);

// Final scaling

double scale_factor2_resp = 1.0 / h3all_resp->Integral();

h3gencstagx_resp->Scale(scale_factor2_resp);
h3genudtagx_resp->Scale(scale_factor2_resp);
h3genxtagx_resp->Scale(scale_factor2_resp);

h3gencstagud_resp->Scale(scale_factor2_resp);
h3genudtagud_resp->Scale(scale_factor2_resp);
h3genxtagud_resp->Scale(scale_factor2_resp);

h3gencstagcs_resp->Scale(scale_factor2_resp);
h3genudtagcs_resp->Scale(scale_factor2_resp);
h3genxtagcs_resp->Scale(scale_factor2_resp);

h3gencsall_resp->Scale(scale_factor2_resp);
h3genudall_resp->Scale(scale_factor2_resp);
h3genxall_resp->Scale(scale_factor2_resp);

h3tagcsall_resp->Scale(scale_factor2_resp);
h3tagudall_resp->Scale(scale_factor2_resp);
h3tagxall_resp->Scale(scale_factor2_resp);

h3all_resp->Scale(scale_factor2_resp);

      
    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

        // Create a canvas
TCanvas *c3 = new TCanvas("c3", "Canvas with 3x3 Grid", 3840, 2160);  // Increase size to 4K resolution

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
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagx_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagx_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagx->GetXaxis()->SetTitle("Mass (GeV)");
    h3gencstagx->GetYaxis()->SetTitle("N frac");
    h3gencstagx->SetTitle("");

    // Set the x-axis range
    h3gencstagx->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3gencstagx->GetYaxis()->SetRangeUser(0, 0.0115);

    TLegend *leg1= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg1->AddEntry(h3gencstagx, "MC", "FPLE");
    leg1->AddEntry(h3gencstagx_scaled, "MC ("")", "FPLE");
    //leg1->AddEntry(h3gencstagx_resp, "MC (RESP)", "FPLE");

     // Update the canvas to reflect the changes
    TLatex *tex1 = new TLatex();
    tex1->SetNDC(); tex1->SetTextSize(0.045);
    tex1->DrawLatex(0.17,0.85,"gencs,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(2);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

        // Draw the histogram
    tdrDraw(h3genudtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagx_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagx_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagx->GetXaxis()->SetTitle("Mass (GeV)");
    h3genudtagx->GetYaxis()->SetTitle("N frac");
    h3genudtagx->SetTitle("");

    // Set the x-axis range
    h3genudtagx->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3genudtagx->GetYaxis()->SetRangeUser(0, 0.0159);

    TLegend *leg2= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg2->AddEntry(h3genudtagx, "MC", "FPLE");
    leg2->AddEntry(h3genudtagx_scaled, "MC ("")", "FPLE");
    //leg2->AddEntry(h3genudtagx_resp, "MC (RESP)", "FPLE");

    TLatex *tex2 = new TLatex();
    tex2->SetNDC(); tex2->SetTextSize(0.045);
    tex2->DrawLatex(0.17,0.85,"genud,tagx");
    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(3);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagx->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagx_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagx_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagx_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagx_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagx->GetXaxis()->SetTitle("Mass (GeV)");
    h3genxtagx->GetYaxis()->SetTitle("N frac");
    h3genxtagx->SetTitle("");

        // Set the x-axis range
    h3genxtagx->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3genxtagx->GetYaxis()->SetRangeUser(0, 0.0105);

    TLegend *leg3= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg3->AddEntry(h3genxtagx, "MC", "FPLE");
    leg3->AddEntry(h3genxtagx_scaled, "MC ("")", "FPLE");
    //leg3->AddEntry(h3genxtagx_resp, "MC (RESP)", "FPLE");


    TLatex *tex3 = new TLatex();
    tex3->SetNDC(); tex3->SetTextSize(0.045);
    tex3->DrawLatex(0.17,0.85,"genx,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(4);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagxall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagxall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagxall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagxall->GetYaxis()->SetTitle("N frac");
    h3tagxall->SetTitle("");

        // Set the x-axis range
    h3tagxall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.038); // Adjust the range as needed

    TLegend *leg4 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg4->AddEntry(h3tagxall, "MC", "FPLE");
    leg4->AddEntry(h3tagxall_scaled, "MC ("")", "FPLE");
    //leg4->AddEntry(h3tagxall_resp, "MC (RESP)", "FPLE");
    leg4->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex4 = new TLatex();
    tex4->SetNDC(); tex4->SetTextSize(0.045);
    tex4->DrawLatex(0.17,0.85,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();



    c3->cd(5);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagud_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagud_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagud->GetXaxis()->SetTitle("Mass (GeV)");
    h3gencstagud->GetYaxis()->SetTitle("N frac");
    h3gencstagud->GetYaxis()->SetTitleOffset(1.3);
    h3gencstagud->SetTitle("");

        // Set the x-axis range
    h3gencstagud->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3gencstagud->GetYaxis()->SetRangeUser(0, 0.0018); // Adjust the range as needed


    TLegend *leg5= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg5->AddEntry(h3gencstagud, "MC", "FPLE");
    leg5->AddEntry(h3gencstagud_scaled, "MC ("")", "FPLE");
    //leg5->AddEntry(h3gencstagud_resp, "MC (RESP)", "FPLE");


    TLatex *tex5 = new TLatex();
    tex5->SetNDC(); tex5->SetTextSize(0.045);
    tex5->DrawLatex(0.17,0.85,"gencs,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(6);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagud_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagud_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagud->GetXaxis()->SetTitle("Mass (GeV)");
    h3genudtagud->GetYaxis()->SetTitle("N frac");
    h3genudtagud->GetYaxis()->SetTitleOffset(1.3);
    h3genudtagud->SetTitle("");

        // Set the x-axis range
    h3genudtagud->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg6= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg6->AddEntry(h3genudtagud, "MC", "FPLE");
    leg6->AddEntry(h3genudtagud_scaled, "MC ("")", "FPLE");
    //leg6->AddEntry(h3genudtagud_resp, "MC (RESP)", "FPLE");

    TLatex *tex6 = new TLatex();
    tex6->SetNDC(); tex6->SetTextSize(0.045);
    tex6->DrawLatex(0.17,0.85,"genud,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(7);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagud->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagud_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagud_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagud_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagud_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagud->GetXaxis()->SetTitle("Mass (GeV)");
    h3genxtagud->GetYaxis()->SetTitle("N frac");
    h3genxtagud->GetYaxis()->SetTitleOffset(1.3);
    h3genxtagud->SetTitle("");

        // Set the x-axis range
    h3genxtagud->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg7 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg7->AddEntry(h3genxtagud, "MC", "FPLE");
    leg7->AddEntry(h3genxtagud_scaled, "MC ("")", "FPLE");
    //leg7->AddEntry(h3genxtagud_resp, "MC (RESP)", "FPLE");

    TLatex *tex7 = new TLatex();
    tex7->SetNDC(); tex7->SetTextSize(0.045);
    tex7->DrawLatex(0.17,0.85,"genx,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(8);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagudall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagudall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagudall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagudall->GetYaxis()->SetTitle("N frac");
    h3tagudall->SetTitle("");

        // Set the x-axis range
    h3tagudall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagudall_data->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg8 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg8->AddEntry(h3tagudall, "MC", "FPLE");
    leg8->AddEntry(h3tagudall_scaled, "MC ("")", "FPLE");
    //leg8->AddEntry(h3tagudall_resp, "MC (RESP)", "FPLE");
    leg8->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex8 = new TLatex();
    tex8->SetNDC(); tex8->SetTextSize(0.045);
    tex8->DrawLatex(0.17,0.85,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(9);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencstagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencstagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencstagcs_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencstagcs_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3gencstagcs->GetXaxis()->SetTitle("Mass (GeV)");
    h3gencstagcs->GetYaxis()->SetTitle("N frac");
    h3gencstagcs->GetYaxis()->SetTitleOffset(1.3);
    h3gencstagcs->SetTitle("");

        // Set the x-axis range
    h3gencstagcs->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    h3gencstagcs->GetYaxis()->SetRangeUser(0, 0.0045); // Adjust the range as needed

    TLegend *leg9 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg9->AddEntry(h3gencstagcs, "MC", "FPLE");
    leg9->AddEntry(h3gencstagcs_scaled, "MC ("")", "FPLE");
    //leg9->AddEntry(h3gencstagcs_resp, "MC (RESP)", "FPLE");

    TLatex *tex9 = new TLatex();
    tex9->SetNDC(); tex9->SetTextSize(0.045);
    tex9->DrawLatex(0.17,0.85,"gencs,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(10);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudtagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudtagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudtagcs_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudtagcs_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genudtagcs->GetXaxis()->SetTitle("Mass (GeV)");
    h3genudtagcs->GetYaxis()->SetTitle("N frac");
    h3genudtagcs->SetTitle("");
    h3genudtagcs->GetYaxis()->SetNoExponent();


    TLegend *leg10 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg10->AddEntry(h3genudtagcs, "MC", "FPLE");
    leg10->AddEntry(h3genudtagcs_scaled, "MC ("")", "FPLE");
    //leg10->AddEntry(h3genudtagcs_resp, "MC (RESP)", "FPLE");

        // Set the x-axis range
    h3genudtagcs->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLatex *tex10 = new TLatex();
    tex10->SetNDC(); tex10->SetTextSize(0.045);
    tex10->DrawLatex(0.17,0.85,"genud,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(11);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagcs->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxtagcs_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxtagcs_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxtagcs_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxtagcs_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genxtagcs->GetXaxis()->SetTitle("Mass (GeV)");
    h3genxtagcs->GetYaxis()->SetTitle("N frac");
    h3genxtagcs->SetTitle("");

        // Set the x-axis range
    h3genxtagcs->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg11 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg11->AddEntry(h3genxtagcs, "MC", "FPLE");
    leg11->AddEntry(h3genxtagcs_scaled, "MC ("")", "FPLE");
    //leg11->AddEntry(h3genxtagcs_resp, "MC (RESP)", "FPLE");

    TLatex *tex11 = new TLatex();
    tex11->SetNDC(); tex11->SetTextSize(0.045);
    tex11->DrawLatex(0.17,0.85,"genx,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(12);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagcsall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3tagcsall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3tagcsall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");
    h3tagcsall->GetYaxis()->SetTitleOffset(1.3);
    h3tagcsall->SetTitle("");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagcsall_data->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg12 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg12->AddEntry(h3tagcsall, "MC", "FPLE");
    leg12->AddEntry(h3tagcsall_scaled, "MC ("")", "FPLE");
    //leg12->AddEntry(h3tagcsall_resp, "MC (RESP)", "FPLE");
    leg12->AddEntry(h3tagcsall_data, "DATA", "PLE");

    TLatex *tex12 = new TLatex();
    tex12->SetNDC(); tex12->SetTextSize(0.045);
    tex12->DrawLatex(0.17,0.85,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(13);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3gencsall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3gencsall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3gencsall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3gencsall_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3gencsall->GetXaxis()->SetTitle("Mass (GeV)");
    h3gencsall->GetYaxis()->SetTitle("N frac");
    h3gencsall->SetTitle("");

        // Set the x-axis range
    h3gencsall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg13= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg13->AddEntry(h3gencsall, "MC", "FPLE");
    leg13->AddEntry(h3gencsall_scaled, "MC ("")", "FPLE");
    //leg13->AddEntry(h3gencsall_resp, "MC (RESP)", "FPLE");

    TLatex *tex13 = new TLatex();
    tex13->SetNDC(); tex13->SetTextSize(0.045);
    tex13->DrawLatex(0.17,0.85,"gencs,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(14);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genudall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genudall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genudall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genudall_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genudall->GetXaxis()->SetTitle("Mass (GeV)");
    h3genudall->GetYaxis()->SetTitle("N frac");
    h3genudall->SetTitle("");

        // Set the x-axis range
    h3genudall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg14= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg14->AddEntry(h3genudall, "MC", "FPLE");
    leg14->AddEntry(h3genudall_scaled, "MC ("")", "FPLE");
    //leg14->AddEntry(h3genudall_resp, "MC (RESP)", "FPLE");

    TLatex *tex14 = new TLatex();
    tex14->SetNDC(); tex14->SetTextSize(0.045);
    tex14->DrawLatex(0.17,0.85,"genud,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(15);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3genxall_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3genxall_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3genxall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3genxall_resp->SetFillColorAlpha(kSpring-5,0.25);

    h3genxall->GetXaxis()->SetTitle("Mass (GeV)");
    h3genxall->GetYaxis()->SetTitle("N frac");
    h3genxall->SetTitle("");

        // Set the x-axis range
    h3genxall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed

    TLegend *leg15= tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg15->AddEntry(h3genxall, "MC", "FPLE");
    leg15->AddEntry(h3genxall_scaled, "MC ("")", "FPLE");
    //leg15->AddEntry(h3genxall_resp, "MC (RESP)", "FPLE");

    TLatex *tex15 = new TLatex();
    tex15->SetNDC(); tex15->SetTextSize(0.045);
    tex15->DrawLatex(0.17,0.85,"genx,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(16);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_scaled,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3all_scaled->SetFillColorAlpha(kSpring-5,0.25);

    //tdrDraw(h3all_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    //h3all_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("Mass (GeV)");
    h3all->GetYaxis()->SetTitle("N frac");
    h3all->SetTitle("");

        // Set the x-axis range
    h3all->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3all->GetYaxis()->SetRangeUser(0, 0.05); // Adjust the range as needed

    TLegend *leg16 = tdrLeg(0.75,0.8-0.05*3,0.9,0.82);
    leg16->AddEntry(h3all, "MC", "FPLE");
    leg16->AddEntry(h3all_scaled, "MC ("")", "FPLE");
    //leg16->AddEntry(h3all_resp, "MC (RESP)", "FPLE");
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
    h3genudtagcs->GetYaxis()->SetTitleOffset(1.1);
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
    c3->SaveAs("pdf/MCvsDATA_Mass_Run2_RESP.pdf");







setTDRStyle();
        // Create a canvas
    TCanvas *c100 = new TCanvas("c100", "Canvas with 3x3 Grid", 2560, 1140);
    

    // Divide the canvas into a 3x3 grid
    c100->Divide(4, 2);

    c100->cd(1);
    
    //TH1D *h101 = tdrHist("h101","Average resp (GeV)",0,0.1,"Jet pair, reco",55,115);
    //TCanvas *c101 = tdrCanvas("c101",h101,8,11,kSquare);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_scaled,"HPz",kNone,kViolet-3,kSolid,-1,1001,kViolet-3);
    h3tagcsall_scaled->SetFillColorAlpha(kViolet-3,0.25);

    tdrDraw(h3tagcsall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagcsall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(-1, 1); // Adjust the range as needed
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.0049);
    h3tagcsall->GetYaxis()->SetTitleSize(0.045);
    h3tagcsall->GetXaxis()->SetTitleSize(0.045);
    h3tagcsall->GetYaxis()->SetLabelSize(0.045);
    h3tagcsall->GetXaxis()->SetLabelSize(0.045);
    h3tagcsall->GetYaxis()->SetNdivisions(505);
    // Adjust the y-axis title location
    h3tagcsall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg103 = tdrLeg(0.68,1.15-0.05*4,0.85,0.7);
    leg103->AddEntry(h3tagcsall, "MC", "FPLE");
    leg103->AddEntry(h3tagcsall_scaled, "MC ("")", "FPLE");
    leg103->AddEntry(h3tagcsall_resp, "MC (RESP)", "FPLE");
    leg103->AddEntry(h3tagcsall_data, "DATA", "PLE");
    //leg103->SetTextSize(0.037);

    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.045);
    tex103->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(2);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_scaled,"HPz",kNone,kViolet-3,kSolid,-1,1001,kViolet-3);
    h3tagudall_scaled->SetFillColorAlpha(kViolet-3,0.25);

    tdrDraw(h3tagudall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagudall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagudall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagudall->GetXaxis()->SetRangeUser(-1, 1); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetRangeUser(0, 0.0086); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetTitleSize(0.045);
    h3tagudall->GetXaxis()->SetTitleSize(0.045);
    h3tagudall->GetYaxis()->SetLabelSize(0.045);
    h3tagudall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagudall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg102 = tdrLeg(0.68,1.15-0.05*4,0.85,0.7);
    leg102->AddEntry(h3tagudall, "MC", "FPLE");
    leg102->AddEntry(h3tagudall_scaled, "MC ("")", "FPLE");
    leg102->AddEntry(h3tagudall_resp, "MC (RESP)", "FPLE");
    leg102->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex102 = new TLatex();
    tex102->SetNDC(); tex102->SetTextSize(0.045);
    tex102->DrawLatex(0.17,0.92,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(3);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_scaled,"HPz",kNone,kViolet-3,kSolid,-1,1001,kViolet-3);
    h3tagxall_scaled->SetFillColorAlpha(kViolet-3,0.25);

    tdrDraw(h3tagxall_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3tagxall_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagxall->GetYaxis()->SetTitle("N frac");
    h3tagxall->GetYaxis()->SetTitleSize(0.045);
    h3tagxall->GetXaxis()->SetTitleSize(0.045);
    h3tagxall->GetYaxis()->SetLabelSize(0.045);
    h3tagxall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagxall->GetYaxis()->SetTitleOffset(1.65);
    

    // Set the x-axis range
    h3tagxall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.044); // Adjust the range as needed

    TLegend *leg101 = tdrLeg(0.68,1.15-0.05*4,0.85,0.7);
    leg101->AddEntry(h3tagxall, "MC", "FPLE");
    leg101->AddEntry(h3tagxall_scaled, "MC ("")", "FPLE");
    leg101->AddEntry(h3tagxall_resp, "MC (RESP)", "FPLE");
    leg101->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex101 = new TLatex();
    tex101->SetNDC(); tex101->SetTextSize(0.045);
    tex101->DrawLatex(0.17,0.92,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(4);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram

    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_scaled,"HPz",kNone,kViolet-3,kSolid,-1,1001,kViolet-3);
    h3all_scaled->SetFillColorAlpha(kViolet-3,0.25);

    tdrDraw(h3all_resp,"HPz",kNone,kSpring-5,kSolid,-1,1001,kSpring-5);
    h3all_resp->SetFillColorAlpha(kSpring-5,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("Mass (GeV)");
    h3all->GetYaxis()->SetTitle("N frac");
    h3all->GetYaxis()->SetTitleSize(0.045);
    h3all->GetXaxis()->SetTitleSize(0.045);
    h3all->GetYaxis()->SetLabelSize(0.045);
    h3all->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3all->GetYaxis()->SetTitleOffset(1.65);

    // Set the x-axis range
    h3all->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3all->GetYaxis()->SetRangeUser(0, 0.052); // Adjust the range as needed

    TLegend *leg104 = tdrLeg(0.68,1.15-0.05*4,0.85,0.7);
    leg104->AddEntry(h3all, "MC", "FPLE");
    leg104->AddEntry(h3all_scaled, "MC ("")", "FPLE");
    leg104->AddEntry(h3all_resp, "MC (RESP)", "FPLE");
    leg104->AddEntry(h3all_data, "DATA", "PLE");
    

    TLatex *tex104 = new TLatex();
    tex104->SetNDC(); tex104->SetTextSize(0.045);
    tex104->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(5);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagcsvsdata = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata");
    h_tagcsvsdata->Divide(h3tagcsall_data);

    TH1D *h_tagcsvsdata_resp = (TH1D*)h3tagcsall_resp->Clone("h_tagcsvsdata_resp");
    h_tagcsvsdata_resp->Divide(h3tagcsall_data);

    TH1D *h_tagcsvsdata_scaled = (TH1D*)h3tagcsall_scaled->Clone("h_tagcsvsdata_scaled");
    h_tagcsvsdata_scaled->Divide(h3tagcsall_data);

    tdrDraw(h_tagcsvsdata,"Pz",kFullCircle,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1);

    tdrDraw(h_tagcsvsdata_resp,"Pz",kFullCircle,kSpring-5);
    h_tagcsvsdata_resp->SetMarkerSize(1);

    tdrDraw(h_tagcsvsdata_scaled,"Pz",kFullCircle,kSpring-5);
    h_tagcsvsdata_scaled->SetMarkerSize(1);

    h_tagcsvsdata->GetYaxis()->SetRangeUser(0, 3.1); // Adjust the range as needed
    h_tagcsvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagcsvsdata->GetYaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetYaxis()->SetLabelSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h_tagcsvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line107 = new TLine(55, 1, 115, 1);
    line107->SetLineColor(kGray);
    line107->SetLineStyle(kDashed);
    line107->Draw();

    tdrDraw(h_tagcsvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagcsvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_tagcsvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_tagcsvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_tagcsvsdata_resp->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f7 = new TF1("f7", "[0]", 70, 100);
    f7->FixParameter(0, 1);
    TF1* f7m = new TF1("f7m", "[0]", 70, 100);
    f7m->FixParameter(0, 1);
    TF1* f7f = new TF1("f7f", "[0]", 70, 100);
    f7f->FixParameter(0, 1);

    h_tagcsvsdata->Fit(f7, "RN");
    double chi27 = f7->GetChisquare();
    int ndf7 = f7->GetNDF();
    h_tagcsvsdata_scaled->Fit(f7m, "RN");
    double chi27m = f7m->GetChisquare();
    int ndf7m = f7m->GetNDF();
    h_tagcsvsdata_resp->Fit(f7f, "RN");
    double chi27f = f7f->GetChisquare();
    int ndf7f = f7f->GetNDF();

    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi27, ndf7));
    tex103->DrawLatex(0.17, 0.78, Form("#chi_{MC ("")}^{2} / NDF = %1.1f / %d", chi27m, ndf7m));
    tex103->DrawLatex(0.17, 0.71, Form("#chi_{MC (RESP)}^{2} / NDF = %1.1f / %d", chi27f, ndf7f));


    TLegend *leg107 = tdrLeg(0.62,1.10-0.05*3,0.8,0.8);
    leg107->AddEntry(h_tagcsvsdata, "MC", "PLE");
    leg107->AddEntry(h_tagcsvsdata_scaled, "MC ("")", "PLE");
    leg107->AddEntry(h_tagcsvsdata_resp, "MC (RESP)", "PLE");
    

    TLatex *tex107 = new TLatex();
    tex107->SetNDC(); tex107->SetTextSize(0.045);
    tex107->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();
    gPad->RedrawAxis();

    c100->cd(6);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagudvsdata = (TH1D*)h3tagudall->Clone("h_tagudvsdata");
    h_tagudvsdata->Divide(h3tagudall_data);

    TH1D *h_tagudvsdata_resp = (TH1D*)h3tagudall_resp->Clone("h_tagudvsdata_resp");
    h_tagudvsdata_resp->Divide(h3tagudall_data);

    TH1D *h_tagudvsdata_scaled = (TH1D*)h3tagudall_scaled->Clone("h_tagudvsdata_scaled");
    h_tagudvsdata_scaled->Divide(h3tagudall_data);


    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagudvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_tagudvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_tagudvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_tagudvsdata_resp->SetMarkerSize(1.3);
    

    h_tagudvsdata->GetYaxis()->SetRangeUser(0, 3.1); // Adjust the range as needed
    h_tagudvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_tagudvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line106 = new TLine(55, 1, 115, 1);
    line106->SetLineColor(kGray);
    line106->SetLineStyle(kDashed);
    line106->Draw();

    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagudvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_tagudvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_tagudvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_tagudvsdata_resp->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f6 = new TF1("f6", "[0]", 70, 100);
    f6->FixParameter(0, 1);
    TF1* f6m = new TF1("f6m", "[0]", 70, 100);
    f6m->FixParameter(0, 1);
    TF1* f6f = new TF1("f6f", "[0]", 70, 100);
    f6f->FixParameter(0, 1);

    h_tagudvsdata->Fit(f6, "RN");
    double chi26 = f6->GetChisquare();
    int ndf6 = f6->GetNDF();
    h_tagudvsdata_scaled->Fit(f6m, "RN");
    double chi26m = f6m->GetChisquare();
    int ndf6m = f6m->GetNDF();
    h_tagudvsdata_resp->Fit(f6f, "RN");
    double chi26f = f6f->GetChisquare();
    int ndf6f = f6f->GetNDF();


    TLegend *leg106 = tdrLeg(0.62,1.10-0.05*3,0.8,0.8);
    leg106->AddEntry(h_tagudvsdata, "MC", "PLE");
    leg106->AddEntry(h_tagudvsdata_scaled, "MC ("")", "PLE");
    leg106->AddEntry(h_tagudvsdata_resp, "MC (RESP)", "PLE");
    

    TLatex *tex106 = new TLatex();
    tex106->SetNDC(); tex106->SetTextSize(0.045);
    tex106->DrawLatex(0.17,0.92,"genall,tagud");

    // Draw the chi2/ndf value on the plot
    tex106->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi26, ndf6));
    tex106->DrawLatex(0.17, 0.78, Form("#chi_{MC ("")}^{2} / NDF = %1.1f / %d", chi26m, ndf6m));
    tex106->DrawLatex(0.17, 0.71, Form("#chi_{MC (RESP)}^{2} / NDF = %1.1f / %d", chi26f, ndf6f));

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(7);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagxvsdata = (TH1D*)h3tagxall->Clone("h_tagxvsdata");
    h_tagxvsdata->Divide(h3tagxall_data);

    TH1D *h_tagxvsdata_resp = (TH1D*)h3tagxall_resp->Clone("h_tagxvsdata_resp");
    h_tagxvsdata_resp->Divide(h3tagxall_data);

    TH1D *h_tagxvsdata_scaled = (TH1D*)h3tagxall_scaled->Clone("h_tagxvsdata_scaled");
    h_tagxvsdata_scaled->Divide(h3tagxall_data);

    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagxvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_tagxvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_tagxvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_tagxvsdata_resp->SetMarkerSize(1.3);

    h_tagxvsdata->GetYaxis()->SetRangeUser(0, 2.3); // Adjust the range as needed
    h_tagxvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_tagxvsdata->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line105 = new TLine(55, 1, 115, 1);
    line105->SetLineColor(kGray);
    line105->SetLineStyle(kDashed);
    line105->Draw();

    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);

    tdrDraw(h_tagxvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_tagxvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_tagxvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_tagxvsdata_resp->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f5 = new TF1("f5", "[0]", 70, 100);
    f5->FixParameter(0, 1);
    TF1* f5m = new TF1("f5m", "[0]", 70, 100);
    f5m->FixParameter(0, 1);
    TF1* f5f = new TF1("f5f", "[0]", 70, 100);
    f5f->FixParameter(0, 1);

    h_tagxvsdata->Fit(f5, "RN");
    double chi25 = f5->GetChisquare();
    int ndf5 = f5->GetNDF();
    h_tagxvsdata_scaled->Fit(f5m, "RN");
    double chi25m = f5m->GetChisquare();
    int ndf5m = f5m->GetNDF();
    h_tagxvsdata_resp->Fit(f5f, "RN");
    double chi25f = f5f->GetChisquare();
    int ndf5f = f5f->GetNDF();

    TLegend *leg105 = tdrLeg(0.6,1.10-0.05*3,0.8,0.8);
    leg105->AddEntry(h_tagxvsdata, "MC", "PLE");
    leg105->AddEntry(h_tagxvsdata_scaled, "MC ("")", "PLE");
    leg105->AddEntry(h_tagxvsdata_resp, "MC (RESP)", "PLE");
    

    TLatex *tex105 = new TLatex();
    tex105->SetNDC(); tex105->SetTextSize(0.045);
    tex105->DrawLatex(0.17,0.92,"genall,tagx");

    // Draw the chi2/ndf value on the plot
    tex105->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi25, ndf5));
    tex105->DrawLatex(0.17, 0.78, Form("#chi_{MC ("")}^{2} / NDF = %1.1f / %d", chi25m, ndf5m));
    tex105->DrawLatex(0.17, 0.71, Form("#chi_{MC (RESP)}^{2} / NDF = %1.1f / %d", chi25f, ndf5f));

    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(8);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_allvsdata = (TH1D*)h3all->Clone("h_allvsdata");
    h_allvsdata->Divide(h3all_data);

    TH1D *h_allvsdata_resp = (TH1D*)h3all_resp->Clone("h_allvsdata_resp");
    h_allvsdata_resp->Divide(h3all_data);

    TH1D *h_allvsdata_scaled = (TH1D*)h3all_scaled->Clone("h_allvsdata_scaled");
    h_allvsdata_scaled->Divide(h3all_data);

    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);

    tdrDraw(h_allvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_allvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_allvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_allvsdata_resp->SetMarkerSize(1.3);

    h_allvsdata->GetYaxis()->SetRangeUser(0, 2.6); // Adjust the range as needed
    h_allvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_allvsdata->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line108 = new TLine(55, 1, 115, 1);
    line108->SetLineColor(kGray);
    line108->SetLineStyle(kDashed);
    line108->Draw();

    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);

    tdrDraw(h_allvsdata_scaled,"Pz",kFullCircle,kViolet-3);
    h_allvsdata_scaled->SetMarkerSize(1.2);

    tdrDraw(h_allvsdata_resp,"Pz",kFullDiamond,kSpring-5);
    h_allvsdata_resp->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f8 = new TF1("f8", "[0]", 70, 100);
    f8->FixParameter(0, 1);
    TF1* f8m = new TF1("f8m", "[0]", 70, 100);
    f8m->FixParameter(0, 1);
    TF1* f8f = new TF1("f8f", "[0]", 70, 100);
    f8f->FixParameter(0, 1);

    h_allvsdata->Fit(f8, "RN");
    h_allvsdata_scaled->Fit(f8m, "RN");
    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();
    double chi28m = f8m->GetChisquare();
    int ndf8m = f8m->GetNDF();
    h_allvsdata_resp->Fit(f8f, "RN");
    double chi28f = f8f->GetChisquare();
    int ndf8f = f8f->GetNDF();

    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi28, ndf8));
    tex103->DrawLatex(0.17, 0.78, Form("#chi_{MC ("")}^{2} / NDF = %1.1f / %d", chi28m, ndf8m));
    tex103->DrawLatex(0.17, 0.71, Form("#chi_{MC (RESP)}^{2} / NDF = %1.1f / %d", chi28f, ndf8f));


    TLegend *leg108 = tdrLeg(0.62,1.10-0.05*3,0.8,0.8);
    leg108->AddEntry(h_allvsdata, "MC", "PLE");
    leg108->AddEntry(h_allvsdata_scaled, "MC ("")", "PLE");
    leg108->AddEntry(h_allvsdata_resp, "MC (RESP)", "PLE");
    

    TLatex *tex108 = new TLatex();
    tex108->SetNDC(); tex108->SetTextSize(0.045);
    tex108->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();



    // Save the canvas as a .pdf file
    c100->SaveAs("pdf/MCvsDATA_Mass2_Run2.pdf");

    // Open a new ROOT file to store histograms
TFile* outputFile = new TFile("output_histograms.root", "RECREATE");
// Write histograms to file after drawing
h3tagcsall->Write();
h3tagcsall_scaled->Write();
h3tagcsall_data->Write();

h3tagudall->Write();
h3tagudall_scaled->Write();
h3tagudall_data->Write();

h3tagxall->Write();
h3tagxall_scaled->Write();
h3tagxall_data->Write();

h3all->Write();
h3all_scaled->Write();
h3all_data->Write();
    
// Finally, close the file to save histograms
outputFile->Close();

}