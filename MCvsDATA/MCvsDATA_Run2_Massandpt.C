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

void MCvsDATA_Run2_Massandpt() {
    TFile *file = new TFile("../output/output_MCSCALEDRun2_tagprobe_test6.root", "READ");

    TH3D* h3PtFlavorPairs_DATAMC_MC = (TH3D*)file->Get("h3PtFlavorPairs_DATAMC");
    TH3D* h3PtFlavorPairs_DATAMC_GEN = (TH3D*)file->Get("h3PtFlavorPairs_DATAMC_gen");

    // Project Z for different indices and draw
    TH1D* h3all_pt_reco = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3all_pt_reco", 1, 3, 1, 3);

    TH1D* h3all_pt_gen = h3PtFlavorPairs_DATAMC_GEN->ProjectionZ("h3all_pt_gen", 1, 3, 1, 3);

    TH3D* h3MassFlavorPairs_DATAMC_MC = (TH3D*)file->Get("h3MassFlavorPairs_DATAMC");
    TH3D* h3MassFlavorPairs_DATAMC_GEN = (TH3D*)file->Get("h3MassFlavorPairs_DATAMC_gen");
    
    TH1D* h3all_mass_reco = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3all_mass_reco", 1, 3, 1, 3);
    TH1D* h3all_mass_gen = h3MassFlavorPairs_DATAMC_GEN->ProjectionZ("h3all_mass_gen", 1, 3, 1, 3);

    TFile *file2 = new TFile("../output/output_DATARun2_tagprobe_test4.root", "READ");
    TH3D* h3PtFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3PtFlavorPairs_DATAMC");
    TH3D* h3MassFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATAMC");

    TH1D* h3all_pt_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_pt_data", 1, 3, 1, 3);
    TH1D* h3all_mass_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_mass_data", 1, 3, 1, 3);

    // Perform scaling for all histograms
    h3all_pt_reco->Scale(1.0 / h3all_pt_reco->Integral());
    h3all_pt_gen->Scale(1.0 / h3all_pt_gen->Integral());
    h3all_mass_reco->Scale(1.0 / h3all_mass_reco->Integral());
    h3all_mass_gen->Scale(1.0 / h3all_mass_gen->Integral());
    h3all_pt_data->Scale(1.0 / h3all_pt_data->Integral());
    h3all_mass_data->Scale(1.0 / h3all_mass_data->Integral());

      
    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

    setTDRStyle();
        // Create a canvas
    TCanvas *c100 = new TCanvas("c100", "Canvas with 3x3 Grid", 2560, 1140);
    

    // Divide the canvas into a 3x3 grid
    c100->Divide(3, 1);

        // Add CMS label to the left corner of the canvas
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextSize(0.04);   // Text size
    cmsLabel.SetTextFont(61);     // CMS font
    cmsLabel.DrawLatex(0.04, 0.95, "CMS");  // Coordinates (0.1, 0.92) for top left corner

    // Add "Preliminary" label next to CMS
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextSize(0.035);
    prelimLabel.SetTextFont(52);  // Italic font for "Preliminary"
    prelimLabel.DrawLatex(0.08, 0.95, "Private");  // Coordinates (0.2, 0.92)

    // Add dataset and luminosity information in the top right corner
    TLatex dataLabel;
    dataLabel.SetNDC();
    dataLabel.SetTextSize(0.04);
    dataLabel.SetTextFont(42);
    dataLabel.DrawLatex(0.75, 0.955, "138 fb^{-1} (Run2 Legacy, 13 TeV)");

    //----------------------------- Pad 1: Mass --------------------------------//
    c100->cd(1);

    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.08);
    gPad->SetTopMargin(0.05);

    // Create histograms for mass
    TH1D *h3all_mass_reco_vs_data = (TH1D*)h3all_mass_reco->Clone("h3all_mass_reco_vs_data");
    //h3all_mass_reco_vs_data->Divide(h3all_mass_data);

    TH1D *h3all_mass_gen_vs_data = (TH1D*)h3all_mass_gen->Clone("h3all_mass_gen_vs_data");
    //h3all_mass_gen_vs_data->Divide(h3all_mass_data);

    // Draw histograms for mass ratio
    tdrDraw(h3all_mass_reco_vs_data, "Pz", kFullSquare, kAzure+7);
    h3all_mass_reco_vs_data->SetMarkerSize(1.3);

    tdrDraw(h3all_mass_gen_vs_data, "Pz", kFullCircle, kSpring-5);
    h3all_mass_gen_vs_data->SetMarkerSize(1.2);

    h3all_mass_reco_vs_data->GetYaxis()->SetRangeUser(0, 0.1); // Adjust the range as needed
    h3all_mass_reco_vs_data->GetXaxis()->SetRangeUser(66, 104);
    h3all_mass_reco_vs_data->GetYaxis()->SetTitle("MC/DATA");
    h3all_mass_reco_vs_data->GetXaxis()->SetTitle("Mass (GeV)");
    // Adjust the y-axis title location
    //h3all_mass_reco_vs_data->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line108 = new TLine(66, 1, 104, 1);
    line108->SetLineColor(kGray);
    line108->SetLineStyle(kDashed);
    line108->Draw("same");

    // Fit the ratio histograms with a constant function
    TF1* f8 = new TF1("f8", "[0]", 70, 100);
    f8->FixParameter(0, 1);
    TF1* f8m = new TF1("f8m", "[0]", 70, 100);
    f8m->FixParameter(0, 1);

    h3all_mass_reco_vs_data->Fit(f8, "RN");
    h3all_mass_gen_vs_data->Fit(f8m, "RN");

    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();
    double chi28m = f8m->GetChisquare();
    int ndf8m = f8m->GetNDF();

    // Draw the chi2/ndf value on the plot
    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.035);
    tex103->DrawLatex(0.17, 0.87, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi28, ndf8));
    tex103->DrawLatex(0.17, 0.83, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi28m, ndf8m));

    // Add legend
    TLegend *leg108 = tdrLeg(0.75, 0.9 - 0.05 * 2, 0.9, 0.9);
    leg108->AddEntry(h3all_mass_reco_vs_data, "MC_{reco}", "PLE");
    leg108->AddEntry(h3all_mass_gen_vs_data, "MC_{gen}", "PLE");

    // Calculate and display RMS for mass histograms
    double rms_mass_reco = h3all_mass_reco->GetRMS();
    double rms_mass_gen = h3all_mass_gen->GetRMS();
    double rms_mass_data = h3all_mass_data->GetRMS();

    tex103->DrawLatex(0.17, 0.78, Form("RMS_{MC_{reco}} = %1.3f", rms_mass_reco));
    tex103->DrawLatex(0.17, 0.74, Form("RMS_{MC_{gen}} = %1.3f", rms_mass_gen));
    tex103->DrawLatex(0.17, 0.7, Form("RMS_{DATA} = %1.3f", rms_mass_data));

    // Update the canvas
    gPad->Update();

    //----------------------------- Pad 2: Pt --------------------------------//
    c100->cd(2);

    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.08);
    gPad->SetTopMargin(0.05);

    // Create histograms for pt
    TH1D *h3all_pt_reco_vs_data = (TH1D*)h3all_pt_reco->Clone("h3all_pt_reco_vs_data");
    //h3all_pt_reco_vs_data->Divide(h3all_pt_data);

    TH1D *h3all_pt_gen_vs_data = (TH1D*)h3all_pt_gen->Clone("h3all_pt_gen_vs_data");
    //h3all_pt_gen_vs_data->Divide(h3all_pt_data);

    // Draw histograms for pt ratio
    tdrDraw(h3all_pt_reco_vs_data, "Pz", kFullSquare, kAzure+7);
    h3all_pt_reco_vs_data->SetMarkerSize(1.3);

    tdrDraw(h3all_pt_gen_vs_data, "Pz", kFullCircle, kSpring-5);
    h3all_pt_gen_vs_data->SetMarkerSize(1.2);

    h3all_pt_reco_vs_data->GetYaxis()->SetRangeUser(0, 0.06); // Adjust the range as needed
    h3all_pt_reco_vs_data->GetXaxis()->SetRangeUser(-0.8, 0.8); // Adjust the range as needed
    h3all_pt_reco_vs_data->GetYaxis()->SetTitle("MC/DATA");
    h3all_pt_reco_vs_data->GetXaxis()->SetTitle("(s-c)/(s+c)");
    // Adjust the y-axis title location
    h3all_pt_reco_vs_data->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line109 = new TLine(-0.8, 1, 0.8, 1);
    line109->SetLineColor(kGray);
    line109->SetLineStyle(kDashed);
    line109->Draw("same");

    // Fit the ratio histograms with a constant function
    TF1* f9 = new TF1("f9", "[0]", -0.7, 0.7);
    f9->FixParameter(0, 1);
    TF1* f9m = new TF1("f9m", "[0]", -0.7, 0.7);
    f9m->FixParameter(0, 1);

    h3all_pt_reco_vs_data->Fit(f9, "RN");
    h3all_pt_gen_vs_data->Fit(f9m, "RN");

    double chi29 = f9->GetChisquare();
    int ndf9 = f9->GetNDF();
    double chi29m = f9m->GetChisquare();
    int ndf9m = f9m->GetNDF();

    // Draw the chi2/ndf value on the plot
    TLatex *tex104 = new TLatex();
    tex104->SetNDC(); tex104->SetTextSize(0.035);
    tex104->DrawLatex(0.17, 0.87, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi29, ndf9));
    tex104->DrawLatex(0.17, 0.83, Form("#chi_{MC (init.)}^{2} / NDF = %1.1f / %d", chi29m, ndf9m));

    // Add legend
    TLegend *leg109 = tdrLeg(0.75, 0.9 - 0.05 * 2, 0.9, 0.9);
    leg109->AddEntry(h3all_pt_reco_vs_data, "MC_{reco}", "PLE");
    leg109->AddEntry(h3all_pt_gen_vs_data, "MC_{gen}", "PLE");

    // Calculate and display RMS for pt histograms
    double rms_pt_reco = h3all_pt_reco->GetRMS();
    double rms_pt_gen = h3all_pt_gen->GetRMS();
    double rms_pt_data = h3all_pt_data->GetRMS();

    tex104->DrawLatex(0.17, 0.78, Form("RMS_{MC_{reco}} = %1.3f", rms_pt_reco));
    tex104->DrawLatex(0.17, 0.74, Form("RMS_{MC_{gen}} = %1.3f", rms_pt_gen));
    tex104->DrawLatex(0.17, 0.7, Form("RMS_{DATA} = %1.3f", rms_pt_data));

    // Update the canvas
    gPad->Update();

   //----------------------------- Pad 3: ptrecovsgen --------------------------------//
    c100->cd(3);

    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.08);
    gPad->SetTopMargin(0.05);

    // Check if 'ptrecovsgen' is already defined
    TH1D* ptrecovsgen = nullptr;

    // Option 1: If 'ptrecovsgen' is loaded from a file or already defined
    ptrecovsgen = (TH1D*)file->Get("ptrecovsgen");

    // Option 2: If 'ptrecovsgen' needs to be created
    if (!ptrecovsgen) {
        // Create 'ptrecovsgen' as the difference between 'h3all_pt_reco' and 'h3all_pt_gen'
        ptrecovsgen = (TH1D*)h3all_pt_reco->Clone("ptrecovsgen");
        ptrecovsgen->Add(h3all_pt_gen, -1); // Subtract h3all_pt_gen from h3all_pt_reco
    }

    // Normalize 'ptrecovsgen' if needed
    // ptrecovsgen->Scale(1.0 / ptrecovsgen->Integral());

    // Set histogram title and axis labels
    ptrecovsgen->GetXaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
    ptrecovsgen->GetYaxis()->SetTitle("N");

    // Use the specified drawing style
    tdrDraw(ptrecovsgen, "HPz", kNone, kAzure+8, kSolid, -1, 0, kViolet+1);
    ptrecovsgen->SetLineWidth(1);

    double rms_ptrecovsgen = ptrecovsgen->GetRMS();
    tex104->DrawLatex(0.7, 0.8, Form("RMS = %1.3f", rms_ptrecovsgen));

    // Update the canvas
    gPad->Update();
    // Save the canvas as a .pdf file
    c100->SaveAs("pdf/MCvsDATA_Run2_massandpt.pdf");

}