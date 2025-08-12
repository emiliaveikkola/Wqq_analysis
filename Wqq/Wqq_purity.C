#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

// Function to compute normalized fraction: bin content / sum of all bins in histogram
void ComputeNormalizedFraction(TH1D *h_input, TH1D *h_frac) {
    // Compute sum of all bins
    double sum = 0;
    double sum_err2 = 0; // For error propagation: sum of squared errors
    for (int k = 1; k <= 7; ++k) {
        double content = h_input->GetBinContent(k);
        double err = h_input->GetBinError(k);
        sum += content;
        sum_err2 += err * err; // Sum of variances
    }
    double sum_err = sqrt(sum_err2); // Error on the sum

    // Compute normalized fraction for each bin
    for (int k = 1; k <= 7; ++k) {
        double content = h_input->GetBinContent(k);
        double content_err = h_input->GetBinError(k);
        
        if (sum > 0) {
            double frac = content / sum;
            double rel_err = 0;
            if (content > 0 && sum > 0) {
                // Relative error: sqrt((σ_content/content)^2 + (σ_sum/sum)^2)
                rel_err = sqrt(pow(content_err / content, 2) + pow(sum_err / sum, 2));
            }
            double frac_err = frac * rel_err;
            
            h_frac->SetBinContent(k, frac);
            h_frac->SetBinError(k, frac_err);
        } else {
            h_frac->SetBinContent(k, 0);
            h_frac->SetBinError(k, 0);
        }
    }
}

void Wqq_purity() {
    setTDRStyle();
    extraText = "Private";

    // Open the ROOT file and get the 1D histograms
    TFile *file = new TFile("output_MCRun2_Wqq_fractions.root", "READ");
    TH1D *h_cs = (TH1D*)file->Get("h_cs");
    TH1D *h_ud = (TH1D*)file->Get("h_ud");
    TH1D *h_xx = (TH1D*)file->Get("h_xx");
    TH1D *h_all = (TH1D*)file->Get("h_all");

    // Create fraction histograms
    TH1D *h_cs_frac = new TH1D("h_cs_frac", ";Flavor Pair;Fraction", 7, 0, 7);
    TH1D *h_ud_frac = new TH1D("h_ud_frac", ";Flavor Pair;Fraction", 7, 0, 7);
    TH1D *h_xx_frac = new TH1D("h_xx_frac", ";Flavor Pair;Fraction", 7, 0, 7);
    TH1D *h_all_frac = new TH1D("h_all_frac", ";Flavor Pair;Fraction", 7, 0, 7);


    // Set bin labels for fraction histograms
    const char *labels[] = {"cs", "ud", "cd", "us", "cb", "ub", "xx"};
    for (int i = 1; i <= 7; ++i) {
        h_cs_frac->GetXaxis()->SetBinLabel(i, labels[i-1]);
        h_ud_frac->GetXaxis()->SetBinLabel(i, labels[i-1]);
        h_xx_frac->GetXaxis()->SetBinLabel(i, labels[i-1]);
    }


    // Compute normalized fractions
    ComputeNormalizedFraction(h_cs, h_cs_frac);
    ComputeNormalizedFraction(h_ud, h_ud_frac);
    ComputeNormalizedFraction(h_xx, h_xx_frac);
    ComputeNormalizedFraction(h_all, h_all_frac);

    // Set up the canvas
    TH1D *h1 = tdrHist("h1", "N fraction", 0, 1, "Flavor pair", 0, 7);
    TCanvas *c1 = tdrCanvas("c1", h1, 4, 0, kSquare);

    // Copy bin labels to template histogram
    for (int i = 1; i <= 7; ++i) {
        h1->GetXaxis()->SetBinLabel(i, labels[i-1]);
    }
    h1->GetXaxis()->LabelsOption("h"); // Horizontal labels
    h1->GetXaxis()->SetLabelSize(0.05);


    
    // Style and draw histograms
    tdrDraw(h_all_frac,"HPz",kNone,kBlack,kSolid,kBlack,kNone,kBlack);
    h_all_frac->SetLineWidth(2);

    tdrDraw(h_cs_frac,"HPz",kNone,kAzure+8,kSolid,kAzure+7,1001,kAzure+8);
    h_cs_frac->SetFillColorAlpha(kAzure+8,0.8);

    tdrDraw(h_ud_frac,"HPz",kNone,kGreen-7,kSolid,kGreen-6,1001,kGreen-7);
    h_ud_frac->SetFillColorAlpha(kGreen-7,0.55);

    tdrDraw(h_xx_frac,"HPz",kNone,kRed-9,kSolid,kRed-7,1001,kRed-9);
    h_xx_frac->SetFillColorAlpha(kRed-9,0.35);




    // Add legend
    TLegend *leg = tdrLeg(0.6, 0.7, 0.85, 0.85);
    leg->AddEntry(h_all_frac, "all", "FL");
    leg->AddEntry(h_cs_frac, "cs-tag", "FPLE");
    leg->AddEntry(h_ud_frac, "ud-tag", "FPLE");
    leg->AddEntry(h_xx_frac, "xx-tag", "FPLE");
    leg->Draw();

    c1->RedrawAxis();
    gPad->Update();

    // Save the canvas as a .pdf file
    c1->SaveAs("pdf/Wqq_purity.pdf");
}