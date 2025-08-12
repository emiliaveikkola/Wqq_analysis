#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

// Define g(x) as a cubic function without a constant term: g(x) = a * (x - 83)^3 + b * (x - 83)^2 + c * (x - 83)
double gFunction(double *x, double *params) {
    double a = params[1]; // coefficient for (x - 83)^3
    double b = params[2]; // coefficient for (x - 83)^2
    double c = params[3]; // coefficient for (x - 83)
    double shifted_x = x[0] - 83;
    return a * pow(shifted_x, 3) + b * pow(shifted_x, 2) + c * shifted_x;
}

// Define f(x) = 1 - alpha * g(x)
double fFunction(double *x, double *params) {
    double alpha = params[0];  // First parameter is alpha
    return 1 - alpha * gFunction(x, params);
}

// Function to fit a, b, and c from the first histogram
void fitInitialParameters(TH1D* ratioHist, double& a, double& b, double& c) {
    TF1 *f1 = new TF1("f1_initial", fFunction, 66, 104, 4);  // Fit range from 66 to 104 GeV
    f1->SetParameter(0, 0.5);   // Initial guess for alpha
    f1->SetParameter(1, 0.1);   // Initial guess for a
    f1->SetParameter(2, 0.1);   // Initial guess for b
    f1->SetParameter(3, 0.1);   // Initial guess for c

    // Fit the function to the histogram and retrieve a, b, and c
    ratioHist->Fit(f1, "QRN");

    a = f1->GetParameter(1);
    b = f1->GetParameter(2);
    c = f1->GetParameter(3);

    std::cout << "Initial fit parameters: a=" << a << ", b=" << b << ", c=" << c << std::endl;
}

// Fit function for subsequent histograms using fixed a, b, and c values, allowing only alpha to vary
void fitFSRWithFixedParams(TH1D* ratioHist, const std::string& suffix, double a, double b, double c, int color) {
    TF1 *f1 = new TF1(("f1_" + suffix).c_str(), fFunction, 66, 104, 4);  // Fit range from 66 to 104 GeV
    f1->FixParameter(1, a);  // Fix a
    f1->FixParameter(2, b);  // Fix b
    f1->FixParameter(3, c);  // Fix c
    f1->SetParameter(0, 0.5);  // Initial guess for alpha

    // Fit only the alpha parameter
    ratioHist->Fit(f1, "QRN");

    double alphaFit = f1->GetParameter(0);
    std::cout << "Fitted alpha for " << suffix << ": " << alphaFit << std::endl;

    // Set line style, color, and width for the fitted function
    f1->SetLineColor(color);
    f1->SetLineStyle(2);  // Dashed line for distinction
    f1->SetLineWidth(2);

    // Draw the fitted function on top of the histogram
    f1->Draw("SAME");
}

void allratio_FSR() {
    setTDRStyle();
    extraText = "Private";

    // Open the file containing the processed histograms
    TFile* file = new TFile("processed_histograms/all.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open processed_histograms.root." << std::endl;
        return;
    }

    // Define suffixes, colors, and histogram name
    std::vector<std::string> suffixes = {"FSR0992", "FSR0995", "FSR0996", "FSR0998"};
    std::vector<int> colors = {kRed-4, kOrange+1, kSpring-5, kAzure+7, kViolet-5};
    std::vector<int> markers = {kOpenCircle, kFullSquare, kOpenSquare, kFullDiamond, kOpenDiamond};
    std::string histName = "h3all";

    TH1D *h2 = tdrHist("h2", "MC/MC_{ref}", 0.8,1.4, "Mass (GeV)", 66,104);
    TCanvas *c2 = tdrCanvas("c2", h2, 4, 11, kSquare);

    h2->GetYaxis()->SetLabelSize(0.045);
    h2->GetXaxis()->SetLabelSize(0.045);
    h2->GetYaxis()->SetTitleSize(0.055);
    h2->GetXaxis()->SetTitleSize(0.055);

    // Add pT and eta region text
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.32, 0.85, "p_{T} > 30 GeV");
    tex.DrawLatex(0.32, 0.8, "|#eta| < 2.5");

    // Retrieve the primary data and unscaled histograms for the ratio calculation
    TH1D* hist_data = (TH1D*)file->Get((histName + "_unscaled").c_str());
    TH1D* hist_unscaled = (TH1D*)file->Get((histName + "_data").c_str());
    if (!hist_data || !hist_unscaled) {
        std::cerr << "Error: Histograms for data or unscaled MC not found." << std::endl;
        file->Close();
        return;
    }

    // Add a legend for the ratio plot
    TLegend* ratioLegend = tdrLeg(0.5, 0.7, 0.85, 0.9);
    ratioLegend->SetTextSize(0.035);

    // Calculate the ratio of data to unscaled MC
    TH1D* ratioHist_data = (TH1D*)hist_data->Clone();
    ratioHist_data->Divide(hist_unscaled);

    // Style the ratio histogram
    tdrDraw(ratioHist_data, "Pz", kFullCircle, kGray+2);
    ratioHist_data->SetMarkerSize(1);
    ratioHist_data->SetLineWidth(2);

    ratioLegend->SetNColumns(2);
    ratioLegend->AddEntry(ratioHist_data, "DATA", "PLE");

    // Fit the ratio histograms with a constant function
    TF1* f8 = new TF1("f8", "[0]", 70, 100);
    f8->FixParameter(0, 1);

    ratioHist_data->Fit(f8, "RN");

    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();

    // Draw the chi2/ndf value on the plot
    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.035);
    tex103->DrawLatex(0.5, 0.65, Form("#chi_{DATA}^{2} / NDF = %1.1f / %d", chi28, ndf8));

    // Draw a dashed line at y=1 for reference
    TLine* line = new TLine(66, 1, 104, 1);
    line->SetLineColor(kGray);
    line->SetLineStyle(kDashed);
    line->Draw("SAME");

    // Draw additional lines at specified y-values
    TLine* line_1_05 = new TLine(66, 1.05, 104, 1.05);
    line_1_05->SetLineColor(kRed);
    line_1_05->SetLineStyle(kDashed);
    line_1_05->Draw("SAME");

    TLine* line_0_995 = new TLine(66, 0.95, 104, 0.95);
    line_0_995->SetLineColor(kRed);
    line_0_995->SetLineStyle(kDashed);
    line_0_995->Draw("SAME");

    TLine* line_1_02 = new TLine(66, 1.02, 104, 1.02);
    line_1_02->SetLineColor(kGreen-2);
    line_1_02->SetLineStyle(kDashed);
    line_1_02->Draw("SAME");

    TLine* line_0_998 = new TLine(66, 0.98, 104, 0.98);
    line_0_998->SetLineColor(kGreen-2);
    line_0_998->SetLineStyle(kDashed);
    line_0_998->Draw("SAME");

    // Fit the initial parameters a, b, and c using the first histogram in suffixes
    double a, b, c;
    std::string firstHistName = histName + "_" + suffixes[0];
    TH1D* firstHist = (TH1D*)file->Get(firstHistName.c_str());
    TH1D* firstRatioHist = (TH1D*)firstHist->Clone();
    firstRatioHist->Divide(hist_unscaled);
    fitInitialParameters(firstRatioHist, a, b, c);  // Fit a, b, and c only once

    // Add each ratio histogram to the legend, fit, and plot
    for (size_t j = 0; j < suffixes.size(); ++j) {
        std::string histFullName = histName + "_" + suffixes[j];
        TH1D* hist = (TH1D*)file->Get(histFullName.c_str());
        TH1D* ratioHist = (TH1D*)hist->Clone();
        ratioHist->Divide(hist_unscaled);

        //tdrDraw(ratioHist, "Pz", markers[j], colors[j]);
        ratioHist->SetMarkerSize(1);
        ratioHist->SetLineWidth(2);

        // Call the fit function with the specific color for this histogram
        //fitFSRWithFixedParams(ratioHist, suffixes[j], a, b, c, colors[j]);
        
        //ratioLegend->AddEntry(ratioHist, suffixes[j].c_str(), "PLE");
        gPad->Update();
        gPad->RedrawAxis();
    }
    ratioLegend->Draw();

    // Save the single ratio histogram as a PDF
    c2->SaveAs("pdf/all_ratio_FSR.pdf");

    // Clean up
    file->Close();
    delete file;
}