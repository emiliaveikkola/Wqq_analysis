#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

// Transform x to the range [-1, 1] for use with Chebyshev polynomials
double transformX(double x) {
    return 2 * (x - 70) / (100 - 70) - 1;
}

// Define g(x) using Chebyshev polynomials T_4, T_3, T_2, and T_1, plus a constant term
double gFunction(double *x, double *params) {
    double a = params[1]; // Coefficient for T_4(x)
    double b = params[2]; // Coefficient for T_3(x)
    double c = params[3]; // Coefficient for T_2(x)
    double d = params[4]; // Coefficient for T_1(x)
    double e = params[5]; // Constant term

    double tx = transformX(x[0]);

    double T1 = tx;
    double T2 = 2 * pow(tx, 2) - 1;
    double T3 = 4 * pow(tx, 3) - 3 * tx;
    double T4 = 8 * pow(tx, 4) - 8 * pow(tx, 2) + 1;

    return a * T4 + b * T3 + c * T2 + d * T1 + e;
}

// Define f(x) = 1 - alpha * g(x)
double fFunction(double *x, double *params) {
    double alpha = params[0];  // First parameter is alpha
    return 1 - alpha * gFunction(x, params);
}

// Function to fit initial parameters a, b, c, d, and e from the first histogram
void fitInitialParameters(TH1D* ratioHist, double& a, double& b, double& c, double& d, double& e) {
    TF1 *f1 = new TF1("f1_initial", fFunction, 66, 104, 6);  // Fit range from 66 to 104 GeV
    f1->SetParameter(0, 0.5);   // Initial guess for alpha
    f1->SetParameter(1, 0.1);   // Initial guess for a
    f1->SetParameter(2, 0.1);   // Initial guess for b
    f1->SetParameter(3, 0.1);   // Initial guess for c
    f1->SetParameter(4, 0.1);   // Initial guess for d
    f1->SetParameter(5, 0.1);   // Initial guess for e

    // Fit the function to the histogram and retrieve a, b, c, d, and e
    ratioHist->Fit(f1, "QRN");

    a = f1->GetParameter(1);
    b = f1->GetParameter(2);
    c = f1->GetParameter(3);
    d = f1->GetParameter(4);
    e = f1->GetParameter(5);

    std::cout << "Initial fit parameters: a=" << a << ", b=" << b << ", c=" << c 
              << ", d=" << d << ", e=" << e << std::endl;
}

// Fit function for subsequent histograms using fixed a, b, c, d, and e values, allowing only alpha to vary
void fitFSRWithFixedParams(TH1D* ratioHist, const std::string& suffix, double a, double b, double c, double d, double e, int color) {
    TF1 *f1 = new TF1(("f1_" + suffix).c_str(), fFunction, 66, 104, 6);  // Fit range from 66 to 104 GeV
    f1->FixParameter(1, a);  // Fix a
    f1->FixParameter(2, b);  // Fix b
    f1->FixParameter(3, c);  // Fix c
    f1->FixParameter(4, d);  // Fix d
    f1->FixParameter(5, e);  // Fix e
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

// Function to plot only the last ratio histogram and save as a PDF
void allratio_ISR() {
    setTDRStyle();
    extraText = "Private";

    // Open the file containing the processed histograms
    TFile* file = new TFile("processed_histograms/all.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open processed_histograms.root." << std::endl;
        return;
    }

    // Create a new ROOT file to store results
    TFile* outputFile = new TFile("processed_histograms/histograms_ISR.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create results.root." << std::endl;
        file->Close();
        return;
    }

    // Define suffixes, colors, and histogram name
    std::vector<std::string> suffixes = {"ISR040", "ISR070", "ISR130", "ISR160"};
    std::vector<int> colors = {kRed-4, kOrange+1, kSpring-5, kAzure+7, kViolet-5};
    std::vector<int> markers = {kOpenCircle, kFullSquare, kOpenSquare, kFullDiamond, kOpenDiamond};
    std::string histName = "h3all";

    // Create a new canvas for the single ratio histogram
    TH1D *h2 = tdrHist("h2", "MC/MC_{ref}", 0.55, 1.6, "Mass (GeV)", 66, 104);
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
    TH1D* hist_data = (TH1D*)file->Get((histName + "_data").c_str());
    TH1D* hist_unscaled = (TH1D*)file->Get((histName + "_unscaled").c_str());
    if (!hist_data || !hist_unscaled) {
        std::cerr << "Error: Histograms for data or unscaled MC not found." << std::endl;
        file->Close();
        outputFile->Close();
        return;
    }

    // Add a legend for the ratio plot
    TLegend* ratioLegend = tdrLeg(0.5, 0.75, 0.85, 0.9);
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

    // Draw a dashed line at y=1 for reference
    TLine* line = new TLine(66, 1, 104, 1);
    line->SetLineColor(kGray);
    line->SetLineStyle(kDashed);
    line->Draw("SAME");

    // Fit the initial parameters a, b, c, d, and e using the first histogram in suffixes
    double a, b, c, d, e;
    std::string firstHistName = histName + "_" + suffixes[3];
    TH1D* firstHist = (TH1D*)file->Get(firstHistName.c_str());
    TH1D* firstRatioHist = (TH1D*)firstHist->Clone();
    firstRatioHist->Divide(hist_unscaled);
    fitInitialParameters(firstRatioHist, a, b, c, d, e);  // Fit a, b, c, d, and e only once

    // Add each ratio histogram to the legend and plot
    for (size_t j = 0; j < suffixes.size(); ++j) {
        std::string histFullName = histName + "_" + suffixes[j];
        TH1D* hist = (TH1D*)file->Get(histFullName.c_str());
        TH1D* ratioHist = (TH1D*)hist->Clone();
        ratioHist->Divide(hist_unscaled);

        tdrDraw(ratioHist, "Pz", markers[j], colors[j]);
        ratioHist->SetMarkerSize(1);
        ratioHist->SetLineWidth(2);
        if (suffixes[j] == "ISR130") {
            ratioHist->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
        }
        if (suffixes[j] == "ISR160") {
            ratioHist->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
        }

        // Call the fit function with the specific color for this histogram
        //fitFSRWithFixedParams(ratioHist, suffixes[j], a, b, c, d, e, colors[j]);

        // Fit the ratio histogram and save the fit function
        TF1* fitFunc = new TF1(("FitFunction_" + suffixes[j]).c_str(), fFunction, 66, 104, 6);
        fitFunc->FixParameter(1, a);
        fitFunc->FixParameter(2, b);
        fitFunc->FixParameter(3, c);
        fitFunc->FixParameter(4, d);
        fitFunc->FixParameter(5, e);
        fitFunc->SetParameter(0, 0.5);

        ratioHist->Fit(fitFunc, "QRN");
        fitFunc->SetLineColor(colors[j]);
        fitFunc->SetLineStyle(2);
        fitFunc->SetLineWidth(2);
        fitFunc->Write();

        ratioLegend->AddEntry(ratioHist, suffixes[j].c_str(), "PLE");

        gPad->Update();
        gPad->RedrawAxis();

        // Save the ratio histogram to the ROOT file
        ratioHist->Write(("RatioHist_" + suffixes[j]).c_str());
    }

    // Save the primary data ratio histogram to the file
    ratioHist_data->Write("RatioHist_Data");

    ratioLegend->Draw();

    // Save the canvas as a PDF
    c2->SaveAs("pdf/all_ratio_ISR.pdf");

    // Close files
    file->Close();
    outputFile->Write();  // Ensure all objects are written
    outputFile->Close();
}