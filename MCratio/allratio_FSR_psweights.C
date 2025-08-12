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

double gFunction(double *x, double *params) {
    // Zero constraints
    /*double zero_factor = (x[0] - 81) * (x[0] - 99);

    // Lower range: Polynomial for x < 82
    if (x[0] < 82) {
        double shifted_x = x[0] - 81; // Shift to stabilize computation
        return zero_factor * (params[0] * pow(shifted_x, 2) + params[1] * shifted_x + params[2]);
    }
    // Upper range: Polynomial for x >= 82
    else {
        double shifted_x = x[0] - 90; // Shift to center around 90
        return zero_factor * (params[3] * pow(shifted_x, 2) + params[4] * shifted_x + params[5]);
    }
    */
   return 1;
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
    TF1 *f1 = new TF1(("f1_" + suffix).c_str(), fFunction, 66, 104, 7);  // Fit range from 66 to 104 GeV
    f1->SetParameters(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5); // Initial guesses

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

void allratio_FSR_psweights() {
    setTDRStyle();
    extraText = "Private";

    // Open the file containing the processed histograms
    TFile* file = new TFile("processed_histograms/psweights.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open processed_histograms.root." << std::endl;
        return;
    }

    // Open the file containing the processed histograms
    TFile* file2 = new TFile("processed_histograms/all.root", "READ");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Could not open processed_histograms.root." << std::endl;
        return;
    }

        // Create a new ROOT file to store results
    TFile* outputFile = new TFile("processed_histograms/histograms_FSR_psweights.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create results.root." << std::endl;
        file->Close();
        return;
    }

    // Define suffixes, legend labels, colors, and histogram name
    std::vector<std::string> suffixes = {"fsr:murfac=0.25", "fsr:murfac=0.5", "fsr:murfac=2.0", "fsr:murfac=4.0"};
    std::vector<std::string> legendNames = {"FSR (0.25)", "FSR (0.5)", "FSR (2.0)", "FSR (4.0)"};
    std::vector<int> colors = {kRed-4, kOrange+1, kSpring-5, kAzure+7, kViolet-5};
    std::vector<int> markers = {kOpenCircle, kFullSquare, kOpenSquare, kFullDiamond, kOpenDiamond};
    std::string histName = "h3all";

    TH1D *h2 = tdrHist("h2", "MC/MC_{ref}", 0.55,1.6, "Mass (GeV)", 66,104);
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
    TH1D* hist_data = (TH1D*)file2->Get((histName + "_data").c_str());
    TH1D* hist_unscaled = (TH1D*)file2->Get((histName + "_unscaled").c_str());
    if (!hist_data || !hist_unscaled) {
        std::cerr << "Error: Histograms for data or unscaled MC not found." << std::endl;
        file2->Close();
        return;
    }

    // Add a legend for the ratio plot
    TLegend* ratioLegend = tdrLeg(0.52, 0.75, 0.91, 0.9);
    ratioLegend->SetTextSize(0.03);

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

        tdrDraw(ratioHist, "Pz", markers[j], colors[j]);
        ratioHist->SetMarkerSize(1);
        ratioHist->SetLineWidth(2);

        // Call the fit function with the specific color for this histogram
        //fitFSRWithFixedParams(ratioHist, suffixes[j], a, b, c, colors[j]);
        // Fit the ratio histogram and save the fit function
        TF1* fitFunc = new TF1(("FitFunction_" + suffixes[j]).c_str(), fFunction, 66, 104, 6);
        fitFunc->FixParameter(1, a);
        fitFunc->FixParameter(2, b);
        fitFunc->FixParameter(3, c);
        fitFunc->SetParameter(0, 0.5);

        ratioHist->Fit(fitFunc, "QRN");
        fitFunc->SetLineColor(colors[j]);
        fitFunc->SetLineStyle(2);
        fitFunc->SetLineWidth(2);
        if (suffixes[j] == "fsr:murfac=2.0") {
            ratioHist->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
        }
        if (suffixes[j] == "fsr:murfac=4.0") {
            ratioHist->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
        }
        fitFunc->Write();

        ratioLegend->AddEntry(ratioHist, legendNames[j].c_str(), "PLE");

        gPad->Update();
        gPad->RedrawAxis();

        // Save the ratio histogram to the ROOT file
        ratioHist->Write(("RatioHist_" + suffixes[j]).c_str());
    }

    // --- Configurable extra histogram block ---
    std::vector<std::string> extraNames       = {"h3all_FSR0995", "h3all_FSR1005"};
    std::vector<std::string> extraLegendNames = {"FSR (0.995)", "FSR (1.005)"};
    std::vector<int>         extraColors      = {kMagenta-4, kViolet+1};
    std::vector<int>         extraMarkers     = {kFullStar, kOpenStar};

    if (file2 && !file2->IsZombie()) {
        for (size_t k = 0; k < extraNames.size(); ++k) {
            TH1D* hist_extra = (TH1D*)file2->Get(extraNames[k].c_str());
            if (hist_extra) {
                std::string ratioName = "RatioHist_" + extraNames[k];
                TH1D* ratioHist_extra = (TH1D*)hist_extra->Clone(ratioName.c_str());
                ratioHist_extra->Divide(hist_unscaled);

                // Style and draw the extra histogram
                tdrDraw(ratioHist_extra, "Pz", extraMarkers[k], extraColors[k]);
                ratioHist_extra->SetMarkerSize(1);
                ratioHist_extra->SetLineWidth(2);
                ratioLegend->AddEntry(ratioHist_extra, extraLegendNames[k].c_str(), "PLE");
                ratioHist_extra->Write(ratioName.c_str());

                if (extraNames[k] == "h3all_FSR0995") {
                    ratioHist_extra->SetMarkerSize(1.5);  // Specific range for "h3tagxall"
                }
                if (extraNames[k] == "h3all_FSR1005") {
                    ratioHist_extra->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
                }
            } else {
                std::cerr << "Error: Histogram '" << extraNames[k]
                          << "' not found in processed_histograms_extra.root." << std::endl;
            }
        }
    }
    // --- End extra block ---

    // Save the primary data ratio histogram to the file
    ratioHist_data->Write("RatioHist_Data");

    ratioLegend->Draw();

    // Save the single ratio histogram as a PDF
    c2->SaveAs("pdf/all_ratio_FSR_psweights.pdf");

    // Close files
    if (file2 && !file2->IsZombie()) {
        file2->Close();
    }
    file->Close();
    outputFile->Write();  // Ensure all objects are written
    outputFile->Close();
}