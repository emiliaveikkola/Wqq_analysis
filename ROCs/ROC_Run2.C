#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <string>
#include "tdrstyle_mod22.C"


// Helper function to extract the base name from a file path
std::string getBaseFileName(const std::string &filePath) {
    size_t lastSlash = filePath.find_last_of("/");
    size_t lastDot = filePath.find_last_of(".");
    std::string baseName = filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    return baseName; 
}



double interpolate(TGraph* graph, double x) {
    int n = graph->GetN();
    double *xPoints = graph->GetX();
    double *yPoints = graph->GetY();

    if (x < xPoints[0] || x > xPoints[n-1]) return 0;

    for (int i = 0; i < n - 1; i++) {
        if (x >= xPoints[i] && x <= xPoints[i+1]) {
            double slope = (yPoints[i+1] - yPoints[i]) / (xPoints[i+1] - xPoints[i]);
            return yPoints[i] + slope * (x - xPoints[i]);
        }
    }
    return 0;
}

void ROC_Run2() {
    setTDRStyle();
    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

    // List of ROOT files to process
    std::vector<std::string> fileNames = {
        "/Users/macbookpro/Downloads/Run2/Muo16_MC.root",
        "/Users/macbookpro/Downloads/Run2/Mikael/Muo16_MC.root",
        "/Users/macbookpro/Downloads/Run2/Muo16APV_MC.root",
        "/Users/macbookpro/Downloads/Run2/Mikael/Muo16APV_MC.root",
        "/Users/macbookpro/Downloads/Run2/Muo17_MC.root",
        "/Users/macbookpro/Downloads/Run2/Mikael/Muo17_MC.root",
        "/Users/macbookpro/Downloads/Run2/Muo18_MC.root",
        "/Users/macbookpro/Downloads/Run2/Mikael/Muo18_MC.root" // Add more as needed
    };
    // Human-friendly labels corresponding to each file in fileNames
    std::vector<std::string> friendlyNames = {
        "Muo16", "Mikael Muo16",
        "Muo16 APV", "Mikael Muo16 APV",
        "Muo17", "Mikael Muo17",
        "Muo18", "Mikael Muo18"
    };

    // Canvas and base histogram
    TH1D *h2 = tdrHist("h2", "Mis-id rate ", 1e-3-5e-4, 1, "Jet efficiency", 0.079, 1);
    TCanvas *c2 = tdrCanvas("c2", h2, 4, 11, kSquare);
    c2->SetLogy();
    
    // Add custom grid lines
    double gridXPositions[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    double gridYPositions[] = {0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005,
                               0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                               0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    // Get the axis ranges for grid lines
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
    double yMin = h2->GetMinimum();
    double yMax = h2->GetMaximum();

    // Draw vertical grid lines
    for (double gridX : gridXPositions) {
        TLine *line = new TLine(gridX, yMin, gridX, yMax);
        line->SetLineColor(kBlack);  // Set grid color
        line->SetLineStyle(3);       // Dashed line for grid
        line->Draw("same");
    }

    // Draw horizontal grid lines
    for (double gridY : gridYPositions) {
        TLine *line = new TLine(xMin, gridY, xMax, gridY);
        line->SetLineColor(kBlack);
        line->SetLineStyle(3);       // Dashed line for grid
        line->Draw("same");
    }

    // Legend for ROC curves
    TLegend *legend = tdrLeg(0.5, 0.55-0.04*fileNames.size(), 0.8, 0.55);

    // Colors for different ROC curves
    std::vector<int> colors = {kViolet-2, kViolet+2, kGreen+1, kOrange-3, kBlue-4, kAzure+7,kOrange+7, kGreen-2};

    for (size_t i = 0; i < fileNames.size(); ++i) {
        TFile *file = new TFile(fileNames[i].c_str(), "READ");
        if (!file->IsOpen()) {
            std::cerr << "Error: could not open file " << fileNames[i] << std::endl;
            continue;
        }

        TTree *tree = nullptr;
        file->GetObject("tree", tree);

        if (!tree) {
            std::cerr << "Error: could not find tree in file " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        // Histograms for ctag1 and ctag2
        TH1D *h_c_ctag1 = new TH1D("h_c_ctag1", "ctag1 for abs(flav1)==4; ctag1; Fraction", 100, 0, 1);
        TH1D *h_uds_ctag1 = new TH1D("h_uds_ctag1", "ctag1 for abs(flav1)==1 && abs(flav1)==2; ctag1; Fraction", 100, 0, 1);

        TH1D *h_c_ctag2 = new TH1D("h_c_ctag2", "ctag1 for abs(flav1)==4; ctag2; Fraction", 100, 0, 1);
        TH1D *h_uds_ctag2 = new TH1D("h_uds_ctag2", "ctag1 for abs(flav1)==1 && abs(flav1)==2; ctag1; Fraction", 100, 0, 1);

        tree->Project("h_c_ctag1", "ctag1/(udstag1+gtag1+ctag1)", "abs(flav1)==4");
        tree->Project("h_uds_ctag1", "ctag1/(udstag1+gtag1+ctag1)", "abs(flav1)==1 || abs(flav1)==2 || abs(flav1)==3");

        tree->Project("h_c_ctag2", "ctag2/(udstag2+gtag2+ctag2)", "abs(flav2)==4");
        tree->Project("h_uds_ctag2", "ctag2/(udstag2+gtag2+ctag2)", "abs(flav2)==1 || abs(flav2)==2 || abs(flav2)==3");

        TH1D *h_c_ctag3 = (TH1D*)h_c_ctag1->Clone("h_c_ctag3");
        TH1D *h_uds_ctag3 = (TH1D*)h_uds_ctag1->Clone("h_uds_ctag3");

        TH1D *h_c_ctag4 = (TH1D*)h_c_ctag2->Clone("h_c_ctag4");
        TH1D *h_uds_ctag4 = (TH1D*)h_uds_ctag2->Clone("h_uds_ctag4");

        h_c_ctag3->Add(h_c_ctag4);
        h_uds_ctag3->Add(h_uds_ctag4);

        // Normalize histograms
        h_c_ctag3->Scale(1. / h_c_ctag3->Integral());
        h_uds_ctag3->Scale(1. / h_uds_ctag3->Integral());

        // ROC curve
        TGraph *rocCurveCUDS = new TGraph(100);

        double dx5_min(1), x5_020(1), x5_cut_020(1), y5_020(1);
        double dx5_min2(1), x5_080(1), x5_cut_080(1), y5_080(1);

        for (int j = 0; j < 100; ++j) {
            double y5 = h_uds_ctag3->Integral(j, 100);
            double x5 = h_c_ctag3->Integral(j, 100);
            rocCurveCUDS->SetPoint(j, x5, y5);

            // For target efficiency 0.20
            if (fabs(x5 - 0.2) < dx5_min) {
                x5_cut_020 = h_c_ctag3->GetBinLowEdge(j);
                x5_020 = x5;
                dx5_min = fabs(x5 - 0.2);
                y5_020 = y5;
            }

            // For target efficiency 0.80
            if (fabs(x5 - 0.8) < dx5_min2) {
                x5_cut_080 = h_c_ctag3->GetBinLowEdge(j);
                x5_080 = x5;
                dx5_min2 = fabs(x5 - 0.8);
                y5_080 = y5;
            }
        }

        // Print the results for x5_020, y5_020, x5_080, y5_080
        std::cout << "File: " << getBaseFileName(fileNames[i]) << std::endl;
        std::cout << "For ctag > " << x5_cut_020 << " ctag_eff = " << x5_020 << " (target = 0.20)"
                  << " mis-tag = " << y5_020 << " (target = 0.01)" << std::endl;

        std::cout << "For ctag < " << x5_cut_080 << " ctag_eff = " << x5_080 << " (target = 0.80)"
                  << " mis-tag = " << y5_080 << " (target = 0.60)" << std::endl;

        // Set style for ROC curve
        rocCurveCUDS->SetLineColor(colors[i % colors.size()]);
        rocCurveCUDS->SetLineWidth(2);
        rocCurveCUDS->Draw("same");

        // Add to legend with human-friendly label
        std::string legendEntry = friendlyNames[i];
        legend->AddEntry(rocCurveCUDS, legendEntry.c_str(), "L");
        legend->SetTextSize(0.035);

        file->Close();
    }

    // Draw diagonal line (ideal ROC curve)
    TF1 *l = new TF1("l", "x", 0.08, 1);
    l->SetLineColor(kBlack);
    l->Draw("same");

    // Draw legend
    legend->Draw();

    // Update canvas
    c2->Modified();
    c2->Update();

    // Save canvas
    c2->SaveAs("pdf/ROC_Curves_Run2.pdf");
}