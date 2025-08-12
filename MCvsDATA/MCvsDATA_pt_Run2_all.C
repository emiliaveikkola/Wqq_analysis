#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

// Main function to plot histograms and their ratios with correct divisions from separate files
void MCvsDATA_pt_Run2_all() {
    setTDRStyle();

    // Open the file containing primary histograms with suffixes
    TFile* file = new TFile("../processed_histograms/processed_histograms_pt.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open processed_histograms.root." << std::endl;
        return;
    }

    // Define histogram names, suffixes, and colors
    std::vector<std::string> histNames = {"h3tagcsall", "h3tagudall", "h3tagxall", "h3all"};
    std::vector<std::string> suffixes = {"JER110"};
    std::vector<int> colors = {kSpring-5};
    std::vector<int> markers = {kOpenSquare};
    //std::vector<std::string> suffixes = {"FSR0990", "JER110", "ISR130"};
    //std::vector<int> colors = {kRed-4, kSpring-5, kViolet-5 };
    //std::vector<int> colors = {kRed-4, kOrange+1, kSpring-5, kAzure+7, kViolet-5 };
    //std::vector<int> markers = {kOpenCircle, kOpenSquare, kOpenDiamond};
    //std::vector<int> markers = {kOpenCircle, kFullSquare, kOpenSquare, kFullDiamond, kOpenDiamond};
    //std::vector<std::string> suffixes = {"JER80", "JER90", "JER105", "JER110", "JER120"};
    //std::vector<std::string> suffixes = {"ISR40", "ISR62", "ISR70", "ISR77", "ISR120", "ISR130", "ISR160"};
    //std::vector<std::string> suffixes = {"ISR40", "ISR70", "ISR120", "ISR130", "ISR160"};
    //std::vector<int> colors2 = {kRed-4, kOrange+1, kSpring-5, kAzure+7, kViolet-5, kMagenta-7, kCyan-7};

    // Create a canvas with a 2x4 grid layout
    setTDRStyle();
    TCanvas* c1 = new TCanvas("c1", "Histograms and Ratios", 2560, 1140);

    c1->Divide(4, 2);  // 4 columns, 2 rows

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
    

    // Loop over histograms to display them in the first row and their ratios in the second row
    for (size_t i = 0; i < histNames.size(); ++i) {
        c1->cd(i + 1);  // Move to the correct pad in the first row
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.01);
        gPad->SetBottomMargin(0.1);
        gPad->SetTopMargin(0.1);

        // Create a standard TLegend for each pad in the first row
        TLegend* legend = new TLegend(0.75, 0.65-(suffixes.size()+1)*0.02, 0.95, 0.8);  // Position in the top right corner
        legend->SetTextSize(0.045);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);  // Transparent legend background

        TH1D* hist_data = (TH1D*)file->Get((histNames[i] + "_data").c_str());
        if (!hist_data) {
            std::cerr << "Error: Could not retrieve data histogram " << histNames[i] + "_data" << std::endl;
            continue;
        }

        TH1D* hist_unscaled = (TH1D*)file->Get((histNames[i] + "_unscaled").c_str());
        if (!hist_unscaled) {
            std::cerr << "Error: Could not retrieve unscaled histogram " << histNames[i] + "_unscaled" << std::endl;
            continue;
        }

        tdrDraw(hist_unscaled, "HPz", kNone, kGray+2, kSolid, -1, 1001, kGray+2);
        hist_unscaled->SetFillColorAlpha(kGray+2, 0.25);

        legend->AddEntry(hist_unscaled, "MC", "f");  // "f" for fill style

        TLatex tex;
        tex.SetNDC();
        // Set font and size if needed (optional)
        tex.SetTextSize(0.04); // Adjust size according to your plot's appearance
        // Add text for pT region (adjust coordinates based on plot needs)
        tex.DrawLatex(0.18, 0.8, "p_{T} > 30 GeV");
        // Add text for eta region (adjust coordinates based on plot needs)
        tex.DrawLatex(0.18, 0.75, "|#eta| < 2.5");

        // Calculate and print RMS for each histogram in the "all" category
        double rms_data = hist_data->GetRMS();
        double rms_unscaled = hist_unscaled->GetRMS();
        std::cout << "RMS of " << histNames[i] << "_data: " << rms_data << std::endl;
        std::cout << "RMS of " << histNames[i] << "_unscaled: " << rms_unscaled << std::endl;

        
        for (size_t j = 0; j < suffixes.size(); ++j) {
            std::string histFullName = histNames[i] + "_" + suffixes[j];
            // Retrieve the primary histogram with suffix from the processed file
            TH1D* hist = (TH1D*)file->Get(histFullName.c_str());
            if (!hist) {
                std::cerr << "Warning: Histogram " << histFullName << " not found in processed_histograms.root." << std::endl;
                continue;
            }

            // Draw the first histogram with "HPz" to set the axis
            tdrDraw(hist, "HPz", kNone, colors[j], kSolid, -1, 1001, colors[j]);
            hist->SetFillColorAlpha(colors[j], 0.25);

            hist_unscaled->GetXaxis()->SetTitle("(s-c)/(s+c)");
            hist_unscaled->GetYaxis()->SetTitle("N frac");
            hist_unscaled->GetYaxis()->SetLabelSize(0.045);
            hist_unscaled->GetXaxis()->SetLabelSize(0.045);
            hist_unscaled->GetXaxis()->SetRangeUser(-1, 1.0);
            hist_unscaled->GetYaxis()->SetTitleSize(0.045);
            hist_unscaled->GetXaxis()->SetTitleSize(0.045);
            hist_unscaled->GetYaxis()->SetTitleOffset(1.65);
            
            // Set specific y-axis range for "h3tagxall" histograms with suffix[j]
            if (histNames[i] == "h3tagxall" && suffixes[j] == suffixes[j]) {
                hist_unscaled->GetYaxis()->SetRangeUser(0, 0.035);  // Specific range for "h3tagxall"
            }

            // Add this histogram to the legend
            legend->AddEntry(hist, suffixes[j].c_str(), "f");  // "f" for fill style
            
             // Calculate and print RMS for each histogram with suffix
            double rms = hist->GetRMS();
            std::cout << "RMS of " << histFullName << ": " << rms << std::endl;
        }
        tdrDraw(hist_data,"Pz",kFullCircle,kPink-9);
        hist_data->SetMarkerSize(1);
        
        legend->AddEntry(hist_data, "DATA", "PLE");  // "f" for fill style

        // Draw the legend on the pad
        legend->Draw();

        // Calculate the ratio
        c1->cd(i + 5);  // Move to the correct pad in the second row
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.01);
        gPad->SetBottomMargin(0.1);
        gPad->SetTopMargin(0);
        
        // Create a legend for the second row
        TLegend* ratioLegend = new TLegend(0.75, 0.78, 0.95, 0.94);  // Position in the top right corner
        ratioLegend->SetTextSize(0.045);
        ratioLegend->SetBorderSize(0);
        ratioLegend->SetFillStyle(0);  // Transparent background

        // Clone the primary histogram and divide by its specific `_data` histogram
        TH1D* ratioHist_unscaled = (TH1D*)hist_unscaled->Clone();  // Fresh clone
        ratioHist_unscaled->Divide(hist_data);

        //TH1D* ratioHist_data = (TH1D*)hist_data->Clone();  // Fresh clone
        //ratioHist_data->Divide(hist_unscaled);

        tdrDraw(ratioHist_unscaled, "Pz", kFullCircle, kGray+2);
        ratioHist_unscaled->SetMarkerSize(1.3);

        //tdrDraw(ratioHist_data, "Pz", kFullCircle, kGray+2);
        //ratioHist_data->SetMarkerSize(1.25);

        // Fit the ratio histogram with a constant function
        TF1* f2 = new TF1("f2", "[0]", -0.7, 0.7);
        f2->FixParameter(0, 1);
        ratioHist_unscaled->Fit(f2, "RN");
        //ratioHist_data->Fit(f2, "RN");
        double chi2_2 = f2->GetChisquare();
        int ndf2 = f2->GetNDF();
        
        tex.SetNDC();
        tex.SetTextSize(0.04);
        tex.DrawLatex(0.4, 0.909, Form("#chi^{2} / NDF = %.1f / %d", chi2_2, ndf2));
        // Add each ratio histogram to the ratio legend
        ratioLegend->AddEntry(ratioHist_unscaled, "MC", "PLE");  // "p" for marker style
        //ratioLegend->SetNColumns(2);
        //ratioLegend->AddEntry(ratioHist_data, "DATA", "PLE");  // "p" for marker style

        // Draw the reference line at y=1 after all histograms to ensure visibility
        TLine* line = new TLine(-0.75, 1, 0.75, 1);
        line->SetLineColor(kGray);
        line->SetLineStyle(kDashed);
        line->Draw("SAME");

        for (size_t j = 0; j < suffixes.size(); ++j) {
            std::string histFullName = histNames[i] + "_" + suffixes[j];

            // Retrieve the histogram again and clone for ratio calculation
            TH1D* hist = (TH1D*)file->Get(histFullName.c_str());

            // Clone the primary histogram and divide by its specific `_data` histogram
            TH1D* ratioHist = (TH1D*)hist->Clone();  // Fresh clone
            ratioHist->Divide(hist_data);

            tdrDraw(ratioHist, "Pz", markers[j], colors[j]);
            ratioHist->SetMarkerSize(1);
            if (histNames[i] == histNames[i] && suffixes[j] == "ISR120") {
                ratioHist->SetMarkerSize(2);  // Specific range for "h3tagxall"
            }
            if (histNames[i] == histNames[i] && suffixes[j] == "JER105") {
                ratioHist->SetMarkerSize(1.3);  // Specific range for "h3tagxall"
            }
            if (histNames[i] == histNames[i] && suffixes[j] == "ISR130") {
                ratioHist->SetMarkerSize(1.25);  // Specific range for "h3tagxall"
            }

            
            ratioHist_unscaled->GetXaxis()->SetTitle("(s-c)/(s+c)");
            ratioHist_unscaled->GetYaxis()->SetTitle("MC/DATA");
            ratioHist_unscaled->GetYaxis()->SetLabelSize(0.045);
            ratioHist_unscaled->GetXaxis()->SetLabelSize(0.045);
            ratioHist_unscaled->GetXaxis()->SetRangeUser(-0.75, 0.75);
            ratioHist_unscaled->GetYaxis()->SetTitleSize(0.045);
            ratioHist_unscaled->GetXaxis()->SetTitleSize(0.045);
            ratioHist_unscaled->GetYaxis()->SetTitleOffset(1.6);
            /*
            ratioHist_data->GetXaxis()->SetTitle("Mass (GeV)");
            ratioHist_data->GetYaxis()->SetTitle("MC/MC_{ref}");
            ratioHist_data->GetXaxis()->SetRangeUser(66, 104);
            //ratioHist_data->GetYaxis()->SetRangeUser(0.7, 1.35);
            ratioHist_data->GetYaxis()->SetLabelSize(0.045);
            ratioHist_data->GetXaxis()->SetLabelSize(0.045);
            ratioHist_data->GetYaxis()->SetTitleSize(0.045);
            ratioHist_data->GetXaxis()->SetTitleSize(0.045);
            ratioHist_data->GetYaxis()->SetTitleOffset(1.5);
            */
            if (histNames[i] == "h3tagcsall" && suffixes[j] == suffixes[j]) {
                ratioHist_unscaled->GetYaxis()->SetRangeUser(0.7, 2.05);  // Specific range for "h3tagxall"
            }

            if (histNames[i] == "h3tagudall" && suffixes[j] == suffixes[j]) {
                ratioHist_unscaled->GetYaxis()->SetRangeUser(0.7, 1.9);  // Specific range for "h3tagxall"
            }

            if (histNames[i] == "h3tagxall" && suffixes[j] == suffixes[j]) {
                ratioHist_unscaled->GetYaxis()->SetRangeUser(0.75, 1.21);  // Specific range for "h3tagxall"
            }

            if (histNames[i] == "h3all" && suffixes[j] == suffixes[j]) {
                ratioHist_unscaled->GetYaxis()->SetRangeUser(0.85, 1.2);  // Specific range for "h3tagxall"
            }

            // Fit the ratio histogram with a constant function
            TF1* f = new TF1("f", "[0]", -0.7, 0.7);
            f->FixParameter(0, 1);
            ratioHist->Fit(f, "RN");
            double chi2 = f->GetChisquare();
            int ndf = f->GetNDF();

            // Draw chi2/NDF on the plot
            TLatex tex;
            tex.SetNDC();
            tex.SetTextSize(0.04);
            double yPosition = 0.87 - j * 0.04;  // Adjust the y position for each suffix to avoid overlap
            tex.DrawLatex(0.4, yPosition, Form("#chi^{2} / NDF = %.1f / %d", chi2, ndf));
            
            // Add each ratio histogram to the ratio legend
            ratioLegend->AddEntry(ratioHist, suffixes[j].c_str(), "PLE");  // "p" for marker style
  
        }
        // Draw the legend on the ratio pad
        ratioLegend->Draw();

        // Set font and size if needed (optional)
        tex.SetTextSize(0.04); // Adjust size according to your plot's appearance
        // Add text for pT region (adjust coordinates based on plot needs)
        tex.DrawLatex(0.18, 0.9, "p_{T} > 30 GeV");
        // Add text for eta region (adjust coordinates based on plot needs)
        tex.DrawLatex(0.18, 0.85, "|#eta| < 2.5");

        gPad->Update();
        gPad->RedrawAxis();
    }

    // Save the canvas to a PDF
    c1->SaveAs("pdf/MCvsDATA_pt_Run2_JER10.pdf");

    // Clean up
    file->Close();
    delete file;
}