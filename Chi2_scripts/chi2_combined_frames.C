#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <TEllipse.h>  // Add the TEllipse header
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <TGraphErrors.h>
#include <TLatex.h>
#include "tdrstyle_mod22.C"

void chi2_combined_frames() {
    setTDRStyle();
    extraText = "Private";

    std::vector<double> udScaleFactors;
    for (double scale_ud = 0.98; scale_ud <= 1.02; scale_ud += 0.002)
        udScaleFactors.push_back(scale_ud);

    for (double scale_ud : udScaleFactors) {
        gROOT->cd(); // Ensure we are in the global directory before deleting
        if (gROOT->FindObject("h2")) delete gROOT->FindObject("h2");
        if (gROOT->FindObject("h2_contours")) delete gROOT->FindObject("h2_contours");
        if (gROOT->FindObject("h")) delete gROOT->FindObject("h");
        if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");

        TFile *file0 = TFile::Open(Form("frames/mass/rootfiles/chi2_3D_mass_cs_ud_%.3f.root", scale_ud));
        TFile *file1 = TFile::Open(Form("frames/pt/rootfiles/chi2_3D_pt_cs_ud_%.3f.root", scale_ud));
        if (!file0 || !file1) continue;

        TH2D *chi2_hist2_saved = (TH2D*)file1->Get(Form("chi2_hist2_saved_%.3f", scale_ud));
        TH2D *chi2_hist_saved = (TH2D*)file0->Get(Form("chi2_hist_saved_%.3f", scale_ud));
        if (!chi2_hist2_saved || !chi2_hist_saved) continue;

        // Clone the first histogram
        TH2D *h2 = (TH2D*)chi2_hist2_saved->Clone("h2");

        // Add the second histogram to the cloned one
        h2->Add(chi2_hist_saved);
        
        // Define effective NDF value for the combined histogram
        const double NDF = 77.0;  // 30 + 47 = 77

        // Loop through all bins to divide the bin content by NDF
        for (int i = 1; i <= h2->GetNbinsX(); ++i) {
            for (int j = 1; j <= h2->GetNbinsY(); ++j) {
                double zValue = h2->GetBinContent(i, j);
                h2->SetBinContent(i, j, zValue / NDF);  // Divide each bin content by NDF
                h2->SetBinError(i, j, sqrt(47. + 30.) / NDF);  // Adjust error accordingly
            }
        }

        // Create a separate histogram for contours only
        TH2D *h2_contours = (TH2D*)h2->Clone("h2_contours");

        // Set contour levels and line properties
        const Int_t NContours = 1;
        Double_t deltaChi2[NContours] = {10};
        Double_t deltaChi2NDF[NContours];
        for (int i = 0; i < NContours; ++i) {
            deltaChi2NDF[i] = deltaChi2[i] / NDF;
        }

        // Set contour levels for the contour-only histogram
        Double_t contours[NContours];
        double minZ = h2->GetMinimum();
        for (int i = 0; i < NContours; ++i) {
            contours[i] = minZ + deltaChi2NDF[i];
            cout << contours[i] << endl;
        }
        h2_contours->SetContour(NContours, contours);

        // Set contour-only histogram properties
        h2_contours->SetLineColor(kWhite);  // Line color for contours
        h2_contours->SetLineWidth(1);
        h2_contours->SetLineStyle(2);     // Set line width
        h2_contours->SetFillColor(0);     // No fill color
        h2_contours->SetFillStyle(0);     // Transparent fill

        // Draw the combined histogram
        TH1D *h = tdrHist("h", "Strange jet response (R_{s}) ", 0.979, 1.021, "Charm jet response (R_{c})",0.979, 1.021);
        if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");
        TCanvas *c1 = tdrCanvas("c1", h, 4, 0, kRectangular);
        gPad->SetRightMargin(0.17);
        gPad->SetLeftMargin(0.13);

        h2->Draw("colz same");
        h2_contours->Draw("CONT3 SAME");  // Overlay the contours without fill
        h2->GetZaxis()->SetRangeUser(3,32);


        // Find the minimum z-value and its corresponding x and y values
        int minBinX = -1, minBinY = -1;
        double minX, minY;

        for (int i = 1; i <= h2->GetNbinsX(); ++i) {
            for (int j = 1; j <= h2->GetNbinsY(); ++j) {
                double z = h2->GetBinContent(i, j);
                if (z == minZ) {
                    minBinX = i;
                    minBinY = j;
                    break;
                }
            }
        }

        double maxZ = h2->GetMaximum();
        if (minBinX != -1 && minBinY != -1) {
            minX = h2->GetXaxis()->GetBinCenter(minBinX);
            minY = h2->GetYaxis()->GetBinCenter(minBinY);
            std::cout << "Maximum Z-value: " << maxZ << std::endl;
            std::cout << "Minimum Z-value: " << minZ << std::endl;
            std::cout << "Corresponding X-value: " << minX << std::endl;
            std::cout << "Corresponding Y-value: " << minY << std::endl;
        } else {
            std::cerr << "Error: Could not find the minimum z-value in the histogram." << std::endl;
        }

        h2->GetZaxis()->SetTitle("#frac{#chi^{2}}{NDF}");

        h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetLabelSize(0.045);
        h2->GetZaxis()->SetLabelSize(0.045);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h2->GetZaxis()->SetTitleSize(0.045);

        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.4);
        h2->GetZaxis()->SetTitleOffset(1.1);

        gPad->RedrawAxis();
        gPad->Update();

        TLine *l = new TLine();
        l->DrawLine(0.98-0.001,0.98-0.001,1.02+0.001,1.02+0.001);
        l->SetLineStyle(kDotted);
        l->DrawLine(1,0.98-0.001,1,1.02+0.001);
        l->DrawLine(0.98-0.001,1.0,1.02+0.001,1.0);;

        TLine *l2 = new TLine();
        l2->DrawLine(0.98-0.001,1.02+0.001,1.02+0.001,0.98-0.001);
        l2->SetLineStyle(kDotted);
        l2->DrawLine(1,0.98-0.001,1,1.02+0.001);
        l2->DrawLine(0.98-0.001,1.0,1.02+0.001,1.0);

        // 2D fit to the full histogram to find the global minimum
        TF2* f2 = new TF2("f2", "[0]*x*x + [1]*y*y + [2]*x*y + [3]*x + [4]*y + [5]", 0.98, 1.02, 0.98, 1.02);
        h2->Fit("f2", "Q0");
        double p0 = f2->GetParameter(0);
        double p1 = f2->GetParameter(1);
        double p2 = f2->GetParameter(2);
        double p3 = f2->GetParameter(3);
        double p4 = f2->GetParameter(4);

        double denom = 4*p0*p1 - p2*p2;
        double x_star = (p2*p4 - 2*p1*p3) / denom;
        double y_star = (p2*p3 - 2*p0*p4) / denom;
        double chi2_min_fit = f2->Eval(x_star, y_star);
        TFile* fitOutFile = new TFile(Form("frames/combined/rootfiles/fit_results_ud_%.3f.root", scale_ud), "RECREATE");
        TH2D* h2_fit_points = (TH2D*)h2->Clone(Form("fit_input_ud_%.3f", scale_ud));
        h2_fit_points->SetDirectory(0);  // Make independent clone
        TF2* fitSaved = (TF2*)f2->Clone(Form("fit2D_ud_%.3f", scale_ud));
        fitSaved->Write();
        h2_fit_points->Write();
        fitOutFile->Close();
        delete fitOutFile;

        // White star under fitted one
        TGraph* g_old = new TGraph();
        g_old->SetPoint(0, minX,minY);
        g_old->SetMarkerStyle(kStar);
        g_old->SetMarkerSize(2);
        g_old->SetMarkerColor(kWhite);
        g_old->Draw("same p");

        // Add a second star to mark the parabola-fit minimum
        TGraph* g_fit = new TGraph();
        g_fit->SetPoint(0, x_star, y_star);
        g_fit->SetMarkerStyle(kFullStar);
        g_fit->SetMarkerColor(kRed);
        g_fit->SetMarkerSize(2.5);
        g_fit->Draw("same p");

         // ---- Begin contour extraction for the 2D fit ----
        // Define a contour level above the minimum; adjust deltaLevel as needed.
        double deltaLevel = 10./NDF; 
        double level = chi2_min_fit + deltaLevel;
        // Set the contour level for f2
        f2->SetContour(1, &level);
        f2->SetLineColor(kRed);
        f2->SetLineStyle(2);
        f2->SetLineWidth(2);
        f2->Draw("CONT3 LIST SAME");
        gPad->Update();

        // Create a new transparent pad on top of the current plot
        TPad *overlayPad = new TPad("overlayPad", "", 0, 0, 1, 1);
        overlayPad->SetFillStyle(0);  // Make it transparent
        overlayPad->SetFrameFillStyle(0);
        overlayPad->Draw();
        overlayPad->cd(); 

        // Create a TLatex object for annotation
        TLatex latex;
        latex.SetTextSize(0.04); // Adjust size according to your plot's appearance
        latex.DrawLatex(0.5, 0.86, "p_{T} > 30 GeV");
        latex.DrawLatex(0.5, 0.81, "|#eta| < 2.5");
        latex.DrawLatex(0.5, 0.755, "#frac{#chi^{2}}{NDF}, NDF = 77");
        latex.DrawLatex(0.32, 0.86, Form("R_{ud} = %.3f", scale_ud));
        latex.DrawLatex(0.25, 0.81, Form("min #chi^{2}/NDF = %.2f", chi2_min_fit));
        TLatex latexOld;
        latexOld.SetTextSize(0.04);
        latexOld.SetTextColor(kWhite);  // Or any other ROOT color code
        latexOld.DrawLatex(0.2, 0.76, Form("old min #chi^{2}/NDF = %.2f", minZ));

        gPad->Modified();
        gPad->Update();

        // Save the canvas
        c1->SaveAs(Form("frames/combined/chi2_combined_frame_ud_%.3f.png", scale_ud));

        delete h2;
        delete h2_contours;
        delete h;
        delete c1;
    }
    gSystem->Exec("magick -delay 50 -loop 0 frames/combined/chi2_combined_frame_ud_*.png frames/combined/chi2_combined_scan.gif");
    std::cout << "[INFO] chi2_combined_scan.gif created in frames/combined/" << std::endl;
}