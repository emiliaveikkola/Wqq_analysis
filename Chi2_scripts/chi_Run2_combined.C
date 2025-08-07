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

void chi_Run2_combined() {
    setTDRStyle();
    extraText = "Private";

  // Load ROOT libraries
  TFile *file0 = TFile::Open("chi2_hist_run2_pt_chi_FSR_after.root");
  TFile *file1 = TFile::Open("chi2_hist_run2_mass_chi_FSR_after.root");

  if (!file0 || !file1) {
    std::cerr << "Error: Could not open one or both files." << std::endl;
    return;
  }

  // Access histograms from the files
  TH2D *chi2_hist2_saved = (TH2D*)file0->Get("chi2_hist2_saved");
  TH2D *chi2_hist_saved = (TH2D*)file1->Get("chi2_hist_saved");

  if (!chi2_hist2_saved || !chi2_hist_saved) {
    std::cerr << "Error: Could not retrieve one or both histograms." << std::endl;
    return;
  }

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
            h2->SetBinError(i, j, sqrt(47 + 30) / NDF);  // Adjust error accordingly
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
    h2_contours->SetLineColor(kRed);  // Line color for contours
    h2_contours->SetLineWidth(1);
    h2_contours->SetLineStyle(2);     // Set line width
    h2_contours->SetFillColor(0);     // No fill color
    h2_contours->SetFillStyle(0);     // Transparent fill

  // Draw the combined histogram
  TH1D *h = tdrHist("h", "Strange jet response (R_{s}) ", 0.979, 1.021, "Charm jet response (R_{c})",0.979, 1.021);
  TCanvas *c1 = tdrCanvas("c1", h, 4, 0, kRectangular);
  gPad->SetRightMargin(0.17);
  gPad->SetLeftMargin(0.13);

  h2->Draw("colz same");
  h2_contours->Draw("CONT3 SAME");  // Overlay the contours without fill


    // Find the minimum z-value and its corresponding x and y values
  int minBinX = -1, minBinY = -1;
  //double minZ = h2->GetMinimum();
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
/*
// Draw the histogram with contour lines
h2->Draw("CONT3 same"); // Or any other suitable option


        // Define the 2D quadratic fit function
TF2 *fitFunc = new TF2("fitFunc",
    "[5] + [0]*(x - [1])*(x - [1]) + [2]*(y - [3])*(y - [3]) + [4]*(x - [1])*(y - [3])",
    h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
    h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
    // Use the previously found minimum values as initial guesses
fitFunc->SetParameters(
    3.405e4,  // a (initial curvature in x)
    minX,    // x0 (from your minimum search)
    1.002,  // b (initial curvature in y)
    minY,    // y0 (from your minimum search)
    0.0,     // c (initial cross-term)
    minZ     // chi^2_min (from your minimum search)
);

// **Calculate contour levels based on the fit**
double chi2_min = fitFunc->GetParameter(5);
const int nContours = 1;
double deltaChi2[nContours] = {0.2};
double contours[nContours];

for (int i = 0; i < nContours; ++i) {
    contours[i] = chi2_min + deltaChi2[i];
    std::cout << "Contour level " << i << ": " << contours[i] << std::endl;
}
// Check the fit status
int fitStatus = h2->Fit("fitFunc", "R");
if (fitStatus != 0) {
    std::cerr << "Fit did not converge properly. Fit status: " << fitStatus << std::endl;
    return;
}

// **Set contour levels on the fit function**
fitFunc->SetContour(nContours, contours);

// **Adjust line properties for the contours**
fitFunc->SetLineColor(kBlack);
fitFunc->SetLineWidth(2);

// **Draw the contours from the fit function**
fitFunc->Draw("CONT1 SAME");
*/

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



TGraph *g = new TGraph();
g->SetPoint(0,minX,minY);
g->SetMarkerStyle(kFullStar);
g->SetMarkerColor(kOrange+10);
g->SetMarkerSize(2.5);
g->Draw("same p");


// Create a new transparent pad on top of the current plot
TPad *overlayPad = new TPad("overlayPad", "", 0, 0, 1, 1);
overlayPad->SetFillStyle(0);  // Make it transparent
overlayPad->SetFrameFillStyle(0);
overlayPad->Draw();
overlayPad->cd(); 

// Create a TLatex object for annotation
TLatex latex;
// Set font and size if needed (optional)
latex.SetTextSize(0.04); // Adjust size according to your plot's appearance
// Add text for pT region (adjust coordinates based on plot needs)
latex.DrawLatex(0.5, 0.86, "p_{T} > 30 GeV");
// Add text for eta region (adjust coordinates based on plot needs)
latex.DrawLatex(0.5, 0.81, "|#eta| < 2.5");
latex.DrawLatex(0.5, 0.755, "#frac{#chi^{2}}{NDF}, NDF = 77");

gPad->Modified();
gPad->Update();

  // Save the canvas
  c1->SaveAs("pdf/chi_run2_combined_colz_chi_FSR_after.pdf");

TFile* outputFile = new TFile("chi2_combinedhist.root", "RECREATE");
h2->Write();

  // Create a TPolyMarker3D to plot the minimum point
TPolyMarker3D *minMarker = new TPolyMarker3D(1);  // '1' for one point
minMarker->SetPoint(0, minX, minY, 300);         // Set to (minX, minY, minZ)
minMarker->SetMarkerStyle(23);                    // Star marker style
minMarker->SetMarkerColor(kRed);                  // Use a contrasting color
minMarker->SetMarkerSize(2.5);                    // Increase size for visibility

// Draw the histogram with LEGO2 option
TCanvas *c2 = new TCanvas("c2", "LEGO2 Histogram", 1000, 600);
h2->Draw("lego2");


// Adjust axis titles
h2->GetXaxis()->SetTitleOffset(1.2);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetZaxis()->SetTitleOffset(1.2);
c2->SetRightMargin(1.5);
h2->GetXaxis()->SetLabelSize(0.04);
h2->GetYaxis()->SetLabelSize(0.04);

// Draw the marker on the same canvas
minMarker->Draw("same");

// Adjust the view angles if necessary
gPad->SetTheta(35);   // Adjust as needed
gPad->SetPhi(-35);    // Adjust as needed

// Update the canvas to ensure everything is drawn
gPad->Modified();
gPad->Update();

// Save the canvas
c2->SaveAs("pdf/chi_run2_combined_lego_chi_FSR_after.pdf");

// Create a new canvas for the logarithmic plot
TCanvas *c3 = new TCanvas("c3", "LEGO2 Histogram with Log Z-axis", 1000, 600);

// Clone the original histogram
TH2D* h2_log = (TH2D*)h2->Clone("h2_log");

// Apply logarithmic transformation to bin contents
for (int i = 1; i <= h2_log->GetNbinsX(); ++i) {
    for (int j = 1; j <= h2_log->GetNbinsY(); ++j) {
        double z = h2_log->GetBinContent(i, j);
        if (z > 0) {
            h2_log->SetBinContent(i, j, log10(z));
        } else {
            // Handle zero or negative values
            h2_log->SetBinContent(i, j, log10(1e-10));
        }
    }
}

// Update the minimum Z-value for the marker
double minZ_log = log10(minZ);

// Create a TPolyMarker3D to plot the minimum point
TPolyMarker3D *minMarker_log = new TPolyMarker3D(1);
minMarker_log->SetPoint(0, minX, minY, minZ_log);
minMarker_log->SetMarkerStyle(23);      // Star marker
minMarker_log->SetMarkerColor(kRed);    // Contrasting color
minMarker_log->SetMarkerSize(2.5);      // Larger size for visibility

// Draw the transformed histogram with LEGO2 option
h2_log->Draw("lego2");

// Adjust the view angles if necessary
gPad->SetTheta(35); // Adjust as needed
gPad->SetPhi(-35);  // Adjust as needed

// Adjust axis titles
h2->GetXaxis()->SetTitleOffset(1.2);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetZaxis()->SetTitleOffset(1.2);
c2->SetRightMargin(1.5);
h2->GetXaxis()->SetLabelSize(0.04);
h2->GetYaxis()->SetLabelSize(0.04);

// Draw the marker on the same canvas
minMarker_log->Draw("same");

// Update and save the canvas
gPad->Modified();
gPad->Update();
c3->SaveAs("pdf/chi_run2_combined_lego_logz_chi_FSR_after.pdf");
}