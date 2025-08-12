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
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <TGraphErrors.h>
#include <TLatex.h>
#include "../minitools/tdrstyle_mod22.C"

void chi_Run2_mass_ud() {
    setTDRStyle();

    // Open the output files
    TFile *file = new TFile("output_MCRun2.root", "READ");
    TFile *file2 = new TFile("output_DATARun2.root", "READ");

    std::vector<double> scaleFactors;
    for (double scale = 0.98; scale <= 1.02; scale += 0.002) {
        scaleFactors.push_back(scale);
    }

    // Create a TH1D histogram to store chi-squared values for scale1
    int nBins = scaleFactors.size();
    double dx = 0.002 / 2;
    TH1D* chi2_hist = new TH1D("chi2_hist", ";R_{ud};#chi^{2}", 
                                nBins, 0.98 - dx, 1.02 + dx);

    for (size_t i = 0; i < scaleFactors.size(); ++i) {
        double scale1 = scaleFactors[i];
        
        // Format the histogram name with 4 decimal places
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << scale1;
        std::string histName = "h3MassFlavorPairs_DATAMC_UD_" + oss.str();

        // Retrieve the histogram
        TH3D* h3MassFlavorPairs_DATAMC_MC = (TH3D*)file->Get(histName.c_str());
        if (!h3MassFlavorPairs_DATAMC_MC) {
            std::cerr << "Histogram " << histName << " not found!" << std::endl;
            continue;
        }

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

    
        TH3D* h3MassFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATAMC");

        TH1D* h3all_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

        TH1D* h3tagcsall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
        TH1D* h3tagudall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
        TH1D* h3tagxall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);


        // Cloning and renaming histograms
        TH1D* h3all_scaled = (TH1D*)h3all->Clone("h3all_scaled");
        TH1D* h3all_scaled2 = (TH1D*)h3all->Clone("h3all_scaled2");

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

        h3tagcsall_data->Scale(1./h3all_data->Integral());
        h3tagudall_data->Scale(1./h3all_data->Integral());
        h3tagxall_data->Scale(1./h3all_data->Integral());

        h3all_data->Scale(1./h3all_data->Integral());

        // Scaling histograms
        double integral = h3all->Integral();
        h3gencstagcs->Scale(1.0 / integral);
        h3gencstagud->Scale(1.0 / integral);
        h3gencstagx->Scale(1.0 / integral);
        
        h3genudtagcs->Scale(1.0 / integral);
        h3genudtagud->Scale(1.0 / integral);
        h3genudtagx->Scale(1.0 / integral);

        h3genxtagcs->Scale(1.0 / integral);
        h3genxtagud->Scale(1.0 / integral);
        h3genxtagx->Scale(1.0 / integral);

        h3tagcsall->Scale(1.0 / integral);
        h3tagudall->Scale(1.0 / integral);
        h3tagxall->Scale(1.0 / integral);

        h3all->Scale(1.0 / integral);

        // Scaling cloned histograms
        //h3gencstagx_scaled->Scale(1.91);
        //h3genudtagx_scaled->Scale(0.55);
        //h3genxtagx_scaled->Scale(0.97);

        //h3gencstagud_scaled->Scale(1.42);
        //h3genudtagud_scaled->Scale(0.71);
        //h3genxtagud_scaled->Scale(0.71);

        //h3gencstagcs_scaled->Scale(0.82);
        //h3genudtagcs_scaled->Scale(0.63);
        //h3genxtagcs_scaled->Scale(0.73);

        // Reset and sum scaled histograms
        h3all_scaled->Reset();
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

        // Final scaling
        double scale_factor2 = 1.0 / h3all_scaled->Integral();

        h3tagcsall_scaled->Scale(scale_factor2);
        h3tagudall_scaled->Scale(scale_factor2);
        h3tagxall_scaled->Scale(scale_factor2);

        h3all_scaled->Scale(scale_factor2);

               // Calculate chi-squared values
    TH1D *h_tagcsvsdata_scaled = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata_scaled");
    h_tagcsvsdata_scaled->Divide(h3tagcsall_data);

    // Fit the ratio histogram with a constant function
    TF1* f7m = new TF1("f7m", "[0]", 70, 100);
    f7m->FixParameter(0, 1);
    h_tagcsvsdata_scaled->Fit(f7m, "RN");
    double chi27m = f7m->GetChisquare();

    // Fill the TH2D histogram for chi-squared (if necessary for plotting)
    chi2_hist->Fill(scale1, chi27m);

    // Apply SetNoExponent only to the X-axis for scale1
    h_tagcsvsdata_scaled->GetXaxis()->SetNoExponent();

    cout << "scale1 = " << scale1 << ", chi27m = " << chi27m << endl << flush;
}


// Draw the 1D chi-squared histogram
setTDRStyle();
lumi_136TeV = "Run3 simulation";
extraText = "Private";
TH1D *h1 = tdrHist("h1","Chi-Squared",3000,3100,"R_{ud}",0.9775,1.0225);
TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
//c1->SetLogy();

tdrDraw(chi2_hist,"P",kFullCircle,kAzure+7);
chi2_hist->SetMarkerSize(1);
// Update the canvas to reflect the changes
// Suppress the error bars
for (int i = 0; i < chi2_hist->GetNbinsX(); ++i) {
    chi2_hist->SetBinError(i+1, 0); // Set bin errors to zero
}
h1->GetXaxis()->SetLabelSize(0.04);
h1->GetYaxis()->SetLabelSize(0.04);
h1->GetXaxis()->SetTitleSize(0.045);
h1->GetYaxis()->SetTitleSize(0.045);

c1->Update();

// Save the canvas
c1->SaveAs("pdf/chi_Run2_mass_ud.pdf");


// Save the chi2_hist to a ROOT file
TFile* outputFile = new TFile("chi2_hist_run2_mass_ud.root", "RECREATE");
chi2_hist->Write();
//g2->Write();
outputFile->Close();

}