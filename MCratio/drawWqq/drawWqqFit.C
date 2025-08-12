// Purpose: Draw Wqq variations for ISR, FSR and JER from Emilia with data
//          Perform TF1 fit to Data/MC ratio given these variations
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <TLatex.h>
#include "../minitools/tdrstyle_mod22.C"


// Define the fit function: Gaussian plus linear background
Double_t fit_func(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t amplitude = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    Double_t amplitude2 = par[3];
    Double_t mean2 = par[4];
    Double_t sigma2 = par[5];
    Double_t slope = par[6];
    Double_t intercept = par[7];

    
    // Gaussian function
    Double_t gauss = amplitude * TMath::Gaus(xx, mean, sigma, kTRUE);
    Double_t gauss2 = amplitude2 * TMath::Gaus(xx, mean2, sigma2, kTRUE);
    Double_t linear = intercept + slope * (xx - 83) / 10. ;
    
    Double_t fitval = (gauss + gauss2) * (1 + linear * 0.01);
    return fitval;
}


void drawWqqFit() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open the output files
  TFile *file1 = new TFile("output_MCSCALEDRun2_tagprobe_after.root", "READ");
  TFile *file2 = new TFile("output_DATARun2_tagprobe_before.root", "READ");

  curdir->cd();
  
  // 40,160 or 70,130
  TH3D* h3MassFlavorPairs_DATAMC_MC = (TH3D*)file1->Get("h3MassFlavorPairs_DATAMC");
  TH3D* h3MassFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATAMC");

  TH1D* h3all_mc = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3all_mc", 1, 3, 1, 3);
  TH1D* h3all_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

  //h3all_data->Rebin(2);
  //h3all_mc->Rebin(2);

    // **New Section: Set contents of specific bins to zero**
    std::vector<double> x_values_to_zero = {86.5, 87.5, 88.5, 89.5};

    // Loop over the x-values and set the corresponding bins to zero
    /*for (size_t i = 0; i < x_values_to_zero.size(); ++i) {
        double x = x_values_to_zero[i];
        int bin = h3all_data->FindBin(x);
        h3all_data->SetBinContent(bin, 0.0);
        h3all_data->SetBinError(bin, 0.0);
        h3all_mc->SetBinContent(bin-1, 0.0);
        h3all_mc->SetBinError(bin-1, 0.0);
    }
*/
  h3all_mc->Scale(1./h3all_mc->Integral());
  h3all_data->Scale(1./h3all_data->Integral());

  // TH1D *h = tdrHist("h","Ratio to MC",0.68,1.36,"Mass (GeV)",65.,100.);
  // [65,100] => [66,97] to avoid edge effects in JER
  TH1D *h1 = tdrHist("h1","N frac",0,0.06,"Mass (GeV)",60.,110.);
  //lumi_13TeV = "Run2";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,4,8,kSquare);

  TH1D *h2 = tdrHist("h2","Data/fit -1 (%)",-6,5,"Mass (GeV)",60.,110.);
  TCanvas *c2 = tdrCanvas("c2",h2,4,8,kSquare);

  TCanvas *c = tdrDiCanvas("c",h1,h2,4,11);

  c->cd(1);

  tdrDraw(h3all_mc,"Pz",kFullSquare,kPink-9);
  tdrDraw(h3all_data,"Pz",kFullCircle,kBlue-9);

  //TF1 *fit = new TF1("fit",fit_func,66,97,3);
  TF1 *fit = new TF1("fit",fit_func,72,96,8); // use co-opted range 97-100 GeV
  // Set initial parameter values
  fit->SetParameters(0.9, 83, 10, 0.1, 90, 10, 0, 0);
  fit->FixParameter(6, 0);
  fit->FixParameter(7, 0);
  h3all_data->Fit(fit,"RNI");
  fit->SetRange(60,110);
  fit->SetLineColor(kBlue-9);
  fit->SetLineWidth(3);
  fit->Draw("SAME");

  TF1 *fit2 = new TF1("fit2",fit_func,72,96,8); // use co-opted range 97-100 GeV
  // Set initial parameter values
  fit2->SetParameters(0.9, 83, 10, 0.1, 90, 10, 0, 3);
  fit2->FixParameter(6, 0);
  fit2->FixParameter(7, 0);
  h3all_mc->Fit(fit2,"RNI");
  fit2->SetRange(60,110);
  fit2->SetLineColor(kPink-9);
  fit2->SetLineWidth(2);
  fit2->Draw("SAME");


double amplitude = fit->GetParameter(0);
double mean = fit->GetParameter(1);
double sigma = fit->GetParameter(2);
double err_amplitude = fit->GetParError(0);
double err_mean = fit->GetParError(1);
double err_sigma = fit->GetParError(2);

double amplitude2 = fit2->GetParameter(0);
double mean2 = fit2->GetParameter(1);
double sigma2 = fit2->GetParameter(2);
double err_amplitude2 = fit2->GetParError(0);
double err_mean2 = fit2->GetParError(1);
double err_sigma2 = fit2->GetParError(2);

TLatex latex;
// Add fit parameters to the canvas
latex.DrawLatex(0.32, 0.80, Form("Amplitude: %.3f #pm %.3f", amplitude, err_amplitude));
latex.DrawLatex(0.32, 0.75, Form("Mean: %.3f #pm %.3f", mean, err_mean));
latex.DrawLatex(0.32, 0.70, Form("Sigma: %.3f #pm %.3f", sigma, err_sigma));

  // Add the legend after all drawings
  TLegend *legend = tdrLeg(0.7, 0.75-0.05*4, 0.95, 0.75); 

  legend->AddEntry(h3all_mc, "MC", "P");
  legend->AddEntry(h3all_data, "DATA", "P");
  legend->AddEntry(fit, "Fit DATA", "L");
  legend->AddEntry(fit2, "Fit MC", "L");
  //legend->AddEntry(fit_mod, "Fit Mod", "L");

  legend->Draw("SAME");

  c->cd(2);

  // Create a histogram for the fit function values
  TH1D *h_fit_data = (TH1D*)h3all_data->Clone("h_fit_data");
  h_fit_data->Reset(); // Clear any existing content

  TH1D *h_fit_mc = (TH1D*)h3all_mc->Clone("h_fit_mc");
  h_fit_mc->Reset(); // Clear any existing content

  int nBins = h_fit_data->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
      double x = h_fit_data->GetBinCenter(i);
      double x_min = h_fit_data->GetBinLowEdge(i);
      double x_max = h_fit_data->GetBinLowEdge(i+1);
      //double fit_value = fit->Eval(x);
      double fit_value = fit->Integral(x_min,x_max)/(x_max-x_min);
      h_fit_data->SetBinContent(i, fit_value);
      h_fit_data->SetBinError(i, 0); // Assuming negligible error on fit values for now
  }

  int nBins2 = h_fit_mc->GetNbinsX();
  for (int i = 1; i <= nBins2; ++i) {
      double x = h_fit_mc->GetBinCenter(i);
      double x_min = h_fit_data->GetBinLowEdge(i);
      double x_max = h_fit_data->GetBinLowEdge(i+1);
      //double fit_value = fit2->Eval(x);
      double fit_value = fit2->Integral(x_min,x_max)/(x_max-x_min);
      h_fit_mc->SetBinContent(i, fit_value);
      h_fit_mc->SetBinError(i, 0); // Assuming negligible error on fit values for now
  }

  // Clone the data histogram to create the ratio histogram
  TH1D *h_ratio_data = (TH1D*)h3all_data->Clone("h_ratio_data");

  // Clone the data histogram to create the ratio histogram
  TH1D *h_ratio_mc = (TH1D*)h3all_mc->Clone("h_ratio_mc");

  // Divide data histogram by fit histogram
  h_ratio_data->Divide(h_fit_data);

  // Divide data histogram by fit histogram
  h_ratio_mc->Divide(h_fit_mc);

  // Transform the ratio histogram to (data/fit - 1) * 100
  for (int i = 1; i <= nBins; ++i) {
      double ratio = h_ratio_data->GetBinContent(i);
      double ratio_error = h_ratio_data->GetBinError(i);

      // Compute percentage difference
      double percent_diff = (ratio - 1.0) * 100.0;
      double percent_error = ratio_error * 100.0;

      h_ratio_data->SetBinContent(i, percent_diff);
      h_ratio_data->SetBinError(i, percent_error);
  }

  // Transform the ratio histogram to (mc/fit - 1) * 100
  for (int i = 1; i <= nBins; ++i) {
      double ratio = h_ratio_mc->GetBinContent(i);
      double ratio_error = h_ratio_mc->GetBinError(i);

      // Compute percentage difference
      double percent_diff = (ratio - 1.0) * 100.0;
      double percent_error = ratio_error * 100.0;

      h_ratio_mc->SetBinContent(i, percent_diff);
      h_ratio_mc->SetBinError(i, percent_error);
  }



    // Add a horizontal line at zero
  TLine *line = new TLine(60., 0.0, 110., 0.0);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2); // Dashed line
  line->Draw("SAME");

  h_ratio_data->Draw("same");
  h_ratio_mc->Draw("same");

  double chi2_data = fit->GetChisquare();
  double chi2_mc = fit2->GetChisquare();

  int ndf_data = fit->GetNDF();
  int ndf_mc = fit2->GetNDF();

  latex.SetNDC();
  latex.SetTextSize(0.045*1.5);
  latex.DrawLatex(0.4, 0.45, Form("#chi^{2}/NDF_{DATA} = %.2f / %d", chi2_data, ndf_data));
  latex.DrawLatex(0.4, 0.38, Form("#chi^{2}/NDF_{MC} = %.2f / %d", chi2_mc, ndf_mc));
  


   // Save the canvas
  c->SaveAs("pdf/drawWqqFit_after.pdf");   // Save as PDF
  
} // drawWqqFromEmilia
