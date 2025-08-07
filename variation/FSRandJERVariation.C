#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include <TLine.h>
#include "tdrstyle_mod22.C"

// Define the custom fit function:
// f(x) = 1 + p0*log(x) + p1*(x-1)
double customLogLin(double *x, double *p) {
    // x[0]: independent variable.
    // p[0]: coefficient of log(x)
    // p[1]: coefficient of (x-1)
    return 1.0 + p[0]*TMath::Log(x[0]) + p[1]*(x[0] - 1);
}


void FSRandJERVariation() {
    setTDRStyle();
    extraText = "Private";
    // Example data: replace these arrays with your actual data points
const int nPoints = 5;
// FSR variations
const char* variations[nPoints] = {"fsr:murfac=0.25", "fsr:murfac=0.5", "fsr:murfac=2.0", "fsr:murfac=4.0"};
// JER values (y-axis)
double JER[nPoints] = {1.00156, 0.999056, 1., 0.999864, 0.999536};
// FSR values (y-axis)
double FSR[nPoints] = {0.994625, 0.99798, 1., 1.00121, 1.00206};
// FSR errors
double FSRerr[nPoints] = {0.000566105, 0.000175924, 0.00000000001, 0.000101863, 0.000172785};
// JER errors
double JERerr[nPoints] = {0.00780942, 0.00233569, 0.00000000001, 0.00134332, 0.00228045};
    // Create a canvas
    TH1D *h2 = tdrHist("h2", "FSR", 0.99,1.01, "JER", 0.99,1.01);
    TCanvas *c2 = tdrCanvas("c2", h2, 4, 11, kSquare);
    gPad->SetLeftMargin(0.15);
    h2->GetYaxis()->SetTitleOffset(1.6);
/*
    TF1* f1 = new TF1("f1","[0]+[1]*(x-1)", 0.98,1.02);
    f1->SetParameters(1,-5);
    f1->SetLineColor(kBlue);
    f1->DrawClone("same");
    f1->SetParameters(1,5);
    f1->SetLineColor(kGreen-2);
    f1->DrawClone("same");


    // Draw vertical reference line at FSR = 1
    double yMin = gPad->GetUymin();
    double yMax = gPad->GetUymax();
    TLine *vline = new TLine(1.0, yMin, 1.0, yMax);
    vline->SetLineColor(kRed);
    vline->SetLineStyle(2);
    vline->Draw();

    double target_fx = 0.996;

    // Solve: f(x) = 1 + a*(x - 1) = target_fx
    double x_value = 1 - (1.0 - target_fx) / 0.5;

    std::cout << "For slope 0.5, x for f(x) = " << target_fx << " is: " << x_value << std::endl;
    
*/

    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    
    h2->GetYaxis()->SetLabelSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.04);

    // Create a TGraph with the data points
    TGraph *gr = new TGraphErrors(nPoints, JER, FSR, JERerr, FSRerr);

    // Set marker style and size
    gr->SetMarkerStyle(kFullSquare);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kAzure+7);
    // Draw the graph with axis and points
    gr->Draw("Same P");

    TF1* frac = new TF1("frac",
        "(1 + [0]*(x-1)) / (1 + [1]*(x-1))",
        0.99, 1.01);
    // enforces f(1)=1 by construction
    frac->SetParameters(0.0, 0.0);
    gr->Fit(frac, "RN");
    frac->SetLineColor(kBlue);
    frac->SetRange(0.999056, 1.01);  // only draw where JER points exist
    frac->Draw("same");

    // --- Alternative fit ignoring the first data point ---
    // Perform an alternative fit ignoring the first data point
    // Clone the graph and remove the first point (index 0)
    TGraphErrors *gr_fit = (TGraphErrors*)gr->Clone("gr_fit");
    gr_fit->RemovePoint(0);

    // Example: power-law fit f(x) = x^[0], anchored at f(1)=1
    TF1 *fpow = new TF1("fpow", "pow(x, [0])", 0.999056, 1.00156);
    fpow->SetParameter(0, 0.0);         // initial exponent guess
    gr_fit->Fit(fpow, "RN");            // fit without default drawing
    fpow->SetLineColor(kRed);
    fpow->Draw("same");

    TLegend *legend = tdrLeg(0.2, 0.8-0.04*4, 0.5, 0.8); 
    legend->SetTextSize(0.04);
    legend->AddEntry(gr, "FSR + JER", "P");

    // Update and save the canvas to a file (e.g. PNG)
    c2->Update();
    c2->SaveAs("pdf/FSR_JER_Variation.pdf");
}