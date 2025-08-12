#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include <TLine.h>
#include "../minitools/tdrstyle_mod22.C"

// Define the custom fit function:
// f(x) = 1 + p0*log(x) + p1*(x-1)
double customLogLin(double *x, double *p) {
    // x[0]: independent variable.
    // p[0]: coefficient of log(x)
    // p[1]: coefficient of (x-1)
    return 1.0 + p[0]*TMath::Log(x[0]) + p[1]*(x[0] - 1);
}


void FSR_JER_ISRVariation() {
    setTDRStyle();
    extraText = "Private";
    // Example data: replace these arrays with your actual data points
const int nPoints = 5;
// FSR variations
const char* variations[nPoints] = {"fsr:murfac=0.25", "fsr:murfac=0.5", "fsr:murfac=2.0", "fsr:murfac=4.0"};
// JER values (y-axis)
double JER[nPoints] = {1.00349, 1.00225, 1.0, 0.997692, 0.996032};
// FSR values (y-axis)
double FSR[nPoints] = {0.994931, 0.998472, 1.0, 1.00087, 1.00152};
// FSR errors
double FSRerr[nPoints] = {0.000684851, 0.00030114, 0.0, 0.000180627, 0.000302437};
// JER errors
double JERerr[nPoints] = {0.00827534, 0.00282685, 0.0, 0.00165145, 0.00278666};

    // Create a canvas
    TH1D *h2 = tdrHist("h2", "FSR", 0.992,1.004, "JER", 0.98,1.02);
    TCanvas *c2 = tdrCanvas("c2", h2, 4, 11, kSquare);
    gPad->SetLeftMargin(0.15);
    h2->GetYaxis()->SetTitleOffset(1.6);
/*
    TF1* f1 = new TF1("f1","[0]+[1]*(x-1)", 0.98,1.02);
    f1->SetParameters(1,-0.75);
    f1->SetLineColor(kBlue);
    f1->DrawClone("same");
    f1->SetParameters(1,-0.5);
    f1->SetLineColor(kGreen-2);
    f1->DrawClone("same");
    f1->SetParameters(1,-1);
    f1->SetLineColor(kRed);
    f1->DrawClone("same");

   // Draw vertical reference line at FSR = 1
    double yMin = gPad->GetUymin();
    double yMax = gPad->GetUymax();
    TLine *vline = new TLine(1.0, yMin, 1.0, yMax);
    vline->SetLineColor(kRed);
    vline->SetLineStyle(2);
    vline->Draw();
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
    gr->SetMarkerColor(kPink-9);
    // Draw the graph with axis and points
    gr->Draw("Same P");


    // Exponential fit: FSR = exp(a*(JER - 1))
    TF1* fexp = new TF1("fexp", "exp([0]*(x-1))", 0.985, 1.015);
    fexp->SetParameter(0, 0.0);        // initial guess for exponent
    gr->Fit(fexp, "RN");              // fit without default drawing
    fexp->SetLineColor(kBlue);
    fexp->Draw("same");

    TLegend *legend = tdrLeg(0.6, 0.8-0.04*4, 0.8, 0.8); 
    legend->SetTextSize(0.04);
    legend->AddEntry(gr, "FSR + JER + ISR", "P");
    // Update and save the canvas to a file (e.g. PNG)
    c2->Update();
    c2->SaveAs("pdf/FSR_JER_ISR_Variation.pdf");
}