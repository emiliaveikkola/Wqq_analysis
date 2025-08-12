#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

// Define the custom fit function:
// f(x) = 1 + p0*log(x) + p1*(x-1)
double customLogLin(double *x, double *p) {
    // x[0]: independent variable.
    // p[0]: coefficient of log(x)
    // p[1]: coefficient of (x-1)
    return 1.0 + p[0]*TMath::Log(x[0]) + p[1]*(x[0] - 1);
}


void JERVariation() {
    setTDRStyle();
    extraText = "Private";
    // Example data: replace these arrays with your actual data points
    const int nPoints = 5;
    // Variation values (x-axis), ranging from 0 to 5
    double variation[nPoints] = {1, 2, 3, 4, 5};
    double variation2[nPoints] = {1.1, 2.1, 3.1, 4.1, 5.1};
    double variation3[nPoints] = {1.3, 2.3, 3.3, 4.3, 5.3};
    double variation4[nPoints] = {1.2, 2.2, 3.2, 4.2, 5.2};
    // JER values (y-axis), expected to lie in the range 0 to 2
    double JER[nPoints] = {1.04292, 1.01584, 1.0, 0.989669, 0.982111}; //only JER
    double JER2[nPoints] = {1.00156, 0.999056, 1., 0.999864, 0.999536}; //with FSR
    double JER3[nPoints] = {1.04045, 1.01311, 1., 0.991527, 0.985237}; //with ISR
    double JER4[nPoints] = {1.00349, 1.00225, 1., 0.997692, 0.996032}; //with ISR and FSR

    double JERerr[nPoints] = {0.0149816, 0.00479335, 0., 0.00277459, 0.00473271};

    double JERerr2[nPoints] = {0.00780942, 0.00233569, 0., 0.00134332, 0.00228045};

    double JERerr3[nPoints] = {0.0133279, 0.00311672, 0., 0.00173287, 0.00298778};

    double JERerr4[nPoints] = {0.00827534, 0.00282685, 0.0, 0.00165145, 0.00278666};

    // Assume negligible errors in variation (x) values.
    double xErr[nPoints] = {0, 0, 0, 0, 0};

    // Create a canvas
    TH1D *h2 = tdrHist("h2", "JER", 0.89,1.14, "FSR PSWeight", 0,6);
    TCanvas *c2 = tdrCanvas("c2", h2, 4, 11, kSquare);
    gPad->SetLeftMargin(0.15);
    h2->GetYaxis()->SetTitleOffset(1.6);
    
    // Now set the custom labels for each bin:
    h2->GetXaxis()->SetBinLabel(1, "");
    h2->GetXaxis()->SetBinLabel(2, "");
    h2->GetXaxis()->SetBinLabel(3, "");
    h2->GetXaxis()->SetBinLabel(4, "");
    h2->GetXaxis()->SetBinLabel(5, "");

    // Now manually draw the custom labels at x = 1, 2, 3, 4, 5.
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42); 
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.245, 0.085, "0.25");
    tex.DrawLatex(0.392, 0.085, "0.5");
    tex.DrawLatex(0.54, 0.085, "1");
    tex.DrawLatex(0.675, 0.085, "2");
    tex.DrawLatex(0.807, 0.085, "4");


    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    
    h2->GetYaxis()->SetLabelSize(0.045);

    // Create a TGraph with the data points
    TGraph *gr = new TGraphErrors(nPoints, variation, JER, xErr, JERerr);
    TGraph *gr2 = new TGraphErrors(nPoints, variation2, JER2, xErr, JERerr2);
    TGraph *gr3 = new TGraphErrors(nPoints, variation3, JER3, xErr, JERerr3);
    TGraph *gr4 = new TGraphErrors(nPoints, variation4, JER4, xErr, JERerr4);

    // Set marker style and size
    gr->SetMarkerStyle(kFullSquare);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kAzure+7);

    gr2->SetMarkerStyle(kFullCircle);
    gr2->SetMarkerSize(0.8);
    gr2->SetMarkerColor(kPink-9);
    
    gr4->SetMarkerStyle(kFullDiamond);
    gr4->SetMarkerSize(1);
    gr4->SetMarkerColor(kSpring-5);

    gr3->SetMarkerStyle(kFullStar);
    gr3->SetMarkerSize(1);
    gr3->SetMarkerColor(kOrange-3);
    // Draw the graph with axis and points
    gr->Draw("Same P");
    gr2->Draw("Same P");
    gr3->Draw("Same P");
    gr4->Draw("Same P");

        // Fit and draw a shifted quadratic for each graph, anchoring f(3)=1
    auto makeFit = [&](TGraph* g, const char* name) {
        TF1* f = new TF1(name, "[0] + [1]*(x-3) + [2]*(x-3)*(x-3)", 0, 6);
        f->FixParameter(0, 1.0);       // enforce f(3)=1
        f->SetParameter(1, 0.0);       // initial slope
        f->SetParameter(2, 0.0);       // initial curvature
        g->Fit(f, "RN");               // fit without drawing default
        f->SetLineColor(g->GetMarkerColor());
        f->Draw("same");
    };
    makeFit(gr,  "fit1");
    makeFit(gr2, "fit2");
    makeFit(gr3, "fit3");
    makeFit(gr4, "fit4");

    TLegend *legend = tdrLeg(0.6, 0.85-0.04*4, 0.8, 0.85); 
    legend->SetTextSize(0.04);
    legend->AddEntry(gr, "only JER", "P");
   legend->AddEntry(gr2, "JER + FSR", "P");
   legend->AddEntry(gr3, "JER + ISR", "P");
   legend->AddEntry(gr4, "JER + FSR + ISR", "P");
    // Update and save the canvas to a file (e.g. PNG)
    c2->Update();
    c2->SaveAs("pdf/JER_Variation.pdf");
}