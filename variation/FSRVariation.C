#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include "tdrstyle_mod22.C"

// Define the custom fit function:
// f(x) = 1 + p0*log(x) + p1*(x-1)
double customLogLin(double *x, double *p) {
    // x[0]: independent variable.
    // p[0]: coefficient of log(x)
    // p[1]: coefficient of (x-1)
    return 1.0 + p[0]*TMath::Log(x[0]) + p[1]*(x[0] - 1);
}


void FSRVariation() {
    setTDRStyle();
    extraText = "Private";
    // Example data: replace these arrays with your actual data points
    const int nPoints = 5;
    // Variation values (x-axis), ranging from 0 to 5
    double variation[nPoints] = {1, 2, 3, 4, 5};
    double variation2[nPoints] = {1.1, 2.1, 3.1, 4.1, 5.1};
    double variation3[nPoints] = {1.2, 2.2, 3.2, 4.2, 5.2};
    double variation4[nPoints] = {1.3, 2.3, 3.3, 4.3, 5.3};
    
    // FSR values (y-axis), expected to lie in the range 0 to 2
    double FSR[nPoints] = {0.994562, 0.998024, 1., 1.00121, 1.00208}; //only FSR
    double FSR2[nPoints] = {0.994625, 0.99798, 1., 1.00121, 1.00206}; //FSR and JER
    double FSR3[nPoints] = {0.994757, 0.998291, 1., 1.00107, 1.00185}; //FSR and ISR
    double FSR4[nPoints] = {0.994931, 0.998472, 1., 1.00087, 1.00152}; //with ISR and JER

    double FSRerr[nPoints] = {0.00054985, 0.000157339, 0., 9.3504e-05, 0.000159752};

    double FSRerr2[nPoints] = {0.000566105, 0.000175924, 0., 0.000101863, 0.000172785};

    double FSRerr3[nPoints] = {0.000639225, 0.000227772, 0., 0.00013805, 0.000234091};

    double FSRerr4[nPoints] = {0.000684851, 0.00030114, 0., 0.000180627, 0.000302437};

    // Assume negligible errors in variation (x) values.
    double xErr[nPoints] = {0, 0, 0, 0, 0};

    // Create a canvas
    TH1D *h2 = tdrHist("h2", "FSR", 0.990,1.006, "FSR PSWeight", 0,6);
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
    TGraph *gr = new TGraphErrors(nPoints, variation, FSR, xErr, FSRerr);
    TGraph *gr2 = new TGraphErrors(nPoints, variation2, FSR2, xErr, FSRerr2);
    TGraph *gr3 = new TGraphErrors(nPoints, variation3, FSR3, xErr, FSRerr3);
    TGraph *gr4 = new TGraphErrors(nPoints, variation4, FSR4, xErr, FSRerr4);

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

    /*TF1* f1 = new TF1("f1","[0]+[1]*(x-3)", 0,6);
    f1->SetParameters(1,0.0025);
    f1->SetLineColor(kBlue);
    f1->DrawClone("same");
    f1->SetParameters(1,0.003);
    f1->SetLineColor(kGreen-2);
    f1->DrawClone("same");
    f1->SetParameters(1,0.002);
    f1->SetLineColor(kRed);
    f1->DrawClone("same");
*/
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

    TLegend *legend = tdrLeg(0.6, 0.5-0.04*4, 0.8, 0.5); 
    legend->SetTextSize(0.04);
    legend->AddEntry(gr, "only FSR", "P");
    legend->AddEntry(gr2, "FSR + JER", "P");
    legend->AddEntry(gr3, "FSR + ISR", "P");
    legend->AddEntry(gr4, "FSR + JER + ISR", "P");
    // Update and save the canvas to a file (e.g. PNG)
    c2->Update();
    c2->SaveAs("pdf/FSR_Variation.pdf");
}