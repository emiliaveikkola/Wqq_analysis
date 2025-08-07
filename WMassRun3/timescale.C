#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMarker.h>
#include "TAxis.h"
#include <iostream>
#include "tdrstyle_mod22.C"  // Assuming you have this for styling

void timescale() {
    // Apply TDR style
    setTDRStyle();
    extraText = "Private";

    // Step 1: Load the TH2D histogram from the specified ROOT file
    TFile* file = TFile::Open("output6.root");  // Use your specified ROOT file
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file 'output3.root'." << std::endl;
        return;
    }

    TProfile* prof_top_improved = (TProfile*)file->Get("prof_top_improved");
    TProfile* prof_W = (TProfile*)file->Get("prof_W_inWindow");
    TProfile* h_Topmass_improved = (TProfile*)file->Get("h_Topmass_improved");

    // First, project TProfile to a TH1D:
    TH1D* h_proj = prof_top_improved->ProjectionX("h_proj");
    TH1D* h_proj2 = prof_W->ProjectionX("h_proj2");

    // Now, loop over the bins and filter out those with relative error > 1%
    int nbins = h_proj->GetNbinsX();
    for (int i = 1; i <= nbins; i++){
        double val = h_proj->GetBinContent(i);
        double err = h_proj->GetBinError(i);
        // If the relative error is greater than 1% (and the bin isn't empty), clear the bin:
        if (val != 0 && (err/val) > 0.01) {
            h_proj->SetBinContent(i, 0);
            h_proj->SetBinError(i, 0);
        }
    }

    // Now, loop over the bins and filter out those with relative error > 1%
    int nbins2 = h_proj2->GetNbinsX();
    for (int i = 1; i <= nbins2; i++){
        double val = h_proj2->GetBinContent(i);
        double err = h_proj2->GetBinError(i);
        // If the relative error is greater than 1% (and the bin isn't empty), clear the bin:
        if (val != 0 && (err/val) > 0.01) {
            h_proj2->SetBinContent(i, 0);
            h_proj2->SetBinError(i, 0);
        }
        //if (val != 0 && (err/val) <= 0.02) cout << "Not Bin: " << i << " Value: " << val  << endl;

    }

    double overallMean = prof_top_improved->GetMean(2);
    cout << "Top: " << overallMean << endl;

    // Now normalize every bin by the overall mean:
    if (overallMean != 0) {
        for (int i = 1; i <= nbins; i++){
            double val = h_proj->GetBinContent(i);
            double err = h_proj->GetBinError(i);
            h_proj->SetBinContent(i, val / overallMean);
            h_proj->SetBinError(i, err / overallMean);
            if (val != 0){
                //cout << "Bin: " << i << " Value: " << val << " Content: "  << val / overallMean << endl;
            }
        }
    }

    double overallMean2 = prof_W->GetMean(2);
        cout << "W: " << overallMean2 << endl;


    // Now normalize every bin by the overall mean:
    if (overallMean2 != 0) {
        for (int i = 1; i <= nbins2; i++){
            double val = h_proj2->GetBinContent(i);
            double err = h_proj2->GetBinError(i);
            h_proj2->SetBinContent(i, val / overallMean2);
            h_proj2->SetBinError(i, err / overallMean2);
            if (val != 0){
                //cout << "Bin: " << i << " Value: " << val << " Content: "  << val / overallMean << endl;
            }
        }
    }

    


    TH1D *h = tdrHist("h", "Relative Mass", 0.93, 1.08, "Run", 380.2e3, 387.1e3);
    TCanvas *canvas = tdrCanvas("canvas", h, 8, 11, kSquare);


    tdrDraw(h_proj2,"Pz",kFullCircle,kSpring-5);
    h_proj2->SetMarkerSize(0.8);

    tdrDraw(h_proj,"Pz",kFullCircle,kBlue-9);
    h_proj->SetMarkerSize(0.8);
    
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetNdivisions(505);

      // Add the legend after all drawings
    TLegend *legend = tdrLeg(0.7, 0.87-0.05*2, 0.95, 0.87); 

    legend->AddEntry(h_proj, "Top", "P");
    legend->AddEntry(h_proj2, "W", "P");
    legend->Draw("SAME");

    gPad->RedrawAxis();
    gPad->Update();

    // Save the canvas to a PDF file as specified
    canvas->SaveAs("pdf/timescale.pdf");

    TH1D *h2 = tdrHist("h2", "N", 0, 15000, "Top Mass (GeV)", 70, 400);
    TCanvas *canvas2 = tdrCanvas("canvas2", h2, 8, 11, kSquare);

    tdrDraw(h_Topmass_improved,"Pz",kFullCircle,kViolet+1);
    h_Topmass_improved->SetMarkerSize(0.8);
    
    h2->GetYaxis()->SetLabelSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.04);
    h2->GetYaxis()->SetTitleSize(0.045);
    h2->GetXaxis()->SetTitleSize(0.045);

    gPad->RedrawAxis();
    gPad->Update();

    // Save the canvas to a PDF file as specified
    canvas2->SaveAs("pdf/top_mass.pdf");

    // Close the file after processing
    file->Close();
}