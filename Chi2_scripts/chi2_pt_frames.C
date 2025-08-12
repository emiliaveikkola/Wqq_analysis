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

void chi2_pt_frames() {
    setTDRStyle();
    extraText = "Private";

    // Open the output files
    TFile *file = new TFile("output_MCRun2_Wqq_test4.root", "READ");
    TFile *file2 = new TFile("output_DATARun2_Wqq_test4.root", "READ");

    std::vector<double> scaleFactors;
    for (double scale = 0.98; scale <= 1.02; scale += 0.002) {
        scaleFactors.push_back(scale);
    }
    std::vector<double> udScaleFactors;
    for (double scale_ud = 0.98; scale_ud <= 1.02; scale_ud += 0.002) {
        udScaleFactors.push_back(scale_ud);
    }
    // Create a TH2D histogram to store chi-squared values
    int nBins = scaleFactors.size();
    double dx = 0.002/2;
    TH2D* chi2_hist2 = new TH2D("chi2_hist", ";R_{c};R_{s};#chi^{2}",
                               nBins, 0.98-dx, 1.02+dx, nBins, 0.98-dx, 1.02+dx);

    TH2D* chi2_NDF_hist = new TH2D("chi2_NDF_hist", ";R_{c};R_{s};#chi^{2}/NDF",
    nBins, 0.98-dx, 1.02+dx, nBins, 0.98-dx, 1.02+dx);

    TH3D* h3PtFlavorPairs_DATA = (TH3D*)file2->Get("h3PtFlavorPairs_DATA");

    TH1D* h3all_data = h3PtFlavorPairs_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

    TH1D* h3tagcsall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
    TH1D* h3tagudall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
    TH1D* h3tagxall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);

    h3tagcsall_data->Scale(1./h3tagcsall_data->Integral());
    h3tagudall_data->Scale(1./h3tagudall_data->Integral());
    h3tagxall_data->Scale(1./h3tagxall_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

    for (double scale_ud : udScaleFactors) {
        chi2_NDF_hist->Reset("ICES");
        chi2_hist2->Reset("ICES");
        for (size_t i = 0; i < scaleFactors.size(); ++i) {
            for (size_t j = 0; j < scaleFactors.size(); ++j) {
                double scale1 = scaleFactors[i];
                double scale2 = scaleFactors[j];
                
                // Format the histogram name with 4 decimal places
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3)  << scale1 << "_" << scale2 << "_" << scale_ud; // c, s & ud
                std::string histName = "h3PtFlavorPairs_MC_" + oss.str();

                // Retrieve the histogram
                TH3D* h3PtFlavorPairs_MC = (TH3D*)file->Get(histName.c_str());
                if (!h3PtFlavorPairs_MC) {
                    std::cerr << "Histogram " << histName << " not found!" << std::endl;
                    continue;
                }

                // Project Z for different indices and draw
                TH1D* h3all = h3PtFlavorPairs_MC->ProjectionZ("h3all", 1, 3, 1, 3);

                TH1D* h3gencsall = h3PtFlavorPairs_MC->ProjectionZ("h3gencsall", 1, 1, 1, 3);
                TH1D* h3genudall = h3PtFlavorPairs_MC->ProjectionZ("h3genudall", 2, 2, 1, 3);
                TH1D* h3genxall = h3PtFlavorPairs_MC->ProjectionZ("h3genxall", 3, 3, 1, 3);

                TH1D* h3tagcsall = h3PtFlavorPairs_MC->ProjectionZ("h3tagcsall", 1, 3, 1, 1);
                TH1D* h3tagudall = h3PtFlavorPairs_MC->ProjectionZ("h3tagudall", 1, 3, 2, 2);
                TH1D* h3tagxall = h3PtFlavorPairs_MC->ProjectionZ("h3tagxall", 1, 3, 3, 3);

                TH1D* h3gencstagcs = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagcs", 1, 1, 1, 1);
                TH1D* h3gencstagud = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagud", 1, 1, 2, 2);
                TH1D* h3gencstagx = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagx", 1, 1, 3, 3);

                TH1D* h3genudtagcs = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagcs", 2, 2, 1, 1);
                TH1D* h3genudtagud = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagud", 2, 2, 2, 2);
                TH1D* h3genudtagx = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagx", 2, 2, 3, 3);

                TH1D* h3genxtagcs = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagcs", 3, 3, 1, 1);
                TH1D* h3genxtagud = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagud", 3, 3, 2, 2);
                TH1D* h3genxtagx = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagx", 3, 3, 3, 3);

                // Scaling histograms
                double integral = h3all->Integral();
                h3gencstagcs->Scale(1.0 / h3gencstagcs->Integral()); //integral);
                h3gencstagud->Scale(1.0 / h3gencstagud->Integral()); //integral);
                h3gencstagx->Scale(1.0 / h3gencstagx->Integral()); //integral);

                h3gencsall->Scale(1.0 / h3gencsall->Integral()); //integral);
                h3genudall->Scale(1.0 / h3genudall->Integral()); //integral);
                h3genxall->Scale(1.0 / h3genxall->Integral()); //integral);
                
                h3genudtagcs->Scale(1.0 / h3genudtagcs->Integral()); //integral);
                h3genudtagud->Scale(1.0 / h3genudtagud->Integral()); //integral);
                h3genudtagx->Scale(1.0 / h3genudtagx->Integral()); //integral);

                h3genxtagcs->Scale(1.0 / h3genxtagcs->Integral()); //integral);
                h3genxtagud->Scale(1.0 / h3genxtagud->Integral()); //integral);
                h3genxtagx->Scale(1.0 / h3genxtagx->Integral()); //integral);

                h3tagcsall->Scale(1.0 / h3tagcsall->Integral()); //integral);
                h3tagudall->Scale(1.0 / h3tagudall->Integral()); //integral);
                h3tagxall->Scale(1.0 / h3tagxall->Integral()); //integral);

                h3all->Scale(1.0 / integral);

                /////////// CS //////////////

                // Calculate chi-squared values
                TH1D *h_tagcsvsdata = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata");
                h_tagcsvsdata->Divide(h3tagcsall_data);

                // Fit the ratio histogram with a constant function
                TF1* f7m_cs = new TF1("f7m_cs", "[0]", -0.7, 0.7);
                f7m_cs->FixParameter(0, 1);
                h_tagcsvsdata->Fit(f7m_cs, "QRN");
                double chi27m_cs = f7m_cs->GetChisquare();
                double NDF_cs = f7m_cs->GetNDF();

                chi2_hist2->Fill(scale1, scale2, chi27m_cs);

                chi2_NDF_hist->Fill(scale1, scale2, (chi27m_cs)/((NDF_cs)));

                // Fill the TH2D histogram
                cout << "scale_ud = " << scale_ud << ", scale1 = " << scale1 << ", scale2 = " << scale2 
                << ", chi27m_cs = " << chi27m_cs  << endl << flush;

                delete h3all;
                delete h3gencsall;
                delete h3genudall;
                delete h3genxall;
                delete h3tagcsall;
                delete h3tagudall;
                delete h3tagxall;
                delete h3gencstagcs;
                delete h3gencstagud;
                delete h3gencstagx;
                delete h3genudtagcs;
                delete h3genudtagud;
                delete h3genudtagx;
                delete h3genxtagcs;
                delete h3genxtagud;
                delete h3genxtagx;
                delete f7m_cs;
                
            }
        }
        // Before creating new hist and canvas
        if (gDirectory->FindObject("h")) delete gDirectory->FindObject("h");
        if (gROOT->FindObject("c2")) delete gROOT->FindObject("c2");
        
        TH1D *h = tdrHist("h", "Strange jet response (R_{s}) ", 0.979, 1.021, "Charm jet response (R_{c})",0.979, 1.021);
        TCanvas *c2 = tdrCanvas("c2", h, 4, 0, kRectangular);
        gPad->SetRightMargin(0.17);
        gPad->SetLeftMargin(0.13);
        
        chi2_NDF_hist->Draw("COLZ same");
        chi2_NDF_hist->GetZaxis()->SetRangeUser(3,9);
        chi2_NDF_hist->GetZaxis()->SetTitle("#frac{#chi^{2}}{NDF}");

        h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetLabelSize(0.045);
        chi2_NDF_hist->GetZaxis()->SetLabelSize(0.045);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        chi2_NDF_hist->GetZaxis()->SetTitleSize(0.045);

        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.4);
        chi2_NDF_hist->GetZaxis()->SetTitleOffset(1.1);


        gPad->RedrawAxis();
        gPad->Update();

        TLine *l = new TLine();
        l->DrawLine(0.98-0.001,1.02+0.001,1.02+0.001,0.98-0.001);
        l->SetLineStyle(kDotted);
        l->DrawLine(1,0.98-0.001,1,1.02+0.001);
        l->DrawLine(0.98-0.001,1.0,1.02+0.001,1.0);

        TGraph *g2 = new TGraph();

/*
        double chi2_min(99999), x_min(-1), y_min(-1);
        for (int i = 1; i != chi2_NDF_hist->GetNbinsX()+1; ++i){
            int j = chi2_NDF_hist->GetNbinsY()+1-i;
            double chi2 = chi2_NDF_hist->GetBinContent(i,j);
            double x = chi2_NDF_hist->GetXaxis()->GetBinCenter(i);
            g2->SetPoint(i-1,x,chi2);
            if (chi2 < chi2_min){
                chi2_min = chi2;
                x_min = chi2_NDF_hist->GetXaxis()->GetBinCenter(i);
                y_min = chi2_NDF_hist->GetYaxis()->GetBinCenter(j);
            }
        }
*/
        // New version: fit off-diagonal chi2 points with a parabola
        int nBinsX = chi2_NDF_hist->GetNbinsX();
        for (int i = 1; i <= nBinsX; ++i){
            int j = nBinsX + 1 - i; // i + j = N + 1
            double chi2 = chi2_NDF_hist->GetBinContent(i, j);
            double x = chi2_NDF_hist->GetXaxis()->GetBinCenter(i);
            g2->SetPoint(i-1, x, chi2);
        }

        g2->Fit("pol2", "Q");
        TF1 *fitFunc = g2->GetFunction("pol2");
        double a = fitFunc->GetParameter(2);
        double b = fitFunc->GetParameter(1);
        double x_min = -b / (2. * a);
        double chi2_min = fitFunc->Eval(x_min);

        // Get the lower x bin center that brackets x_min
        int i_bin = chi2_NDF_hist->GetXaxis()->FindBin(x_min);
        double x_lower = chi2_NDF_hist->GetXaxis()->GetBinCenter(i_bin);
        double x_upper = chi2_NDF_hist->GetXaxis()->GetBinCenter(i_bin+1);

        // Compute fraction of how far x_min is between these two bin centers
        double frac = (x_min - x_lower) / (x_upper - x_lower);

        // For the off-diagonal, if the corresponding y bin for x_lower is given by 
        // j_lower = (nBinsX + 1 - i_bin), then the next bin is j_upper = (nBinsX + 1 - (i_bin+1))
        int j_lower = chi2_NDF_hist->GetNbinsY() + 1 - i_bin;
        int j_upper = chi2_NDF_hist->GetNbinsY() + 1 - (i_bin + 1);
        double y_lower = chi2_NDF_hist->GetYaxis()->GetBinCenter(j_lower);
        double y_upper = chi2_NDF_hist->GetYaxis()->GetBinCenter(j_upper);

        // Interpolate linearly to get y_min corresponding to x_min
        double y_min = y_lower + frac * (y_upper - y_lower);

        // Use the exact interpolated values for marking the minimum
        double x_star = x_min;
        double y_star = y_min;

        cout << "chi2_min: " << chi2_min << ", x_min: " << x_min << ", y_min: " << y_min << endl << flush;
        cout << "chi2_NDF_hist->GetMinimum(): " << chi2_NDF_hist->GetMinimum() << endl << flush;

        // Old method: directly find minimum along off-diagonal
        double chi2_min_old(99999), x_min_old(-1), y_min_old(-1);
        for (int i = 1; i <= chi2_NDF_hist->GetNbinsX(); ++i){
            int j = chi2_NDF_hist->GetNbinsY() + 1 - i;
            double chi2 = chi2_NDF_hist->GetBinContent(i, j);
            if (chi2 < chi2_min_old){
                chi2_min_old = chi2;
                x_min_old = chi2_NDF_hist->GetXaxis()->GetBinCenter(i);
                y_min_old = chi2_NDF_hist->GetYaxis()->GetBinCenter(j);
            }
        }
    
        // Draw old method minimum as a hollow star
        TGraph* g_old = new TGraph();
        g_old->SetPoint(0, x_min_old, y_min_old);
        g_old->SetMarkerStyle(kStar);
        g_old->SetMarkerSize(2);
        g_old->SetMarkerColor(kWhite);
        g_old->Draw("same p");

        // Create a TGraph to mark the minimum point with a star
        TGraph *g = new TGraph();
        g->SetPoint(0, x_star, y_star);  // Mark the minimum point
        g->SetMarkerStyle(kFullStar);
        g->SetMarkerSize(2.5);
        g->Draw("same p");  // Draw the graph with the star on the same plot

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
        latex.DrawLatex(0.5, 0.755, "#frac{#chi^{2}}{NDF}, NDF = 47");
        latex.DrawLatex(0.32, 0.86, Form("R_{ud} = %.3f", scale_ud));
        latex.DrawLatex(0.25, 0.81, Form("min #chi^{2}/NDF = %.2f", chi2_min));
        TLatex latexOld;
        latexOld.SetTextSize(0.04);
        latexOld.SetTextColor(kWhite);  // Or any other ROOT color code
        latexOld.DrawLatex(0.2, 0.76, Form("old min #chi^{2}/NDF = %.2f", chi2_min_old));

        gPad->Modified();
        gPad->Update();
        c2->SaveAs(Form("frames/pt/chi2_pt_frame_ud_%.3f.png", scale_ud));
        c2->SaveAs(Form("frames/pt/chi2_pt_frame_ud_%.3f.pdf", scale_ud));
        TFile* outFile = new TFile(Form("frames/pt/rootfiles/chi2_3D_pt_cs_ud_%.3f.root", scale_ud), "RECREATE");
        chi2_hist2->SetName(Form("chi2_hist2_saved_%.3f", scale_ud));  // distinguish from mass
        chi2_hist2->Write();
        TF1* savedFitFunc = (TF1*)fitFunc->Clone(Form("fit_offdiag_ud_%.3f", scale_ud));
        savedFitFunc->Write();
        TGraph* savedGraph = (TGraph*)g2->Clone(Form("offdiag_graph_ud_%.3f", scale_ud));
        savedGraph->Write();
        outFile->Close();
    }
        // Create animated GIF from generated PNGs
        gSystem->Exec("magick -delay 50 -loop 0 frames/pt/chi2_pt_frame_ud_*.png frames/pt/chi2_pt_scan.gif");
        std::cout << "[INFO] chi2_pt_scan.gif created in frames/pt/" << std::endl;
}