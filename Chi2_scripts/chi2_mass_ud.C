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

void chi2_mass_ud() {
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

    TProfile2D* chi2_prof = new TProfile2D("chi2_prof", ";R_{ud};R_{cs};#chi^{2}",
                                nBins, 0.98-dx, 1.02+dx, nBins, 0.98-dx, 1.02+dx);

    TH3D* h3MassFlavorPairs_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATA");

    TH1D* h3all_data = h3MassFlavorPairs_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

    TH1D* h3tagcsall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
    TH1D* h3tagudall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
    TH1D* h3tagxall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);

    h3tagcsall_data->Scale(1./h3tagcsall_data->Integral());
    h3tagudall_data->Scale(1./h3tagudall_data->Integral());
    h3tagxall_data->Scale(1./h3tagxall_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

    for (double scale_ud : udScaleFactors) {
        for (size_t i = 0; i < scaleFactors.size(); ++i) {
            for (size_t j = 0; j < scaleFactors.size(); ++j) {
                double scale1 = scaleFactors[i];
                double scale2 = scaleFactors[j];
                
                // Format the histogram name with 4 decimal places
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3)  << scale1 << "_" << scale2 << "_" << scale_ud; // c, s & ud
                std::string histName = "h3MassFlavorPairs_MC_" + oss.str();

                // Retrieve the histogram
                TH3D* h3MassFlavorPairs_MC = (TH3D*)file->Get(histName.c_str());
                if (!h3MassFlavorPairs_MC) {
                    std::cerr << "Histogram " << histName << " not found!" << std::endl;
                    continue;
                }

                // Project Z for different indices and draw
                TH1D* h3all = h3MassFlavorPairs_MC->ProjectionZ("h3all", 1, 3, 1, 3);

                TH1D* h3gencsall = h3MassFlavorPairs_MC->ProjectionZ("h3gencsall", 1, 1, 1, 3);
                TH1D* h3genudall = h3MassFlavorPairs_MC->ProjectionZ("h3genudall", 2, 2, 1, 3);
                TH1D* h3genxall = h3MassFlavorPairs_MC->ProjectionZ("h3genxall", 3, 3, 1, 3);

                TH1D* h3tagcsall = h3MassFlavorPairs_MC->ProjectionZ("h3tagcsall", 1, 3, 1, 1);
                TH1D* h3tagudall = h3MassFlavorPairs_MC->ProjectionZ("h3tagudall", 1, 3, 2, 2);
                TH1D* h3tagxall = h3MassFlavorPairs_MC->ProjectionZ("h3tagxall", 1, 3, 3, 3);

                TH1D* h3gencstagcs = h3MassFlavorPairs_MC->ProjectionZ("h3gencstagcs", 1, 1, 1, 1);
                TH1D* h3gencstagud = h3MassFlavorPairs_MC->ProjectionZ("h3gencstagud", 1, 1, 2, 2);
                TH1D* h3gencstagx = h3MassFlavorPairs_MC->ProjectionZ("h3gencstagx", 1, 1, 3, 3);

                TH1D* h3genudtagcs = h3MassFlavorPairs_MC->ProjectionZ("h3genudtagcs", 2, 2, 1, 1);
                TH1D* h3genudtagud = h3MassFlavorPairs_MC->ProjectionZ("h3genudtagud", 2, 2, 2, 2);
                TH1D* h3genudtagx = h3MassFlavorPairs_MC->ProjectionZ("h3genudtagx", 2, 2, 3, 3);

                TH1D* h3genxtagcs = h3MassFlavorPairs_MC->ProjectionZ("h3genxtagcs", 3, 3, 1, 1);
                TH1D* h3genxtagud = h3MassFlavorPairs_MC->ProjectionZ("h3genxtagud", 3, 3, 2, 2);
                TH1D* h3genxtagx = h3MassFlavorPairs_MC->ProjectionZ("h3genxtagx", 3, 3, 3, 3);

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
                TF1* f7m_cs = new TF1("f7m_cs", "[0]", 70, 100);
                f7m_cs->FixParameter(0, 1);
                h_tagcsvsdata->Fit(f7m_cs, "QRN");
                double chi27m_cs = f7m_cs->GetChisquare();
                double NDF_cs = f7m_cs->GetNDF();

                ///////////// UD /////////////

                // Calculate chi-squared values
                TH1D *h_tagudvsdata = (TH1D*)h3tagudall->Clone("h_tagudvsdata");
                h_tagudvsdata->Divide(h3tagudall_data);

                // Fit the ratio histogram with a constant function
                TF1* f7m_ud = new TF1("f7m_ud", "[0]", 70, 100);
                f7m_ud->FixParameter(0, 1);
                h_tagudvsdata->Fit(f7m_ud, "QRN");
                double chi27m_ud = f7m_ud->GetChisquare();
                double NDF_ud = f7m_ud->GetNDF();
/*
                ////////// X ////////////

                // Calculate chi-squared values
                TH1D *h_tagxvsdata = (TH1D*)h3tagxall->Clone("h_tagxvsdata");
                h_tagxvsdata->Divide(h3tagxall_data);

                // Fit the ratio histogram with a constant function
                TF1* f7m_x = new TF1("f7m_x", "[0]", 70, 100);
                f7m_x->FixParameter(0, 1);
                h_tagxvsdata->Fit(f7m_x, "QRN");
                double chi27m_x = f7m_x->GetChisquare();
                double NDF_x = f7m_x->GetNDF();

                ////////// ALL ////////////

                // Calculate chi-squared values
                TH1D *h_allvsdata = (TH1D*)h3all->Clone("h_allvsdata");
                h_allvsdata->Divide(h3all_data);

                // Fit the ratio histogram with a constant function
                TF1* f7m_all = new TF1("f7m_all", "[0]", 70, 100);
                f7m_all->FixParameter(0, 1);
                h_allvsdata->Fit(f7m_all, "QRN");
                double chi27m_all = f7m_all->GetChisquare();
                double NDF_all = f7m_all->GetNDF();


*/
                chi2_prof->Fill(scale_ud, sqrt(scale1 * scale2), (chi27m_cs + chi27m_ud)/(NDF_cs+ NDF_ud));

                // Fill the TH2D histogram
                cout << "scale_ud = " << scale_ud << ", scale1 = " << scale1 << ", scale2 = " << scale2 
                << ", chi27m_cs = " << chi27m_cs << endl << flush;

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
                delete f7m_ud;
                //delete f7m_x;
                //delete f7m_all;
                
            }
        }
    }



        TH1D *h2 = tdrHist("h2", "CS jet response (R_{cs}) ", 0.979, 1.021, "UD jet response (R_{ud})",0.979, 1.021);
        TCanvas *c3 = tdrCanvas("c3", h2, 4, 0, kRectangular);
        gPad->SetRightMargin(0.17);
        gPad->SetLeftMargin(0.13);
        
        chi2_prof->Draw("COLZ same");
        chi2_prof->GetZaxis()->SetTitle("#frac{#chi^{2}}{NDF}");

        h2->GetXaxis()->SetLabelSize(0.045);
        h2->GetYaxis()->SetLabelSize(0.045);
        chi2_prof->GetZaxis()->SetLabelSize(0.045);
        h2->GetXaxis()->SetTitleSize(0.045);
        h2->GetYaxis()->SetTitleSize(0.045);
        chi2_prof->GetZaxis()->SetTitleSize(0.045);

        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.4);
        chi2_prof->GetZaxis()->SetTitleOffset(1.1);


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


        // Old method: directly find minimum along diagonal
        double chi2_min_old2 = 99999, x_min_old2 = -1, y_min_old2 = -1;
        for (int i = 1; i <= chi2_prof->GetNbinsX(); ++i){
            double chi2 = chi2_prof->GetBinContent(i, i);
            if (chi2 < chi2_min_old2){
                chi2_min_old2 = chi2;
                x_min_old2 = chi2_prof->GetXaxis()->GetBinCenter(i);
                y_min_old2 = chi2_prof->GetYaxis()->GetBinCenter(i);
            }
        }
        // Old method: directly find minimum along off-diagonal
        double chi2_min_old(99999), x_min_old(-1), y_min_old(-1);
        for (int i = 1; i <= chi2_prof->GetNbinsX(); ++i){
            int j = chi2_prof->GetNbinsY() + 1 - i;
            double chi2 = chi2_prof->GetBinContent(i, j);
            if (chi2 < chi2_min_old){
                chi2_min_old = chi2;
                x_min_old = chi2_prof->GetXaxis()->GetBinCenter(i);
                y_min_old = chi2_prof->GetYaxis()->GetBinCenter(j);
            }
        }
    
        TGraph *g2 = new TGraph();
        int nBinsX = chi2_prof->GetNbinsX();
        for (int i = 1; i <= nBinsX; ++i){
            int j = nBinsX + 1 - i; // i + j = N + 1
            double chi2 = chi2_prof->GetBinContent(i, j);
            double x = chi2_prof->GetXaxis()->GetBinCenter(i);
            g2->SetPoint(i-1, x, chi2);
        }

        g2->Fit("pol2", "Q");
        TF1 *fitFunc = g2->GetFunction("pol2");
        double a = fitFunc->GetParameter(2);
        double b = fitFunc->GetParameter(1);
        double x_min = -b / (2. * a);
        double chi2_min = fitFunc->Eval(x_min);

        // Get the lower x bin center that brackets x_min
        int i_bin = chi2_prof->GetXaxis()->FindBin(x_min);
        double x_lower = chi2_prof->GetXaxis()->GetBinCenter(i_bin);
        double x_upper = chi2_prof->GetXaxis()->GetBinCenter(i_bin+1);

        // Compute fraction of how far x_min is between these two bin centers
        double frac = (x_min - x_lower) / (x_upper - x_lower);

        // For the off-diagonal, if the corresponding y bin for x_lower is given by 
        // j_lower = (nBinsX + 1 - i_bin), then the next bin is j_upper = (nBinsX + 1 - (i_bin+1))
        int j_lower = chi2_prof->GetNbinsY() + 1 - i_bin;
        int j_upper = chi2_prof->GetNbinsY() + 1 - (i_bin + 1);
        double y_lower = chi2_prof->GetYaxis()->GetBinCenter(j_lower);
        double y_upper = chi2_prof->GetYaxis()->GetBinCenter(j_upper);

        // Interpolate linearly to get y_min corresponding to x_min
        double y_min = y_lower + frac * (y_upper - y_lower);

        // Use the exact interpolated values for marking the minimum
        double x_star = x_min;
        double y_star = y_min;

        TGraph *g3 = new TGraph();
        
        // New version: fit diagonal chi2 points with a parabola
        int nBinsX2 = chi2_prof->GetNbinsX();
        for (int i = 1; i != chi2_prof->GetNbinsX()+1; ++i){
            double chi2 = chi2_prof->GetBinContent(i,i);
            double x = chi2_prof->GetXaxis()->GetBinCenter(i);
            g3->SetPoint(i-1, x, chi2);
        }
        
        g3->Fit("pol2", "Q");
        TF1 *fitFunc2 = g3->GetFunction("pol2");
        double a2 = fitFunc2->GetParameter(2);
        double b2 = fitFunc2->GetParameter(1);
        double x_min2 = -b2 / (2. * a2);
        double chi2_min2 = fitFunc2->Eval(x_min2);
        double x_star2 = x_min2;  // Use the exact value from the fit
        double y_star2 = x_min2;  // For a diagonal, x = y
        cout << "chi2_min2: " <<chi2_min2 << ", x_min2: " << x_min2 << ", y_star2: " << x_min2 << endl << flush;
        cout << "chi2_prof->GetMinimum(): " << chi2_prof->GetMinimum() << endl << flush;

        // Draw old method minimum as a hollow star
        TGraph* g_old = new TGraph();
        g_old->SetPoint(0, x_min_old, y_min_old);
        g_old->SetMarkerStyle(kStar);
        g_old->SetMarkerSize(2);
        g_old->SetMarkerColor(kWhite);
        g_old->Draw("same p");
        
        TGraph* g_old2 = new TGraph();
        g_old2->SetPoint(0, x_min_old2, y_min_old2);
        g_old2->SetMarkerStyle(kStar);
        g_old2->SetMarkerSize(2);
        g_old2->SetMarkerColor(kWhite);
        g_old2->Draw("same p");

        // Create a TGraph to mark the minimum point with a star
        TGraph *g4 = new TGraph();
        g4->SetPoint(0, x_star2, y_star2);  // Mark the minimum point
        g4->SetMarkerStyle(kFullStar);
        g4->SetMarkerSize(2.5);
        g4->Draw("same p");  // Draw the graph with the star on the same plot

        TGraph *g5 = new TGraph();
        g5->SetPoint(0, x_star, y_star);  // Mark the minimum point
        g5->SetMarkerStyle(kFullStar);
        g5->SetMarkerSize(2.5);
        g5->Draw("same p"); 


        // Find the minimum z-value and its corresponding x and y values
        int minBinX = -1, minBinY = -1;
        double minX, minY;
        double minZ = chi2_prof->GetMinimum();

        for (int i = 1; i <= chi2_prof->GetNbinsX(); ++i) {
            for (int j = 1; j <= chi2_prof->GetNbinsY(); ++j) {
                double z = chi2_prof->GetBinContent(i, j);
                if (z == minZ) {
                    minBinX = i;
                    minBinY = j;
                    break;
                }
            }
        }
        minX = chi2_prof->GetXaxis()->GetBinCenter(minBinX);
        minY = chi2_prof->GetYaxis()->GetBinCenter(minBinY);

        // 2D fit to the full histogram to find the global minimum
        TF2* f2 = new TF2("f2", "[0]*x*x + [1]*y*y + [2]*x*y + [3]*x + [4]*y + [5]", 0.98, 1.02, 0.98, 1.02);
        TGraph2D* g2d = new TGraph2D();
        int idx = 0;
        for (int i = 1; i <= chi2_prof->GetNbinsX(); ++i) {
            for (int j = 1; j <= chi2_prof->GetNbinsY(); ++j) {
                double x = chi2_prof->GetXaxis()->GetBinCenter(i);
                double y = chi2_prof->GetYaxis()->GetBinCenter(j);
                double z = chi2_prof->GetBinContent(i, j);
                if (z > 0) {
                    g2d->SetPoint(idx++, x, y, z);
                }
            }
        }

        // Now fit the graph with f2
        g2d->Fit(f2, "Q0");
        //f2->Draw("CONT3 SAME");
        double p0 = f2->GetParameter(0);
        double p1 = f2->GetParameter(1);
        double p2 = f2->GetParameter(2);
        double p3 = f2->GetParameter(3);
        double p4 = f2->GetParameter(4);

        double denom = 4*p0*p1 - p2*p2;
        double x_star3 = (p2*p4 - 2*p1*p3) / denom;
        double y_star3 = (p2*p3 - 2*p0*p4) / denom;
        double chi2_min_fit = f2->Eval(x_star3, y_star3);

        // White star under fitted one
        TGraph* g_old3 = new TGraph();
        g_old3->SetPoint(3, minX, minY);
        g_old3->SetMarkerStyle(kStar);
        g_old3->SetMarkerSize(2);
        g_old3->SetMarkerColor(kWhite);
        g_old3->Draw("same p");

        cout << "minX = " << minX << ", minY = " << minY << endl;

        // Add a second star to mark the parabola-fit minimum
        TGraph* g_fit = new TGraph();
        g_fit->SetPoint(3, x_star3, y_star3);
        g_fit->SetMarkerStyle(kFullStar);
        g_fit->SetMarkerColor(kRed);
        g_fit->SetMarkerSize(2.5);
        g_fit->Draw("same p");

        // Create a new transparent pad on top of the current plot
        TPad *overlayPad2 = new TPad("overlayPad2", "", 0, 0, 1, 1);
        overlayPad2->SetFillStyle(0);  // Make it transparent
        overlayPad2->SetFrameFillStyle(0);
        overlayPad2->Draw();
        overlayPad2->cd(); 

        // Create a TLatex object for annotation
        TLatex latex2;
        // Set font and size if needed (optional)
        latex2.SetTextSize(0.04); // Adjust size according to your plot's appearance
        // Add text for pT region (adjust coordinates based on plot needs)
        latex2.DrawLatex(0.56, 0.86, "p_{T} > 30 GeV");
        // Add text for eta region (adjust coordinates based on plot needs)
        latex2.DrawLatex(0.56, 0.81, "|#eta| < 2.5");
        latex2.DrawLatex(0.56, 0.755, "#frac{#chi^{2}}{NDF}, NDF = 30");
        latex2.DrawLatex(0.16, 0.86, Form("dia: min #chi^{2}/NDF = %.2f, ", chi2_min2));
        latex2.DrawLatex(0.16, 0.81, Form("off: min #chi^{2}/NDF = %.2f, ", chi2_min));
        latex2.DrawLatex(0.16, 0.75, Form("glb: min #chi^{2}/NDF = %.2f, ", chi2_min_fit));
        TLatex latexOld;
        latexOld.SetTextSize(0.04);
        latexOld.SetTextColor(kWhite);  // Or any other ROOT color code
        latexOld.DrawLatex(0.448, 0.86, Form("%.2f", chi2_min_old2));
        latexOld.DrawLatex(0.441, 0.81, Form("%.2f", chi2_min_old));
        latexOld.DrawLatex(0.448, 0.75, Form("%.2f", minZ));

        gPad->Modified();
        gPad->Update();
        c3->SaveAs("frames/csvsud.pdf");
        TFile* outFileDiag = new TFile("frames/chi2_prof.root", "RECREATE");
        chi2_prof->Write("chi2_prof");
        f2->Write("f2_function");
        g2->Write("g2_function");
        g3->Write("g3_function");
        outFileDiag->Close();
}