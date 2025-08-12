#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>
#include "../minitools/tdrstyle_mod22.C"

TH1D *hisr_up(0), *hisr_dw(0);
TH1D *hfsr_up(0), *hfsr_dw(0);
TH1D *hjer_up(0), *hjer_dw(0);
TH1D *hist_up(0), *hist_dw(0);

TH1D *hisr_up_scaled(0), *hisr_dw_scaled(0);
TH1D *hfsr_up_scaled(0), *hfsr_dw_scaled(0);
TH1D *hjer_up_scaled(0), *hjer_dw_scaled(0);

Double_t hist_func(Double_t *x, Double_t *p) {
    assert(hist_up);
    assert(hist_dw);
    double mass = x[0];
    int i = hisr_up->FindBin(mass);
    double up = hist_up->GetBinContent(i);
    double dw = hist_dw->GetBinContent(i);
    double slope = (up - dw) / 2.;
    return p[0] * slope;
}

Double_t fit_func(Double_t *x, Double_t *p) {
    hist_up = hisr_up; hist_dw = hisr_dw;
    double isr = hist_func(x, &p[0]);

    hist_up = hfsr_up; hist_dw = hfsr_dw;
    double fsr = hist_func(x, &p[1]);

    hist_up = hjer_up; hist_dw = hjer_dw;
    double jer = hist_func(x, &p[2]);

    double mass = x[0];
    int i = hist_up->FindBin(mass);
    if (i == 98) { isr = fsr = jer = 0; isr = p[0]; }
    if (i == 99) { isr = fsr = jer = 0; fsr = p[1]; }
    if (i == 100) { isr = fsr = jer = 0; jer = p[2]; }

    return (1 + isr + fsr + jer);
}

Double_t fit_func_scaled(Double_t *x, Double_t *p) {
    hist_up = hisr_up_scaled; hist_dw = hisr_dw_scaled;
    double isr = hist_func(x, &p[0]);

    hist_up = hfsr_up_scaled; hist_dw = hfsr_dw_scaled;
    double fsr = hist_func(x, &p[1]);

    hist_up = hjer_up_scaled; hist_dw = hjer_dw_scaled;
    double jer = hist_func(x, &p[2]);

    double mass = x[0];
    int i = hist_up->FindBin(mass);
    if (i == 98) { isr = fsr = jer = 0; isr = p[0]; }
    if (i == 99) { isr = fsr = jer = 0; fsr = p[1]; }
    if (i == 100) { isr = fsr = jer = 0; jer = p[2]; }

    return (1 + isr + fsr + jer);
}
            // Vectors to store FSR and JER values and errors
            std::vector<double> fsr_values, fsr_errors, jer_values, jer_errors;
            std::vector<std::string> valid_suffixes;

void drawWqqFSRscales() {
    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    TFile *fisr = new TFile("processed_histograms/histograms_ISR_psweights.root", "READ");
    if (!fisr || fisr->IsZombie()) {
        std::cerr << "Error: Could not open histograms_ISR.root." << std::endl;
        return;
    }

    TFile *ffsr = new TFile("processed_histograms/histograms_FSR.root", "READ");
    if (!ffsr || ffsr->IsZombie()) {
        std::cerr << "Error: Could not open histograms_FSR.root." << std::endl;
        fisr->Close();
        delete fisr;
        return;
    }

    TFile *fjer = new TFile("processed_histograms/histograms_JER.root", "READ");
    if (!fjer || fjer->IsZombie()) {
        std::cerr << "Error: Could not open histograms_JER.root." << std::endl;
        fisr->Close();
        ffsr->Close();
        delete fisr;
        delete ffsr;
        return;
    }

    TFile *fdata = new TFile("processed_histograms/histograms_FSR_psweights.root", "READ");
    if (!fdata || fdata->IsZombie()) {
        std::cerr << "Error: Could not open histograms_FSR_psweights.root." << std::endl;
        fisr->Close();
        ffsr->Close();
        fjer->Close();
        delete fisr;
        delete ffsr;
        delete fjer;
        return;
    }

    // Load ISR and JER histograms
    //hisr_dw = (TH1D*)fisr->Get("RatioHist_ISR070");
    //hisr_up = (TH1D*)fisr->Get("RatioHist_ISR130");
    hisr_dw = (TH1D*)fisr->Get("RatioHist_isr:murfac=0.25");
    hisr_up = (TH1D*)fisr->Get("RatioHist_isr:murfac=4.0");
    hfsr_dw = (TH1D*)ffsr->Get("RatioHist_FSR0995");
    hfsr_up = (TH1D*)ffsr->Get("RatioHist_FSR1005");
    hjer_dw = (TH1D*)fjer->Get("RatioHist_JER090");
    hjer_up = (TH1D*)fjer->Get("RatioHist_JER110");

    if (!hisr_dw || !hisr_up || !hfsr_dw || !hfsr_up || !hjer_dw || !hjer_up) {
        std::cerr << "Error: One or more ISR/FSR/JER histograms not found." << std::endl;
        fisr->Close();
        ffsr->Close();
        fjer->Close();
        fdata->Close();
        delete fisr;
        delete ffsr;
        delete fjer;
        delete fdata;
        return;
    }

    // Clone scaled histograms
    hisr_up_scaled = (TH1D*)hisr_up->Clone("hisr_up_scaled");
    hisr_dw_scaled = (TH1D*)hisr_dw->Clone("hisr_dw_scaled");
    hfsr_up_scaled = (TH1D*)hfsr_up->Clone("hfsr_up_scaled");
    hfsr_dw_scaled = (TH1D*)hfsr_dw->Clone("hfsr_dw_scaled");
    hjer_up_scaled = (TH1D*)hjer_up->Clone("hjer_up_scaled");
    hjer_dw_scaled = (TH1D*)hjer_dw->Clone("hjer_dw_scaled");

    // Define FSR suffixes
    std::vector<std::string> fsr_suffixes = {
        "fsr:murfac=0.25",  "fsr:murfac=0.5", "fsr:murfac=2.0", "fsr:murfac=4.0"
    };

    // Loop over each FSR suffix
    for (const auto& suffix : fsr_suffixes) {
        // Load data histogram for the current suffix
        std::string histName = "RatioHist_" + suffix;
        TH1D *hdata = (TH1D*)fdata->Get(histName.c_str());
        if (!hdata) {
            std::cerr << "Error: Histogram " << histName << " not found." << std::endl;
            continue;
        }

        // Co-opt bins 98,99,100 for parameter constraints
        hdata->SetBinContent(98, 1); hdata->SetBinError(98, 1);
        hdata->SetBinContent(99, 1); hdata->SetBinError(99, 1);
        hdata->SetBinContent(100, 1); hdata->SetBinError(100, 1);

        // Create histograms and canvas
        TH1D *h1 = tdrHist("h1", "Ratio to MC", 0.86, 1.21, "Mass (GeV)", 66., 97.);
        TH1D *h2 = tdrHist("h2", "Data/fit -1 (%)", -2, 2, "Mass (GeV)", 66., 97.);
        extraText = "Private";
        TCanvas *c = tdrDiCanvas(("c_" + suffix).c_str(), h1, h2, 4, 11);

        c->cd(1);

        // First fit
        TF1 *fit1 = new TF1("fit1", fit_func, 70, 100, 3);
        fit1->SetParName(0, "isr_up (+30%)");
        fit1->SetParName(1, "fsr_up (+0.5%)");
        fit1->SetParName(2, "jer_up (+10%)");
        //fit1->FixParameter(0,0);
        //fit1->FixParameter(1,0);
        //fit1->FixParameter(2,0);
        hdata->Fit(fit1, "RN");
        fit1->SetRange(66, 100);

        double chi2_1 = fit1->GetChisquare();
        int ndf_1 = fit1->GetNDF();
        double scaleFactor = (ndf_1 > 0) ? TMath::Sqrt(chi2_1 / ndf_1) : 1.0;

        // Scale errors
        auto scaleBinErrors = [&](TH1D* h, double sf) {
            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                double oldErr = h->GetBinError(i);
                double newErr = oldErr * sf;
                h->SetBinError(i, newErr);
            }
        };

        // Clone histograms for this iteration
        hisr_up_scaled = (TH1D*)hisr_up->Clone(("hisr_up_scaled_" + suffix).c_str());
        hisr_dw_scaled = (TH1D*)hisr_dw->Clone(("hisr_dw_scaled_" + suffix).c_str());
        hfsr_up_scaled = (TH1D*)hfsr_up->Clone(("hfsr_up_scaled_" + suffix).c_str());
        hfsr_dw_scaled = (TH1D*)hfsr_dw->Clone(("hfsr_dw_scaled_" + suffix).c_str());
        hjer_up_scaled = (TH1D*)hjer_up->Clone(("hjer_up_scaled_" + suffix).c_str());
        hjer_dw_scaled = (TH1D*)hjer_dw->Clone(("hjer_dw_scaled_" + suffix).c_str());

        // Scale errors
        scaleBinErrors(hisr_up_scaled, scaleFactor);
        scaleBinErrors(hisr_dw_scaled, scaleFactor);
        scaleBinErrors(hfsr_up_scaled, scaleFactor);
        scaleBinErrors(hfsr_dw_scaled, scaleFactor);
        scaleBinErrors(hjer_up_scaled, scaleFactor);
        scaleBinErrors(hjer_dw_scaled, scaleFactor);

         // Scale factor for errors
        double errorScale = 1.0;
        if (ndf_1 > 0) {
            errorScale = TMath::Sqrt(chi2_1 / (double)ndf_1);
        }

        TH1D *hdata_scaled = (TH1D*)hdata->Clone("hdata_scaled");
        for (int i = 1; i <= hdata_scaled->GetNbinsX(); ++i) {
            double oldErr = hdata_scaled->GetBinError(i);
            double newErr = oldErr * errorScale;
            hdata_scaled->SetBinError(i, newErr);
        }

        // Draw variations and data
        TLine *l = new TLine();
        l->SetLineStyle(kDashed);
        l->SetLineColor(kGray+2);
        l->DrawLine(66, 1, 97, 1);

        tdrDraw(hisr_up_scaled, "Pz", kOpenSquare, kBlue-9);
        tdrDraw(hisr_dw_scaled, "Pz", kFullSquare, kBlue);
        tdrDraw(hfsr_up_scaled, "Pz", kOpenCircle, kRed-9);
        tdrDraw(hfsr_dw_scaled, "Pz", kFullCircle, kRed);
        tdrDraw(hjer_up_scaled, "Pz", kFullDiamond, kGreen+2);
        tdrDraw(hjer_dw_scaled, "Pz", kOpenDiamond, kGreen-9);
        tdrDraw(hdata_scaled, "Pz", kFullCircle, kBlack);

        // Second fit with scaled errors
        TF1 *fit2 = new TF1("fit2", fit_func_scaled, 70, 100, 3);
        fit2->SetParName(0, "isr_up (+30%)");
        fit2->SetParName(1, "fsr_up (+0.5%)");
        fit2->SetParName(2, "jer_up (+10%)");
        //fit2->FixParameter(0,0);
        //fit2->FixParameter(1,0);
        //fit2->FixParameter(2,0);
        //if(suffix == "fsr:murfac=0.25") fit2->SetParLimits(2,-1,0);
        fit2->SetRange(66,100);
        fit2->SetLineColor(kBlack);
        fit2->SetLineWidth(3);
        fit2->Draw("SAME");
        hdata_scaled->Fit(fit2, "RN");

        // Add fit results
        double chi2 = fit2->GetChisquare();
        int ndf = fit2->GetNDF();
        double par0 = fit2->GetParameter(0);
        double par1 = fit2->GetParameter(1);
        double par2 = fit2->GetParameter(2);
        double err0 = fit2->GetParError(0);
        double err1 = fit2->GetParError(1);
        double err2 = fit2->GetParError(2);

        // Inside the loop, after fit2
        // Store FSR and JER values and errors
        fsr_values.push_back(1.+par1*0.005); // Store as 1 + par1 to match example
        fsr_errors.push_back(err1*0.005);
        jer_values.push_back(1.+par2*0.1); // Store as 1 + par2 to match example
        jer_errors.push_back(err2*0.1);
        valid_suffixes.push_back(suffix);

        TLatex latex;
        latex.SetTextSize(0.045);
        latex.SetNDC();
        latex.DrawLatex(0.32, 0.85, Form("#chi^{2}/NDF = %.2f / %d", chi2, ndf));
        latex.DrawLatex(0.32, 0.80, Form("ISR: %.6f #pm %.6f", par0, err0));
        latex.DrawLatex(0.32, 0.75, Form("FSR: %.6f #pm %.6f", par1, err1));
        latex.DrawLatex(0.32, 0.70, Form("JER: %.6f #pm %.6f", par2, err2));

        // Add legend
        TLegend *legend = tdrLeg(0.65, 0.88 - 0.04 * 8, 0.95, 0.88);
        legend->AddEntry(hisr_up_scaled, "ISR Up", "P");
        legend->AddEntry(hisr_dw_scaled, "ISR Down", "P");
        legend->AddEntry(hfsr_up_scaled, "FSR Up", "P");
        legend->AddEntry(hfsr_dw_scaled, "FSR Down", "P");
        legend->AddEntry(hjer_up_scaled, "JER Up", "P");
        legend->AddEntry(hjer_dw_scaled, "JER Down", "P");
        legend->AddEntry(hdata_scaled, ("Data (" + suffix + ")").c_str(), "P");
        legend->AddEntry(fit2, "Fit", "L");
        legend->Draw("SAME");

        c->cd(2);

        // Create fit histogram
        TH1D *h_fit = (TH1D*)hdata_scaled->Clone("h_fit");
        h_fit->Reset();
        int nBins = h_fit->GetNbinsX();
        for (int i = 1; i <= nBins; ++i) {
            double x = h_fit->GetBinCenter(i);
            double fit_value = fit2->Eval(x);
            h_fit->SetBinContent(i, fit_value);
            h_fit->SetBinError(i, 0);
        }

        // Create ratio histogram
        TH1D *h_ratio = (TH1D*)hdata_scaled->Clone("h_ratio");
        // ---- Get the fit results ----
  Double_t chi2_fit = fit2->GetChisquare(); // total chi-squared
  Int_t    ndf_fit  = fit2->GetNDF();       // number of degrees of freedom
        h_ratio->Divide(h_fit);
        for (int i = 1; i <= nBins; ++i) {
            double ratio = h_ratio->GetBinContent(i);
            double ratio_error = h_ratio->GetBinError(i) * TMath::Sqrt(chi2_fit / ndf_fit);
            double percent_diff = (ratio - 1.0) * 100.0;
            double percent_error = ratio_error * 100.0;
            h_ratio->SetBinContent(i, percent_diff);
            h_ratio->SetBinError(i, percent_error);
        }

        TLine *line = new TLine(66., 0.0, 97., 0.0);
        line->SetLineColor(kGray+1);
        line->SetLineStyle(2);
        line->Draw("SAME");

        h_ratio->Draw("SAME");
        h_ratio->SetMarkerColor(kBlack);
        h_ratio->SetMarkerStyle(kFullCircle);
        h_ratio->SetLineWidth(2);
        h_ratio->SetLineColor(kBlack);

        // Save the canvas
        std::string safe_suffix = suffix;
        std::replace(safe_suffix.begin(), safe_suffix.end(), ':', '_'); // Replace ':' with '_' for filename
        c->SaveAs(("FSRscales/FSR_" + safe_suffix + ".pdf").c_str());

        // Clean up
        delete c;
        delete h1;
        delete h2;
        delete hdata_scaled;
        delete h_fit;
        delete h_ratio;
        delete fit1;
        delete fit2;
        delete hisr_up_scaled;
        delete hisr_dw_scaled;
        delete hfsr_up_scaled;
        delete hfsr_dw_scaled;
        delete hjer_up_scaled;
        delete hjer_dw_scaled;

        std::cout << "Suffix: " << suffix << ", chi2_1/ndf_1 = " << chi2_1 << "/" << ndf_1 
        << ", scaleFactor = " << scaleFactor << std::endl;
    }

    // After the loop, before closing files
    // Print stored values as arrays
    std::cout << "const int nPoints = " << fsr_values.size() << ";" << std::endl;
    std::cout << "// FSR variations" << std::endl;
    std::cout << "const char* variations[nPoints] = {";
    for (size_t i = 0; i < valid_suffixes.size(); ++i) {
        std::cout << "\"" << valid_suffixes[i] << "\"";
        if (i < valid_suffixes.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;

    std::cout << "// JER values (y-axis)" << std::endl;
    std::cout << "double JER[nPoints] = {";
    for (size_t i = 0; i < jer_values.size(); ++i) {
        std::cout << jer_values[i];
        if (i < jer_values.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;

    std::cout << "// FSR values (y-axis)" << std::endl;
    std::cout << "double FSR[nPoints] = {";
    for (size_t i = 0; i < fsr_values.size(); ++i) {
        std::cout << fsr_values[i];
        if (i < fsr_values.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;

    std::cout << "// FSR errors" << std::endl;
    std::cout << "double FSRerr[nPoints] = {";
    for (size_t i = 0; i < fsr_errors.size(); ++i) {
        std::cout << fsr_errors[i];
        if (i < fsr_errors.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;

    std::cout << "// JER errors" << std::endl;
    std::cout << "double JERerr[nPoints] = {";
    for (size_t i = 0; i < jer_errors.size(); ++i) {
        std::cout << jer_errors[i];
        if (i < jer_errors.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;
        

    // Close files
    fisr->Close();
    ffsr->Close();
    fjer->Close();
    fdata->Close();
    delete fisr;
    delete ffsr;
    delete fjer;
    delete fdata;
}