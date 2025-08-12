#include <TArrow.h>
#include <TArrow.h>

// Helper function to draw overflow arrows for bins above yMax
void drawOverflowArrow(TH1* h, double xMin, double xMax, double yMax) {
    int nBins = h->GetNbinsX();
    double arrowHeadSize = 0.02 * (xMax - xMin);

    for (int i = 1; i <= nBins; ++i) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);

        if (x < xMin || x > xMax) continue;
        if (y > yMax) {
            TArrow* arrow = new TArrow(x, yMax * 0.98, x, yMax * 0.85, arrowHeadSize, "|>");
            arrow->SetLineColor(kRed + 2);
            arrow->SetFillColor(kRed + 2);
            arrow->SetLineWidth(2);
            arrow->Draw("same");
        }
    }
}
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <algorithm>
#include "../minitools/tdrstyle_mod22.C"


void Wqq_MassPt() {
    TFile *file = new TFile("output_MCRun2_Wqq_MassPt.root", "READ");

    TH3D* h3MassFlavorPairs_MC = (TH3D*)file->Get("h3MassFlavorPairs_MC");
    
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

      TFile *file2 = new TFile("output_DATARun2_Wqq_test4.root", "READ");
      TH3D* h3MassFlavorPairs_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATA");

      TH1D* h3all_data = h3MassFlavorPairs_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

      TH1D* h3tagcsall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
      TH1D* h3tagudall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
      TH1D* h3tagxall_data = h3MassFlavorPairs_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);
      TH3D* h3PtFlavorPairs_MC = (TH3D*)file->Get("h3PtFlavorPairs_MC");
    
      // Project Z for different indices and draw
      TH1D* h3all_pt = h3PtFlavorPairs_MC->ProjectionZ("h3all_pt", 1, 3, 1, 3);

      TH1D* h3gencsall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3gencsall_pt", 1, 1, 1, 3);
      TH1D* h3genudall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genudall_pt", 2, 2, 1, 3);
      TH1D* h3genxall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genxall_pt", 3, 3, 1, 3);

      TH1D* h3tagcsall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3tagcsall_pt", 1, 3, 1, 1);
      TH1D* h3tagudall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3tagudall_pt", 1, 3, 2, 2);
      TH1D* h3tagxall_pt = h3PtFlavorPairs_MC->ProjectionZ("h3tagxall_pt", 1, 3, 3, 3);

      TH1D* h3gencstagcs_pt = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagcs_pt", 1, 1, 1, 1);
      TH1D* h3gencstagud_pt = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagud_pt", 1, 1, 2, 2);
      TH1D* h3gencstagx_pt = h3PtFlavorPairs_MC->ProjectionZ("h3gencstagx_pt", 1, 1, 3, 3);
      
      TH1D* h3genudtagcs_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagcs_pt", 2, 2, 1, 1);
      TH1D* h3genudtagud_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagud_pt", 2, 2, 2, 2);
      TH1D* h3genudtagx_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genudtagx_pt", 2, 2, 3, 3);
      
      TH1D* h3genxtagcs_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagcs_pt", 3, 3, 1, 1);
      TH1D* h3genxtagud_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagud_pt", 3, 3, 2, 2);
      TH1D* h3genxtagx_pt = h3PtFlavorPairs_MC->ProjectionZ("h3genxtagx_pt", 3, 3, 3, 3);

      TH3D* h3PtFlavorPairs_DATA = (TH3D*)file2->Get("h3PtFlavorPairs_DATA");

      TH1D* h3all_data_pt = h3PtFlavorPairs_DATA->ProjectionZ("h3all_data_pt", 1, 3, 1, 3);

      TH1D* h3tagcsall_data_pt = h3PtFlavorPairs_DATA->ProjectionZ("h3tagcsall_data_pt", 1, 3, 1, 1);
      TH1D* h3tagudall_data_pt = h3PtFlavorPairs_DATA->ProjectionZ("h3tagudall_data_pt", 1, 3, 2, 2);
      TH1D* h3tagxall_data_pt = h3PtFlavorPairs_DATA->ProjectionZ("h3tagxall_data_pt", 1, 3, 3, 3);


      h3gencstagcs->Scale(1./h3gencstagcs->Integral());
    h3gencstagud->Scale(1./h3gencstagud->Integral());
    h3gencstagx->Scale(1./h3gencstagx->Integral());
    
    h3genudtagcs->Scale(1./h3genudtagcs->Integral());
    h3genudtagud->Scale(1./h3genudtagud->Integral());
    h3genudtagx->Scale(1./h3genudtagx->Integral());

    h3genxtagcs->Scale(1./h3genxtagcs->Integral());
    h3genxtagud->Scale(1./h3genxtagud->Integral());
    h3genxtagx->Scale(1./h3genxtagx->Integral());

    h3gencsall->Scale(1./h3gencsall->Integral());
    h3genudall->Scale(1./h3genudall->Integral());
    h3genxall->Scale(1./h3genxall->Integral());

    h3tagcsall->Scale(1./h3tagcsall->Integral());
    h3tagudall->Scale(1./h3tagudall->Integral());
    h3tagxall->Scale(1./h3tagxall->Integral());

    h3all->Scale(1./h3all->Integral());

    h3tagcsall_data->Scale(1./h3tagcsall_data->Integral());
    h3tagudall_data->Scale(1./h3tagudall_data->Integral());
    h3tagxall_data->Scale(1./h3tagxall_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

    h3gencstagcs_pt->Scale(1./h3gencstagcs_pt->Integral());
    h3gencstagud_pt->Scale(1./h3gencstagud_pt->Integral());
    h3gencstagx_pt->Scale(1./h3gencstagx_pt->Integral());
    
    h3genudtagcs_pt->Scale(1./h3genudtagcs_pt->Integral());
    h3genudtagud_pt->Scale(1./h3genudtagud_pt->Integral());
    h3genudtagx_pt->Scale(1./h3genudtagx_pt->Integral());

    h3genxtagcs_pt->Scale(1./h3genxtagcs_pt->Integral());
    h3genxtagud_pt->Scale(1./h3genxtagud_pt->Integral());
    h3genxtagx_pt->Scale(1./h3genxtagx_pt->Integral());

    h3gencsall_pt->Scale(1./h3gencsall_pt->Integral());
    h3genudall_pt->Scale(1./h3genudall_pt->Integral());
    h3genxall_pt->Scale(1./h3genxall_pt->Integral());

    h3tagcsall_pt->Scale(1./h3tagcsall_pt->Integral());
    h3tagudall_pt->Scale(1./h3tagudall_pt->Integral());
    h3tagxall_pt->Scale(1./h3tagxall_pt->Integral());

    h3all_pt->Scale(1./h3all_pt->Integral());

    h3tagcsall_data_pt->Scale(1./h3tagcsall_data_pt->Integral());
    h3tagudall_data_pt->Scale(1./h3tagudall_data_pt->Integral());
    h3tagxall_data_pt->Scale(1./h3tagxall_data_pt->Integral());

    h3all_data_pt->Scale(1./h3all_data_pt->Integral());

    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

    setTDRStyle();
        // Create a canvas
    TCanvas *c100 = new TCanvas("c100", "Canvas with 3x3 Grid", 1860, 1440);
    

    // Divide the canvas into a 3x3 grid
    c100->Divide(2, 2);

    // Add CMS label to the left corner of the canvas
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextSize(0.04);   // Text size
    cmsLabel.SetTextFont(61);     // CMS font
    cmsLabel.DrawLatex(0.07, 0.95, "CMS");  // Coordinates (0.1, 0.92) for top left corner

    // Add "Preliminary" label next to CMS
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextSize(0.035);
    prelimLabel.SetTextFont(52);  // Italic font for "Preliminary"
    prelimLabel.DrawLatex(0.15, 0.95, "Private");  // Coordinates (0.2, 0.92)

    // Add dataset and luminosity information in the top right corner
    TLatex dataLabel;
    dataLabel.SetNDC();
    dataLabel.SetTextSize(0.038);
    dataLabel.SetTextFont(42);
    dataLabel.DrawLatex(0.62, 0.955, "138 fb^{-1} (Run2 Legacy, 13 TeV)");
    


    c100->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.09);
    gPad->SetTopMargin(0.1);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("Mass (GeV)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(55, 115); // Adjust the range as needed
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.051);
    h3tagcsall->GetYaxis()->SetTitleSize(0.045);
    h3tagcsall->GetXaxis()->SetTitleSize(0.045);
    h3tagcsall->GetYaxis()->SetLabelSize(0.045);
    h3tagcsall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagcsall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg101 = tdrLeg(0.68,0.8-0.05*2,0.85,0.8);
    leg101->AddEntry(h3tagcsall, "MC", "FPLE");
    leg101->AddEntry(h3tagcsall_data, "DATA", "PLE");
    //leg103->SetTextSize(0.037);

    TLatex *tex101 = new TLatex();
    tex101->SetNDC(); tex101->SetTextSize(0.045);
    tex101->DrawLatex(0.17,0.82,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    // Draw overflow arrows (toggle with if (false))
    if (false) {
        drawOverflowArrow(h3tagcsall_data, 55, 115, 0.051);
    }

    c100->cd(2);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.09);
    gPad->SetTopMargin(0.1);

    // Draw the histogram
    tdrDraw(h3tagcsall_pt,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall_pt->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_data_pt,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data_pt->SetMarkerSize(1);

    h3tagcsall_pt->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagcsall_pt->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall_pt->GetXaxis()->SetRangeUser(-1, 1); // Adjust the range as needed
    h3tagcsall_pt->GetYaxis()->SetRangeUser(0, 0.051);
    h3tagcsall_pt->GetYaxis()->SetTitleSize(0.045);
    h3tagcsall_pt->GetXaxis()->SetTitleSize(0.045);
    h3tagcsall_pt->GetYaxis()->SetLabelSize(0.045);
    h3tagcsall_pt->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagcsall_pt->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg102 = tdrLeg(0.68,0.8-0.05*2,0.85,0.8);
    leg102->AddEntry(h3tagcsall, "MC", "FPLE");
    leg102->AddEntry(h3tagcsall_data, "DATA", "PLE");
    //leg103->SetTextSize(0.037);

    TLatex *tex102 = new TLatex();
    tex102->SetNDC(); tex102->SetTextSize(0.045);
    tex102->DrawLatex(0.17,0.82,"cs-tagged");

    // Update the canvas to reflect the changes
    gPad->Update();

    // Draw overflow arrows (toggle with if (false))
    if (false) {
        drawOverflowArrow(h3tagcsall_data_pt, -1, 1, 0.051);
    }

    c100->cd(3);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagcsvsdata = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata");
    h_tagcsvsdata->Divide(h3tagcsall_data);

    tdrDraw(h_tagcsvsdata,"Pz",kFullCircle,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1);

    h_tagcsvsdata->GetYaxis()->SetRangeUser(0.8, 1.2); // Adjust the range as needed
    h_tagcsvsdata->GetXaxis()->SetRangeUser(55, 115);
    h_tagcsvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagcsvsdata->GetYaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetYaxis()->SetLabelSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h_tagcsvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line103 = new TLine(55, 1, 115, 1);
    line103->SetLineColor(kGray);
    line103->SetLineStyle(kDashed);
    line103->Draw();

    tdrDraw(h_tagcsvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f7 = new TF1("f7", "[0]", 70, 100);
    f7->FixParameter(0, 1);

    h_tagcsvsdata->Fit(f7, "RN");
    double chi27 = f7->GetChisquare();
    int ndf7 = f7->GetNDF();

    TLegend *leg103 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg103->AddEntry(h_tagcsvsdata, "MC", "PLE");    

    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.045);
    tex103->DrawLatex(0.17,0.92,"cs-tagged");

    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi27, ndf7));

    // Update the canvas to reflect the changes
    gPad->Update();
    gPad->RedrawAxis();

    // Draw overflow arrows (toggle with if (false)) for ratio plot
    if (false) {
        drawOverflowArrow(h_tagcsvsdata, 55, 115, 1.2);
    }

    c100->cd(4);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagcsvsdata_pt = (TH1D*)h3tagcsall_pt->Clone("h_tagcsvsdata_pt");
    h_tagcsvsdata_pt->Divide(h3tagcsall_data_pt);

    tdrDraw(h_tagcsvsdata_pt,"Pz",kFullCircle,kAzure+7);
    h_tagcsvsdata_pt->SetMarkerSize(1);

    h_tagcsvsdata_pt->GetYaxis()->SetRangeUser(0.8, 1.2); // Adjust the range as needed
    h_tagcsvsdata_pt->GetXaxis()->SetRangeUser(-1, 1);
    h_tagcsvsdata_pt->GetYaxis()->SetTitle("MC/DATA");
    h_tagcsvsdata_pt->GetYaxis()->SetTitleSize(0.045);
    h_tagcsvsdata_pt->GetXaxis()->SetTitleSize(0.045);
    h_tagcsvsdata_pt->GetYaxis()->SetLabelSize(0.045);
    h_tagcsvsdata_pt->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h_tagcsvsdata_pt->GetYaxis()->SetTitleOffset(1.5);

    TLine *line104 = new TLine(-1, 1, 1, 1);
    line104->SetLineColor(kGray);
    line104->SetLineStyle(kDashed);
    line104->Draw();

    tdrDraw(h_tagcsvsdata_pt,"Pz",kFullSquare,kAzure+7);
    h_tagcsvsdata_pt->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f8 = new TF1("f8", "[0]", -0.7, 0.7);
    f8->FixParameter(0, 1);

    h_tagcsvsdata_pt->Fit(f8, "RN");
    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();

    TLegend *leg104 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg104->AddEntry(h_tagcsvsdata_pt, "MC", "PLE");    

    TLatex *tex104 = new TLatex();
    tex104->SetNDC(); tex104->SetTextSize(0.045);
    tex104->DrawLatex(0.17,0.92,"genall,tagcs");

    // Draw the chi2/ndf value on the plot
    tex104->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi28, ndf8));

    // Update the canvas to reflect the changes
    gPad->Update();
    gPad->RedrawAxis();

    // Draw overflow arrows (toggle with if (false)) for ratio plot
    if (false) {
        drawOverflowArrow(h_tagcsvsdata_pt, -1, 1, 1.2);
    }


    std::cout << "[Ratio mass] Underflow bin: " << h_tagcsvsdata->GetBinContent(0) << std::endl;
std::cout << "[Ratio mass] Overflow bin: " << h_tagcsvsdata->GetBinContent(h_tagcsvsdata->GetNbinsX() + 1) << std::endl;
std::cout << "[Ratio pt] Underflow bin: " << h_tagcsvsdata_pt->GetBinContent(0) << std::endl;
std::cout << "[Ratio pt] Overflow bin: " << h_tagcsvsdata_pt->GetBinContent(h_tagcsvsdata_pt->GetNbinsX() + 1) << std::endl;
    // Save the canvas as a .pdf file
    c100->SaveAs("pdf/Wqq_masspt.pdf");
}