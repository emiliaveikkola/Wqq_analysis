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
#include "tdrstyle_mod22.C"

void Wqq_Pt() {
    TFile *file = new TFile("output_MCRun2_Wqq_MassPt.root", "READ");
    TH3D* h3PtFlavorPairs_MC = (TH3D*)file->Get("h3PtFlavorPairs_MC");
    
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

      TFile *file2 = new TFile("output_DATARun2_Wqq_test4.root", "READ");
      TH3D* h3PtFlavorPairs_DATA = (TH3D*)file2->Get("h3PtFlavorPairs_DATA");

      TH1D* h3all_data = h3PtFlavorPairs_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

      TH1D* h3tagcsall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
      TH1D* h3tagudall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
      TH1D* h3tagxall_data = h3PtFlavorPairs_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);

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

    lumi_136TeV = "Run3 simulation";
    extraText = "Private";

        // Create a canvas
    TCanvas *c3 = new TCanvas("c3", "Canvas with 3x3 Grid", 3840, 2160);  // Increase size to 4K resolution

    // Divide the canvas into a 3x3 grid
    c3->Divide(4, 4);
    gStyle->SetOptStat(0);

    for (int i = 1; i <= 16; ++i) {
    c3->cd(i);
    
    // Set minimal margins inside each pad
    gPad->SetLeftMargin(0.12);  // Leave enough space for the y-axis title
    gPad->SetRightMargin(0.02);  // Minimal right margin
    gPad->SetTopMargin(0.02);  // Minimal top margin
    gPad->SetBottomMargin(0.12);
    gPad->Update();
    }
    c3->cd(1);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagx->SetFillColorAlpha(kAzure+7,0.25);

    h3gencstagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagx->GetYaxis()->SetTitle("N frac");
    h3gencstagx->SetTitle("");

    // Set the x-axis range
    h3gencstagx->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3gencstagx->GetYaxis()->SetRangeUser(0, 0.051);


    TLegend *leg1= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg1->AddEntry(h3gencstagx, "MC", "FPLE");

     // Update the canvas to reflect the changes
    TLatex *tex1 = new TLatex();
    tex1->SetNDC(); tex1->SetTextSize(0.045);
    tex1->DrawLatex(0.17,0.85,"gencs,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(2);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

        // Draw the histogram
    tdrDraw(h3genudtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagx->SetFillColorAlpha(kAzure+7,0.25);

    h3genudtagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagx->GetYaxis()->SetTitle("N frac");
    h3genudtagx->SetTitle("");

    // Set the x-axis range
    h3genudtagx->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3genudtagx->GetYaxis()->SetRangeUser(0, 0.051);

    TLegend *leg2= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg2->AddEntry(h3genudtagx, "MC", "FPLE");

    TLatex *tex2 = new TLatex();
    tex2->SetNDC(); tex2->SetTextSize(0.045);
    tex2->DrawLatex(0.17,0.85,"genud,tagx");
    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(3);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagx,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagx->SetFillColorAlpha(kAzure+7,0.25);

    h3genxtagx->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagx->GetYaxis()->SetTitle("N frac");
    h3genxtagx->SetTitle("");

        // Set the x-axis range
    h3genxtagx->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3genxtagx->GetYaxis()->SetRangeUser(0, 0.051);

    TLegend *leg3= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg3->AddEntry(h3genxtagx, "MC", "FPLE");


    TLatex *tex3 = new TLatex();
    tex3->SetNDC(); tex3->SetTextSize(0.045);
    tex3->DrawLatex(0.17,0.85,"genx,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(4);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagxall->GetYaxis()->SetTitle("N frac");
    h3tagxall->SetTitle("");

        // Set the x-axis range
    h3tagxall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg4 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg4->AddEntry(h3tagxall, "MC", "FPLE");
    leg4->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex4 = new TLatex();
    tex4->SetNDC(); tex4->SetTextSize(0.045);
    tex4->DrawLatex(0.17,0.85,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();



    c3->cd(5);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagud->SetFillColorAlpha(kAzure+7,0.25);

    h3gencstagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagud->GetYaxis()->SetTitle("N frac");
    h3gencstagud->GetYaxis()->SetTitleOffset(1.3);
    h3gencstagud->SetTitle("");

        // Set the x-axis range
    h3gencstagud->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3gencstagud->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed


    TLegend *leg5= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg5->AddEntry(h3gencstagud, "MC", "FPLE");


    TLatex *tex5 = new TLatex();
    tex5->SetNDC(); tex5->SetTextSize(0.045);
    tex5->DrawLatex(0.17,0.85,"gencs,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(6);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagud->SetFillColorAlpha(kAzure+7,0.25);

    h3genudtagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagud->GetYaxis()->SetTitle("N frac");
    h3genudtagud->GetYaxis()->SetTitleOffset(1.3);
    h3genudtagud->SetTitle("");

        // Set the x-axis range
    h3genudtagud->GetXaxis()->SetRangeUser(-1,1);
    h3genudtagud->GetYaxis()->SetRangeUser(0, 0.051);  // Adjust the range as needed

    TLegend *leg6= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg6->AddEntry(h3genudtagud, "MC", "FPLE");

    TLatex *tex6 = new TLatex();
    tex6->SetNDC(); tex6->SetTextSize(0.045);
    tex6->DrawLatex(0.17,0.85,"genud,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(7);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagud,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagud->SetFillColorAlpha(kAzure+7,0.25);

    h3genxtagud->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagud->GetYaxis()->SetTitle("N frac");
    h3genxtagud->GetYaxis()->SetTitleOffset(1.3);
    h3genxtagud->SetTitle("");

        // Set the x-axis range
    h3genxtagud->GetXaxis()->SetRangeUser(-1,1);
    h3genxtagud->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg7 = tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg7->AddEntry(h3genxtagud, "MC", "FPLE");

    TLatex *tex7 = new TLatex();
    tex7->SetNDC(); tex7->SetTextSize(0.045);
    tex7->DrawLatex(0.17,0.85,"genx,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(8);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagudall->GetYaxis()->SetTitle("N frac");
    h3tagudall->SetTitle("");

        // Set the x-axis range
    h3tagudall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg8 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg8->AddEntry(h3tagudall, "MC", "FPLE");
    leg8->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex8 = new TLatex();
    tex8->SetNDC(); tex8->SetTextSize(0.045);
    tex8->DrawLatex(0.17,0.85,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(9);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencstagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencstagcs->SetFillColorAlpha(kAzure+7,0.25);

    h3gencstagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencstagcs->GetYaxis()->SetTitle("N frac");
    h3gencstagcs->GetYaxis()->SetTitleOffset(1.3);
    h3gencstagcs->SetTitle("");

        // Set the x-axis range
    h3gencstagcs->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3gencstagcs->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg9 = tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg9->AddEntry(h3gencstagcs, "MC", "FPLE");


    TLatex *tex9 = new TLatex();
    tex9->SetNDC(); tex9->SetTextSize(0.045);
    tex9->DrawLatex(0.17,0.85,"gencs,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(10);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudtagcs->SetFillColorAlpha(kAzure+7,0.25);



    h3genudtagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudtagcs->GetYaxis()->SetTitle("N frac");
    h3genudtagcs->SetTitle("");
    h3genudtagcs->GetYaxis()->SetNoExponent();


    TLegend *leg10 = tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg10->AddEntry(h3genudtagcs, "MC", "FPLE");


        // Set the x-axis range
    h3genudtagcs->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3genudtagcs->GetYaxis()->SetRangeUser(0, 0.051);
    
    TLatex *tex10 = new TLatex();
    tex10->SetNDC(); tex10->SetTextSize(0.045);
    tex10->DrawLatex(0.17,0.85,"genud,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(11);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxtagcs,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxtagcs->SetFillColorAlpha(kAzure+7,0.25);


    h3genxtagcs->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxtagcs->GetYaxis()->SetTitle("N frac");
    h3genxtagcs->SetTitle("");

        // Set the x-axis range
    h3genxtagcs->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3genxtagcs->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed


    TLegend *leg11 = tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg11->AddEntry(h3genxtagcs, "MC", "FPLE");


    TLatex *tex11 = new TLatex();
    tex11->SetNDC(); tex11->SetTextSize(0.045);
    tex11->DrawLatex(0.17,0.85,"genx,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(12);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);


    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");
    h3tagcsall->GetYaxis()->SetTitleOffset(1.3);
    h3tagcsall->SetTitle("");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg12 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg12->AddEntry(h3tagcsall, "MC", "FPLE");
    leg12->AddEntry(h3tagcsall_data, "DATA", "PLE");

    TLatex *tex12 = new TLatex();
    tex12->SetNDC(); tex12->SetTextSize(0.045);
    tex12->DrawLatex(0.17,0.85,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(13);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3gencsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3gencsall->SetFillColorAlpha(kAzure+7,0.25);


    h3gencsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3gencsall->GetYaxis()->SetTitle("N frac");
    h3gencsall->SetTitle("");

        // Set the x-axis range
    h3gencsall->GetXaxis()->SetRangeUser(-1,1);
    h3gencsall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg13= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg13->AddEntry(h3gencsall, "MC", "FPLE");

    TLatex *tex13 = new TLatex();
    tex13->SetNDC(); tex13->SetTextSize(0.045);
    tex13->DrawLatex(0.17,0.85,"gencs,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();

    c3->cd(14);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genudall->SetFillColorAlpha(kAzure+7,0.25);

    h3genudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genudall->GetYaxis()->SetTitle("N frac");
    h3genudall->SetTitle("");

        // Set the x-axis range
    h3genudall->GetXaxis()->SetRangeUser(-1,1);
    h3genudall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg14= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg14->AddEntry(h3genudall, "MC", "FPLE");

    TLatex *tex14 = new TLatex();
    tex14->SetNDC(); tex14->SetTextSize(0.045);
    tex14->DrawLatex(0.17,0.85,"genud,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(15);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3genxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3genxall->SetFillColorAlpha(kAzure+7,0.25);

    h3genxall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3genxall->GetYaxis()->SetTitle("N frac");
    h3genxall->SetTitle("");

        // Set the x-axis range
    h3genxall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3genxall->GetYaxis()->SetRangeUser(0, 0.051);

    TLegend *leg15= tdrLeg(0.75,0.8-0.05,0.9,0.82);
    leg15->AddEntry(h3genxall, "MC", "FPLE");

    TLatex *tex15 = new TLatex();
    tex15->SetNDC(); tex15->SetTextSize(0.045);
    tex15->DrawLatex(0.17,0.85,"genx,tagall");

    // Update the canvas to reflect the changes
    gPad->Update();


    c3->cd(16);

    //gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.08);
    //gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3all->GetYaxis()->SetTitle("N frac");
    h3all->SetTitle("");

        // Set the x-axis range
    h3all->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3all->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg16 = tdrLeg(0.75,0.8-0.05*2,0.9,0.82);
    leg16->AddEntry(h3all, "MC", "FPLE");
    leg16->AddEntry(h3all_data, "DATA", "PLE");
    

    TLatex *tex16 = new TLatex();
    tex16->SetNDC(); tex16->SetTextSize(0.045);
    tex16->DrawLatex(0.17,0.85,"all");



 h3gencstagx->GetXaxis()->SetTitleSize(0.06);  
h3gencstagx->GetXaxis()->SetLabelSize(0.05);  
h3gencstagx->GetYaxis()->SetTitleSize(0.06);  
h3gencstagx->GetYaxis()->SetLabelSize(0.05);  

h3genudtagx->GetXaxis()->SetTitleSize(0.06);  
h3genudtagx->GetXaxis()->SetLabelSize(0.05);  
h3genudtagx->GetYaxis()->SetTitleSize(0.06);  
h3genudtagx->GetYaxis()->SetLabelSize(0.05);  

h3genxtagx->GetXaxis()->SetTitleSize(0.06);  
h3genxtagx->GetXaxis()->SetLabelSize(0.05);  
h3genxtagx->GetYaxis()->SetTitleSize(0.06);  
h3genxtagx->GetYaxis()->SetLabelSize(0.05);  

h3tagxall->GetXaxis()->SetTitleSize(0.06);  
h3tagxall->GetXaxis()->SetLabelSize(0.05);  
h3tagxall->GetYaxis()->SetTitleSize(0.06);  
h3tagxall->GetYaxis()->SetLabelSize(0.05);  

h3gencstagud->GetXaxis()->SetTitleSize(0.06);  
h3gencstagud->GetXaxis()->SetLabelSize(0.05);  
h3gencstagud->GetYaxis()->SetTitleSize(0.06);  
h3gencstagud->GetYaxis()->SetLabelSize(0.05);  

h3genudtagud->GetXaxis()->SetTitleSize(0.06);  
h3genudtagud->GetXaxis()->SetLabelSize(0.05);  
h3genudtagud->GetYaxis()->SetTitleSize(0.06);  
h3genudtagud->GetYaxis()->SetLabelSize(0.05);  

h3genxtagud->GetXaxis()->SetTitleSize(0.06);  
h3genxtagud->GetXaxis()->SetLabelSize(0.05);  
h3genxtagud->GetYaxis()->SetTitleSize(0.06);  
h3genxtagud->GetYaxis()->SetLabelSize(0.05);  

h3tagudall->GetXaxis()->SetTitleSize(0.06);  
h3tagudall->GetXaxis()->SetLabelSize(0.05);  
h3tagudall->GetYaxis()->SetTitleSize(0.06);  
h3tagudall->GetYaxis()->SetLabelSize(0.05);  

h3gencstagcs->GetXaxis()->SetTitleSize(0.06);  
h3gencstagcs->GetXaxis()->SetLabelSize(0.05);  
h3gencstagcs->GetYaxis()->SetTitleSize(0.06);  
h3gencstagcs->GetYaxis()->SetLabelSize(0.05);  

h3genudtagcs->GetXaxis()->SetTitleSize(0.06);  
h3genudtagcs->GetXaxis()->SetLabelSize(0.05);  
h3genudtagcs->GetYaxis()->SetTitleSize(0.06);  
h3genudtagcs->GetYaxis()->SetLabelSize(0.05);  

h3genxtagcs->GetXaxis()->SetTitleSize(0.06);  
h3genxtagcs->GetXaxis()->SetLabelSize(0.05);  
h3genxtagcs->GetYaxis()->SetTitleSize(0.06);  
h3genxtagcs->GetYaxis()->SetLabelSize(0.05);  

h3tagcsall->GetXaxis()->SetTitleSize(0.06);  
h3tagcsall->GetXaxis()->SetLabelSize(0.05);  
h3tagcsall->GetYaxis()->SetTitleSize(0.06);  
h3tagcsall->GetYaxis()->SetLabelSize(0.05);  

h3gencsall->GetXaxis()->SetTitleSize(0.06);  
h3gencsall->GetXaxis()->SetLabelSize(0.05);  
h3gencsall->GetYaxis()->SetTitleSize(0.06);  
h3gencsall->GetYaxis()->SetLabelSize(0.05);  

h3genudall->GetXaxis()->SetTitleSize(0.06);  
h3genudall->GetXaxis()->SetLabelSize(0.05);  
h3genudall->GetYaxis()->SetTitleSize(0.06);  
h3genudall->GetYaxis()->SetLabelSize(0.05);  

h3genxall->GetXaxis()->SetTitleSize(0.06);  
h3genxall->GetXaxis()->SetLabelSize(0.05);  
h3genxall->GetYaxis()->SetTitleSize(0.06);  
h3genxall->GetYaxis()->SetLabelSize(0.05);  

h3all->GetXaxis()->SetTitleSize(0.06);  
h3all->GetXaxis()->SetLabelSize(0.05);  
h3all->GetYaxis()->SetTitleSize(0.06);  
h3all->GetYaxis()->SetLabelSize(0.05);


    
    h3gencstagcs->GetYaxis()->SetTitleOffset(1.1);
    h3gencstagud->GetYaxis()->SetTitleOffset(1.1);
    h3gencstagx->GetYaxis()->SetTitleOffset(1.1);
    h3all->GetYaxis()->SetTitleOffset(1.1);
    h3genudtagcs->GetYaxis()->SetTitleOffset(1.1);
    h3genudtagud->GetYaxis()->SetTitleOffset(1.1);
    h3genudtagx->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagcs->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagud->GetYaxis()->SetTitleOffset(1.1);
    h3genxtagx->GetYaxis()->SetTitleOffset(1.1);

    h3gencsall->GetYaxis()->SetTitleOffset(1.1);
    h3genudall->GetYaxis()->SetTitleOffset(1.1);
    h3genxall->GetYaxis()->SetTitleOffset(1.1);
    h3tagcsall->GetYaxis()->SetTitleOffset(1.1);
    h3tagudall->GetYaxis()->SetTitleOffset(1.1);
    h3tagxall->GetYaxis()->SetTitleOffset(1.1);

    // Update the canvas to reflect the changes
    gPad->Update();


    // Save the canvas as a .pdf file
    c3->SaveAs("pdf/Wqq_pt1.pdf");







setTDRStyle();
        // Create a canvas
    TCanvas *c100 = new TCanvas("c100", "Canvas with 3x3 Grid", 2560, 1140);
    

    // Divide the canvas into a 3x3 grid
    c100->Divide(4, 2);

    c100->cd(1);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagcsall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagcsall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagcsall_data,"Pz",kFullCircle,kPink-9);
    h3tagcsall_data->SetMarkerSize(1);

    h3tagcsall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagcsall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagcsall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagcsall->GetYaxis()->SetRangeUser(0, 0.051);
    h3tagcsall->GetYaxis()->SetTitleSize(0.045);
    h3tagcsall->GetXaxis()->SetTitleSize(0.045);
    h3tagcsall->GetYaxis()->SetLabelSize(0.045);
    h3tagcsall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagcsall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg103 = tdrLeg(0.68,0.9-0.05*2,0.85,0.9);
    leg103->AddEntry(h3tagcsall, "MC", "FPLE");
    leg103->AddEntry(h3tagcsall_data, "DATA", "PLE");
    //leg103->SetTextSize(0.037);

    TLatex *tex103 = new TLatex();
    tex103->SetNDC(); tex103->SetTextSize(0.045);
    tex103->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(2);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagudall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagudall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagudall_data,"Pz",kFullCircle,kPink-9);
    h3tagudall_data->SetMarkerSize(1);

    h3tagudall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagudall->GetYaxis()->SetTitle("N frac");

        // Set the x-axis range
    h3tagudall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed
    h3tagudall->GetYaxis()->SetTitleSize(0.045);
    h3tagudall->GetXaxis()->SetTitleSize(0.045);
    h3tagudall->GetYaxis()->SetLabelSize(0.045);
    h3tagudall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagudall->GetYaxis()->SetTitleOffset(1.65);

    TLegend *leg102 = tdrLeg(0.68,0.9-0.05*2,0.85,0.9);
    leg102->AddEntry(h3tagudall, "MC", "FPLE");
    leg102->AddEntry(h3tagudall_data, "DATA", "PLE");


    TLatex *tex102 = new TLatex();
    tex102->SetNDC(); tex102->SetTextSize(0.045);
    tex102->DrawLatex(0.17,0.92,"genall,tagud");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(3);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram
    tdrDraw(h3tagxall,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3tagxall->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3tagxall_data,"Pz",kFullCircle,kPink-9);
    h3tagxall_data->SetMarkerSize(1);

    h3tagxall->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3tagxall->GetYaxis()->SetTitle("N frac");
    h3tagxall->GetYaxis()->SetTitleSize(0.045);
    h3tagxall->GetXaxis()->SetTitleSize(0.045);
    h3tagxall->GetYaxis()->SetLabelSize(0.045);
    h3tagxall->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3tagxall->GetYaxis()->SetTitleOffset(1.65);
    

    // Set the x-axis range
    h3tagxall->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3tagxall->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg101 = tdrLeg(0.68,0.9-0.05*2,0.85,0.9);
    leg101->AddEntry(h3tagxall, "MC", "FPLE");
    leg101->AddEntry(h3tagxall_data, "DATA", "PLE");

    TLatex *tex101 = new TLatex();
    tex101->SetNDC(); tex101->SetTextSize(0.045);
    tex101->DrawLatex(0.17,0.92,"genall,tagx");

    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(4);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    // Draw the histogram

    tdrDraw(h3all,"HPz",kNone,kAzure+7,kSolid,-1,1001,kAzure+7);
    h3all->SetFillColorAlpha(kAzure+7,0.25);

    tdrDraw(h3all_data,"Pz",kFullCircle,kPink-9);
    h3all_data->SetMarkerSize(1);

    h3all->GetXaxis()->SetTitle("(s-c)/(s+c)");
    h3all->GetYaxis()->SetTitle("N frac");
    h3all->GetYaxis()->SetTitleSize(0.045);
    h3all->GetXaxis()->SetTitleSize(0.045);
    h3all->GetYaxis()->SetLabelSize(0.045);
    h3all->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h3all->GetYaxis()->SetTitleOffset(1.65);

    // Set the x-axis range
    h3all->GetXaxis()->SetRangeUser(-1,1); // Adjust the range as needed
    h3all->GetYaxis()->SetRangeUser(0, 0.051); // Adjust the range as needed

    TLegend *leg104 = tdrLeg(0.68,0.9-0.05*2,0.85,0.9);
    leg104->AddEntry(h3all, "MC", "FPLE");
    leg104->AddEntry(h3all_data, "DATA", "PLE");
    

    TLatex *tex104 = new TLatex();
    tex104->SetNDC(); tex104->SetTextSize(0.045);
    tex104->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(5);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagcsvsdata = (TH1D*)h3tagcsall->Clone("h_tagcsvsdata");
    h_tagcsvsdata->Divide(h3tagcsall_data);

    tdrDraw(h_tagcsvsdata,"Pz",kFullCircle,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1);


    h_tagcsvsdata->GetYaxis()->SetRangeUser(0, 2.1); // Adjust the range as needed
    h_tagcsvsdata->GetYaxis()->SetTitle("MC/DATA");
    h_tagcsvsdata->GetYaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetTitleSize(0.045);
    h_tagcsvsdata->GetYaxis()->SetLabelSize(0.045);
    h_tagcsvsdata->GetXaxis()->SetLabelSize(0.045);
    // Adjust the y-axis title location
    h_tagcsvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line107 = new TLine(-1, 1, 1, 1);
    line107->SetLineColor(kGray);
    line107->SetLineStyle(kDashed);
    line107->Draw();

    tdrDraw(h_tagcsvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagcsvsdata->SetMarkerSize(1.3);


    // Fit the ratio histogram with a constant function
    TF1* f7 = new TF1("f7", "[0]", -0.7, 0.7);
    f7->FixParameter(0, 1);

    h_tagcsvsdata->Fit(f7, "RN");
    double chi27 = f7->GetChisquare();
    int ndf7 = f7->GetNDF();


    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi27, ndf7));



    TLegend *leg107 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg107->AddEntry(h_tagcsvsdata, "MC", "PLE");    

    TLatex *tex107 = new TLatex();
    tex107->SetNDC(); tex107->SetTextSize(0.045);
    tex107->DrawLatex(0.17,0.92,"genall,tagcs");

    // Update the canvas to reflect the changes
    gPad->Update();
    gPad->RedrawAxis();

    c100->cd(6);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagudvsdata = (TH1D*)h3tagudall->Clone("h_tagudvsdata");
    h_tagudvsdata->Divide(h3tagudall_data);


    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);

    h_tagudvsdata->GetYaxis()->SetRangeUser(0, 2.1); // Adjust the range as needed
    h_tagudvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_tagudvsdata->GetYaxis()->SetTitleOffset(1.5);

    TLine *line106 = new TLine(-1, 1, 1, 1);
    line106->SetLineColor(kGray);
    line106->SetLineStyle(kDashed);
    line106->Draw();

    tdrDraw(h_tagudvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagudvsdata->SetMarkerSize(1.3);


    // Fit the ratio histogram with a constant function
    TF1* f6 = new TF1("f6", "[0]", -0.7, 0.7);
    f6->FixParameter(0, 1);

    h_tagudvsdata->Fit(f6, "RN");
    double chi26 = f6->GetChisquare();
    int ndf6 = f6->GetNDF();


    TLegend *leg106 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg106->AddEntry(h_tagudvsdata, "MC", "PLE");


    TLatex *tex106 = new TLatex();
    tex106->SetNDC(); tex106->SetTextSize(0.045);
    tex106->DrawLatex(0.17,0.92,"genall,tagud");

    // Draw the chi2/ndf value on the plot
    tex106->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi26, ndf6));


    // Update the canvas to reflect the changes
    gPad->Update();

    c100->cd(7);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_tagxvsdata = (TH1D*)h3tagxall->Clone("h_tagxvsdata");
    h_tagxvsdata->Divide(h3tagxall_data);


    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);

    h_tagxvsdata->GetYaxis()->SetRangeUser(0, 2.1); // Adjust the range as needed
    h_tagxvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_tagxvsdata->GetYaxis()->SetTitleOffset(1.5);

    // Add reference line at y = 1
    TLine *line105 = new TLine(-1, 1, 1, 1);
    line105->SetLineColor(kGray);
    line105->SetLineStyle(kDashed);
    line105->Draw();

    tdrDraw(h_tagxvsdata,"Pz",kFullSquare,kAzure+7);
    h_tagxvsdata->SetMarkerSize(1.3);


    // Fit the ratio histogram with a constant function
    TF1* f5 = new TF1("f5", "[0]", -0.7, 0.7);
    f5->FixParameter(0, 1);

    h_tagxvsdata->Fit(f5, "RN");
    double chi25 = f5->GetChisquare();
    int ndf5 = f5->GetNDF();
  

    TLegend *leg105 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg105->AddEntry(h_tagxvsdata, "MC", "PLE");

    

    TLatex *tex105 = new TLatex();
    tex105->SetNDC(); tex105->SetTextSize(0.045);
    tex105->DrawLatex(0.17,0.92,"genall,tagx");

    // Draw the chi2/ndf value on the plot
    tex105->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi25, ndf5));


    // Update the canvas to reflect the changes
    gPad->Update();


    c100->cd(8);

    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0);

    TH1D *h_allvsdata = (TH1D*)h3all->Clone("h_allvsdata");
    h_allvsdata->Divide(h3all_data);



    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);


    h_allvsdata->GetYaxis()->SetRangeUser(0, 2.1); // Adjust the range as needed
    h_allvsdata->GetYaxis()->SetTitle("MC/DATA");
    // Adjust the y-axis title location
    h_allvsdata->GetYaxis()->SetTitleOffset(1.5);
    h_allvsdata->GetXaxis()->SetTitleOffset(1.);

    // Add reference line at y = 1
    TLine *line108 = new TLine(-1, 1, 1, 1);
    line108->SetLineColor(kGray);
    line108->SetLineStyle(kDashed);
    line108->Draw();

    tdrDraw(h_allvsdata,"Pz",kFullSquare,kAzure+7);
    h_allvsdata->SetMarkerSize(1.3);

    // Fit the ratio histogram with a constant function
    TF1* f8 = new TF1("f8", "[0]", -0.7, 0.7);
    f8->FixParameter(0, 1);
    h_allvsdata->Fit(f8, "RN");
    double chi28 = f8->GetChisquare();
    int ndf8 = f8->GetNDF();
  
    // Draw the chi2/ndf value on the plot
    tex103->DrawLatex(0.17, 0.85, Form("#chi_{MC}^{2} / NDF = %1.1f / %d", chi28, ndf8));



    TLegend *leg108 = tdrLeg(0.68,0.9-0.1,0.85,0.9);
    leg108->AddEntry(h_allvsdata, "MC", "PLE");

    

    TLatex *tex108 = new TLatex();
    tex108->SetNDC(); tex108->SetTextSize(0.045);
    tex108->DrawLatex(0.17,0.92,"all");

    // Update the canvas to reflect the changes
    gPad->Update();



    // Save the canvas as a .pdf file
    c100->SaveAs("pdf/Wqq_pt2.pdf");

}