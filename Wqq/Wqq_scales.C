#define TagandProbe_cxx
#include "TagandProbe.h"
#include <TH2.h>
#include <TH3.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h>

#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <TLegend.h>
#include <TColor.h>
#include <TStopwatch.h>
#include "tdrstyle_mod22.C"
#include "JERHandler.h"

int Tag(double x) {
   if (x > 0.43) return 1;
   if (x < 0.06) return 2;
   return 3;
}

int Flav(double x) {
   if (x == 4 || x == -4) return 1; //c
   if (x == 1 || x == -1 || x == 2 || x == -2 || x == 3 || x  == -3) return 2; //uds
   return 3;                        //other
}

int DetermineJetPairDATA(double x, double y) {
   if (x > 0.43 || y > 0.43) return 1;
   if (x < 0.06 && y < 0.06) return 2;
   return 3;
}

int DetermineJetPairDATA_test(double x, double y) {
   if (x > 0.37 || y > 0.37) return 1;
   if (x > 0.38 || y > 0.38) return 2;
   if (x > 0.43 || y > 0.43) return 3;
   if (x > 0.60 || y > 0.60) return 4;
   if (x > 0.67 || y > 0.67) return 5;
   if (x < 0.06 && y < 0.06) return 6;
   return 7;
}

int DetermineJetPairMC(int x, int y) {
   if (x == 3 && y == -4 || x == -3 && y == 4 ||
   x == 4 && y == -3 || x == -4 && y == 3) return 1;
   if (x == 1 && y == -2 || x == -1 && y == 2 ||
   x == 2 && y == -1 || x == -2 && y == 1) return 2;
   return 3;
}

bool isExcludedPair(int a, int b) {
    // a should be 4, -4, 5, or -5 and b should not be from the excluded values (4, -4, 5, -5, 21)
    return (a == 4 || a == -4 || a == 5 || a == -5) &&
           (b != 4 && b != -4 && b != 5 && b != -5 && b != 21);
}

int DetermineJetPairMC_test(int x, int y) {
    // Check for category 1 (cX, anticX, bX, antibX)
    if (isExcludedPair(x, y) || isExcludedPair(y, x)) {
        return 1;
    }

    // Check for category 2 (dQ, antidQ, uQ, antiuQ)
    if ((x == 1 || x == -1 || x == 2 || x == -2) &&
        (y == 1 || y == -1 || y == 2 || y == -2)) {
        return 2;
    }

    // Default case
    return 3;
}


// Define the bin labels
const char* binLabels3[3] = {"c", "l", "x"};
const char* binLabels12[12] = {"antib", "antic", "antis", "antiu", "antid", "0", "d", "u", "s", "c", "b", "g"};
const char* binLabels7[7] = {"0", "d", "u", "s", "c", "b", "g"};

// Helper function to set bin labels for 2D histograms (TH2, TProfile2D)
void SetBinLabels2D(TAxis* xAxis, TAxis* yAxis, const char** labels, int n) {
    for (int i = 0; i < n; ++i) {
        xAxis->SetBinLabel(i + 1, labels[i]);
        yAxis->SetBinLabel(i + 1, labels[i]);
    }
}

// Helper function to set bin labels for 2D histograms (TH2, TProfile2D)
void SetBinLabels1D(TAxis* xAxis, const char** labels, int n) {
    for (int i = 0; i < n; ++i) {
        xAxis->SetBinLabel(i + 1, labels[i]);
    }
}

// Helper function to set bin labels for 3D histograms (TH3D)
void SetBinLabels3D(TH3* hist, const char** labels, int n) {
    for (int i = 0; i < n; ++i) {
        hist->GetXaxis()->SetBinLabel(i + 1, labels[i]);
        hist->GetYaxis()->SetBinLabel(i + 1, labels[i]);
    }
}

std::vector<std::vector<TH3D*>> initializeHistograms(const std::vector<double>& scaleFactors, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> innerVec;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 200, 0, 200);
            innerVec.push_back(hist);
        }
        histograms.push_back(innerVec);
    }

    return histograms;
}


std::vector<std::vector<TH3D*>> initializeHistograms2(const std::vector<double>& scaleFactors, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> innerVec;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 67, -1, 1);
            innerVec.push_back(hist);
        }
        histograms.push_back(innerVec);
    }

    return histograms;
}

std::vector<std::vector<TH3D*>> initializeHistograms1(const std::vector<double>& scales, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale : scales) {
        std::vector<TH3D*> innerHistograms;
        TString histName = TString::Format("%s_%.3f", namePrefix.c_str(), scale);
        TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 200, 0, 200);
        innerHistograms.push_back(hist);
        histograms.push_back(innerHistograms);
    }

    return histograms;
}

std::vector<std::vector<TH3D*>> initializeHistograms3(const std::vector<double>& scales, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale : scales) {
        std::vector<TH3D*> innerHistograms;
        TString histName = TString::Format("%s_%.3f", namePrefix.c_str(), scale);
        TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 67, -1, 1);
        innerHistograms.push_back(hist);
        histograms.push_back(innerHistograms);
    }

    return histograms;
}

std::vector<std::vector<TH3D*>> initializeHistograms4(const std::vector<double>& scaleFactors, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> innerVec;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 200, 0, 200);
            innerVec.push_back(hist);
        }
        histograms.push_back(innerVec);
    }

    return histograms;
}


std::vector<std::vector<TH3D*>> initializeHistograms5(const std::vector<double>& scaleFactors, const std::string& namePrefix, const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;

    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> innerVec;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 67, -1, 1);
            innerVec.push_back(hist);
        }
        histograms.push_back(innerVec);
    }

    return histograms;
}

// ---------------------------------------------------------------------------
// Per‑jet scale factors derived from the response fits
//  - Jet‑1:   cubic polynomial
//  - Jet‑2:   exponential + power‑law
// ---------------------------------------------------------------------------
inline double ScaleJet1Func(double pt)
{
    // c * (1 + a * pt^{n}) from fit_lin1
    constexpr double c = 1.01395406;
    constexpr double a = 3543.28304352;
    constexpr double n = -2.75716765;

    return c * (1.0 + a * std::pow(pt, n));
}

inline double ScaleJet2Func(double pt)
{
    const double inv  = 1.0 / pt;
    const double inv2 = inv * inv;
    const double inv3 = inv2 * inv;
    const double inv4 = inv3 * inv;

    constexpr double p0 = 0.80000000;
    constexpr double p1 = 37.51785306;
    constexpr double p2 = -2759.71029660;
    constexpr double p3 = 83772.06072387;
    constexpr double p4 = -737684.87977565;

    return p0 + p1*inv + p2*inv2 + p3*inv3 + p4*inv4;
}

void TagandProbe::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TagandProbe.C
//      root> TagandProbe t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   JERHandler jer;
   
      // Activate only necessary branches
   fChain->SetBranchStatus("*", 0); // Disable all branches

   fChain->SetBranchStatus("ctag1", 1);
   fChain->SetBranchStatus("ctag2", 1);
   fChain->SetBranchStatus("flav1", 1);
   fChain->SetBranchStatus("flav2", 1);
   fChain->SetBranchStatus("gen_pt1", 1);
   fChain->SetBranchStatus("gen_eta1", 1);
   fChain->SetBranchStatus("gen_phi1", 1);
   fChain->SetBranchStatus("gen_m1", 1);
   fChain->SetBranchStatus("gen_pt2", 1);
   fChain->SetBranchStatus("gen_eta2", 1);
   fChain->SetBranchStatus("gen_phi2", 1);
   fChain->SetBranchStatus("gen_m2", 1);
   fChain->SetBranchStatus("pt1", 1);
   fChain->SetBranchStatus("eta1", 1);
   fChain->SetBranchStatus("phi1", 1);
   fChain->SetBranchStatus("m1", 1);
   fChain->SetBranchStatus("pt2", 1);
   fChain->SetBranchStatus("eta2", 1);
   fChain->SetBranchStatus("phi2", 1);
   fChain->SetBranchStatus("m2", 1);
   fChain->SetBranchStatus("weight", 1);
   fChain->SetBranchStatus("fitProb", 1);
   fChain->SetBranchStatus("PSWgts", 1); // Needed for ISR/FSR systematics
   fChain->SetBranchStatus("pfRho", 1);

   TDirectory *curdir = gDirectory;
   // Create the output file based on a condition
   TFile* fout;

   bool DATA = false;
   bool SCALED = false;
   bool NOTSCALED = false;
   bool ITERATION = false;
   bool CHI2 = false;
   bool PSWEIGHTS = false;
   bool SMEAR = false;
   if (DATA) {fout = new TFile("output_DATARun2_tagprobe.root", "RECREATE");}
   else if (SCALED) {fout = new TFile("output_MCSCALEDRun2_tagprobe_ptgen1.2.root", "RECREATE");}
   else if (NOTSCALED) {fout = new TFile("output_MCRun2_tagprobe.root", "RECREATE");}
   else if (ITERATION) {fout = new TFile("output_MCRun2_tagprobe_iteration.root", "RECREATE");}
   else if (CHI2) {fout = new TFile("output_MCSCALEDRun2_tagprobe_chi2_FSR_up2_before.root", "RECREATE");}
   else if (PSWEIGHTS) {fout = new TFile("output_MCSCALEDRun2_tagprobe_scales.root", "RECREATE");}
   else if (SMEAR) {fout = new TFile("output_MCSCALEDRun2_tagprobe_JER_smearing2.root", "RECREATE");}



   else {fout = new TFile("false.root", "RECREATE");}

   TH1::SetDefaultSumw2();

   TH2D *h2tag_c = new TH2D("h2tag_c", ";reco probe;gen probe; N", 3, 1, 4, 3, 1, 4);
   TH2D *h2tag_u = new TH2D("h2tag_u", ";reco probe;gen probe; N", 3, 1, 4, 3, 1, 4);
   TH2D *h2tag_x = new TH2D("h2tag_x", ";reco probe;gen probe; N", 3, 1, 4, 3, 1, 4);

   TH2D* h2matrixN = new TH2D("h2matrixN", ";flav1;flav2; N", 12,-5.5,6.5,12,-5.5,6.5);
   TProfile2D* p2matrixMass = new TProfile2D("p2matrixMass", ";flav1;flav2; Mass", 12,-5.5,6.5,12,-5.5,6.5);
   TProfile2D* p2matrixPt = new TProfile2D("p2matrixPt", ";flav1;flav2; Pt", 12,-5.5,6.5,12,-5.5,6.5);
   TH3D* h3matrixMass = new TH3D("h3matrixMass", ";flav1;flav2; Mass", 12,-5.5,6.5,12,-5.5,6.5,200,0,200);
   TH3D *h3matrixPt = new TH3D("h3matrixPt", ";flav1;flav2; Pt", 12,-5.5,6.5,12,-5.5,6.5,67,-1,1);

   TH1D *h_gen_nom = new TH1D("h_gen_nom", ";N;Mass", 200,0,200);
   TH1D *h_gen_fsr = new TH1D("h_gen_fsr", ";N;Mass", 200,0,200);
   TH1D *h_gen_jer = new TH1D("h_gen_jer", ";N;Mass", 200,0,200);
   TH1D *h_gen_isr = new TH1D("h_gen_isr", ";N;Mass", 200,0,200);

   TH3D *h3MassFlavorPairs_DATAMC  = new TH3D("h3MassFlavorPairs_DATAMC", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC = new TH3D("h3PtFlavorPairs_DATAMC", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_up  = new TH3D("h3MassFlavorPairs_DATAMC_up", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_up = new TH3D("h3PtFlavorPairs_DATAMC_up", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_down  = new TH3D("h3MassFlavorPairs_DATAMC_down", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_down = new TH3D("h3PtFlavorPairs_DATAMC_down", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_up2  = new TH3D("h3MassFlavorPairs_DATAMC_up2", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_up2 = new TH3D("h3PtFlavorPairs_DATAMC_up2", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_down2  = new TH3D("h3MassFlavorPairs_DATAMC_down2", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_down2 = new TH3D("h3PtFlavorPairs_DATAMC-down2", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_gen  = new TH3D("h3MassFlavorPairs_DATAMC_gen", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_gen = new TH3D("h3PtFlavorPairs_DATAMC_gen", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_RESP = new TH3D("h3MassFlavorPairs_DATAMC_RESP", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_RESP = new TH3D("h3PtFlavorPairs_DATAMC_RESP", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH1D *h1gen_vs_reco_mass = new TH1D("h1gen_vs_reco_mass", ";reco/gen;N", 250, 0, 2.5);

   TH1D *ptrecovsgen = new TH1D("ptrecovsgen", "; reco/gen; N", 100, 0, 3);

   TH1D *x = new TH1D("x", ";x;N", 200, -10, 10);
   TH1D *xx = new TH1D("xx", ";xx;N", 200, -10, 10);
   TH1D *y = new TH1D("y", ";y;N", 200, -10, 10);
   TH1D *yy = new TH1D("yy", ";yy;N", 200, -10, 10);
   TH1D *xy = new TH1D("xy", ";xy;N", 200, -10, 10);

   TH2D* h2_70 = new TH2D("h2_70", ";flav1;flav2; N", 7,-0.5,6.5,7,-0.5,6.5);
   TH2D* h2_97 = new TH2D("h2_97", ";flav1;flav2; N", 7,-0.5,6.5,7,-0.5,6.5);
   TH2D* h2_77_87 = new TH2D("h2_77_87", ";flav1;flav2; N", 7,-0.5,6.5,7,-0.5,6.5);

   TH1D* h_smear_mass = new TH1D("h_smear_mass", ";Mass;N", 200,0,200);
   TH1D* h_smear_pt = new TH1D("h_smear_pt", ";pt diff;N", 67,-1,1);

   TH3D *h3MassFlavorPairs_DATAMC_JER_smear  = new TH3D("h3MassFlavorPairs_DATAMC_JER_smear", ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATAMC_JER_smear = new TH3D("h3PtFlavorPairs_DATAMC_JER_smear", ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

      // Response histograms: pt_reco/pt_gen vs. pt_gen for each jet
   TProfile* prof_resp_jet1 = new TProfile("prof_resp_jet1", ";gen_pt1 (GeV); reco_pt1/gen_pt1", 100, 0, 200);
   TProfile* prof_resp_jet2 = new TProfile("prof_resp_jet2", ";gen_pt2 (GeV); reco_pt2/gen_pt2", 100, 0, 200);
   TH2D*    h2_resp_jet1    = new TH2D("h2_resp_jet1",    ";gen_pt1 (GeV); reco_pt1/gen_pt1", 100, 0, 200, 100, 0.5, 1.5);
   TH2D*    h2_resp_jet2    = new TH2D("h2_resp_jet2",    ";gen_pt2 (GeV); reco_pt2/gen_pt2", 100, 0, 200, 100, 0.5, 1.5);

      // Scaled MC response histograms: pt_reco_scaled/pt_gen vs pt_gen
   TProfile* prof_resp_jet1_scaled = new TProfile(
       "prof_resp_jet1_scaled",
       ";gen_pt1 (GeV); reco_pt1_scaled/gen_pt1",
       100, 0, 200
   );
   TProfile* prof_resp_jet2_scaled = new TProfile(
       "prof_resp_jet2_scaled",
       ";gen_pt2 (GeV); reco_pt2_scaled/gen_pt2",
       100, 0, 200
   );
   TH2D*    h2_resp_jet1_scaled    = new TH2D(
       "h2_resp_jet1_scaled",
       ";gen_pt1 (GeV); reco_pt1_scaled/gen_pt1",
       100, 0, 200, 100, 0.5, 1.5
   );
   TH2D*    h2_resp_jet2_scaled    = new TH2D(
       "h2_resp_jet2_scaled",
       ";gen_pt2 (GeV); reco_pt2_scaled/gen_pt2",
       100, 0, 200, 100, 0.5, 1.5
   );

   // Full-4vector scaled response histograms
   TProfile* prof_resp4v_jet1 = new TProfile(
       "prof_resp4v_jet1",
       ";gen_pt1 (GeV); reco_pt1_full4v/gen_pt1",
       100, 0, 200
   );
   TProfile* prof_resp4v_jet2 = new TProfile(
       "prof_resp4v_jet2",
       ";gen_pt2 (GeV); reco_pt2_full4v/gen_pt2",
       100, 0, 200
   );
   TH2D*    h2_resp4v_jet1    = new TH2D(
       "h2_resp4v_jet1",
       ";gen_pt1 (GeV); reco_pt1_full4v/gen_pt1",
       100, 0, 200, 100, 0.5, 1.5
   );
   TH2D*    h2_resp4v_jet2    = new TH2D(
       "h2_resp4v_jet2",
       ";gen_pt2 (GeV); reco_pt2_full4v/gen_pt2",
       100, 0, 200, 100, 0.5, 1.5
   );



   // Define pt scaling factors
    std::vector<double> ptscales = {1.05, 1.1, 1.2, 0.8, 0.9};
        //std::vector<double> ptscales;
    //for (double ptscale = 0.8; ptscale <= 1.3; ptscale += 0.05) {
        //ptscales.push_back(ptscale);
    //}
    std::vector<TH3D*> h3MassFlavorPairs_DATAMC_JER;
    std::vector<TH3D*> h3PtFlavorPairs_DATAMC_JER;

    std::vector<TH3D*> h3MassFlavorPairs_DATAMC_JER2;
    std::vector<TH3D*> h3PtFlavorPairs_DATAMC_JER2;

    // Initialize histograms for each pt scale factor with specific binning
    for (double ptscale : ptscales) {
        std::string histNameMass2 = Form("h3MassFlavorPairs_DATAMC_JER%.0f", ptscale * 100);
        std::string histNamePt2 = Form("h3PtFlavorPairs_DATAMC_JER%.0f", ptscale * 100);

        // Create the mass histogram for the current pt scale
        h3MassFlavorPairs_DATAMC_JER.push_back(new TH3D(histNameMass2.c_str(), ";MC;DATA; Mass", 
                                                    3, 1, 4, 3, 1, 4, 200, 0, 200));
        // Create the pt histogram for the current pt scale
        h3PtFlavorPairs_DATAMC_JER.push_back(new TH3D(histNamePt2.c_str(), ";MC;DATA; Pt", 
                                                3, 1, 4, 3, 1, 4, 67, -1, 1));
    }

    for (double ptscale : ptscales) {
        std::string histNameMass2 = Form("h3MassFlavorPairs_DATAMC_JER2%.0f", ptscale * 100);
        std::string histNamePt2 = Form("h3PtFlavorPairs_DATAMC_JER2%.0f", ptscale * 100);

        // Create the mass histogram for the current pt scale
        h3MassFlavorPairs_DATAMC_JER2.push_back(new TH3D(histNameMass2.c_str(), ";MC;DATA; Mass", 
                                                    3, 1, 4, 3, 1, 4, 200, 0, 200));
        // Create the pt histogram for the current pt scale
        h3PtFlavorPairs_DATAMC_JER2.push_back(new TH3D(histNamePt2.c_str(), ";MC;DATA; Pt", 
                                                3, 1, 4, 3, 1, 4, 67, -1, 1));
    }

    // Define weight factors
    std::vector<double> weights = {1.2, 1.3, 0.7, 1.6, 0.4, 1/1.3, 1/1.6};
    //std::vector<double> weights;
    //for (double weight = 0.4; weight <= 1.7; weight += 0.01) {
        //weights.push_back(weight);
    //}
    std::vector<TH3D*> h3MassFlavorPairs_DATAMC_ISR;
    std::vector<TH3D*> h3PtFlavorPairs_DATAMC_ISR;

    // Constants for the weighting formula
    double totalEvents = 10044793.0; //10174932.0;
    double gXPairEvents = 2329762.0; //2387032.0;
    double otherPairFactor = (totalEvents - gXPairEvents) / totalEvents;
    double gXPairFactor = gXPairEvents / totalEvents;

    double totalEvents2 = 0.0; //10174932.0;
    double gXPairEvents2 = 0.0; //2387032.0;
    double otherPairFactor2 = (totalEvents2 - gXPairEvents2) / totalEvents2;
    double gXPairFactor2 = gXPairEvents2 / totalEvents2;

    // Initialize histograms for each weight
    for (double weight : weights) {
        std::string histNameMass = Form("h3MassFlavorPairs_DATAMC_ISR%.0f", weight * 100);
        std::string histNamePt = Form("h3PtFlavorPairs_DATAMC_ISR%.0f", weight * 100);

        h3MassFlavorPairs_DATAMC_ISR.push_back(new TH3D(histNameMass.c_str(), ";MC;DATA; Mass", 
                                                    3, 1, 4, 3, 1, 4, 200, 0, 200));
        h3PtFlavorPairs_DATAMC_ISR.push_back(new TH3D(histNamePt.c_str(), ";MC;DATA; Pt", 
                                                3, 1, 4, 3, 1, 4, 67, -1, 1));
    }

   TH2D* h2ctag = new TH2D("h2ctag", ";flavours;ctag; N", 7,-0.5,6.5,100,0,1);


   TH2D* h1cstag = new TH2D("h1cstag", ";flav1;flav2; N",  12,-5.5,6.5,12,-5.5,6.5);
   TH2D* h1udtag = new TH2D("h1udtag", ";flav1;flav2; N",  12,-5.5,6.5,12,-5.5,6.5);
   TH2D* h1xxtag = new TH2D("h1xxtag", ";flav1;flav2; N",  12,-5.5,6.5,12,-5.5,6.5);

   // Apply the helper function for 2D histograms and TProfile2D histograms
   SetBinLabels2D(h2tag_c->GetXaxis(), h2tag_c->GetYaxis(), binLabels3, 3);
   SetBinLabels2D(h2tag_u->GetXaxis(), h2tag_u->GetYaxis(), binLabels3, 3);
   SetBinLabels2D(h2tag_x->GetXaxis(), h2tag_x->GetYaxis(), binLabels3, 3);

   // Apply the helper function for TProfile2D histograms
   SetBinLabels2D(p2matrixMass->GetXaxis(), p2matrixMass->GetYaxis(), binLabels12, 12);
   SetBinLabels2D(p2matrixPt->GetXaxis(), p2matrixPt->GetYaxis(), binLabels12, 12);

   // Apply the helper function for h2matrixN (2D histogram)
   SetBinLabels2D(h2matrixN->GetXaxis(), h2matrixN->GetYaxis(), binLabels12, 12);
   SetBinLabels2D(h1cstag->GetXaxis(), h1cstag->GetYaxis(), binLabels12, 12);
   SetBinLabels2D(h1udtag->GetXaxis(), h1udtag->GetYaxis(), binLabels12, 12);
   SetBinLabels2D(h1xxtag->GetXaxis(), h1xxtag->GetYaxis(), binLabels12, 12);


   // Apply the helper function for 3D histograms
   SetBinLabels3D(h3matrixMass, binLabels12, 12);
   SetBinLabels3D(h3matrixPt, binLabels12, 12);

   SetBinLabels1D(h2ctag->GetXaxis(), binLabels7, 7);
   
   SetBinLabels2D(h2_70->GetXaxis(), h2_70->GetYaxis(), binLabels7, 7);
   SetBinLabels2D(h2_97->GetXaxis(), h2_97->GetYaxis(), binLabels7, 7);
   SetBinLabels2D(h2_77_87->GetXaxis(), h2_77_87->GetYaxis(), binLabels7, 7);   

   std::vector<double> scales;
    for (double scale = 0.95; scale <= 1.05; scale += 0.001) {
        scales.push_back(scale);
    }

    std::vector<double> scaleFactors;
    for (double scale = 0.98; scale <= 1.02; scale += 0.002) {
        scaleFactors.push_back(scale);
    }

   std::vector<std::vector<TH3D*>> h3MassFlavorPairs_DATAMC_Vector1 = initializeHistograms1(scales, "h3MassFlavorPairs_DATAMC_FSR", ";MC;DATA; Mass");
   std::vector<std::vector<TH3D*>> h3PtFlavorPairs_DATAMC_Vector1 = initializeHistograms3(scales, "h3PtFlavorPairs_DATAMC_FSR", ";MC;DATA; Pt");

    // Initialize the histograms
    std::vector<std::vector<TH3D*>> h3MassFlavorPairs_DATAMC_Vector = initializeHistograms(scaleFactors, "h3MassFlavorPairs_DATAMC", ";MC;DATA; Mass");
    std::vector<std::vector<TH3D*>> h3PtFlavorPairs_DATAMC_Vector = initializeHistograms2(scaleFactors, "h3PtFlavorPairs_DATAMC", ";MC;DATA; Pt");

    std::vector<std::vector<TH3D*>> h3MassFlavorPairs_DATAMC_Vector_CHI = initializeHistograms4(scaleFactors, "h3MassFlavorPairs_DATAMC_CHI", ";MC;DATA; Mass");
    std::vector<std::vector<TH3D*>> h3PtFlavorPairs_DATAMC_Vector_CHI = initializeHistograms5(scaleFactors, "h3PtFlavorPairs_DATAMC_CHI", ";MC;DATA; Pt");



   TLorentzVector p4genjet1, p4genjet2, p4recojet1, p4recojet2, p4recojet1_scaled, p4recojet2_scaled, p4recojet1_scaled2, p4recojet2_scaled2, p4recojet_tag, p4recojet_probe;
   TLorentzVector p4genjet1_scaled, p4genjet2_scaled, p4genjet1_scaled2, p4genjet2_scaled2;
   TLorentzVector p4recojet_C, p4recojet_S, p4recojet1_C, p4recojet1_S, p4recojet2_C, p4recojet2_S;
   TLorentzVector p4genjet_C, p4genjet_S, p4genjet1_C, p4genjet1_S, p4genjet2_C, p4genjet2_S;
   TLorentzVector p4recojet_C_JER, p4recojet_S_JER, p4recojet_C_JER2, p4recojet_S_JER2, p4recojet_C_FSR, p4recojet_S_FSR, p4genjet_C_JER2, p4genjet_S_JER2;
   TLorentzVector p4recojet1_FSR, p4recojet2_FSR, p4recojet1_JER, p4recojet2_JER;
   TLorentzVector p4w_scaled, p4w_scaled_RESP, p4_scaled_C, p4_scaled_S;

   // Initialize counters
int count_cs = 0;
int count_cs_true = 0;
int count_cs_false = 0;
int count_sc = 0;
int count_sc_true = 0;
int count_sc_false = 0;

int countAll = 0;
int countXg = 0;
int countOther = 0;

   curdir->cd();
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
   std::cout << std::ctime(&now) << std::endl<< flush;  
   cout << "Processing " << nentries << " events" << endl << flush;
   TStopwatch t;
   t.Start();
   const int nlap = 1000;
   const int nlap2 = 80000;
   //nentries = 100000;

   double totalWeightSum = 0.0;
   double XgWeightSum = 0.0;
   double otherWeightSum = 0.0;

   double totalWeightSum2 = 0.0;
   double XgWeightSum2 = 0.0;
   double otherWeightSum2 = 0.0;

if (SCALED) { //|| CHI2
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        TChain* chain = dynamic_cast<TChain*>(fChain);
        TFile* currentFile = chain->GetFile();

        if (jentry%nlap==0) {
            cout << "." << flush;
        }
        if (jentry%nlap2==0 && jentry!=0) {
            double time = t.RealTime();
        if (time>0) cout << Form("\n\%1.0f ev/s\n",nlap2/time) << flush;
            t.Reset();
            t.Start();
        }

        if (fitProb > 0.2){
            totalWeightSum += weight;

            if (flav1 == 21 || flav2 == 21) {
                XgWeightSum += weight;
            }
            else{
                otherWeightSum += weight;
            }
        }   
    }
    std::cout << "Total Weight Sum: " << totalWeightSum << std::endl;
    std::cout << "Xg Weight Sum: " << XgWeightSum << std::endl;
    std::cout << "Other Weight Sum: " << otherWeightSum << std::endl;
}

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        TChain* chain = dynamic_cast<TChain*>(fChain);
        TFile* currentFile = chain->GetFile();
        jer.update(currentFile->GetName());


        // if (Cut(ientry) < 0) continue;

        if (jentry%nlap==0) {
            cout << "." << flush;
        }
        if (jentry%nlap2==0 && jentry!=0) {
            double time = t.RealTime();
        if (time>0) cout << Form("\n\%1.0f ev/s\n",nlap2/time) << flush;
            t.Reset();
            t.Start();
        }
        p4genjet1.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
        p4genjet2.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);
        p4recojet1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
        p4recojet2.SetPtEtaPhiM(pt2, eta2, phi2, m2);

        double genmass = (p4genjet1 + p4genjet2).M();
        double recomass = (p4recojet1 + p4recojet2).M();
        double reco_pt = (p4recojet1 + p4recojet2).Pt();
        double gen_pt = (p4genjet1 + p4genjet2).Pt(); 

        double recomass2 = 0;

        double reco_dpt2 = (p4recojet1.Pt()-p4recojet2.Pt())/(p4recojet1.Pt()+p4recojet2.Pt());

        int pairIndexDATA = DetermineJetPairDATA(ctag1,ctag2);
        int pairIndexMC_test = DetermineJetPairMC_test(flav1,flav2);
        int pairIndexMC = DetermineJetPairMC(flav1,flav2);

        if (fitProb > 0.2){ //fitProb_scaled
            if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
                ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
                ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry%2 == 0)){
                p4recojet_C.SetPtEtaPhiM(pt1, eta1, phi1, m1);
                p4recojet_S.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                p4genjet_C.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
                p4genjet_S.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);

                    count_cs++;
                        if ((flav1 == 4 || flav1 == -4) && (flav2 == 3 || flav2 == -3)) count_cs_true++;
                        if (( (flav1 == 3 || flav1 == -3) && (flav2 == 4 || flav2 == -4) )) count_cs_false++;
                }
                else {
                p4recojet_S.SetPtEtaPhiM(pt1, eta1, phi1, m1);
                p4recojet_C.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                p4genjet_S.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
                p4genjet_C.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);

                    count_sc++;
                        if ((flav2 == 4 || flav2 == -4) && (flav1 == 3 || flav1 == -3)) count_sc_true++;
                        if (( (flav2 == 3 || flav2 == -3) && (flav1 == 4 || flav1 == -4) )) count_sc_false++;
                }

                double reco_dpt = (p4recojet_S.Pt()-p4recojet_C.Pt())/(p4recojet_S.Pt()+p4recojet_C.Pt());
                double gen_dpt = (p4genjet_S.Pt() - p4genjet_C.Pt()) / (p4genjet_S.Pt() + p4genjet_C.Pt());

                if (SMEAR){
                
                // Apply JER smearing to jets and fill histograms

                    // Smear the pts stochastically
                    double pt1_smeared = jer.smearPt(pt1, eta1, pfRho);
                    double pt2_smeared = jer.smearPt(pt2, eta2, pfRho);

                    // Update the reco jet 4-vectors
                    p4recojet1_scaled2.SetPtEtaPhiM(pt1_smeared, eta1, phi1, m1);
                    p4recojet2_scaled2.SetPtEtaPhiM(pt2_smeared, eta2, phi2, m2);

                    // Set C and S jets based on event topology
                    if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
                        ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
                        ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry % 2 == 0)) {
                        p4recojet_C_JER2 = p4recojet1_scaled2;
                        p4recojet_S_JER2 = p4recojet2_scaled2;
                    } else {
                        p4recojet_C_JER2 = p4recojet2_scaled2;
                        p4recojet_S_JER2 = p4recojet1_scaled2;
                    }

                    // Calculate mass and dpt after smearing
                    double scaledMass_JER = (p4recojet1_scaled2 + p4recojet2_scaled2).M();
                    double scaledpt_JER = (p4recojet_S_JER2.Pt() - p4recojet_C_JER2.Pt()) / (p4recojet_S_JER2.Pt() + p4recojet_C_JER2.Pt());
                    // Fill histograms
                    if (scaledpt_JER <= -0.99 || scaledpt_JER >= 0.99) {
                        std::cout << "Suspicious scaledpt_JER = " << scaledpt_JER
                                << ", pt1_smeared = " << pt1_smeared
                                << ", pt2_smeared = " << pt2_smeared << std::endl;
                    }

                    h_smear_mass->Fill(scaledMass_JER, weight);
                    h_smear_pt->Fill(scaledpt_JER, weight);
                    h3MassFlavorPairs_DATAMC_JER_smear->Fill(pairIndexMC, pairIndexDATA, scaledMass_JER, weight);
                    h3PtFlavorPairs_DATAMC_JER_smear->Fill(pairIndexMC, pairIndexDATA, scaledpt_JER, weight);

                    h3MassFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                    h3PtFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);
                
                }

                if (PSWEIGHTS) {
                    std::vector<std::string> psLabels = {
                        "Baseline", "fsr:murfac=0.707", "fsr:murfac=1.414", "fsr:murfac=0.5",
                        "fsr:murfac=2.0", "fsr:murfac=0.25", "fsr:murfac=4.0", "fsr:g2gg:murfac=0.5",
                        "fsr:g2gg:murfac=2.0", "fsr:g2qq:murfac=0.5", "fsr:g2qq:murfac=2.0",
                        "fsr:q2qg:murfac=0.5", "fsr:q2qg:murfac=2.0", "fsr:x2xg:murfac=0.5",
                        "fsr:x2xg:murfac=2.0", "fsr:g2gg:cns=-2.0", "fsr:g2gg:cns=2.0",
                        "fsr:g2qq:cns=-2.0", "fsr:g2qq:cns=2.0", "fsr:q2qg:cns=-2.0",
                        "fsr:q2qg:cns=2.0", "fsr:x2xg:cns=-2.0", "fsr:x2xg:cns=2.0",
                        "isr:murfac=0.707", "isr:murfac=1.414", "isr:murfac=0.5", "isr:murfac=2.0",
                        "isr:murfac=0.25", "isr:murfac=4.0", "isr:g2gg:murfac=0.5", "isr:g2gg:murfac=2.0",
                        "isr:g2qq:murfac=0.5", "isr:g2qq:murfac=2.0", "isr:q2qg:murfac=0.5",
                        "isr:q2qg:murfac=2.0", "isr:x2xg:murfac=0.5", "isr:x2xg:murfac=2.0",
                        "isr:g2gg:cns=-2.0", "isr:g2gg:cns=2.0", "isr:g2qq:cns=-2.0", "isr:g2qq:cns=2.0",
                        "isr:q2qg:cns=-2.0", "isr:q2qg:cns=2.0", "isr:x2xg:cns=-2.0", "isr:x2xg:cns=2.0"
                    };
                    std::vector<std::pair<std::string, double>> psWeights;
                    for (int i = 0; i < 46; i++) {
                        std::string label = psLabels[i];
                        double value = (*PSWgts)[i + 1] / (*PSWgts)[0];
                        psWeights.push_back(std::make_pair(label, value));
                    }

                    // Fill the nominal histograms (without PS weight variations)
                    h3MassFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                    h3PtFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);

                    // Set the output file as the current directory
                    fout->cd();

                    // Loop over each PS weight scale and fill the corresponding histograms
                    for (const auto &ps : psWeights) {
                        std::string scaleName = ps.first;
                        double scaleFactor = ps.second;

                        // Construct histogram names
                        TString hMassName = TString::Format("h3MassFlavorPairs_DATAMC_%s", scaleName.c_str());
                        TString hPtName = TString::Format("h3PtFlavorPairs_DATAMC_%s", scaleName.c_str());

                        // Retrieve or create histograms
                        TH3D* hMass = dynamic_cast<TH3D*>(fout->Get(hMassName));
                        if (!hMass) {
                            hMass = new TH3D(hMassName, ";MC;DATA; Mass", 3, 1, 4, 3, 1, 4, 200, 0, 200);
                        }

                        TH3D* hPt = dynamic_cast<TH3D*>(fout->Get(hPtName));
                        if (!hPt) {
                            hPt = new TH3D(hPtName, ";MC;DATA; Pt", 3, 1, 4, 3, 1, 4, 67, -1, 1);
                        }

                        // Fill histograms using the PS weight scale
                        hMass->Fill(pairIndexMC, pairIndexDATA, recomass, scaleFactor * weight);
                        hPt->Fill(pairIndexMC, pairIndexDATA, reco_dpt, scaleFactor * weight);
                    }
                }

                if (DATA) {
                    h3MassFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, recomass);
                    h3PtFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, reco_dpt);
                } //DATA

                if (NOTSCALED){
                    for (size_t i = 0; i < scaleFactors.size(); ++i) {
                        for (size_t j = 0; j < scaleFactors.size(); ++j) {
                            double scale_C = scaleFactors[i];
                            double scale_S = scaleFactors[j];
                            // Calculate scaled components separately
                            double scaledPx = p4recojet_C.Px() * scale_C + p4recojet_S.Px() * scale_S;
                            double scaledPy = p4recojet_C.Py() * scale_C + p4recojet_S.Py() * scale_S;
                            double scaledPz = p4recojet_C.Pz() * scale_C + p4recojet_S.Pz() * scale_S;
                            double scaledE  = p4recojet_C.E() * scale_C + p4recojet_S.E() * scale_S;
                            // Set the scaled 4-vector
                            p4w_scaled.SetPxPyPzE(scaledPx, scaledPy, scaledPz, scaledE);
                            // The scaled mass
                            double scaledMass = p4w_scaled.M();
                            // The scaled delta pt
                            double scaledpt = (scale_S * p4recojet_S.Pt()-scale_C * p4recojet_C.Pt())/(scale_S * p4recojet_S.Pt()+scale_C * p4recojet_C.Pt());

                            if (pairIndexMC == 1){ //flav == 3 or flav == 4
                                h3MassFlavorPairs_DATAMC_Vector[i][j]->Fill(pairIndexMC, pairIndexDATA, scaledMass, weight);
                                h3PtFlavorPairs_DATAMC_Vector[i][j]->Fill(pairIndexMC, pairIndexDATA, scaledpt, weight);
                            }
                            else {
                                h3MassFlavorPairs_DATAMC_Vector[i][j]->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                                h3PtFlavorPairs_DATAMC_Vector[i][j]->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);
                            }
                        }
                    }
                    h3MassFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                    h3PtFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);

                    h2ctag->Fill(abs(min(flav1,6)),ctag1);
                    h2ctag->Fill(abs(min(flav2,6)),ctag2);
                    h2matrixN->Fill(min(flav1,6), min(flav2,6));
                    p2matrixMass->Fill(min(flav1,6), min(flav2,6),recomass,weight);
                    p2matrixPt->Fill(min(flav1,6), min(flav2,6),reco_dpt2,weight);
                    h3matrixMass->Fill(min(flav1,6), min(flav2,6),recomass,weight);
                    h3matrixPt->Fill(min(flav1,6), min(flav2,6),reco_dpt2,weight);

                    // Alternating jets for tag and probe based on event number
                    double tag, probe, flav_tag, ctag_tag, flav_probe, ctag_probe;

                    if (jentry % 2 == 0) {
                        // Even events: jet1 is tag, jet2 is probe
                        p4recojet_tag = p4recojet1;
                        p4recojet_probe = p4recojet2;

                        tag = Tag(ctag1);  // Use ctag1 for the tag jet
                        flav_tag = Flav(flav1);
                        ctag_tag = ctag1;

                        probe = Tag(ctag2); // Use ctag2 for the probe jet
                        flav_probe = Flav(flav2);
                        ctag_probe = ctag2;
                    } else {
                        // Odd events: jet2 is tag, jet1 is probe
                        p4recojet_tag = p4recojet2;
                        p4recojet_probe = p4recojet1;

                        tag = Tag(ctag2);  // Use ctag2 for the tag jet
                        flav_tag = Flav(flav2);
                        ctag_tag = ctag2;

                        probe = Tag(ctag1); // Use ctag1 for the probe jet
                        flav_probe = Flav(flav1);
                        ctag_probe = ctag1;
                    }

                    // Perform your tagging logic for the tag jet
                    if (ctag_tag > 0.43) {
                        h2tag_c->Fill(probe, flav_probe);  // c-tagged
                    } 
                    if (ctag_tag < 0.06) {
                        h2tag_u->Fill(probe, flav_probe);  // u-tagged
                    } 
                    if (ctag_tag < 0.43 && ctag_tag > 0.06) {
                        h2tag_x->Fill(probe, flav_probe);  // intermediate region
                    }

                    if (pairIndexDATA == 1) h1cstag->Fill(min(flav1,6), min(flav2,6));
                    if (pairIndexDATA == 2) h1udtag->Fill(min(flav1,6), min(flav2,6));
                    if (pairIndexDATA == 3) h1xxtag->Fill(min(flav1,6), min(flav2,6));
                }
        
        
            if (SCALED){
                double FSR_UP = (*PSWgts)[5] / (*PSWgts)[0]; //2
                double FSR_UP2 = (*PSWgts)[7] / (*PSWgts)[0]; //4
                double FSR_DOWN = (*PSWgts)[4] / (*PSWgts)[0]; //0.5
                double FSR_DOWN2 = (*PSWgts)[6] / (*PSWgts)[0]; //0.25
                
                double ISR_UP = (*PSWgts)[27] / (*PSWgts)[0]; //2
                double ISR_DOWN = (*PSWgts)[26] / (*PSWgts)[0]; //0.5
                if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
                ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
                ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry%2 == 0)){
                    p4recojet_C.SetPtEtaPhiM(pt1, eta1, phi1, m1);
                    p4recojet_S.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                    p4genjet_C.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
                    p4genjet_S.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);
                }
                else {
                    p4recojet_S.SetPtEtaPhiM(pt1, eta1, phi1, m1);
                    p4recojet_C.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                    p4genjet_S.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
                    p4genjet_C.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);
                }

                double reco_dpt = (p4recojet_S.Pt()-p4recojet_C.Pt())/(p4recojet_S.Pt()+p4recojet_C.Pt());
                double gen_dpt = (p4genjet_S.Pt() - p4genjet_C.Pt()) / (p4genjet_S.Pt() + p4genjet_C.Pt());

                h1gen_vs_reco_mass->Fill(recomass/genmass, weight);
                double scaledPx_RESP = p4recojet_C.Px() * 0.994 + p4recojet_S.Px() * 1.002;
                double scaledPy_RESP = p4recojet_C.Py() * 0.994 + p4recojet_S.Py() * 1.002;
                double scaledPz_RESP = p4recojet_C.Pz() * 0.994 + p4recojet_S.Pz() * 1.002;
                double scaledE_RESP  = p4recojet_C.E() * 0.994 + p4recojet_S.E() * 1.002;
                p4w_scaled_RESP.SetPxPyPzE(scaledPx_RESP, scaledPy_RESP, scaledPz_RESP, scaledE_RESP);
                double scaledMass_RESP = p4w_scaled_RESP.M();
                double scaledpt_RESP = (1.002 * p4recojet_S.Pt() - 0.994 * p4recojet_C.Pt())/(1.002 * p4recojet_S.Pt() + 0.994 * p4recojet_C.Pt());
                if (pairIndexMC == 1){
                    h3MassFlavorPairs_DATAMC_RESP->Fill(pairIndexMC, pairIndexDATA, scaledMass_RESP, weight);
                    h3PtFlavorPairs_DATAMC_RESP->Fill(pairIndexMC, pairIndexDATA, scaledpt_RESP, weight);
                }
                else {
                    h3MassFlavorPairs_DATAMC_RESP->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                    h3PtFlavorPairs_DATAMC_RESP->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);
                }

                h3MassFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, recomass, weight);
                h3PtFlavorPairs_DATAMC->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weight);

                // ISR
                countAll++;
                if (flav1 == 21 || flav2 == 21) {
                    countXg++;
                    // For gX pairs, use each weight directly
                    for (size_t i = 0; i < weights.size(); ++i) {
                        double weightXg = weights[i];
                        h3MassFlavorPairs_DATAMC_ISR[i]->Fill(pairIndexMC, pairIndexDATA, recomass, weightXg * weight);
                        h3PtFlavorPairs_DATAMC_ISR[i]->Fill(pairIndexMC, pairIndexDATA, reco_dpt, weightXg * weight);
                    }
                } else {
                    countOther++;
                    // For other pairs, calculate the adjusted weight
                    for (size_t i = 0; i < weights.size(); ++i) {
                        double weightXg = weights[i];
                        double otherWeight = (1 - (XgWeightSum/totalWeightSum * weightXg)) / (otherWeightSum/totalWeightSum);
                        //(1 - (gXPairFactor * weightXg)) / otherPairFactor;
                        h3MassFlavorPairs_DATAMC_ISR[i]->Fill(pairIndexMC, pairIndexDATA, recomass, otherWeight * weight);
                        h3PtFlavorPairs_DATAMC_ISR[i]->Fill(pairIndexMC, pairIndexDATA, reco_dpt, otherWeight * weight);
                    }
                }

                // Loop over pt scaling factors to calculate and fill histograms
                for (size_t i = 0; i < ptscales.size(); ++i) {  
                    double ptscale = ptscales[i];

                    // Jet-1
                    double f1   = ScaleJet1Func(p4genjet1.Pt());
                    double pt1_scaled = ((p4recojet1.Pt() - p4genjet1.Pt()) * ptscale + p4genjet1.Pt());
                    // Jet-2
                    double f2   = ScaleJet2Func(p4genjet2.Pt());
                    double pt2_scaled = ((p4recojet2.Pt() - p4genjet2.Pt()) * ptscale + p4genjet2.Pt());

                    p4recojet1_scaled.SetPtEtaPhiM(pt1_scaled, eta1, phi1, m1);
                    p4recojet2_scaled.SetPtEtaPhiM(pt2_scaled, eta2, phi2, m2);

                    p4recojet1_scaled2 = ((p4recojet1 - p4genjet1) * ptscale + p4genjet1);
                    p4recojet2_scaled2 = ((p4recojet2 - p4genjet2) * ptscale + p4genjet2);

                    if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
                        ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
                        ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry%2 == 0)){
                        p4recojet_C_JER = p4recojet1_scaled;
                        p4recojet_S_JER = p4recojet2_scaled;

                        p4recojet_C_JER2 = p4recojet1_scaled2;
                        p4recojet_S_JER2 = p4recojet2_scaled2;
                    }
                    else {
                        p4recojet_C_JER = p4recojet2_scaled;
                        p4recojet_S_JER = p4recojet1_scaled;

                        p4recojet_C_JER2 = p4recojet2_scaled2;
                        p4recojet_S_JER2 = p4recojet1_scaled2;
                    }

                    // Calculate the scaled mass and pt difference ratio
                    double scaledMass_JER2 = (p4recojet1_scaled2 + p4recojet2_scaled2).M();
                    double scaledpt_JER2 = (p4recojet_S_JER2.Pt() - p4recojet_C_JER2.Pt()) / (p4recojet_S_JER2.Pt() + p4recojet_C_JER2.Pt());

                    double scaledMass_JER = (p4recojet1_scaled + p4recojet2_scaled).M();
                    double scaledpt_JER = (p4recojet_S_JER.Pt() - p4recojet_C_JER.Pt()) / (p4recojet_S_JER.Pt() + p4recojet_C_JER.Pt());

                    h3MassFlavorPairs_DATAMC_JER[i]->Fill(pairIndexMC, pairIndexDATA, scaledMass_JER, weight);
                    h3PtFlavorPairs_DATAMC_JER[i]->Fill(pairIndexMC, pairIndexDATA, scaledpt_JER, weight);

                    h3MassFlavorPairs_DATAMC_JER2[i]->Fill(pairIndexMC, pairIndexDATA, scaledMass_JER2, weight);
                    h3PtFlavorPairs_DATAMC_JER2[i]->Fill(pairIndexMC, pairIndexDATA, scaledpt_JER2, weight);
                    
                } //JER
                // Fill scaled response: reco_pt_scaled / gen_pt for each jet

                if (p4genjet1.Pt() > 0){

                    prof_resp_jet1->Fill(p4genjet1.Pt(), p4recojet1.Pt()/p4genjet1.Pt());
                    h2_resp_jet1->Fill(p4genjet1.Pt(), p4recojet1.Pt()/p4genjet1.Pt());

                    prof_resp_jet1_scaled->Fill(p4genjet1.Pt(), ((p4recojet1.Pt() - p4genjet1.Pt()) * 1.10 + p4genjet1.Pt())/p4genjet1.Pt());
                    h2_resp_jet1_scaled->Fill(p4genjet1.Pt(), ((p4recojet1.Pt() - p4genjet1.Pt()) * 1.10 + p4genjet1.Pt())/p4genjet1.Pt());
                }
                if (p4genjet2.Pt() > 0){
                    prof_resp_jet2->Fill(p4genjet2.Pt(), p4recojet2.Pt()/p4genjet2.Pt());
                    h2_resp_jet2->Fill(p4genjet2.Pt(), p4recojet2.Pt()/p4genjet2.Pt());

                    prof_resp_jet2_scaled->Fill(p4genjet2.Pt(), ((p4recojet2.Pt() - p4genjet2.Pt()) * 1.10 + p4genjet2.Pt())/p4genjet2.Pt());
                    h2_resp_jet2_scaled->Fill(p4genjet2.Pt(), ((p4recojet2.Pt() - p4genjet2.Pt()) * 1.10 + p4genjet2.Pt())/p4genjet2.Pt());
                }     
                // FSR
                for (size_t i = 0; i < scales.size(); ++i) {
                    double scale_FSR = scales[i];
                    double scale = scales[i];
                    p4recojet1_FSR = scale * p4recojet1;
                    p4recojet2_FSR = scale * p4recojet2;
                    p4recojet_C_FSR = scale * p4recojet_C;
                    p4recojet_S_FSR = scale * p4recojet_S;
                    double recomass_FSR = (p4recojet1_FSR + p4recojet2_FSR).M();
                    double reco_dpt_FSR = (p4recojet_S_FSR.Pt()-p4recojet_C_FSR.Pt())/(p4recojet_S_FSR.Pt()+p4recojet_C_FSR.Pt());
                    h3MassFlavorPairs_DATAMC_Vector1[i][0]->Fill(pairIndexMC, pairIndexDATA, recomass_FSR, weight);  
                    h3PtFlavorPairs_DATAMC_Vector1[i][0]->Fill(pairIndexMC, pairIndexDATA, reco_dpt_FSR, weight);
                } // FSR
            } //SCALED
        } //fitprob
    } // jentry

   // Write full-4vector response histograms
   fout->Write();
   fout->Close();
   exit(0);
}
