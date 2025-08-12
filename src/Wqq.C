#define Wqq_cxx
#include "../interface/Wqq.h"
#include <TH2.h>
#include <TH3.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>

#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <TLegend.h>
#include <TColor.h>
#include <TStopwatch.h>
#include "../minitools/tdrstyle_mod22.C"

int CTAG(double x, double y) {
   if (x > 0.43 || y > 0.43) return 1;
   if (x < 0.06 && y < 0.06) return 2;
   return 3;
}

int FLAV(int x, int y) {
   if (x == 3 && y == -4 || x == -3 && y == 4 ||
   x == 4 && y == -3 || x == -4 && y == 3) return 1;
   if (x == 1 && y == -2 || x == -1 && y == 2 ||
   x == 2 && y == -1 || x == -2 && y == 1) return 2;
   return 3;
}


// Function to initialize 3D histograms with mass binning
std::vector<std::vector<std::vector<TH3D*>>> initializeHistograms(const std::vector<double>& scaleFactors,
                                                                    const std::string& namePrefix,
                                                                    const std::string& axisTitle) {
    std::vector<std::vector<std::vector<TH3D*>>> histograms;
    // Loop over first scale factor (scale1)
    for (double scale1 : scaleFactors) {
        std::vector<std::vector<TH3D*>> vec2D;
        // Loop over second scale factor (scale2)
        for (double scale2 : scaleFactors) {
            std::vector<TH3D*> innerVec;
            // Loop over third scale factor (scale3)
            for (double scale3 : scaleFactors) {
                TString histName = TString::Format("%s_%.3f_%.3f_%.3f", namePrefix.c_str(), scale1, scale2, scale3);
                TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 200, 0, 200);
                innerVec.push_back(hist);
            }
            vec2D.push_back(innerVec);
        }
        histograms.push_back(vec2D);
    }
    return histograms;
}

// Function to initialize 3D histograms with Pt binning
std::vector<std::vector<std::vector<TH3D*>>> initializeHistograms2(const std::vector<double>& scaleFactors,
                                                                       const std::string& namePrefix,
                                                                       const std::string& axisTitle) {
    std::vector<std::vector<std::vector<TH3D*>>> histograms;
    // Loop over first scale factor (scale1)
    for (double scale1 : scaleFactors) {
        std::vector<std::vector<TH3D*>> vec2D;
        // Loop over second scale factor (scale2)
        for (double scale2 : scaleFactors) {
            std::vector<TH3D*> innerVec;
            // Loop over third scale factor (scale3)
            for (double scale3 : scaleFactors) {
                TString histName = TString::Format("%s_%.3f_%.3f_%.3f", namePrefix.c_str(), scale1, scale2, scale3);
                TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 67, -1, 1);
                innerVec.push_back(hist);
            }
            vec2D.push_back(innerVec);
        }
        histograms.push_back(vec2D);
    }
    return histograms;
}

std::vector<std::vector<TH3D*>> initializeHistograms3(const std::vector<double>& scaleFactors,
                                                        const std::string& namePrefix,
                                                        const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;
    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> row;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 200, 0, 200);
            row.push_back(hist);
        }
        histograms.push_back(row);
    }
    return histograms;
}

std::vector<std::vector<TH3D*>> initializeHistograms4(const std::vector<double>& scaleFactors,
                                                        const std::string& namePrefix,
                                                        const std::string& axisTitle) {
    std::vector<std::vector<TH3D*>> histograms;
    for (double scale1 : scaleFactors) {
        std::vector<TH3D*> row;
        for (double scale2 : scaleFactors) {
            TString histName = TString::Format("%s_%.3f_%.3f", namePrefix.c_str(), scale1, scale2);
            TH3D* hist = new TH3D(histName, axisTitle.c_str(), 3, 1, 4, 3, 1, 4, 67, -1, 1);
            row.push_back(hist);
        }
        histograms.push_back(row);
        }
    return histograms;
}

// Function to determine the category bin (1 to 7) based on flav1, flav2
int GetFlavorPairBin(int flav1, int flav2) {
    // Clamp flavors to 0–6
    flav1 = std::min(std::abs(flav1), 6);
    flav2 = std::min(std::abs(flav2), 6);
    
    // Check flavor pairs (symmetric: (i,j) or (j,i))
    if ((flav1 == 4 && flav2 == 3) || (flav1 == 3 && flav2 == 4)) return 0; // cs
    if ((flav1 == 2 && flav2 == 1) || (flav1 == 1 && flav2 == 2)) return 1; // ud
    if ((flav1 == 4 && flav2 == 1) || (flav1 == 1 && flav2 == 4)) return 2; // cd
    if ((flav1 == 2 && flav2 == 3) || (flav1 == 3 && flav2 == 2)) return 3; // us
    if ((flav1 == 4 && flav2 == 5) || (flav1 == 5 && flav2 == 4)) return 4; // cb
    if ((flav1 == 2 && flav2 == 5) || (flav1 == 5 && flav2 == 2)) return 5; // ub
    return 6; // x (all other pairs)
}

void Wqq::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Wqq.C
//      root> Wqq t
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

   TDirectory *curdir = gDirectory;
   // Create the output file based on a condition
   TFile* fout;

   bool DATA = false;
   bool SCALED_MC = false;
   bool MC_fractions = false;
   bool MC = false;
   if (DATA) {fout = new TFile("output_DATARun2_Wqq_test4.root", "RECREATE");}
   else if (SCALED_MC) {fout = new TFile("output_MCRun2_Wqq_test3.root", "RECREATE");}
   else if (MC_fractions) {fout = new TFile("output_MCRun2_Wqq_fractions.root", "RECREATE");}
   else if (MC) {fout = new TFile("output_MCRun2_Wqq.root", "RECREATE");}
   else {fout = new TFile("false.root", "RECREATE");}
   fout->cd();

   TH1::SetDefaultSumw2();

   TH3D *h3MassFlavorPairs_DATA  = new TH3D("h3MassFlavorPairs_DATA", ";TrueFlav;TagFlav; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_DATA = new TH3D("h3PtFlavorPairs_DATA", ";TrueFlav;TagFlav; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   TH3D *h3MassFlavorPairs_MC  = new TH3D("h3MassFlavorPairs_MC", ";TrueFlav;TagFlav; Mass", 3, 1, 4, 3, 1, 4,200,0,200);
   TH3D *h3PtFlavorPairs_MC = new TH3D("h3PtFlavorPairs_MC", ";TrueFlav;TagFlav; Pt", 3, 1, 4, 3, 1, 4,67,-1,1);

   // Response histograms: pt_reco/pt_gen vs. pt_gen for each jet
   TProfile* prof_resp_jet1 = new TProfile("prof_resp_jet1", ";gen_pt1 (GeV); reco_pt1/gen_pt1", 100, 0, 200);
   TProfile* prof_resp_jet2 = new TProfile("prof_resp_jet2", ";gen_pt2 (GeV); reco_pt2/gen_pt2", 100, 0, 200);
   TH2D*    h2_resp_jet1    = new TH2D("h2_resp_jet1",    ";gen_pt1 (GeV); reco_pt1/gen_pt1", 100, 0, 200, 100, 0.5, 1.5);
   TH2D*    h2_resp_jet2    = new TH2D("h2_resp_jet2",    ";gen_pt2 (GeV); reco_pt2/gen_pt2", 100, 0, 200, 100, 0.5, 1.5);


   TH1D *h_cs = new TH1D("h_cs", ";Flav; N",7,0,7);
   TH1D *h_ud = new TH1D("h_ud", ";Flav; N",7,0,7);
   TH1D *h_xx = new TH1D("h_xx", ";Flav; N",7,0,7);
   TH1D *h_all = new TH1D("h_all", ";Flav; N",7,0,7);

    std::vector<double> scaleFactors;
    for (double scale = 0.98; scale <= 1.02; scale += 0.002) { 
        scaleFactors.push_back(scale);
    }

    std::vector<std::vector<std::vector<TH3D*>>>  h3MassFlavorPairs_MC_Vector = initializeHistograms(scaleFactors, "h3MassFlavorPairs_MC", ";MC;DATA; Mass");
    std::vector<std::vector<std::vector<TH3D*>>>  h3PtFlavorPairs_MC_Vector = initializeHistograms2(scaleFactors, "h3PtFlavorPairs_MC", ";MC;DATA; Pt");

    std::vector<std::vector<TH3D*>>  h3MassFlavorPairs_MC_Vector_CS = initializeHistograms3(scaleFactors, "h3MassFlavorPairs_MC", ";MC;DATA; Mass");
    std::vector<std::vector<TH3D*>>  h3PtFlavorPairs_MC_Vector_CS = initializeHistograms4(scaleFactors, "h3PtFlavorPairs_MC", ";MC;DATA; Pt");
    
    TLorentzVector p4genjet1, p4genjet2, p4recojet1, p4recojet2, p4recojet_C, p4recojet_S, p4w_scaled, p4w_reco;

   //curdir->cd();
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
   std::cout << std::ctime(&now) << std::endl<< flush;  
   cout << "Processing " << nentries << " events" << endl << flush;
   TStopwatch t;
   t.Start();
   const int nlap = 1000;
   const int nlap2 = 80000;
   //nentries = 100;

for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      TChain* chain = dynamic_cast<TChain*>(fChain);
      TFile* currentFile = chain->GetFile();

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

      p4recojet1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
      p4recojet2.SetPtEtaPhiM(pt2, eta2, phi2, m2);

      p4genjet1.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
      p4genjet2.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);

      double Px = p4recojet1.Px() + p4recojet2.Px();
      double Py = p4recojet1.Py() + p4recojet2.Py();
      double Pz = p4recojet1.Pz() + p4recojet2.Pz();
      double E  = p4recojet1.E() + p4recojet2.E();

      p4w_reco.SetPxPyPzE(Px, Py, Pz, E);

      double recoMass = p4w_reco.M();

      int trueFlav = FLAV(flav1,flav2);
      int tagFlav = CTAG(ctag1,ctag2);

      bool isUD1 = (abs(flav1) == 1 || abs(flav1) == 2);
      bool isUD2 = (abs(flav2) == 1 || abs(flav2) == 2);
      bool isS1 = (abs(flav1) == 3);
      bool isS2 = (abs(flav2) == 3);
      bool isC1 = (abs(flav1) == 4);
      bool isC2 = (abs(flav2) == 4);

      if (fitProb > 0.2){
        if (SCALED_MC){
            // c and s
            for (size_t i = 0; i < scaleFactors.size(); ++i) {
                for (size_t j = 0; j < scaleFactors.size(); ++j) {
                    for (size_t k = 0; k < scaleFactors.size(); ++k) {
                        double scale_C = scaleFactors[i];
                        double scale_S = scaleFactors[j];
                        double scale_UD = scaleFactors[k];

                        double scale1 = isUD1 ? scale_UD : (isS1 ? scale_S : (isC1 ? scale_C : 1.0)); // 1.0 -> ud
                        double scale2 = isUD2 ? scale_UD : (isS2 ? scale_S : (isC2 ? scale_C : 1.0));

                        if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
                                ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
                                ((ctag1 <= 0.43 && ctag2 <= 0.43) && (jentry % 2 == 0))) {
                                p4recojet_C.SetPxPyPzE(p4recojet1.Px() * scale1, p4recojet1.Py() * scale1, p4recojet1.Pz() * scale1, p4recojet1.E() * scale1);

                                p4recojet_S.SetPxPyPzE(p4recojet2.Px() * scale2, p4recojet2.Py() * scale2, p4recojet2.Pz() * scale2, p4recojet2.E() * scale2);
                            } else {
                                p4recojet_S.SetPxPyPzE(p4recojet1.Px() * scale1, p4recojet1.Py() * scale1, p4recojet1.Pz() * scale1, p4recojet1.E() * scale1);

                                p4recojet_C.SetPxPyPzE(p4recojet2.Px() * scale2, p4recojet2.Py() * scale2, p4recojet2.Pz() * scale2, p4recojet2.E() * scale2);
                            }
                    
                        double scaledPx = p4recojet1.Px() * scale1 + p4recojet2.Px() * scale2;
                        double scaledPy = p4recojet1.Py() * scale1 + p4recojet2.Py() * scale2;
                        double scaledPz = p4recojet1.Pz() * scale1 + p4recojet2.Pz() * scale2;
                        double scaledE  = p4recojet1.E() * scale1 + p4recojet2.E() * scale2;
                        // Set the scaled 4-vector
                        p4w_scaled.SetPxPyPzE(scaledPx, scaledPy, scaledPz, scaledE);
                        // The scaled mass
                        double scaledMass = p4w_scaled.M();
                        double scaleddPt = (p4recojet_S.Pt() - p4recojet_C.Pt()) / (p4recojet_S.Pt() + p4recojet_C.Pt());
                        
                        // Fill the histograms.
                        h3MassFlavorPairs_MC_Vector[i][j][k]->Fill(trueFlav, tagFlav, scaledMass, weight);
                        h3PtFlavorPairs_MC_Vector[i][j][k]->Fill(trueFlav, tagFlav, scaleddPt, weight);
                    }
                }
            }
        } // MC
        if (MC_fractions) {
            int bin = GetFlavorPairBin(flav1, flav2);
            h_all->Fill(bin, weight);
            if (tagFlav == 1) h_cs->Fill(bin, weight);
            if (tagFlav == 2) h_ud->Fill(bin, weight);
            if (tagFlav == 3) h_xx->Fill(bin, weight);
        }
        if (MC) {
            if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
            ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
            ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry%2 == 0)){
                p4recojet_C.SetPxPyPzE(p4recojet1.Px(), p4recojet1.Py(), p4recojet1.Pz(), p4recojet1.E());
                p4recojet_S.SetPxPyPzE(p4recojet2.Px(), p4recojet2.Py(), p4recojet2.Pz(), p4recojet2.E());
            }
            else {
                p4recojet_S.SetPxPyPzE(p4recojet1.Px(), p4recojet1.Py(), p4recojet1.Pz(), p4recojet1.E());

                p4recojet_C.SetPxPyPzE(p4recojet2.Px(), p4recojet2.Py(), p4recojet2.Pz(), p4recojet2.E());
            }

            double recodPt = (p4recojet_S.Pt()-p4recojet_C.Pt())/(p4recojet_S.Pt()+p4recojet_C.Pt());

            h3MassFlavorPairs_MC->Fill(trueFlav, tagFlav, recoMass);
            h3PtFlavorPairs_MC->Fill(trueFlav, tagFlav, recodPt);
        }
        if (DATA) {
            if ((ctag1 > 0.43 && ctag2 <= 0.43) ||
               ((ctag1 > 0.43 && ctag2 > 0.43) && ctag1 > ctag2) ||
               ((ctag1 <= 0.43 && ctag2 <= 0.43) && jentry%2 == 0)){
                p4recojet_C.SetPxPyPzE(p4recojet1.Px(), p4recojet1.Py(), p4recojet1.Pz(), p4recojet1.E());
                p4recojet_S.SetPxPyPzE(p4recojet2.Px(), p4recojet2.Py(), p4recojet2.Pz(), p4recojet2.E());
            }
            else {
                p4recojet_S.SetPxPyPzE(p4recojet1.Px(), p4recojet1.Py(), p4recojet1.Pz(), p4recojet1.E());

                p4recojet_C.SetPxPyPzE(p4recojet2.Px(), p4recojet2.Py(), p4recojet2.Pz(), p4recojet2.E());
            }
            double recodPt = (p4recojet_S.Pt()-p4recojet_C.Pt())/(p4recojet_S.Pt()+p4recojet_C.Pt());
            h3MassFlavorPairs_DATA->Fill(trueFlav, tagFlav, recoMass);
            h3PtFlavorPairs_DATA->Fill(trueFlav, tagFlav, recodPt);
        } // DATA
      } // fitProb
   } // jentry

   // Write response histograms
   prof_resp_jet1->Write();
   prof_resp_jet2->Write();
   h2_resp_jet1->Write();
   h2_resp_jet2->Write();

   std::cout << "\n[Wqq] Output written to: " << fout->GetName() << std::endl;
   fout->Write();
   fout->Close();
   curdir->cd();
   exit(0);
}
