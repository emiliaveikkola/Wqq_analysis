#define StrangeJet_cxx
#include "StrangeJet.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>

#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <TLegend.h>
#include <TColor.h>
#include <TStopwatch.h>
#include "tdrstyle_mod22.C"

bool    debug = false;


void StrangeJet::Loop(){
  //   In a ROOT session, you can do:
  //      root> .L StrangeJet.C
  //      root> StrangeJet t
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

  fChain->SetBranchStatus("*",0);  // disable all branches

  std::vector<std::string> branchNames = {"genWeight",
  "GenJet_partonFlavour", "GenJetGenPartCand_genJetIdx", "GenJetGenPartCand_GenPartCandIdx","GenPartCand_pdgId",
  "nGenJet", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_mass",
  "nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass",
  "Jet_rawFactor",
  "Jet_btagDeepFlavCvL", "Jet_btagDeepFlavCvB","Jet_btagDeepFlavB",
  "Jet_btagPNetCvL", "Jet_btagPNetCvB", "Jet_btagPNetB",
  "Jet_partonFlavour","Jet_hadronFlavour",
  "GenPartCand_pt", "GenPartCand_eta", "GenPartCand_phi",  "GenPartCand_mass", "GenPartCand_charge"};

  for (const auto& branchName : branchNames) {
      if (fChain->GetBranch(branchName.c_str())) {
          fChain->SetBranchStatus(branchName.c_str(), 1);
      }
  }
  
  /*
  fChain->SetBranchStatus("genWeight",1);  // activate branchname
  
  fChain->SetBranchStatus("GenJet_partonFlavour",1);  // activate branchname
  fChain->SetBranchStatus("nGenJetGenPartCand",1);  // activate branchname
  fChain->SetBranchStatus("GenJetGenPartCand_genJetIdx",1);  // activate branchname
  fChain->SetBranchStatus("GenJetGenPartCand_GenPartCandIdx",1);  // activate branchname
  fChain->SetBranchStatus("GenPartCand_pdgId",1);  // activate branchname
  
  fChain->SetBranchStatus("nGenJet",1);  // activate branchname
  fChain->SetBranchStatus("GenJet_pt",1);  // activate branchname
  fChain->SetBranchStatus("GenJet_eta",1);  // activate branchname
  fChain->SetBranchStatus("GenJet_phi",1);  // activate branchname
  fChain->SetBranchStatus("GenJet_mass",1);  // activate branchname
  
  fChain->SetBranchStatus("GenPartCand_pt",1);  // activate branchname
  fChain->SetBranchStatus("GenPartCand_eta",1);
  fChain->SetBranchStatus("GenPartCand_phi",1);  // activate branchname
  fChain->SetBranchStatus("GenPartCand_mass",1);  // activate branchname
  fChain->SetBranchStatus("GenPartCand_charge",1);  // activate branchname
  
  */
  
  
  cout << "Hello World!" << endl;
  
  
  // Bins
  Double_t xbins[] = {0.11633, 0.17847, 0.24341, 0.31129, 0.38224, 0.45639, 0.53389, 0.61489, 0.69955, 0.78804, 
		      0.88053, 0.97719, 1.07823, 1.18382, 1.29419, 1.40955, 1.53011, 1.65613, 1.78784, 1.92549, 
		      2.06937, 2.21975, 2.37693, 2.5412, 2.7129, 2.89235, 3.07992, 3.27595, 3.48085, 3.695, 
		      3.91883, 4.15277, 4.39728, 4.65284, 4.91995, 5.19912, 5.49091, 5.79588, 6.11463, 6.44778, 
		      6.79598, 7.15992, 7.5403, 7.93787, 8.3534, 8.7877, 9.24163, 9.71606, 10.21193, 10.73021, 
		      11.2719, 11.83807, 12.42982, 13.0483, 13.69473, 14.37037, 15.07653, 15.8146, 16.58602, 
		      17.39229, 18.23499, 19.11576, 20.03633, 20.99849, 22.00413, 23.0552, 24.15376, 25.30195, 
		      26.50203, 27.75632, 29.06729, 30.43749, 31.86959, 33.36641, 34.93085, 36.56597, 38.27498, 
		      40.0612, 41.92812, 43.8794, 45.91884, 48.05043, 50.27832, 52.60688, 55.04064, 57.58437, 
		      60.24303, 63.02181, 65.92614, 68.96169, 72.1344, 75.45046, 78.91634, 82.53882, 86.32497, 
		      90.28218, 94.41819, 98.74107, 100}; // Adjust the bin edges as needed
  int nxbins = sizeof(xbins) / sizeof(xbins[0]) - 1;
  
  Double_t xbins_max[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22.0, 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23.0, 23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9, 24.0, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6, 24.7, 24.8, 24.9, 25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9, 29.0, 29.1, 29.2, 29.3, 29.4, 29.5, 29.6, 29.7, 29.8, 29.9, 30.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 30.8, 30.9, 31.0, 31.1, 31.2, 31.3, 31.4, 31.5, 31.6, 31.7, 31.8, 31.9, 32.0, 32.1, 32.2, 32.3, 32.4, 32.5, 32.6, 32.7, 32.8, 32.9, 33.0, 33.1, 33.2, 33.3, 33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34.0, 34.1, 34.2, 34.3, 34.4, 34.5, 34.6, 34.7, 34.8, 34.9, 35.0, 35.1, 35.2, 35.3, 35.4, 35.5, 35.6, 35.7, 35.8, 35.9, 36.0, 36.1, 36.2, 36.3, 36.4, 36.5, 36.6, 36.7, 36.8, 36.9, 37.0, 37.1, 37.2, 37.3, 37.4, 37.5, 37.6, 37.7, 37.8, 37.9, 38.0, 38.1, 38.2, 38.3, 38.4, 38.5, 38.6, 38.7, 38.8, 38.9, 39.0, 39.1, 39.2, 39.3, 39.4, 39.5, 39.6, 39.7, 39.8, 39.9, 40.0, 40.1, 40.2, 40.3, 40.4, 40.5, 40.6, 40.7, 40.8, 40.9, 41.0, 41.1, 41.2, 41.3, 41.4, 41.5, 41.6, 41.7, 41.8, 41.9, 42.0, 42.1, 42.2, 42.3, 42.4, 42.5, 42.6, 42.7, 42.8, 42.9, 43.0, 43.1, 43.2, 43.3, 43.4, 43.5, 43.6, 43.7, 43.8, 43.9, 44.0, 44.1, 44.2, 44.3, 44.4, 44.5, 44.6, 44.7, 44.8, 44.9, 45.0, 45.1, 45.2, 45.3, 45.4, 45.5, 45.6, 45.7, 45.8, 45.9, 46.0, 46.1, 46.2, 46.3, 46.4, 46.5, 46.6, 46.7, 46.8, 46.9, 47.0, 47.1, 47.2, 47.3, 47.4, 47.5, 47.6, 47.7, 47.8, 47.9, 48.0, 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 48.7, 48.8, 48.9, 49.0, 49.1, 49.2, 49.3, 49.4, 49.5, 49.6, 49.7, 49.8, 49.9, 50.0, 50.1, 50.2, 50.3, 50.4, 50.5, 50.6, 50.7, 50.8, 50.9, 51.0, 51.1, 51.2, 51.3, 51.4, 51.5, 51.6, 51.7, 51.8, 51.9, 52.0, 52.1, 52.2, 52.3, 52.4, 52.5, 52.6, 52.7, 52.8, 52.9, 53.0, 53.1, 53.2, 53.3, 53.4, 53.5, 53.6, 53.7, 53.8, 53.9, 54.0, 54.1, 54.2, 54.3, 54.4, 54.5, 54.6, 54.7, 54.8, 54.9, 55.0, 55.1, 55.2, 55.3, 55.4, 55.5, 55.6, 55.7, 55.8, 55.9, 56.0, 56.1, 56.2, 56.3, 56.4, 56.5, 56.6, 56.7, 56.8, 56.9, 57.0, 57.1, 57.2, 57.3, 57.4, 57.5, 57.6, 57.7, 57.8, 57.9, 58.0, 58.1, 58.2, 58.3, 58.4, 58.5, 58.6, 58.7, 58.8, 58.9, 59.0, 59.1, 59.2, 59.3, 59.4, 59.5, 59.6, 59.7, 59.8, 59.9, 60.0, 60.1, 60.2, 60.3, 60.4, 60.5, 60.6, 60.7, 60.8, 60.9, 61.0, 61.1, 61.2, 61.3, 61.4, 61.5, 61.6, 61.7, 61.8, 61.9, 62.0, 62.1, 62.2, 62.3, 62.4, 62.5, 62.6, 62.7, 62.8, 62.9, 63.0, 63.1, 63.2, 63.3, 63.4, 63.5, 63.6, 63.7, 63.8, 63.9, 64.0, 64.1, 64.2, 64.3, 64.4, 64.5, 64.6, 64.7, 64.8, 64.9, 65.0, 65.1, 65.2, 65.3, 65.4, 65.5, 65.6, 65.7, 65.8, 65.9, 66.0, 66.1, 66.2, 66.3, 66.4, 66.5, 66.6, 66.7, 66.8, 66.9, 67.0, 67.1, 67.2, 67.3, 67.4, 67.5, 67.6, 67.7, 67.8, 67.9, 68.0, 68.1, 68.2, 68.3, 68.4, 68.5, 68.6, 68.7, 68.8, 68.9, 69.0, 69.1, 69.2, 69.3, 69.4, 69.5, 69.6, 69.7, 69.8, 69.9, 70.0, 70.1, 70.2, 70.3, 70.4, 70.5, 70.6, 70.7, 70.8, 70.9, 71.0, 71.1, 71.2, 71.3, 71.4, 71.5, 71.6, 71.7, 71.8, 71.9, 72.0, 72.1, 72.2, 72.3, 72.4, 72.5, 72.6, 72.7, 72.8, 72.9, 73.0, 73.1, 73.2, 73.3, 73.4, 73.5, 73.6, 73.7, 73.8, 73.9, 74.0, 74.1, 74.2, 74.3, 74.4, 74.5, 74.6, 74.7, 74.8, 74.9, 75.0, 75.1, 75.2, 75.3, 75.4, 75.5, 75.6, 75.7, 75.8, 75.9, 76.0, 76.1, 76.2, 76.3, 76.4, 76.5, 76.6, 76.7, 76.8, 76.9, 77.0, 77.1, 77.2, 77.3, 77.4, 77.5, 77.6, 77.7, 77.8, 77.9, 78.0, 78.1, 78.2, 78.3, 78.4, 78.5, 78.6, 78.7, 78.8, 78.9, 79.0, 79.1, 79.2, 79.3, 79.4, 79.5, 79.6, 79.7, 79.8, 79.9, 80.0, 80.1, 80.2, 80.3, 80.4, 80.5, 80.6, 80.7, 80.8, 80.9, 81.0, 81.1, 81.2, 81.3, 81.4, 81.5, 81.6, 81.7, 81.8, 81.9, 82.0, 82.1, 82.2, 82.3, 82.4, 82.5, 82.6, 82.7, 82.8, 82.9, 83.0, 83.1, 83.2, 83.3, 83.4, 83.5, 83.6, 83.7, 83.8, 83.9, 84.0, 84.1, 84.2, 84.3, 84.4, 84.5, 84.6, 84.7, 84.8, 84.9, 85.0, 85.1, 85.2, 85.3, 85.4, 85.5, 85.6, 85.7, 85.8, 85.9, 86.0, 86.1, 86.2, 86.3, 86.4, 86.5, 86.6, 86.7, 86.8, 86.9, 87.0, 87.1, 87.2, 87.3, 87.4, 87.5, 87.6, 87.7, 87.8, 87.9, 88.0, 88.1, 88.2, 88.3, 88.4, 88.5, 88.6, 88.7, 88.8, 88.9, 89.0, 89.1, 89.2, 89.3, 89.4, 89.5, 89.6, 89.7, 89.8, 89.9, 90.0, 90.1, 90.2, 90.3, 90.4, 90.5, 90.6, 90.7, 90.8, 90.9, 91.0, 91.1, 91.2, 91.3, 91.4, 91.5, 91.6, 91.7, 91.8, 91.9, 92.0, 92.1, 92.2, 92.3, 92.4, 92.5, 92.6, 92.7, 92.8, 92.9, 93.0, 93.1, 93.2, 93.3, 93.4, 93.5, 93.6, 93.7, 93.8, 93.9, 94.0, 94.1, 94.2, 94.3, 94.4, 94.5, 94.6, 94.7, 94.8, 94.9, 95.0, 95.1, 95.2, 95.3, 95.4, 95.5, 95.6, 95.7, 95.8, 95.9, 96.0, 96.1, 96.2, 96.3, 96.4, 96.5, 96.6, 96.7, 96.8, 96.9, 97.0, 97.1, 97.2, 97.3, 97.4, 97.5, 97.6, 97.7, 97.8, 97.9, 98.0, 98.1, 98.2, 98.3, 98.4, 98.5, 98.6, 98.7, 98.8, 98.9, 99.0, 99.1, 99.2, 99.3, 99.4, 99.5, 99.6, 99.7, 99.8, 99.9, 100.0
  };
  int nxbins_max = sizeof(xbins_max) / sizeof(xbins_max[0]) - 1;
  
  // L2Res pT binning (central+forward hybrid)
  double vptd[] =
  //{59.,85.,104.,170.,236., 302., 370., 460., 575.}; // central
  //{86., 110., 132., 204., 279., 373.} // forward
    {10, 15, 21, 28, 37, 49,
     59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
     638, 737, 846, 967, 1101, 1248,
     1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103, 3450, 4037, 5220};
  double nptd = sizeof(vptd) / sizeof(vptd[0]) - 1;

     // Inclusive jets pT binning
   double vx[] =
     {15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
      2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
      4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
   double nx = sizeof(vx) / sizeof(vx[0]) - 1;
   // Wide version
   double vxw[] =
     {15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 272, 330, 395, 468, 548, 686, 846, 1032, 1248, 1588, 2000, 2500, 3103, 3832, 4713, 5777, 7000};
   double nxw = sizeof(vxw) / sizeof(vxw[0]) - 1;
   // eXtra Wide version
   double vxx[] =
     {15, 28, 49, 84, 153, 272, 395, 548, 846, 1248, 2000, 3103, 4713, 7000};
   double nxx = sizeof(vxx) / sizeof(vxx[0]) - 1;
  
  //Bins
  Double_t xbins_three[] = {0.1,1.05,3.55,100};
  int nxbins_three = sizeof(xbins_three) / sizeof(xbins_three[0]) - 1;
  
  Double_t xbins_lead[] = {0.1,1.3,2.5,4,6,8,10,12,14,16,18,20,22,24,26,
  28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,95,100};
  int nxbins_lead = sizeof(xbins_lead) / sizeof(xbins_lead[0]) - 1;
  
  Double_t xbins_cand[] = {0.1,0.15,0.3,0.6,1,1.5,2,3,4,6,10,20,35,60,100};
  int nxbins_cand = sizeof(xbins_cand) / sizeof(xbins_cand[0]) - 1;
  
  Double_t xbins_jet[] = {0.1,0.15,0.3,0.6,1,1.5,2,3,4,6,10,20,35,60,100};
  int nxbins_jet = sizeof(xbins_jet) / sizeof(xbins_jet[0]) - 1;
  
  // Bins
  Double_t xbins_one[] = {0.1,100}; 
  int nxbins_one = sizeof(xbins_one) / sizeof(xbins_one[0]) - 1;
  
  TDirectory *curdir = gDirectory;
  TFile *fout = new TFile("output_z2.root","recreate");
  
  TProfile *pu0 = new TProfile("pu0",";Ptjet;has non-zero energysum",nptd,vptd);
  
  string vpid[] = {"kaon0l", "kaon0s", "kaonp", "kaonm", "lambdapm",//"lambdap", "lambdam","kaon0"
                   "sigmapm",//"sigma0", "sigmap", "sigmam", "sigmaantip", "sigmaantim",
                   "pionp", "pionm", "xiomega",//"pion0","xi0","xim","antixi0","antixim","omegam","antiomegam",
                   "proton", "antiproton", "neutron", "antineutron",
                   "electron", "muon","photon"}; //,"positron", "antimuon"
  int npid = sizeof(vpid) / sizeof(vpid[0]);

  string vpidS[] = {"kaon0l", "kaon0s", "kaonp", "kaonm", "lambdapm",//"lambdap", "lambdam","kaon0"
                   "sigmapm",//"sigma0", "sigmap", "sigmam", "sigmaantip", "sigmaantim",
                  "xiomega",//"pion0","xi0","xim","antixi0","antixim","omegam","antiomegam",
                   //"proton", "antiproton", "neutron", "antineutron",
                   //"electron", "muon","photon"
                   }; //,"positron", "antimuon"
  int npidS = sizeof(vpidS) / sizeof(vpidS[0]);
  
  //cout << npid << endl;
  
  string vpid2[] = {"K^{0}_{L}", "K^{0}_{S}", "K^{+}", "K^{-}", "#Lambda^{#pm}",//"#Lambda^{+}", "#Lambda^{-}", //"K^{0}",
		    "#Sigma^{#pm}",//"#Sigma^{0}", "#Sigma^{+}", "#Sigma^{-}","#bar{#Sigma}^{u}", "#bar{#Sigma}^{d}",
		    "#Pi^{+}", "#Pi^{-}","#Xi,#Omega", //"#Pi^{0}", "#Xi^{0}","#Xi^{-}","#bar{Xi^{0}}","#bar{#Xi^{-}}","#Omega^{-}","#bar{#Omega^{-}}",
		    "p", "#bar{p}", "n", "#bar{n}", 
		    "e", "#mu","#gamma"}; //"e^{+}", "#mu^{+}"
  int npid2 = sizeof(vpid2) / sizeof(vpid2[0]);

  string vpid2S[] = {"K^{0}_{L}", "K^{0}_{S}", "K^{+}", "K^{-}", "#Lambda^{#pm}",//"#Lambda^{+}", "#Lambda^{-}", //"K^{0}",
		    "#Sigma^{#pm}",//"#Sigma^{0}", "#Sigma^{+}", "#Sigma^{-}","#bar{#Sigma}^{u}", "#bar{#Sigma}^{d}",
		    "#Xi,#Omega", //"#Pi^{0}", "#Xi^{0}","#Xi^{-}","#bar{Xi^{0}}","#bar{#Xi^{-}}","#Omega^{-}","#bar{#Omega^{-}}",
		    //"p", "#bar{p}", "n", "#bar{n}", 
		    //"e", "#mu","#gamma"
        }; //"e^{+}", "#mu^{+}"
  int npid2S = sizeof(vpid2S) / sizeof(vpid2S[0]);
  
  /*int vpid3[] = {130, 310, 321, -321, {3122, -3122}, //311
                   {3222, 3112}, //3212, 3222, 3112, -3222, -3112,
                   211, -211, {3322,-3322,3312,-3312,3334,-3334}, //111, 3322, 3312, -3322, -3312, 3334, -3334,
                   2212, -2212, 2112, -2112,
                   11, 13, -11, -13, 22};
		   int npid3 = sizeof(vpid3) / sizeof(vpid3[0]);*/
  
  std::vector<std::vector<int>> vpid3 = {
    {130}, {310}, {321}, {-321}, {3122, -3122}, //311
    {3222, 3112, -3222, -3112}, //3212, 3222, 3112, -3222, -3112,
    {211}, {-211}, {3322, -3322, 3312, -3312, 3334, -3334}, //111, 3322, 3312, -3322, -3312, 3334, -3334,
    {2212}, {-2212}, {2112}, {-2112},
    {11, -11}, {13, -13}, {22} //, {-11}, {-13}
  };

   std::vector<std::vector<int>> vpid3S = {
    {130}, {310}, {321}, {-321}, {3122, -3122}, //311
    {3222, 3112, -3222, -3112}, //3212, 3222, 3112, -3222, -3112,
    {3322, -3322, 3312, -3312, 3334, -3334}, //111, 3322, 3312, -3322, -3312, 3334, -3334,
    //{11, -11}, {13, -13}, {22} //, {-11}, {-13}
  };

    std::vector<std::vector<int>> vpid3UD = {
    //{130}, {310}, {321}, {-321}, //{3122, -3122}, //311
    //{3222, 3112, -3222, -3112}, //3212, 3222, 3112, -3222, -3112,
    {211}, {-211}, //{3322, -3322, 3312, -3312, 3334, -3334}, //111, 3322, 3312, -3322, -3312, 3334, -3334,
    {2212}, {-2212}, {2112}, {-2112}
    //{11, -11}, {13, -13}, {22} //, {-11}, {-13}
  };
  
  
  std::map < string, int > mid;
  mid["kaon0l"] = 130;
  mid["kaon0s"] = 310;
  mid["kaonp"] = 321;
  mid["kaonm"] = -321;
  mid["kaon0"] = 311;
  mid["lambdap"] = 3122;
  mid["lambdam"] = -3122;
  mid["sigma0"] = 3212;
  mid["sigmap"] = 3222;
  mid["sigmam"] = 3112;
  mid["sigmaantip"] = -3222;
  mid["sigmaantim"] = -3112;
  mid["pion0"] = 111;
  mid["pionp"] = 211;
  mid["pionm"] = -211;
  mid["proton"] = 2212;
  mid["antiproton"] = -2212;
  mid["neutron"] = 2112;
  mid["antineutron"] = -2112;
  mid["electron"] = 11;
  mid["muon"] = 13;
  mid["positron"] = -11;
  mid["antimuon"] = -13;
  
  std::map < string, int > mq;
  mq["s"] = 3;
  mq["u"] = 2;
  mq["d"] = 1;
  
  
  string vq[] = {"d", "u", "s", "c", "b", "g"};
  int nq = sizeof(vq) / sizeof(vq[0]);
  
  int vq2[] = {1, 2, 3, 4, 5, 21};
  int nq2 = sizeof(vq2) / sizeof(vq2[0]);
  
  string vc[] = {"ch", "nh", "ne"};
  int nc = sizeof(vc) / sizeof(vc[0]);
  
  bool visc[] = {false, false, false};
  int nisc = sizeof(visc) / sizeof(visc[0]);
  
  string vxvar[] = {"ptcand", "ptjet", "ptlead"};
  int nxvar = sizeof(vxvar) / sizeof(vxvar[0]);
  
  string vyvar[] = {"flch", "flnh", "flne", "fhch", "fhnh", "fhne"};
  int nyvar = sizeof(vyvar) / sizeof(vyvar[0]);
  
  assert(npid == npid2); //+g
  
  std::map < string, TH1D* > mh;
  
  std::map < string, TProfile* > mp;
  
  
  for (int iq = 0; iq != nq; ++ iq) {
    for (int ix = 0; ix != nxvar; ++ ix) {
      const char *cx = vxvar[ix].c_str();
      const char *cq = vq[iq].c_str();
      Double_t *x = xbins_cand; // &xbins_cand[0]
      int nx = nxbins_cand;
      if (vxvar[ix] == "ptjet") { //swap if doesn't work
	      x = vptd;
	      nx = nptd;
      }
      const char *hname_all = Form("h_all_%s_vs_%s", cq, cx);
      const char *htitle_all = Form(";p_{T,%s}, all;%s N", cx, cq);
      mh[hname_all] = new TH1D(hname_all, htitle_all, nx, x);
      for (int id = 0; id != npid; ++ id) {
	      const char *pid = vpid[id].c_str();
	      const char *pid2 = vpid2[id].c_str();
	      if (vxvar[ix] == "ptjet") { //swap if doesn't work
	        x = vptd;
	        nx = nptd;
	      }
        const char *hname = Form("h_%s_%s_vs_%s", pid, cq, cx);
        const char *htitle = Form(";p_{T,%s}, %s;%s N", cx, pid2, cq);
        mh[hname] = new TH1D(hname, htitle, nx, x);
        for (int iy = 0; iy != nyvar && id == 0; ++ iy) {
          const char *cy = vyvar[iy].c_str();
          const char *pname = Form("h_%s_%s_%s_vs_%s", "all", cq, cy, cx);
          const char *ptitle = Form(";E_{%s}, %s;%s %s", cx, "all", cq, cy);
          mp[pname] = new TProfile(pname, ptitle, nx, x); 
        } // for iy
      } // for id
    } // for ix
  } // for iq
  const char *hsname_all = Form("h_all_%s_vs_ptslead", "s");
  const char *hstitle_all = Form(";p_{T,ptslead}, all;%s N", "s");
  mh[hsname_all] = new TH1D(hsname_all, hstitle_all, nxbins_cand, xbins_cand);
  for (int idS = 0; idS != npidS; ++ idS) {
    const char *pidS = vpidS[idS].c_str();
    const char *pid2S = vpid2S[idS].c_str();
    const char *hsname = Form("h_%s_s_vs_ptslead", pidS);
    const char *hstitle = Form(";p_{T,ptslead}, %s;s N", pid2S);
    mh[hsname] = new TH1D(hsname, hstitle, nxbins_cand, xbins_cand);
  }
/*
for (int iq = 0; iq != nq; ++ iq) {
    for (int ix = 0; ix != nxvar; ++ ix) {
      const char *cx = vxvar[ix].c_str();
      const char *cq = vq[iq].c_str();
      Double_t *x = xbins_cand; // &xbins_cand[0]
      int nx = nxbins_cand;
      if (vxvar[ix] == "ptjet") { //swap if doesn't work
	      x = vptd;
	      nx = nptd;
      }
      const char *hname_all = Form("h_all_sanslead_%s_vs_%s", cq, cx);
      const char *htitle_all = Form(";E_{T,%s}, all sans lead;%s N", cx, cq);
      mh[hname_all] = new TH1D(hname_all, htitle_all, nx, x);
      for (int id = 0; id != npid; ++ id) {
	      const char *pid = vpid[id].c_str();
	      const char *pid2 = vpid2[id].c_str();
	      if (vxvar[ix] == "ptjet") { //swap if doesn't work
	        x = vptd;
	        nx = nptd;
	      }
        const char *hname = Form("h_sanslead_%s_%s_vs_%s", pid, cq, cx);
        const char *htitle = Form(";p_{T,%s}, %s sans lead;%s N", cx, pid2, cq);
        mh[hname] = new TH1D(hname, htitle, nx, x);
      } // for id
    } // for ix
  } // for iq
  */
  std::map < string, TH1D* > mhc;
  
  for (int iq = 0; iq != nq; ++ iq) {
    for (int ic = 0; ic != nc; ++ ic) {
      for (int ix = 0; ix != nxvar; ++ ix) {
        const char *cx = vxvar[ix].c_str();
        const char *cq = vq[iq].c_str();
        const char *cc = vc[ic].c_str();
        Double_t *x = xbins_cand; // &xbins_cand[0]
        int nx = nxbins_cand;
        if (vxvar[ix] == "ptjet") { //swap if doesn't work
          x = vptd;
          nx = nptd;
        }
        const char *hname = Form("h_%s_%s_vs_%s", cq, cc ,cx);
        const char *htitle = Form(";p_{T,%s}, %s;%s N", cx, cq, cc);
        //if (debug) {cout << hname << endl << flush;}
        mhc[hname] = new TH1D(hname, htitle, nx, x);
      } // for ix
    } // for ic
  } // for iq

  for (int ic = 0; ic != nc; ++ ic) {
    const char *cc = vc[ic].c_str();
    const char *hsname = Form("h_s_%s_vs_ptslead", cc);
    const char *hstitle = Form(";p_{T,ptslead}, s; %s N", cc);
    mhc[hsname] = new TH1D(hsname, hstitle, nxbins_cand, xbins_cand);
  }
   std::map < string, TH1D* > mhe;

  for (int iq = 0; iq != nq; ++ iq) {
      for (int ix = 0; ix != nxvar; ++ ix) {
        const char *cx = vxvar[ix].c_str();
        const char *cq = vq[iq].c_str();
        Double_t *x = xbins_cand; // &xbins_cand[0]
        int nx = nxbins_cand;
        if (vxvar[ix] == "ptjet") { //swap if doesn't work
          x = vptd;
          nx = nptd;
        }
        const char *hname = Form("h_elead/ejet_%s_vs_%s", cq,cx);
        const char *htitle = Form(";E_{%s}, %s; N", cx, cq);
        //if (debug) {cout << hname << endl << flush;}
        mhe[hname] = new TH1D(hname, htitle, nx, x);
      } // for ix
  } // for iq
  
  // All jets, all candidates
  TH1D *h_all = new TH1D("h_all",";PtCand, in all jets all candidates;N",nxbins_lead, xbins_lead);
  // All jets, different candidates
  TH1D *h_all_s = new TH1D("h_all_s",";PtCand, in all jets s candidates;N",nxbins_lead, xbins_lead);
  TH1D *h_all_ud = new TH1D("h_all_ud",";PtCand, in all jets ud candidates;N",nxbins_lead, xbins_lead);

  // S jets, all candidates
  TH1D *h_s = new TH1D("h_s",";PtCand, in s jet all candidates;N",nxbins_lead, xbins_lead);
  // S jets, different candidates
  TH1D *h_s_s = new TH1D("h_s_s",";PtCand, in s jets s candidates;N",nxbins_lead, xbins_lead);
  TH1D *h_s_ud = new TH1D("h_s_ud",";PtCand, in s jets ud candidates;N",nxbins_lead, xbins_lead);

  // UD jets, all candidates
  TH1D *h_ud = new TH1D("h_ud",";PtCand, in ud jet all candidates;N",nxbins_lead, xbins_lead);
  // UD jets, different candidates
  TH1D *h_ud_s = new TH1D("h_ud_s",";PtCand, in ud jets s candidates;N",nxbins_lead, xbins_lead);
  TH1D *h_ud_ud = new TH1D("h_ud_ud",";PtCand, in ud jets ud candidates;N",nxbins_lead, xbins_lead);

  // G jets, all candidates
  TH1D *h_g = new TH1D("h_g",";PtCand, in g jet all candidates;N",nxbins_lead, xbins_lead);
  // G jets, different candidates
  TH1D *h_g_s = new TH1D("h_g_s",";PtCand, in g jets s candidates;N",nxbins_lead, xbins_lead);
  TH1D *h_g_ud = new TH1D("h_g_ud",";PtCand, in g jets ud candidates;N",nxbins_lead, xbins_lead);

  // All jets, all leading candidates
  TH1D *h_all_lead = new TH1D("h_all_lead",";PtCand, in all jets all leading candidates;N", nxbins_lead, xbins_lead);
  // All jets, different leading candidates
  TH1D *h_all_slead = new TH1D("h_all_slead",";PtCand, in all jets leading s candidates;N", nxbins_lead, xbins_lead);
  TH1D *h_all_udlead = new TH1D("h_all_udlead",";PtCand, in all jets leading ud candidates;N", nxbins_lead, xbins_lead);
 
  // S jets, all leading candidates
  TH1D *h_s_lead = new TH1D("h_s_lead",";PtCand, in s jet leading candidates;N", nxbins_lead, xbins_lead);
  // S jets, different leading candidates
  TH1D *h_s_slead = new TH1D("h_s_slead",";PtCand, in s jets leading s candidates;N", nxbins_lead, xbins_lead);
  TH1D *h_s_udlead = new TH1D("h_s_udlead",";PtCand, in s jets leading ud candidates;N", nxbins_lead, xbins_lead);

  // S jets, leading is s candidate
  TH1D *h_s_lead_is_slead = new TH1D("h_s_lead_is_slead",";PtCand, in s jets leading ptcl is s candidate;N", nxbins_lead, xbins_lead);
  // S jets, leading is ud candidate
  //TH1D *h_s_lead_is_udlead = new TH1D("h_s_lead_is_udlead",";PtCand, in s jets leading ptcl is ud candidate;N", nxbins_lead, xbins_lead);
  // All jets, leading is s candidate
  TH1D *h_all_lead_is_slead = new TH1D("h_all_lead_is_slead",";PtCand, in all jets leading ptcl is s candidate;N", nxbins_lead, xbins_lead);
  // All jets, leading is ud candidate
  //TH1D *h_all_lead_is_udlead = new TH1D("h_all_lead_is_udlead",";PtCand, in all jets leading ptcl is ud candidate;N", nxbins_lead, xbins_lead);

  // UD jets, all leading candidates
  TH1D *h_ud_lead = new TH1D("h_ud_lead",";PtCand, in ud jets leading candidates;N", nxbins_lead, xbins_lead);
  // UD jets, different leading candidates
  TH1D *h_ud_slead = new TH1D("h_ud_slead",";PtCand, in ud jets leading s candidates;N", nxbins_lead, xbins_lead);
  TH1D *h_ud_udlead = new TH1D("h_ud_udlead",";PtCand, in ud jets leading ud candidates;N", nxbins_lead, xbins_lead);

  // G jets, all leading candidates
  TH1D *h_g_lead = new TH1D("h_g_lead",";PtCand, in g jets leading candidates;N", nxbins_lead, xbins_lead);
  // G jets, different leading candidates
  TH1D *h_g_slead = new TH1D("h_g_slead",";PtCand, in g jets leading s candidates;N", nxbins_lead, xbins_lead);
  TH1D *h_g_udlead = new TH1D("h_g_udlead",";PtCand, in g jets leading ud candidates;N", nxbins_lead, xbins_lead);

  TH1D *hperp = new TH1D("hperp",";PtCand, No jets;N", nxbins_max, xbins_max);
  TH1D *hjet = new TH1D("hjet",";PtCand, All jets;N", nxbins_max, xbins_max);
  //TH1D *h_elead_vs_ejet = new TH1D("h_elead_vs_ejet",";p_{T,jet}; lead/jet energyfraction",100,0,1);
  //TH1D *h_deltaR_lead = new TH1D("h_deltaR_lead",";DeltaR lead/jet ;N_{jet}", 115, 0,0.50178);
  //TH1D *h_deltaR_sec = new TH1D("h_deltaR_sec",";DeltaR sec/jet ;N_{jet}", 115, 0,0.50178);
  //TH1D *h_deltaR_sec_e = new TH1D("h_deltaR_sec_e",";DeltaR sec/jet ;N_{jet} #cdot elead/ejet", 115, 0,0.50178);
  //TH1D *h_ejet = new TH1D("h_ejet",";E_{jet};N",100,0,1);

  TProfile *prq = new TProfile("prq",";p_{T,genjet};Response (ud)",nx,vx);
  TProfile *pru = new TProfile("pru",";p_{T,genjet};Response (u)",nx,vx);
  TProfile *prd = new TProfile("prd",";p_{T,genjet};Response (d)",nx,vx);
  TProfile *prs = new TProfile("prs",";p_{T,genjet};Response (s)",nx,vx);
  TProfile *prc = new TProfile("prc",";p_{T,genjet};Response (c)",nx,vx);
  TProfile *prb = new TProfile("prb",";p_{T,genjet};Response (b)",nx,vx);
  TProfile *prg = new TProfile("prg",";p_{T,genjet};Response (g)",nx,vx);

  TProfile *pcq = new TProfile("pcq",";p_{T,genjet};Response (ud)",nx,vx);
  TProfile *pcu = new TProfile("pcu",";p_{T,genjet};Response (u)",nx,vx);
  TProfile *pcd = new TProfile("pcd",";p_{T,genjet};Response (d)",nx,vx);
  TProfile *pcs = new TProfile("pcs",";p_{T,genjet};Response (s)",nx,vx);
  TProfile *pcc = new TProfile("pcc",";p_{T,genjet};Response (c)",nx,vx);
  TProfile *pcb = new TProfile("pcb",";p_{T,genjet};Response (b)",nx,vx);
  TProfile *pcg = new TProfile("pcg",";p_{T,genjet};Response (g)",nx,vx);

  // Bin width 1./4. of ECAL crystal (2pi/(72.*5))
  TH1D *hdrq = new TH1D("hdrq",";#DeltaR(ud-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdru = new TH1D("hdru",";#DeltaR(u-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdrd = new TH1D("hdrd",";#DeltaR(d-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdrs = new TH1D("hdrs",";#DeltaR(s-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdrc = new TH1D("hdrc",";#DeltaR(c-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdrb = new TH1D("hdrb",";#DeltaR(b-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hdrg = new TH1D("hdrg",";#DeltaR(g-jet,cand);N_{cand}",200,0,0.8727);

  TH1D *hrwq = new TH1D("hrwq",";#DeltaR(ud-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrwu = new TH1D("hrwu",";#DeltaR(u-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrwd = new TH1D("hrwd",";#DeltaR(d-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrws = new TH1D("hrws",";#DeltaR(s-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrwc = new TH1D("hrwc",";#DeltaR(c-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrwb = new TH1D("hrwb",";#DeltaR(b-jet,cand);N_{cand}",200,0,0.8727);
  TH1D *hrwg = new TH1D("hrwg",";#DeltaR(g-jet,cand);N_{cand}",200,0,0.8727);

  // Fraction of energy carried by the leading particle
  TH1D *hfq = new TH1D("hfq",";p_{T,lead}/p_{T,jet};N_{ud-jet}",100,0,1);
  TH1D *hfu = new TH1D("hfu",";p_{T,lead}/p_{T,jet};N_{u-jet}",100,0,1);
  TH1D *hfd = new TH1D("hfd",";p_{T,lead}/p_{T,jet};N_{d-jet}",100,0,1);
  TH1D *hfs = new TH1D("hfs",";p_{T,lead}/p_{T,jet};N_{s-jet}",100,0,1);
  TH1D *hfc = new TH1D("hfc",";p_{T,lead}/p_{T,jet};N_{c-jet}",100,0,1);
  TH1D *hfb = new TH1D("hfb",";p_{T,lead}/p_{T,jet};N_{b-jet}",100,0,1);
  TH1D *hfg = new TH1D("hfg",";p_{T,lead}/p_{T,jet};N_{g-jet}",100,0,1);

  // Fragmentation function
  TH1D *hffq = new TH1D("hffq",";p_{T,cand}/p_{T,jet};N_{ud-jet}",100,0,1);
  TH1D *hffu = new TH1D("hffu",";p_{T,cand}/p_{T,jet};N_{u-jet}",100,0,1);
  TH1D *hffd = new TH1D("hffd",";p_{T,cand}/p_{T,jet};N_{d-jet}",100,0,1);
  TH1D *hffs = new TH1D("hffs",";p_{T,cand}/p_{T,jet};N_{s-jet}",100,0,1);
  TH1D *hffc = new TH1D("hffc",";p_{T,cand}/p_{T,jet};N_{c-jet}",100,0,1);
  TH1D *hffb = new TH1D("hffb",";p_{T,cand}/p_{T,jet};N_{b-jet}",100,0,1);
  TH1D *hffg = new TH1D("hffg",";p_{T,cand}/p_{T,jet};N_{g-jet}",100,0,1);

  // Fragmentation function sans leading particle
  TH1D *hf2q = new TH1D("hf2q",";p_{T,cand}/p_{T,jet};N_{ud-jet}",100,0,1);
  TH1D *hf2u = new TH1D("hf2u",";p_{T,cand}/p_{T,jet};N_{u-jet}",100,0,1);
  TH1D *hf2d = new TH1D("hf2d",";p_{T,cand}/p_{T,jet};N_{d-jet}",100,0,1);
  TH1D *hf2s = new TH1D("hf2s",";p_{T,cand}/p_{T,jet};N_{s-jet}",100,0,1);
  TH1D *hf2c = new TH1D("hf2c",";p_{T,cand}/p_{T,jet};N_{c-jet}",100,0,1);
  TH1D *hf2b = new TH1D("hf2b",";p_{T,cand}/p_{T,jet};N_{b-jet}",100,0,1);
  TH1D *hf2g = new TH1D("hf2g",";p_{T,cand}/p_{T,jet};N_{g-jet}",100,0,1);




  // Histograms for Jet pt > 30 GeV
TH1D *h_c_ctag_cvl_df30 = new TH1D("h_c_ctag_cvl_df30", "CvL tag for c-jets at pt>30; CvL; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvl_pn30 = new TH1D("h_c_ctag_cvl_pn30", "CvL PNet for c-jets at pt>30; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_df30 = new TH1D("h_c_ctag_cvb_df30", "CvB tag for c-jets at pt>30; CvB; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_pn30 = new TH1D("h_c_ctag_cvb_pn30", "CvB PNet for c-jets at pt>30; CvB PNet; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_df30 = new TH1D("h_b_ctag_b_df30", "DeepFlavB for b-jets at pt>30; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_pn30 = new TH1D("h_b_ctag_b_pn30", "PNetB for b-jets at pt>30; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_df30 = new TH1D("h_c_udstag_cvl_df30", "CvL for uds-jets at pt>30; CvL; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_pn30 = new TH1D("h_c_udstag_cvl_pn30", "CvL PNet for uds-jets at pt>30; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_df30 = new TH1D("h_b_udstag_b_df30", "DeepFlavB for uds-jets at pt>30; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_pn30 = new TH1D("h_b_udstag_b_pn30", "PNetB for uds-jets at pt>30; PNetB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_df30 = new TH1D("h_b_btag_b_df30", "DeepFlavB for b-jets in b-jets at pt>30; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_pn30 = new TH1D("h_b_btag_b_pn30", "PNetB for b-jets in b-jets at pt>30; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_df30 = new TH1D("h_c_btag_cvb_df30", "CvB for c-jets in b-jets at pt>30; CvB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_pn30 = new TH1D("h_c_btag_cvb_pn30", "CvB PNet for c-jets in b-jets at pt>30; CvB PNet; Fraction", 100, 0, 1);

// Histograms for Jet pt > 60 GeV
TH1D *h_c_ctag_cvl_df60 = new TH1D("h_c_ctag_cvl_df60", "CvL tag for c-jets at pt>60; CvL; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvl_pn60 = new TH1D("h_c_ctag_cvl_pn60", "CvL PNet for c-jets at pt>60; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_df60 = new TH1D("h_c_ctag_cvb_df60", "CvB tag for c-jets at pt>60; CvB; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_pn60 = new TH1D("h_c_ctag_cvb_pn60", "CvB PNet for c-jets at pt>60; CvB PNet; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_df60 = new TH1D("h_b_ctag_b_df60", "DeepFlavB for b-jets at pt>60; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_pn60 = new TH1D("h_b_ctag_b_pn60", "PNetB for b-jets at pt>60; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_df60 = new TH1D("h_c_udstag_cvl_df60", "CvL for uds-jets at pt>60; CvL; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_pn60 = new TH1D("h_c_udstag_cvl_pn60", "CvL PNet for uds-jets at pt>60; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_df60 = new TH1D("h_b_udstag_b_df60", "DeepFlavB for uds-jets at pt>60; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_pn60 = new TH1D("h_b_udstag_b_pn60", "PNetB for uds-jets at pt>60; PNetB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_df60 = new TH1D("h_b_btag_b_df60", "DeepFlavB for b-jets in b-jets at pt>60; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_pn60 = new TH1D("h_b_btag_b_pn60", "PNetB for b-jets in b-jets at pt>60; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_df60 = new TH1D("h_c_btag_cvb_df60", "CvB for c-jets in b-jets at pt>60; CvB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_pn60 = new TH1D("h_c_btag_cvb_pn60", "CvB PNet for c-jets in b-jets at pt>60; CvB PNet; Fraction", 100, 0, 1);

// Histograms for Jet pt > 80 GeV
TH1D *h_c_ctag_cvl_df80 = new TH1D("h_c_ctag_cvl_df80", "CvL tag for c-jets at pt>80; CvL; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvl_pn80 = new TH1D("h_c_ctag_cvl_pn80", "CvL PNet for c-jets at pt>80; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_df80 = new TH1D("h_c_ctag_cvb_df80", "CvB tag for c-jets at pt>80; CvB; Fraction", 100, 0, 1);
TH1D *h_c_ctag_cvb_pn80 = new TH1D("h_c_ctag_cvb_pn80", "CvB PNet for c-jets at pt>80; CvB PNet; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_df80 = new TH1D("h_b_ctag_b_df80", "DeepFlavB for b-jets at pt>80; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_ctag_b_pn80 = new TH1D("h_b_ctag_b_pn80", "PNetB for b-jets at pt>80; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_df80 = new TH1D("h_c_udstag_cvl_df80", "CvL for uds-jets at pt>80; CvL; Fraction", 100, 0, 1);
TH1D *h_c_udstag_cvl_pn80 = new TH1D("h_c_udstag_cvl_pn80", "CvL PNet for uds-jets at pt>80; CvL PNet; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_df80 = new TH1D("h_b_udstag_b_df80", "DeepFlavB for uds-jets at pt>80; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_udstag_b_pn80 = new TH1D("h_b_udstag_b_pn80", "PNetB for uds-jets at pt>80; PNetB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_df80 = new TH1D("h_b_btag_b_df80", "DeepFlavB for b-jets in b-jets at pt>80; DeepFlavB; Fraction", 100, 0, 1);
TH1D *h_b_btag_b_pn80 = new TH1D("h_b_btag_b_pn80", "PNetB for b-jets in b-jets at pt>80; PNetB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_df80 = new TH1D("h_c_btag_cvb_df80", "CvB for c-jets in b-jets at pt>80; CvB; Fraction", 100, 0, 1);
TH1D *h_c_btag_cvb_pn80 = new TH1D("h_c_btag_cvb_pn80", "CvB PNet for c-jets in b-jets at pt>80; CvB PNet; Fraction", 100, 0, 1);
   
   
  
  
  TLorentzVector p4jet, p4cand, p4perp, p4sum, p4lead, p4reco;
  
  
  curdir->cd();
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  cout << "Processing " << nentries << " events" << endl << flush;
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::cout << std::ctime(&now) << std::endl<< flush;  
  TStopwatch t;
  t.Start();
  const int nlap = 1000;
  const int nlap2 = 80000;
  //nentries = 100000;
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
    //      cout << "Processing Event " << event << endl;
    double w = genWeight;

      ////////////////
      // Particle loop
      ////////////////

    // Loop over two leading jets
    for (int i = 0; i != min(2,nGenJet); ++i) {
      //cout << "  GenJet " << i << ":" << GenJet_pt[i] << ", " << endl;
      
      //if (debug){cout << "here1" << endl;}
      bool isSjet = (abs(GenJet_partonFlavour[i]) == 3);
      bool isDjet = (abs(GenJet_partonFlavour[i]) == 1);
      bool isUjet = (abs(GenJet_partonFlavour[i]) == 2);
      bool isGjet = (abs(GenJet_partonFlavour[i]) == 21);
      
      p4jet.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
      p4perp.SetPtEtaPhiM(GenJet_pt[i], -GenJet_eta[i], GenJet_phi[i]+TMath::Pi()/2, GenJet_mass[i]);
      p4sum.SetPtEtaPhiM(0,0,0,0);
      p4lead.SetPtEtaPhiM(0,0,0,0);

      double ptjet = GenJet_pt[i];
      double ejet = p4jet.E();

      //bool pass = (p4j.Pt()>84 && p4j.Pt()<153 && fabs(p4j.Eta())<1.3)
      bool pass = (p4jet.Pt()>80 && p4jet.Pt()<100 && fabs(p4jet.Eta())<1.3);

      double esum(0), esumch(0), esumnh(0), esumne(0);
      double elead(0), eleadch(0), eleadnh(0), eleadne(0);
      double ptlead(0.);
      double ptslead(0.);
      double ptudlead(0.);
      int iLeadGenCand(-1);
      int iLeadSCand(-1);
      int iLeadUDCand(-1);
      
        //continue;

      // Pre-loop to find leading candidate in jet
      for (int k = 0; k != nGenJetGenPartCand; ++k) {

        if (GenJetGenPartCand_genJetIdx[k]==i) {
          int j = GenJetGenPartCand_GenPartCandIdx[k];
          double pt = GenPartCand_pt[j];
          int pdgId = GenPartCand_pdgId[j];

          // Update ptlead for the highest pt candidate
          if (pt>ptlead) {
            p4lead.SetPtEtaPhiM(GenPartCand_pt[j],GenPartCand_eta[j],
                GenPartCand_phi[j],GenPartCand_mass[j]);
            ptlead = pt;
            elead = p4lead.E(); //leading pt particle energy != leading energy particle
            p4lead = p4cand; 
            iLeadGenCand = j;
            
            // Check if PDG ID is in vpid3S and update ptslead for the highest pt candidate in vpid3S
            bool isSCand = false;
            for (const auto& idList : vpid3S) {
              if (std::find(idList.begin(), idList.end(), pdgId) != idList.end()) {
                isSCand = true;
                break; // PDG ID is in the list
              }
            }

            if (isSCand && pt > ptslead) {
              ptslead = pt;
              iLeadSCand = j;
            }
            
            bool isUDCand = false;
            for (const auto& idList : vpid3UD) {
              if (std::find(idList.begin(), idList.end(), pdgId) != idList.end()) {
                isUDCand = true;
                break; // PDG ID is in the list
              }
            }
            if (isUDCand && pt > ptudlead) {
              ptudlead = pt;
              iLeadUDCand = j;
            }
          } //pt > ptlead
        } // if GenJetGenPartCand_genJetIdx[k]==i
      } // for k

      if (pass) {
        Short_t f = GenJet_partonFlavour[i];
        double fi = ptlead / p4jet.Pt();
        if (abs(f)==1|| abs(f)==2)
          hfq->Fill(fi);
        if (abs(f)==1) hfd->Fill(fi);
        if (abs(f)==2) hfu->Fill(fi);
        if (abs(f)==3) hfs->Fill(fi);
        if (abs(f)==4) hfc->Fill(fi);
        if (abs(f)==5) hfb->Fill(fi);
        if (f==21)
        hfg->Fill(fi);
      
      } // if pass

      // Loop over candidates in these jets
      for (int k = 0; k != nGenJetGenPartCand; ++k) {
        int iGenJet = GenJetGenPartCand_genJetIdx[k];
        int iGenCand = GenJetGenPartCand_GenPartCandIdx[k];
        int pdgId = GenPartCand_pdgId[iGenCand];

        bool isSCand = false;
          for (const auto& idList : vpid3S) {
            if (std::find(idList.begin(), idList.end(), pdgId) != idList.end()) {
              isSCand = true;
              break; // PDG ID is in the list
            }
          }

        bool isUDCand = false;
        for (const auto& idList : vpid3UD) {
          if (std::find(idList.begin(), idList.end(), pdgId) != idList.end()) {
            isUDCand = true;
            break; // PDG ID is in the list
          }
        }
	
        bool isPhoton = (GenPartCand_pdgId[iGenCand] == 22);
	
        if (p4perp.DeltaR(p4cand) < 0.4){
          if (fabs(GenJet_eta[iGenJet]) < 1.3 && ptjet > 80 && ptjet < 100) {
            hperp->Fill(GenPartCand_pt[iGenCand], w);
          }
        }
        if (p4jet.DeltaR(p4cand) < 0.4){
          if (fabs(GenJet_eta[iGenJet]) < 1.3 && ptjet > 80 && ptjet < 100) {
            hjet->Fill(GenPartCand_pt[iGenCand], w);
          }
        }
	
        if (i == iGenJet) {
          // cout << "Cand" << iGenCand << ":" << GenPartCand_pt[iGenCand] << ", " ;
          p4cand.SetPtEtaPhiM(GenPartCand_pt[iGenCand], GenPartCand_eta[iGenCand], GenPartCand_phi[iGenCand], GenPartCand_mass[iGenCand]);
          if (pass) {
            double dr = p4jet.DeltaR(p4cand);
            //double dr = p4lead.DeltaR(p4cand);
            Short_t f = GenJet_partonFlavour[i];
            if (abs(f)==1|| abs(f)==2) hdrq->Fill(dr);
            if (abs(f)==1) hdrd->Fill(dr);
            if (abs(f)==2) hdru->Fill(dr);
            if (abs(f)==3) hdrs->Fill(dr);
            if (abs(f)==4) hdrc->Fill(dr);
            if (abs(f)==5) hdrb->Fill(dr);
            if (f==21) hdrg->Fill(dr);

            double w = p4cand.Pt() / p4jet.Pt();
            if (abs(f)==1|| abs(f)==2) hrwq->Fill(dr, w);
            if (abs(f)==1) hrwd->Fill(dr, w);
            if (abs(f)==2) hrwu->Fill(dr, w);
            if (abs(f)==3) hrws->Fill(dr, w);
            if (abs(f)==4) hrwc->Fill(dr, w);
            if (abs(f)==5) hrwb->Fill(dr, w);
            if (f==21) hrwg->Fill(dr, w);

            double fi = p4cand.Pt() / p4jet.Pt();
            if (abs(f)==1|| abs(f)==2) hffq->Fill(fi);
            if (abs(f)==1) hffd->Fill(fi);
            if (abs(f)==2) hffu->Fill(fi);
            if (abs(f)==3) hffs->Fill(fi);
            if (abs(f)==4) hffc->Fill(fi);
            if (abs(f)==5) hffb->Fill(fi);
            if (f==21) hffg->Fill(fi);

            if (iGenCand != iLeadGenCand) {
              if (abs(f)==1|| abs(f)==2) hf2q->Fill(fi);
              if (abs(f)==1) hf2d->Fill(fi);
              if (abs(f)==2) hf2u->Fill(fi);
              if (abs(f)==3) hf2s->Fill(fi);
              if (abs(f)==4) hf2c->Fill(fi);
              if (abs(f)==5) hf2b->Fill(fi);
              if (f==21) hf2g->Fill(fi);
            } // lead cand
          } // if pass    
          //if (debug){cout << "jetloop" << endl;}
          
          bool isCH = (GenPartCand_charge[iGenCand] != 0);
          bool isNH = (GenPartCand_charge[iGenCand] == 0 && !isPhoton);
          bool isNE = (GenPartCand_charge[iGenCand] == 0 && isPhoton);
	  
          visc[0] = isCH;
          visc[1] = isNH;
          visc[2] = isNE;

	        double ecand = p4cand.E();
          if (iGenCand == iLeadGenCand) { 
            if (isCH) {
              eleadch = ecand;
             }
            if (isNH) {
              eleadnh = ecand;
            }
            if (isNE) { 
              eleadne = ecand;
            }
          } // if iGenCand == iLeadGenCand
	  
            p4sum += p4cand;
            esum += ecand;
            if (isCH) { esumch += ecand; }
            if (isNH) { esumnh += ecand; }
            if (isNE) { esumne += ecand; }
            
            if (fabs(GenJet_eta[iGenJet]) < 1.3) {
              
            if (ptjet > 80 && ptjet < 100) {
              //in all jets all candidates
              h_all->Fill(GenPartCand_pt[iGenCand], w);
              //in all jets s candidates
              if (isSCand) h_all_s->Fill(GenPartCand_pt[iGenCand], w);
              //in all jets ud candidates
              if (isUDCand) h_all_ud->Fill(GenPartCand_pt[iGenCand], w);

              //in s jets all candidates
              if (isSjet) h_s->Fill(GenPartCand_pt[iGenCand], w);
              //in s jets s candidates
              if (isSjet && isSCand) h_s_s->Fill(GenPartCand_pt[iGenCand], w);
              //in s jets ud candidates
              if (isSjet && isUDCand) h_s_ud->Fill(GenPartCand_pt[iGenCand], w);

              //in ud jets all candidates
              if (isUjet || isDjet) h_ud->Fill(GenPartCand_pt[iGenCand], w);
              //in ud jets s candidates
              if ((isUjet || isDjet) && isSCand) h_ud_s->Fill(GenPartCand_pt[iGenCand], w);
              //in ud jets ud candidates
              if ((isUjet || isDjet) && isUDCand) h_ud_ud->Fill(GenPartCand_pt[iGenCand], w);

              // in g jets all candidates
              if (isGjet) h_g->Fill(GenPartCand_pt[iGenCand], w);
              //in g jets s candidates
              if (isGjet && isSCand) h_g_s->Fill(GenPartCand_pt[iGenCand], w);
              //in g jets ud candidates
              if (isGjet && isUDCand) h_g_ud->Fill(GenPartCand_pt[iGenCand], w);

              // in all jets all leading candidates
              if (iLeadGenCand == iGenCand) {h_all_lead->Fill(ptlead, w);}
              //in all jets leading s candidates
              if (iLeadSCand == iGenCand) {h_all_slead->Fill(ptslead, w);}
              //in all jets leading ud candidates
              if (iLeadUDCand == iGenCand) {h_all_udlead->Fill(ptudlead, w);}
              //in all jets leading is s candidate
              if (iLeadGenCand == iGenCand && isSCand) {h_all_lead_is_slead->Fill(ptlead, w);}
              //in all jets leading is ud candidate
              //if (iLeadGenCand == iGenCand && isUDCand) {h_all_lead_is_udlead->Fill(ptlead, w);}


              // in s jet all leading candidates
              if (iLeadGenCand == iGenCand && isSjet) {h_s_lead->Fill(ptlead, w);}
              //in s jet leading s candidates
              if (iLeadSCand == iGenCand && isSjet) {h_s_slead->Fill(ptslead, w);}
              //in s jet leading ud candidates
              if (iLeadUDCand == iGenCand && isSjet) {h_s_udlead->Fill(ptudlead, w);}
              // in s jet leading is s candidate
              if (iLeadGenCand == iGenCand && isSjet && isSCand) {h_s_lead_is_slead->Fill(ptlead, w);}
              // in s jet leading is ud candidate
              //if (iLeadGenCand == iGenCand && isSjet && isUDCand) {h_s_lead_is_udlead->Fill(ptlead, w);}

              // in ud jet all leading candidates
              if (iLeadGenCand == iGenCand && (isUjet || isDjet)) {h_ud_lead->Fill(ptlead, w);}
              //in ud jet leading s candidates
              if (iLeadSCand == iGenCand && (isUjet || isDjet)) {h_ud_slead->Fill(ptslead, w);}   
              //in ud jet leading ud candidates
              if (iLeadUDCand == iGenCand && (isUjet || isDjet)) {h_ud_udlead->Fill(ptudlead, w);}

              // in g jet all leading candidates
              if (iLeadGenCand == iGenCand && isGjet) {h_g_lead->Fill(ptlead, w);}
              // in g jet leading s candidates
              if (iLeadSCand == iGenCand && isGjet) {h_g_slead->Fill(ptslead, w);}
              // in g jet leading udcandidates
              if (iLeadUDCand == iGenCand && isGjet) {h_g_udlead->Fill(ptudlead, w);}
            }
            
            //if (debug) { cout << "isId=" << isId << ", GenPartCand_pdgId[iGenCand]=" << GenPartCand_pdgId[iGenCand] << ", vpid3[id]=" << vpid3[id] << endl; }
            //bool isId = (GenPartCand_pdgId[iGenCand] == mid[pid]);

            bool isId2(false);
            // Check if GenPartCand_pdgId[iGenCand] is in any of the vpid3 vectors
            for (std::size_t id = 0; id < vpid3.size() && !isId2; ++id) { // Iterate over vpid3 vectors
                for (std::size_t elem = 0; elem < vpid3[id].size(); ++elem) { // Iterate over elements within a vector
                    if (GenPartCand_pdgId[iGenCand] == vpid3[id][elem]) {
                        isId2 = true; // Found a match, set isId2 to true
                        break; // No need to continue checking this vector
                    }
                }
                if (isId2) break; // Found a match, no need to check further vectors
            }

            // If GenPartCand_pdgId[iGenCand] is not found in any vector, print the details
            if (!isId2) {
              if (currentFile != nullptr) {
                std::string fullPath(currentFile->GetName());
                std::size_t lastSlashPos = fullPath.find_last_of("/\\"); // Handles both forward and backward slashes
                std::string fileName = fullPath.substr(lastSlashPos + 1); // Extracts the filename
                cout << "File: " << fileName << "," << endl
                << "run: " << run // The run number
                << ", event: " << event // The event number
                << ", jentry: " << jentry // The jentry number
                << ", ientry: " << ientry // The jentry number
                << ", GenJet index: " << iGenJet // Print the GenJet index
                << ", GenPartCand index: " << iGenCand // Print the GenPartCand index
                << ", GenPartCand_pdgId: " << GenPartCand_pdgId[iGenCand]
                << ", GenPartCand_pt: " << GenPartCand_pt[iGenCand]
                << ", GenJet_pt: " << GenJet_pt[iGenJet]
                << endl << flush;
              }
            }
            // Candidate type loop for each flavour
            for (int iq = 0; iq != nq; ++ iq) {
              //if (debug){cout << "iqloop" << endl;}
              bool isQ = (abs(GenJet_partonFlavour[iGenJet]) == vq2[iq]);
              //if (debug && isQ) { cout << "isQ=" << isQ << ", GenJet_partonFlavour[i]=" << GenJet_partonFlavour[iGenJet] << ", vq2[iq]=" << vq2[iq] << endl; }
              //bool isQ = (abs(GenJet_partonFlavour[i]) == mq[cq]);
              for (int ix = 0; ix != nxvar && isQ; ++ ix) {
                const char *cx = vxvar[ix].c_str();
                const char *cq = vq[iq].c_str();
                double x(0);
                const char *hname_all = Form("h_all_%s_vs_%s", cq, cx);
                if (vxvar[ix] == "ptjet") { //swap if doesn't work
                  x = GenJet_pt[iGenJet]; 
                  mh[hname_all]->Fill(x, w);
                }
                if (vxvar[ix] == "ptcand" && ptjet > 80 && ptjet < 100) {
                  x = GenPartCand_pt[iGenCand];
                  mh[hname_all]->Fill(x, w);
                }
                if (iLeadGenCand == iGenCand && vxvar[ix] == "ptlead" && p4jet.Pt()>80 && p4jet.Pt()<100) { //swap if doesn't work
                  x = ptlead; 
                  mh[hname_all]->Fill(x, w);
                }
                bool isId(false);
                //if (debug){cout << "ixloop" << endl;}
                for (std::size_t id = 0; id < vpid3.size() && !isId; ++id) { // Iterate over vectors
                  bool isMatchFound = false;
                  for (std::size_t elem = 0; elem < vpid3[id].size() && !isMatchFound; ++elem) { // Iterate over elements in a vector
                    //if (vpid3[id][elem] == 0) continue; // Skip zeros
                    if (GenPartCand_pdgId[iGenCand] == vpid3[id][elem]) {
                      isMatchFound = true;
                    }
                  }
                  if (isMatchFound) {
                    isId = true;
                    const char *pid = vpid[id].c_str();
                    const char *hname = Form("h_%s_%s_vs_%s", pid, cq, cx);
                    if (vxvar[ix] == "ptjet") {
                      x = GenJet_pt[iGenJet];
                      mh[hname]->Fill(x, w);
                    }
                    if (vxvar[ix] == "ptcand" && ptjet > 80 && ptjet < 100) {
                      x = GenPartCand_pt[iGenCand];
                      mh[hname]->Fill(x, w);
                    } // No need to check other vectors if a match is already found for this particle
                    if (iLeadGenCand == iGenCand && vxvar[ix] == "ptlead" && p4jet.Pt()>80 && p4jet.Pt()<100) {
                      x = ptlead;
                      mh[hname]->Fill(x, w);
                    }
                  } // if match found
                } // for id
              } // for ix
            } // for iq

            if (iLeadSCand == iGenCand && ptslead > 0 && ptjet>80 && ptjet<100) {
              double x(0);
              const char *hsname_all = Form("h_all_%s_vs_ptslead", "s");
              x = ptslead; 
              mh[hsname_all]->Fill(x, w);
            
              bool isSId(false);
              for (std::size_t idS = 0; idS < vpid3S.size() && !isSId; ++idS) { // Iterate over vectors
                bool isSMatchFound = false;
                for (std::size_t elem = 0; elem < vpid3S[idS].size() && !isSMatchFound; ++elem) { // Iterate over elements in a vector
                  //if (vpid3[id][elem] == 0) continue; // Skip zeros
                  if (GenPartCand_pdgId[iGenCand] == vpid3S[idS][elem]) {
                    isSMatchFound = true;
                  }
                }
                if (isSMatchFound) {
                  isSId = true;
                  const char *pidS = vpidS[idS].c_str();
                  const char *hsname = Form("h_%s_s_vs_ptslead", pidS);
                  x = ptslead;
                  mh[hsname]->Fill(x, w);
                }
              }
            }

            // Particle type loop for each flavour
            for (int ic = 0; ic != nc ; ++ ic) {
              bool isT = visc[ic];
              for (int iq = 0; iq != nq && isT; ++ iq) {
                if (debug){cout << "iqloop" << iq << endl << flush;}
                bool isQ = (abs(GenJet_partonFlavour[iGenJet]) == vq2[iq]);
                //bool isQ = (abs(GenJet_partonFlavour[i]) == mq[cq]);
                for (int ix = 0; ix != nxvar && isQ; ++ ix) {
                  if (debug){cout << "ixloop" << ix << endl << flush;}
                  const char *cx = vxvar[ix].c_str();
                  const char *cq = vq[iq].c_str();
                  const char *cc = vc[ic].c_str();
                  double x(0);
                  const char *hname = Form("h_%s_%s_vs_%s", cq, cc, cx);
                  if (vxvar[ix] == "ptjet") { //swap if doesn't work
                    if (debug){cout << "ptjet " << hname << endl << flush;}
                    x = GenJet_pt[iGenJet]; 
                    TH1D*h = mhc[hname]; assert(h);
                    h->Fill(x, w);
                  }
                  if (vxvar[ix] == "ptcand" && ptjet > 80 && ptjet < 100) {
                    if (debug){cout << "ptcand" << endl;}
                    x = GenPartCand_pt[iGenCand];
                    mhc[hname]->Fill(x, w);
                  }
                  if (iLeadGenCand == iGenCand && vxvar[ix] == "ptlead" && p4jet.Pt()>80 && p4jet.Pt()<100) {
                      x = ptlead;
                      mhc[hname]->Fill(x, w);
                    }
                } // for ix
              } // for iq
            } // for ic

            if (iLeadSCand == iGenCand && ptslead > 0 && ptjet > 80 && ptjet < 100) {
              double x(0);
              for (int ic = 0; ic != nc ; ++ ic) {
                bool isT = visc[ic];
                if (isT){
                  const char *cc = vc[ic].c_str();
                  const char *hsname = Form("h_s_%s_vs_ptslead",cc);
                  x = ptslead;
                  mhc[hsname]->Fill(x, w);
                }
              }
            }  
          } //eta loop
        } // i = iGenJet
      } // candloop for j



      // Leading vs jet energy for each flavour
      if (fabs(GenJet_eta[i]) < 1.3 && ptjet > 80 && ptjet < 100) {
        for (int iq = 0; iq != nq; ++ iq) {
          if (debug){cout << "iqloop" << iq << endl << flush;}
          bool isQ = (abs(GenJet_partonFlavour[i]) == vq2[iq]);
          //bool isQ = (abs(GenJet_partonFlavour[i]) == mq[cq]);
          for (int ix = 0; ix != nxvar && isQ; ++ ix) {
            if (debug){cout << "ixloop" << ix << endl << flush;}
            const char *cx = vxvar[ix].c_str();
            const char *cq = vq[iq].c_str();
            double x(0);
            const char *hname = Form("h_elead/ejet_%s_vs_%s", cq, cx);
            if (vxvar[ix] == "ptjet") { //swap if doesn't work
              if (debug){cout << "ptjet " << hname << endl << flush;}
              x = elead/ejet; 
              TH1D*h = mhe[hname]; assert(h);
              h->Fill(x, w);
            } //ptjet
          } // ix
        } // iq
      } // barrel

      // Tprofiles for particle types leading vs jet
      if (fabs(GenJet_eta[i] < 1.3) && esum > 0){
	      pu0->Fill(GenJet_pt[i], esum ? 1 : 0);
	      //for (int id = 0; id != npid; ++ id) {
        //bool isId = (GenPartCand_pdgId[iGenCand] == vpid3[i]);
        //bool isId = (GenPartCand_pdgId[iGenCand] == mid[pid]);
        for (int iq = 0; iq != nq ; ++ iq) { // &&isId
          bool isQ = (abs(GenJet_partonFlavour[i]) == vq2[iq]);
          //bool isQ = (abs(GenJet_partonFlavour[i]) == mq[cq]);
          for (int ix = 0; ix != nxvar && isQ; ++ ix) {
            if (vxvar[ix] == "ptjet"){
              //const char *pid = vpid[id].c_str();
              const char *cx = vxvar[ix].c_str();
              const char *cq = vq[iq].c_str(); 
              double x = GenJet_pt[i]; 
              for (int iy = 0; iy != nyvar; ++ iy) {
                double y(0);
                if (iy == 3){y = esumch/esum;} //for fhc 
                if (iy == 4){y = esumnh/esum;} //for fhn
                if (iy == 5){y = esumne/esum;} //for fhe
                if (iy == 0){y = eleadch/elead;} //for flc 
                if (iy == 1){y = eleadnh/elead;} //for fln
                if (iy == 2){y = eleadne/elead;} //for fle
                const char *cy = vyvar[iy].c_str();
                const char *pname = Form("h_%s_%s_%s_vs_%s", "all", cq, cy, cx); //pid
                mp[pname]->Fill(x, y, w);
              } // iy
            } // ptjet
          } // for ix
        } // for iq
      } // for i
      //} //jetloop for i

      // Loop over two leading jets
      //for (int i = 0; i != min(nGenJet,2); ++i) {

      p4jet.SetPtEtaPhiM(GenJet_pt[i],GenJet_eta[i],GenJet_phi[i],
        GenJet_mass[i]);
    
      bool pass2 = (fabs(p4jet.Eta())<1.3);
      if (!pass2) continue;
    
      // Loop over matching reco jets
      for (int j = 0; j != nJet; ++j) {

      if (fabs(Jet_eta[j]) < 2.5){
        if (Jet_pt[j] > 30) {
          if(abs(Jet_partonFlavour[j]) == 4 && Jet_hadronFlavour[j] == 4){
            // C v UDSG
            h_c_ctag_cvl_df30->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_ctag_cvl_pn30->Fill(Jet_btagPNetCvL[j]);

            // C v B
            h_c_ctag_cvb_df30->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_ctag_cvb_pn30->Fill(Jet_btagPNetCvB[j]);

            // B v C
            h_b_ctag_b_df30->Fill(Jet_btagDeepFlavB[j]);
            h_b_ctag_b_pn30->Fill(Jet_btagPNetB[j]);

          }
          if((abs(Jet_partonFlavour[j]) == 1 || abs(Jet_partonFlavour[j]) == 2 
          || abs(Jet_partonFlavour[j]) == 3 || abs(Jet_partonFlavour[j]) == 21)){
            // C v UDSG
            h_c_udstag_cvl_df30->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_udstag_cvl_pn30->Fill(Jet_btagPNetCvL[j]);

            // B v UDSG
            h_b_udstag_b_df30->Fill(Jet_btagDeepFlavB[j]);
            h_b_udstag_b_pn30->Fill(Jet_btagPNetB[j]);
          }


          if(abs(Jet_partonFlavour[j]) == 5 && Jet_hadronFlavour[j] == 5){
            // B v UDSG & C
            h_b_btag_b_df30->Fill(Jet_btagDeepFlavB[j]);
            h_b_btag_b_pn30->Fill(Jet_btagPNetB[j]);

            // C v B
            h_c_btag_cvb_df30->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_btag_cvb_pn30->Fill(Jet_btagPNetCvB[j]);
          }

        }
        if (Jet_pt[j] > 60) {
          if(abs(Jet_partonFlavour[j]) == 4 && Jet_hadronFlavour[j] == 4){
            // C v UDSG
            h_c_ctag_cvl_df60->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_ctag_cvl_pn60->Fill(Jet_btagPNetCvL[j]);

            // C v B
            h_c_ctag_cvb_df60->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_ctag_cvb_pn60->Fill(Jet_btagPNetCvB[j]);

            // B v C
            h_b_ctag_b_df60->Fill(Jet_btagDeepFlavB[j]);
            h_b_ctag_b_pn60->Fill(Jet_btagPNetB[j]);

          }
          if((abs(Jet_partonFlavour[j]) == 1 || abs(Jet_partonFlavour[j]) == 2 
          || abs(Jet_partonFlavour[j]) == 3 || abs(Jet_partonFlavour[j]) == 21)){
            // C v UDSG
            h_c_udstag_cvl_df60->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_udstag_cvl_pn60->Fill(Jet_btagPNetCvL[j]);

            // B v UDSG
            h_b_udstag_b_df60->Fill(Jet_btagDeepFlavB[j]);
            h_b_udstag_b_pn60->Fill(Jet_btagPNetB[j]);
          }


          // B v UDSG & B v C
          if(abs(Jet_partonFlavour[j]) == 5 && Jet_hadronFlavour[j] == 5){
            // B v UDSG & C
            h_b_btag_b_df60->Fill(Jet_btagDeepFlavB[j]);
            h_b_btag_b_pn60->Fill(Jet_btagPNetB[j]);

            // C v B
            h_c_btag_cvb_df60->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_btag_cvb_pn60->Fill(Jet_btagPNetCvB[j]);
          }
        }
        if (Jet_pt[j] > 80) {
                    if(abs(Jet_partonFlavour[j]) == 4 && Jet_hadronFlavour[j] == 4){
            // C v UDSG
            h_c_ctag_cvl_df80->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_ctag_cvl_pn80->Fill(Jet_btagPNetCvL[j]);

            // C v B
            h_c_ctag_cvb_df80->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_ctag_cvb_pn80->Fill(Jet_btagPNetCvB[j]);

            // B v C
            h_b_ctag_b_df80->Fill(Jet_btagDeepFlavB[j]);
            h_b_ctag_b_pn80->Fill(Jet_btagPNetB[j]);

          }
          if((abs(Jet_partonFlavour[j]) == 1 || abs(Jet_partonFlavour[j]) == 2 
          || abs(Jet_partonFlavour[j]) == 3 || abs(Jet_partonFlavour[j]) == 21)){
            // C v UDSG
            h_c_udstag_cvl_df80->Fill(Jet_btagDeepFlavCvL[j]);
            h_c_udstag_cvl_pn80->Fill(Jet_btagPNetCvL[j]);

            // B v UDSG
            h_b_udstag_b_df80->Fill(Jet_btagDeepFlavB[j]);
            h_b_udstag_b_pn80->Fill(Jet_btagPNetB[j]);
          }


          if(abs(Jet_partonFlavour[j]) == 5 && Jet_hadronFlavour[j] == 5){
            // B v UDSG & C
            h_b_btag_b_df80->Fill(Jet_btagDeepFlavB[j]);
            h_b_btag_b_pn80->Fill(Jet_btagPNetB[j]);

            // C v B
            h_c_btag_cvb_df80->Fill(Jet_btagDeepFlavCvB[j]);
            h_c_btag_cvb_pn80->Fill(Jet_btagPNetCvB[j]);
          }
        }
      }

      
        p4reco.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],
            Jet_mass[j]);

        if (p4jet.DeltaR(p4reco)<0.2) {
          Short_t f = GenJet_partonFlavour[i];
          double pt = p4jet.Pt();

          // Raw response
          double r = p4reco.Pt()*(1-Jet_rawFactor[j]) / pt;
          if (abs(f)==1|| abs(f)==2)
            prq->Fill(pt, r);
          if (abs(f)==1) pru->Fill(pt, r);
          if (abs(f)==2) prd->Fill(pt, r);
          if (abs(f)==3) prs->Fill(pt, r);
          if (abs(f)==4) prc->Fill(pt, r);
          if (abs(f)==5) prb->Fill(pt, r);
          if (f==21) prg->Fill(pt, r);

          // Corrected response
          double c = p4reco.Pt() / pt;
          if (abs(f)==1|| abs(f)==2)
            pcq->Fill(pt, c);
          if (abs(f)==1) pcu->Fill(pt, c);
          if (abs(f)==2) pcd->Fill(pt, c);
          if (abs(f)==3) pcs->Fill(pt, c);
          if (abs(f)==4) pcc->Fill(pt, c);
          if (abs(f)==5) pcb->Fill(pt, c);
          if (f==21) pcg->Fill(pt, c);  

        } // dR<0.2
      } // for j
    } // for jetloop i
  } // for jentry
  fout->Write();
  fout->Close();
  exit(0);
}
