//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 14 12:23:55 2024 by ROOT version 6.32.06
// from TTree tree/tree
// found on file: Muo16_MC.root
//////////////////////////////////////////////////////////

#ifndef TagandProbe_h
#define TagandProbe_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class TagandProbe {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         fitProb;
   Float_t         fitTopMass;
   Float_t         recoRbq;
   Float_t         recoWMass;
   Float_t         recoWEta;
   Float_t         recoWPhi;
   Float_t         pt1;
   Float_t         pt2;
   Float_t         eta1;
   Float_t         eta2;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         m1;
   Float_t         m2;
   Float_t         sc1;
   Float_t         sc2;
   Float_t         muf1;
   Float_t         muf2;
   Float_t         elf1;
   Float_t         elf2;
   Float_t         nhf1;
   Float_t         nhf2;
   Float_t         qgl1;
   Float_t         qgl2;
   Float_t         qgax1;
   Float_t         qgax2;
   Float_t         qgptd1;
   Float_t         qgptd2;
   Int_t           qgmult1;
   Int_t           qgmult2;
   Float_t         btag1;
   Float_t         btag2;
   Float_t         bleptag1;
   Float_t         bleptag2;
   Float_t         ctag1;
   Float_t         ctag2;
   Float_t         udstag1;
   Float_t         udstag2;
   Float_t         gtag1;
   Float_t         gtag2;
   Int_t           flav1;
   Int_t           flav2;
   Int_t           ntrk1;
   Int_t           ntrk2;
   Int_t           nputrk1;
   Int_t           nputrk2;
   Float_t         gen_pt1;
   Float_t         gen_pt2;
   Float_t         gen_eta1;
   Float_t         gen_eta2;
   Float_t         gen_phi1;
   Float_t         gen_phi2;
   Float_t         gen_m1;
   Float_t         gen_m2;
   Float_t         bpt1;
   Float_t         bpt2;
   Float_t         beta1;
   Float_t         beta2;
   Float_t         bphi1;
   Float_t         bphi2;
   Float_t         bm1;
   Float_t         bm2;
   Float_t         bsc1;
   Float_t         bsc2;
   Float_t         bmuf1;
   Float_t         bmuf2;
   Float_t         belf1;
   Float_t         belf2;
   Float_t         bnhf1;
   Float_t         bnhf2;
   Float_t         bqgl1;
   Float_t         bqgl2;
   Float_t         bqgax1;
   Float_t         bqgax2;
   Float_t         bqgptd1;
   Float_t         bqgptd2;
   Int_t           bqgmult1;
   Int_t           bqgmult2;
   Float_t         bbtag1;
   Float_t         bbtag2;
   Float_t         bbleptag1;
   Float_t         bbleptag2;
   Float_t         bctag1;
   Float_t         bctag2;
   Float_t         budstag1;
   Float_t         budstag2;
   Float_t         bgtag1;
   Float_t         bgtag2;
   Int_t           bflav1;
   Int_t           bflav2;
   Int_t           bntrk1;
   Int_t           bntrk2;
   Int_t           bnputrk1;
   Int_t           bnputrk2;
   Float_t         gen_bpt1;
   Float_t         gen_bpt2;
   Float_t         gen_beta1;
   Float_t         gen_beta2;
   Float_t         gen_bphi1;
   Float_t         gen_bphi2;
   Float_t         gen_bm1;
   Float_t         gen_bm2;
   Int_t           gen_bLeadId1;
   Int_t           gen_bLeadId2;
   Int_t           gen_bFlags1;
   Int_t           gen_bFlags2;
   Float_t         gen_bXB1;
   Float_t         gen_bXB2;
   Float_t         gen_bNuCorr1;
   Float_t         gen_bNuCorr2;
   Float_t         gpt1;
   Float_t         geta1;
   Float_t         gphi1;
   Float_t         gm1;
   Float_t         gsc1;
   Float_t         gmuf1;
   Float_t         gelf1;
   Float_t         gnhf1;
   Float_t         gqgl1;
   Float_t         gqgax1;
   Float_t         gqgptd1;
   Int_t           gqgmult1;
   Float_t         gbtag1;
   Float_t         gbleptag1;
   Float_t         gctag1;
   Float_t         gudstag1;
   Float_t         ggtag1;
   Int_t           gflav1;
   Int_t           gntrk1;
   Int_t           gnputrk1;
   Float_t         gen_gpt1;
   Float_t         gen_geta1;
   Float_t         gen_gphi1;
   Float_t         gen_gm1;
   Float_t         lep_pt;
   Float_t         lep_eta;
   Float_t         lep_phi;
   Float_t         lep_m;
   Float_t         hadTopPt;
   Float_t         hadTopPtG;
   Float_t         lepTopPtG;
   Float_t         weight;
   Float_t         tptweight;
   Float_t         mcweight;
   Float_t         puweight;
   Float_t         btagweight;
   Float_t         idweight;
   Float_t         isoweight;
   Float_t         recoweight;
   Float_t         trgweight;
   Int_t           run;
   Int_t           lumiBlock;
   Int_t           event;
   Float_t         pfRho;
   Int_t           TruePU;
   Int_t           NPrVtx;
   Int_t           NPrVtxGood;
   Float_t         METx;
   Float_t         METy;
   Int_t           SIdx;
   Int_t           PIdx;
   Int_t           ComboType;
   vector<float>   *PSWgts;

   // List of branches
   TBranch        *b_fitProb;   //!
   TBranch        *b_fitTopMass;   //!
   TBranch        *b_recoRbq;   //!
   TBranch        *b_recoWMass;   //!
   TBranch        *b_recoWEta;   //!
   TBranch        *b_recoWPhi;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_m1;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_sc1;   //!
   TBranch        *b_sc2;   //!
   TBranch        *b_muf1;   //!
   TBranch        *b_muf2;   //!
   TBranch        *b_elf1;   //!
   TBranch        *b_elf2;   //!
   TBranch        *b_nhf1;   //!
   TBranch        *b_nhf2;   //!
   TBranch        *b_qgl1;   //!
   TBranch        *b_qgl2;   //!
   TBranch        *b_qgax1;   //!
   TBranch        *b_qgax2;   //!
   TBranch        *b_qgptd1;   //!
   TBranch        *b_qgptd2;   //!
   TBranch        *b_qgmult1;   //!
   TBranch        *b_qgmult2;   //!
   TBranch        *b_btag1;   //!
   TBranch        *b_btag2;   //!
   TBranch        *b_bleptag1;   //!
   TBranch        *b_bleptag2;   //!
   TBranch        *b_ctag1;   //!
   TBranch        *b_ctag2;   //!
   TBranch        *b_udstag1;   //!
   TBranch        *b_udstag2;   //!
   TBranch        *b_gtag1;   //!
   TBranch        *b_gtag2;   //!
   TBranch        *b_flav1;   //!
   TBranch        *b_flav2;   //!
   TBranch        *b_ntrk1;   //!
   TBranch        *b_ntrk2;   //!
   TBranch        *b_nputrk1;   //!
   TBranch        *b_nputrk2;   //!
   TBranch        *b_gen_pt1;   //!
   TBranch        *b_gen_pt2;   //!
   TBranch        *b_gen_eta1;   //!
   TBranch        *b_gen_eta2;   //!
   TBranch        *b_gen_phi1;   //!
   TBranch        *b_gen_phi2;   //!
   TBranch        *b_gen_m1;   //!
   TBranch        *b_gen_m2;   //!
   TBranch        *b_bpt1;   //!
   TBranch        *b_bpt2;   //!
   TBranch        *b_beta1;   //!
   TBranch        *b_beta2;   //!
   TBranch        *b_bphi1;   //!
   TBranch        *b_bphi2;   //!
   TBranch        *b_bm1;   //!
   TBranch        *b_bm2;   //!
   TBranch        *b_bsc1;   //!
   TBranch        *b_bsc2;   //!
   TBranch        *b_bmuf1;   //!
   TBranch        *b_bmuf2;   //!
   TBranch        *b_belf1;   //!
   TBranch        *b_belf2;   //!
   TBranch        *b_bnhf1;   //!
   TBranch        *b_bnhf2;   //!
   TBranch        *b_bqgl1;   //!
   TBranch        *b_bqgl2;   //!
   TBranch        *b_bqgax1;   //!
   TBranch        *b_bqgax2;   //!
   TBranch        *b_bqgptd1;   //!
   TBranch        *b_bqgptd2;   //!
   TBranch        *b_bqgmult1;   //!
   TBranch        *b_bqgmult2;   //!
   TBranch        *b_bbtag1;   //!
   TBranch        *b_bbtag2;   //!
   TBranch        *b_bbleptag1;   //!
   TBranch        *b_bbleptag2;   //!
   TBranch        *b_bctag1;   //!
   TBranch        *b_bctag2;   //!
   TBranch        *b_budstag1;   //!
   TBranch        *b_budstag2;   //!
   TBranch        *b_bgtag1;   //!
   TBranch        *b_bgtag2;   //!
   TBranch        *b_bflav1;   //!
   TBranch        *b_bflav2;   //!
   TBranch        *b_bntrk1;   //!
   TBranch        *b_bntrk2;   //!
   TBranch        *b_bnputrk1;   //!
   TBranch        *b_bnputrk2;   //!
   TBranch        *b_gen_bpt1;   //!
   TBranch        *b_gen_bpt2;   //!
   TBranch        *b_gen_beta1;   //!
   TBranch        *b_gen_beta2;   //!
   TBranch        *b_gen_bphi1;   //!
   TBranch        *b_gen_bphi2;   //!
   TBranch        *b_gen_bm1;   //!
   TBranch        *b_gen_bm2;   //!
   TBranch        *b_gen_bLeadId1;   //!
   TBranch        *b_gen_bLeadId2;   //!
   TBranch        *b_gen_bFlags1;   //!
   TBranch        *b_gen_bFlags2;   //!
   TBranch        *b_gen_bXB1;   //!
   TBranch        *b_gen_bXB2;   //!
   TBranch        *b_gen_bnucorr1;   //!
   TBranch        *b_gen_bnucorr2;   //!
   TBranch        *b_gpt1;   //!
   TBranch        *b_geta1;   //!
   TBranch        *b_gphi1;   //!
   TBranch        *b_gm1;   //!
   TBranch        *b_gsc1;   //!
   TBranch        *b_gmuf1;   //!
   TBranch        *b_gelf1;   //!
   TBranch        *b_gnhf1;   //!
   TBranch        *b_gqgl1;   //!
   TBranch        *b_gqgax1;   //!
   TBranch        *b_gqgptd1;   //!
   TBranch        *b_gqgmult1;   //!
   TBranch        *b_gbtag1;   //!
   TBranch        *b_gbleptag1;   //!
   TBranch        *b_gctag1;   //!
   TBranch        *b_gudstag1;   //!
   TBranch        *b_ggtag1;   //!
   TBranch        *b_gflav1;   //!
   TBranch        *b_gntrk1;   //!
   TBranch        *b_gnputrk1;   //!
   TBranch        *b_gen_gpt1;   //!
   TBranch        *b_gen_geta1;   //!
   TBranch        *b_gen_gphi1;   //!
   TBranch        *b_gen_gm1;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_m;   //!
   TBranch        *b_hadTopPt;   //!
   TBranch        *b_hadTopPtGen;   //!
   TBranch        *b_lepTopPtGen;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_tptweight;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_btagweight;   //!
   TBranch        *b_idweight;   //!
   TBranch        *b_isoweight;   //!
   TBranch        *b_recoweight;   //!
   TBranch        *b_trgweight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_pfRho;   //!
   TBranch        *b_TruePU;   //!
   TBranch        *b_NPrVtx;   //!
   TBranch        *b_NPrVtxGood;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_SampleIdx;   //!
   TBranch        *b_PermIdx;   //!
   TBranch        *b_ComboType;   //!
   TBranch        *b_PSWgts;   //!

   TagandProbe(TTree *tree=0);
   virtual ~TagandProbe();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TagandProbe_cxx
TagandProbe::TagandProbe(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Muo16_MC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Muo16_MC.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

TagandProbe::~TagandProbe()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TagandProbe::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TagandProbe::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TagandProbe::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PSWgts = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fitProb", &fitProb, &b_fitProb);
   fChain->SetBranchAddress("fitTopMass", &fitTopMass, &b_fitTopMass);
   fChain->SetBranchAddress("recoRbq", &recoRbq, &b_recoRbq);
   fChain->SetBranchAddress("recoWMass", &recoWMass, &b_recoWMass);
   fChain->SetBranchAddress("recoWEta", &recoWEta, &b_recoWEta);
   fChain->SetBranchAddress("recoWPhi", &recoWPhi, &b_recoWPhi);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("m1", &m1, &b_m1);
   fChain->SetBranchAddress("m2", &m2, &b_m2);
   fChain->SetBranchAddress("sc1", &sc1, &b_sc1);
   fChain->SetBranchAddress("sc2", &sc2, &b_sc2);
   fChain->SetBranchAddress("muf1", &muf1, &b_muf1);
   fChain->SetBranchAddress("muf2", &muf2, &b_muf2);
   fChain->SetBranchAddress("elf1", &elf1, &b_elf1);
   fChain->SetBranchAddress("elf2", &elf2, &b_elf2);
   fChain->SetBranchAddress("nhf1", &nhf1, &b_nhf1);
   fChain->SetBranchAddress("nhf2", &nhf2, &b_nhf2);
   fChain->SetBranchAddress("qgl1", &qgl1, &b_qgl1);
   fChain->SetBranchAddress("qgl2", &qgl2, &b_qgl2);
   fChain->SetBranchAddress("qgax1", &qgax1, &b_qgax1);
   fChain->SetBranchAddress("qgax2", &qgax2, &b_qgax2);
   fChain->SetBranchAddress("qgptd1", &qgptd1, &b_qgptd1);
   fChain->SetBranchAddress("qgptd2", &qgptd2, &b_qgptd2);
   fChain->SetBranchAddress("qgmult1", &qgmult1, &b_qgmult1);
   fChain->SetBranchAddress("qgmult2", &qgmult2, &b_qgmult2);
   fChain->SetBranchAddress("btag1", &btag1, &b_btag1);
   fChain->SetBranchAddress("btag2", &btag2, &b_btag2);
   fChain->SetBranchAddress("bleptag1", &bleptag1, &b_bleptag1);
   fChain->SetBranchAddress("bleptag2", &bleptag2, &b_bleptag2);
   fChain->SetBranchAddress("ctag1", &ctag1, &b_ctag1);
   fChain->SetBranchAddress("ctag2", &ctag2, &b_ctag2);
   fChain->SetBranchAddress("udstag1", &udstag1, &b_udstag1);
   fChain->SetBranchAddress("udstag2", &udstag2, &b_udstag2);
   fChain->SetBranchAddress("gtag1", &gtag1, &b_gtag1);
   fChain->SetBranchAddress("gtag2", &gtag2, &b_gtag2);
   fChain->SetBranchAddress("flav1", &flav1, &b_flav1);
   fChain->SetBranchAddress("flav2", &flav2, &b_flav2);
   fChain->SetBranchAddress("ntrk1", &ntrk1, &b_ntrk1);
   fChain->SetBranchAddress("ntrk2", &ntrk2, &b_ntrk2);
   fChain->SetBranchAddress("nputrk1", &nputrk1, &b_nputrk1);
   fChain->SetBranchAddress("nputrk2", &nputrk2, &b_nputrk2);
   fChain->SetBranchAddress("gen_pt1", &gen_pt1, &b_gen_pt1);
   fChain->SetBranchAddress("gen_pt2", &gen_pt2, &b_gen_pt2);
   fChain->SetBranchAddress("gen_eta1", &gen_eta1, &b_gen_eta1);
   fChain->SetBranchAddress("gen_eta2", &gen_eta2, &b_gen_eta2);
   fChain->SetBranchAddress("gen_phi1", &gen_phi1, &b_gen_phi1);
   fChain->SetBranchAddress("gen_phi2", &gen_phi2, &b_gen_phi2);
   fChain->SetBranchAddress("gen_m1", &gen_m1, &b_gen_m1);
   fChain->SetBranchAddress("gen_m2", &gen_m2, &b_gen_m2);
   fChain->SetBranchAddress("bpt1", &bpt1, &b_bpt1);
   fChain->SetBranchAddress("bpt2", &bpt2, &b_bpt2);
   fChain->SetBranchAddress("beta1", &beta1, &b_beta1);
   fChain->SetBranchAddress("beta2", &beta2, &b_beta2);
   fChain->SetBranchAddress("bphi1", &bphi1, &b_bphi1);
   fChain->SetBranchAddress("bphi2", &bphi2, &b_bphi2);
   fChain->SetBranchAddress("bm1", &bm1, &b_bm1);
   fChain->SetBranchAddress("bm2", &bm2, &b_bm2);
   fChain->SetBranchAddress("bsc1", &bsc1, &b_bsc1);
   fChain->SetBranchAddress("bsc2", &bsc2, &b_bsc2);
   fChain->SetBranchAddress("bmuf1", &bmuf1, &b_bmuf1);
   fChain->SetBranchAddress("bmuf2", &bmuf2, &b_bmuf2);
   fChain->SetBranchAddress("belf1", &belf1, &b_belf1);
   fChain->SetBranchAddress("belf2", &belf2, &b_belf2);
   fChain->SetBranchAddress("bnhf1", &bnhf1, &b_bnhf1);
   fChain->SetBranchAddress("bnhf2", &bnhf2, &b_bnhf2);
   fChain->SetBranchAddress("bqgl1", &bqgl1, &b_bqgl1);
   fChain->SetBranchAddress("bqgl2", &bqgl2, &b_bqgl2);
   fChain->SetBranchAddress("bqgax1", &bqgax1, &b_bqgax1);
   fChain->SetBranchAddress("bqgax2", &bqgax2, &b_bqgax2);
   fChain->SetBranchAddress("bqgptd1", &bqgptd1, &b_bqgptd1);
   fChain->SetBranchAddress("bqgptd2", &bqgptd2, &b_bqgptd2);
   fChain->SetBranchAddress("bqgmult1", &bqgmult1, &b_bqgmult1);
   fChain->SetBranchAddress("bqgmult2", &bqgmult2, &b_bqgmult2);
   fChain->SetBranchAddress("bbtag1", &bbtag1, &b_bbtag1);
   fChain->SetBranchAddress("bbtag2", &bbtag2, &b_bbtag2);
   fChain->SetBranchAddress("bbleptag1", &bbleptag1, &b_bbleptag1);
   fChain->SetBranchAddress("bbleptag2", &bbleptag2, &b_bbleptag2);
   fChain->SetBranchAddress("bctag1", &bctag1, &b_bctag1);
   fChain->SetBranchAddress("bctag2", &bctag2, &b_bctag2);
   fChain->SetBranchAddress("budstag1", &budstag1, &b_budstag1);
   fChain->SetBranchAddress("budstag2", &budstag2, &b_budstag2);
   fChain->SetBranchAddress("bgtag1", &bgtag1, &b_bgtag1);
   fChain->SetBranchAddress("bgtag2", &bgtag2, &b_bgtag2);
   fChain->SetBranchAddress("bflav1", &bflav1, &b_bflav1);
   fChain->SetBranchAddress("bflav2", &bflav2, &b_bflav2);
   fChain->SetBranchAddress("bntrk1", &bntrk1, &b_bntrk1);
   fChain->SetBranchAddress("bntrk2", &bntrk2, &b_bntrk2);
   fChain->SetBranchAddress("bnputrk1", &bnputrk1, &b_bnputrk1);
   fChain->SetBranchAddress("bnputrk2", &bnputrk2, &b_bnputrk2);
   fChain->SetBranchAddress("gen_bpt1", &gen_bpt1, &b_gen_bpt1);
   fChain->SetBranchAddress("gen_bpt2", &gen_bpt2, &b_gen_bpt2);
   fChain->SetBranchAddress("gen_beta1", &gen_beta1, &b_gen_beta1);
   fChain->SetBranchAddress("gen_beta2", &gen_beta2, &b_gen_beta2);
   fChain->SetBranchAddress("gen_bphi1", &gen_bphi1, &b_gen_bphi1);
   fChain->SetBranchAddress("gen_bphi2", &gen_bphi2, &b_gen_bphi2);
   fChain->SetBranchAddress("gen_bm1", &gen_bm1, &b_gen_bm1);
   fChain->SetBranchAddress("gen_bm2", &gen_bm2, &b_gen_bm2);
   fChain->SetBranchAddress("gen_bLeadId1", &gen_bLeadId1, &b_gen_bLeadId1);
   fChain->SetBranchAddress("gen_bLeadId2", &gen_bLeadId2, &b_gen_bLeadId2);
   fChain->SetBranchAddress("gen_bFlags1", &gen_bFlags1, &b_gen_bFlags1);
   fChain->SetBranchAddress("gen_bFlags2", &gen_bFlags2, &b_gen_bFlags2);
   fChain->SetBranchAddress("gen_bXB1", &gen_bXB1, &b_gen_bXB1);
   fChain->SetBranchAddress("gen_bXB2", &gen_bXB2, &b_gen_bXB2);
   fChain->SetBranchAddress("gen_bNuCorr1", &gen_bNuCorr1, &b_gen_bnucorr1);
   fChain->SetBranchAddress("gen_bNuCorr2", &gen_bNuCorr2, &b_gen_bnucorr2);
   fChain->SetBranchAddress("gpt1", &gpt1, &b_gpt1);
   fChain->SetBranchAddress("geta1", &geta1, &b_geta1);
   fChain->SetBranchAddress("gphi1", &gphi1, &b_gphi1);
   fChain->SetBranchAddress("gm1", &gm1, &b_gm1);
   fChain->SetBranchAddress("gsc1", &gsc1, &b_gsc1);
   fChain->SetBranchAddress("gmuf1", &gmuf1, &b_gmuf1);
   fChain->SetBranchAddress("gelf1", &gelf1, &b_gelf1);
   fChain->SetBranchAddress("gnhf1", &gnhf1, &b_gnhf1);
   fChain->SetBranchAddress("gqgl1", &gqgl1, &b_gqgl1);
   fChain->SetBranchAddress("gqgax1", &gqgax1, &b_gqgax1);
   fChain->SetBranchAddress("gqgptd1", &gqgptd1, &b_gqgptd1);
   fChain->SetBranchAddress("gqgmult1", &gqgmult1, &b_gqgmult1);
   fChain->SetBranchAddress("gbtag1", &gbtag1, &b_gbtag1);
   fChain->SetBranchAddress("gbleptag1", &gbleptag1, &b_gbleptag1);
   fChain->SetBranchAddress("gctag1", &gctag1, &b_gctag1);
   fChain->SetBranchAddress("gudstag1", &gudstag1, &b_gudstag1);
   fChain->SetBranchAddress("ggtag1", &ggtag1, &b_ggtag1);
   fChain->SetBranchAddress("gflav1", &gflav1, &b_gflav1);
   fChain->SetBranchAddress("gntrk1", &gntrk1, &b_gntrk1);
   fChain->SetBranchAddress("gnputrk1", &gnputrk1, &b_gnputrk1);
   fChain->SetBranchAddress("gen_gpt1", &gen_gpt1, &b_gen_gpt1);
   fChain->SetBranchAddress("gen_geta1", &gen_geta1, &b_gen_geta1);
   fChain->SetBranchAddress("gen_gphi1", &gen_gphi1, &b_gen_gphi1);
   fChain->SetBranchAddress("gen_gm1", &gen_gm1, &b_gen_gm1);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_m", &lep_m, &b_lep_m);
   fChain->SetBranchAddress("hadTopPt", &hadTopPt, &b_hadTopPt);
   fChain->SetBranchAddress("hadTopPtG", &hadTopPtG, &b_hadTopPtGen);
   fChain->SetBranchAddress("lepTopPtG", &lepTopPtG, &b_lepTopPtGen);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("tptweight", &tptweight, &b_tptweight);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("btagweight", &btagweight, &b_btagweight);
   fChain->SetBranchAddress("idweight", &idweight, &b_idweight);
   fChain->SetBranchAddress("isoweight", &isoweight, &b_isoweight);
   fChain->SetBranchAddress("recoweight", &recoweight, &b_recoweight);
   fChain->SetBranchAddress("trgweight", &trgweight, &b_trgweight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("pfRho", &pfRho, &b_pfRho);
   fChain->SetBranchAddress("TruePU", &TruePU, &b_TruePU);
   fChain->SetBranchAddress("NPrVtx", &NPrVtx, &b_NPrVtx);
   fChain->SetBranchAddress("NPrVtxGood", &NPrVtxGood, &b_NPrVtxGood);
   fChain->SetBranchAddress("METx", &METx, &b_METx);
   fChain->SetBranchAddress("METy", &METy, &b_METy);
   fChain->SetBranchAddress("SIdx", &SIdx, &b_SampleIdx);
   fChain->SetBranchAddress("PIdx", &PIdx, &b_PermIdx);
   fChain->SetBranchAddress("ComboType", &ComboType, &b_ComboType);
   fChain->SetBranchAddress("PSWgts", &PSWgts, &b_PSWgts);
   Notify();
}

bool TagandProbe::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void TagandProbe::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TagandProbe::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
if (entry) {}
   return 1;
}
#endif // #ifdef TagandProbe_cxx
