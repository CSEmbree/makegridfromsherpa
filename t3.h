//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 17 12:23:28 2013 by ROOT version 5.34/03
// from TTree t3/Reconst ntuple
// found on file: sherpatopinputs/NTuple_B-like.root
//////////////////////////////////////////////////////////

#ifndef t3_h
#define t3_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class t3 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           id;
   Int_t           nparticle;
   Float_t         px[2];   //[nparticle]
   Float_t         py[2];   //[nparticle]
   Float_t         pz[2];   //[nparticle]
   Float_t         E[2];   //[nparticle]
   Double_t        alphas;
   Int_t           kf[2];   //[nparticle]
   Double_t        weight;
   Double_t        weight2;
   Double_t        me_wgt;
   Double_t        me_wgt2;
   Double_t        x1;
   Double_t        x2;
   Double_t        x1p;
   Double_t        x2p;
   Int_t           id1;
   Int_t           id2;
   Double_t        fac_scale;
   Double_t        ren_scale;
   Int_t           nuwgt;
   Double_t        usr_wgts[18];   //[nuwgt]
   Char_t          alphasPower;
   Char_t          part[2];

   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_nparticle;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_alphas;   //!
   TBranch        *b_kf;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_weight2;   //!
   TBranch        *b_me_wtg;   //!
   TBranch        *b_me_wtg2;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_x1p;   //!
   TBranch        *b_x2p;   //!
   TBranch        *b_id1;   //!
   TBranch        *b_id2;   //!
   TBranch        *b_fac_scale;   //!
   TBranch        *b_ren_scale;   //!
   TBranch        *b_nuwgt;   //!
   TBranch        *b_usr_wgts;   //!
   TBranch        *b_alphasPower;   //!
   TBranch        *b_part;   //!

   t3(TTree *tree=0);
   virtual ~t3();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t3_cxx
t3::t3(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sherpatopinputs/NTuple_B-like.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sherpatopinputs/NTuple_B-like.root");
      }
      f->GetObject("t3",tree);

   }
   Init(tree);
}

t3::~t3()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t3::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t3::LoadTree(Long64_t entry)
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

void t3::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("nparticle", &nparticle, &b_nparticle);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("alphas", &alphas, &b_alphas);
   fChain->SetBranchAddress("kf", kf, &b_kf);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("weight2", &weight2, &b_weight2);
   fChain->SetBranchAddress("me_wgt", &me_wgt, &b_me_wtg);
   fChain->SetBranchAddress("me_wgt2", &me_wgt2, &b_me_wtg2);
   fChain->SetBranchAddress("x1", &x1, &b_x1);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("x1p", &x1p, &b_x1p);
   fChain->SetBranchAddress("x2p", &x2p, &b_x2p);
   fChain->SetBranchAddress("id1", &id1, &b_id1);
   fChain->SetBranchAddress("id2", &id2, &b_id2);
   fChain->SetBranchAddress("fac_scale", &fac_scale, &b_fac_scale);
   fChain->SetBranchAddress("ren_scale", &ren_scale, &b_ren_scale);
   fChain->SetBranchAddress("nuwgt", &nuwgt, &b_nuwgt);
   fChain->SetBranchAddress("usr_wgts", usr_wgts, &b_usr_wgts);
   fChain->SetBranchAddress("alphasPower", &alphasPower, &b_alphasPower);
   fChain->SetBranchAddress("part", part, &b_part);
   Notify();
}

Bool_t t3::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t3::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t3::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t3_cxx
