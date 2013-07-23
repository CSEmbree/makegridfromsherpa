//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec 16 13:27:36 2012 by ROOT version 5.34/00
// from TTree t3/Reconst ntuple
// found on file: NTuple_R-like.root
//////////////////////////////////////////////////////////

#ifndef sherpatree_h
#define sherpatree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class sherpatree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           id;
   Int_t           nparticle;
   Float_t         px[4];   //[nparticle]
   Float_t         py[4];   //[nparticle]
   Float_t         pz[4];   //[nparticle]
   Float_t         E[4];   //[nparticle]
   Double_t        alphas;
   Int_t           kf[4];   //[nparticle]
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
   Double_t        usr_wgts[1];   //[nuwgt]
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

   sherpatree(TTree *tree=0);
};

#endif

