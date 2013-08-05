//#include "root.h"


void test(){
    TFile *file1 = TFile::Open("outputtop/atlas2012_5fb_top_mtt_ljet_B-ghistos.root");
    TFile *file2 = TFile::Open("atlas2012_5fb_top_mtt_ljet_B_histos.root");

    TH1D *h1 = file1->Get("hreftmpLOsub0atlas2012_5fb_top_mtt_ljet_B");
    TH1D *h2 = file2->Get("LOsubProc-0-convolute_for_B");
    
    TCanvas *c1 = new TCanvas();
    c1->Divide(2,1);
    c1->cd(1);
    
    h1->Draw();
    h2->Draw("same");
    
    TH1D *htmp = (TH1D*) h1->Clone("htmp");
    htmp->Divide(h2);
    
    c1->cd(2);
    htmp->Draw();
}
