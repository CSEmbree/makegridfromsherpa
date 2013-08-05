//#include "root.h"


void test() {
    int numSubProcs = 25;
    string SubProcID[] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"};


    TFile *file1 = TFile::Open("outputtop/atlas2012_5fb_top_mtt_ljet_B-ghistos.root");
    TFile *file2 = TFile::Open("atlas2012_5fb_top_mtt_ljet_B_histos.root");

    for(int iproc=0; iproc<numSubProcs; iproc++) 
    {
        std::cout<<"Checking iproc: "<<iproc<<std::endl;
        if( iproc==0 || iproc==14 || iproc==15 || iproc==16 || iproc==17 ) //sub procs allowed through
        {
            string f1HistoName = "hrefSubProcHistosDiv_"+SubProcID[iproc];
            TH1D *h1 = (TH1D*)file1->Get(f1HistoName.c_str());

            string f2HistoName = "LOsubProc-"+SubProcID[iproc]+"-convolute_for_B";
            TH1D *h2 = (TH1D*)file2->Get(f2HistoName.c_str());

            TCanvas *c1 = new TCanvas();
            c1->Divide( 2 , 1 );
            c1->cd( 1 );

            h1->Draw();
            h2->Draw("same");


            string ratioName = "href/LOconvolute--"+SubProcID[iproc];
            TH1D *htmp = (TH1D*) h1->Clone(ratioName.c_str());
            htmp->Divide(h2); //href/LOconv

            c1->cd(2);
            htmp->Draw();
            
            std::cout<<"Showing subProc: "<<iproc<<" histos "<<f1HistoName<<"/"<<f2HistoName<<std::endl;
        }
    }
}
