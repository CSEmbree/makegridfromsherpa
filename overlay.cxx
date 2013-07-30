//overlay histos
//Cameron Embree
#include "root.h"

// ./overlay atlas2012_5fb_top_mtt_ljet_R_histos.root atlas2012_5fb_top_mtt_ljet_R_histos.root convolute_for_R htest1_R outputtop/

TFile* OpenFile(string fileName) {
    TFile *openedFile = TFile::Open(fileName.c_str(),"READ");
    if(!openedFile) { std::cout<<"Overlay::main: ERROR: Unable to open '"<<fileName<<"'"<<std::endl; exit(0); }
    else { std::cout<<"Overlay::main: Opened '"<<fileName<<"'"<<std::endl; return openedFile; }
}

int main(int argc, char** argv) {
    if(argc<5) {
        std::cout<<"Overlay::main: USAGE: Need 4 parameters. 1-rootFile 2-RootFile 3-h1 4-h2 5-outputPath(Optional)"<<std::endl; return -1; }
    std::cout<<"Overlay::main: Overlaying '"<<argv[3]<<"' ONTO '"<<argv[4]<<"'"<<std::endl;

    string inFileName1 = argv[1];
    string inFileName2 = argv[2];
    string hName1 = argv[3];
    string hName2 = argv[4];

    //allow change of output file path, default is local directory
    string outFilePath="";
    if(argc>=6) outFilePath = argv[5];

    //get histogram file names
    string outFileName = "overlay_of_"+hName1+"_onto_"+hName2+".root";
    string outPath = outFilePath+outFileName;
    
    //Read in histograms from file names
    TFile* inFile1 = OpenFile(inFileName1);
    TFile* inFile2 = OpenFile(inFileName2);
    
    //prep output file for overlay of both histograms
    TFile* outFile = new TFile(outPath.c_str(),"recreate"); //overwrite
    TCanvas *c1 = new TCanvas(outFileName.c_str(),outFileName.c_str(),600,400);
    
    
    //get histo one from inFile1 and plot it
    TH1D* h1 = (TH1D*)inFile1->Get(hName1.c_str());
    if(!h1) { std::cout<<"Overlay::main: ERROR: Unable to open file1: '"<<argv[3]<<"'(from: "<<argv[1]<<")"<<std::endl; exit(0); }
    h1->Draw(); c1->Update();

    //get histo two from inFile2 and overlay it with histo one plot
    TH1D *h2 = (TH1D*)inFile2->Get(hName2.c_str());
    if(!h2) { std::cout<<"Overlay::main: ERROR: Unable to open file2: '"<<argv[4]<<"'(from: "<<argv[2]<<")"<<std::endl; exit(0); };
    h2->Draw("same"); c1->Write();

    outFile->Write(); outFile->Close();
    std::cout<<"Overlay::main: Wrote output to file: "<<outPath<<std::endl;
    return 0;
}
