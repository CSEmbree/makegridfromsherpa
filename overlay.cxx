//overlay histos
//Cameron Embree
#include "root.h"

bool debug;

int main(int argc, char** argv) {

    if(argc<5) {
        if(debug) std::cout<<"Overlay::main: USAGE: Need at least 4 parameters. 1-File1Name, 2-File2Name, 3-h1Name, 4-h2Name"<<std::endl; return -1; }
    if(debug) std::cout<<"Overlay::main: Overlaying '"<<argv[3]<<"'(from: "<<argv[1]<<") ONTO '"<<argv[4]<<"'(from: "<<argv[2]<<")"<<std::endl;

    //allow change of output file path, default is local directory
    string outFilePath="";
    if(argc>=6) {
        outFilePath = argv[5];
        if(debug) std::cout<<"Overlay::main: Overlay path set to: '"<<outFilePath<<"'"<<std::endl;
    }

    //get histogram file names
    string inFileName1 = argv[1];
    string inFileName2 = argv[2];
    string outFileName = "overlay_of_"+string(argv[3])+"_onto_"+string(argv[4])+".root";
    string canvasName = outFilePath+outFileName;
    
    
    //Read in histograms from file names
    TFile* inFile1 = TFile::Open(inFileName1.c_str(),"READ");
    if(!inFile1) { std::cout<<"Overlay::main: ERROR: Unable to open '"<<argv[1]<<"'"<<std::endl; exit(0); }
    else { if(debug)std::cout<<"Overlay::main: Opened '"<<argv[1]<<"'"<<std::endl;}
    //inFile1->Print("all");
    
    TFile* inFile2 = TFile::Open(inFileName2.c_str(),"READ");
    if(!inFile2) { std::cout<<"Overlay::main: ERROR: Unable to open '"<<argv[2]<<"'"<<std::endl; exit(0); }
    else { if(debug)std::cout<<"Overlay::main: Opened '"<<argv[2]<<"'"<<std::endl;}
    //inFile2->Print("all");
    
    //prep output file for overlay of both histograms
    TFile* outFile = new TFile(outFileName.c_str(),"recreate"); //overwrite
    TCanvas *c1 = new TCanvas(canvasName.c_str(),canvasName.c_str(),600,400);
    
    
    //get histo one from inFile1 and plot it
    TH1D* h1 = (TH1D*)inFile1->Get(argv[3]);
    if(!h1) { std::cout<<"Overlay::main: ERROR: Unable to open file1: '"<<argv[3]<<"'(from: "<<argv[1]<<")"<<std::endl; exit(0); }
    //h1->Print("all");
    h1->Draw();
    c1->Update();


    //get histo two from inFile2 and overlay it with histo one plot
    TH1D *h2 = (TH1D*)inFile2->Get(argv[4]);
    if(!h2) { std::cout<<"Overlay::main: ERROR: Unable to open file2: '"<<argv[4]<<"'(from: "<<argv[2]<<")"<<std::endl; exit(0); };
    //h2->Print("all");
    h2->Draw("same");
    c1->Write();
    

    //write and close everything
    //inFile1->Close();
    //inFile2->Close();

    outFile->Write();
    outFile->Close();
    
    std::cout<<"Overlay::main: Wrote output to file: "<<outFileName<<std::endl;
    exit(0);
}
