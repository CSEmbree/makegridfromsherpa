#include <iostream>
#include <string>
#include <vector>
#include <climits>

#include "root.h"
#include "t3.h"
//#include "utils.h"
#include "normtmp.h"
#include "appl_grid/appl_grid.h"
#include "fjClustering.h"

#include "MyEvent.h"
#include "MyData.h"
#include "MyFrameData.h"
#include "MyGrid.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>



//#include "LHAPDF.h"                 //added for convolute
#include "LHAPDF/LHAPDF.h"
//#include <LHAPDF/LHAPDF.h>
#include "VariableDefinitions.h"    //added for convolute

#define PI 3.141592653589793238462

/*
 * EXAMPLE Execution:
 *                   ./makegridfromsherpa
 *                   ./makegridfromsherpa <filename> <numevents>
 *                   ./makegridfromsherpa atlas2012_top.txt 1000000
 */

//enum ntup_types {i_R=0,i_B,i_RB};           //indexes for different histograms depending on "NTuple_*-like" in htest1
//string ntup_names[]= {"_R","_B","_RthenB"}; //names of each htest1 index type, which comes from the NTuples they run over

//changed to do B-like Ntups first
enum ntup_types {i_B=0,i_R,i_RB};           //indexes for different histograms depending on "NTuple_*-like" in htest1
string ntup_names[]= {"_B","_R","_RthenB"}; //names of each htest1 index type, which comes from the NTuples they run over


bool debug=true;

extern "C" void evolvepdf_(const double& , const double& , double* );
extern "C" double alphaspdf_(const double& Q);


string GetEnv( const string & var ) {

    const char* res= getenv( var.c_str() );

    std::string s = res!=NULL? res:"";
    cout<<"s: "<<s<<endl;
    return s;
}


TH1D* divide( const TH1D* h1, const TH1D* h2 ) {

    bool DBG=true;
    if ( h1==NULL || h2==NULL ) return NULL;

    TH1D* h = (TH1D*)h1->Clone();


    if ( DBG ) std::cout << "histograms h1: " << h1->GetTitle() << ", h2: " << h2->GetTitle() << std::endl;



    for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) {
        double b  = h2->GetBinContent(i);
        double be = h2->GetBinError(i);
        double t  = h1->GetBinContent(i);
        double te = h1->GetBinError(i);

        double r  = ( b!=0 ? t/b : 0 );
        //    double re = ( b!=0 ? sqrt((r+1)*r/b) : 0 );
        double re = 0;

        h->SetBinContent( i, r );
        h->SetBinError( i, re ) ;

        //    if ( debug ) std::cout << "\tx=" << h->GetBinCenter(i) << "\tratio=" << r << std::endl;
    }

    double hmin = h->GetBinContent(1);
    double hmax = h->GetBinContent(1);

    for ( int i=2 ; i<=h->GetNbinsX() ; i++ ) {
        double d = h->GetBinContent(i);
        if ( hmin>d ) hmin=d;
        if ( hmax<d ) hmax=d;
    }

    if ( DBG ) std::cout << "\tmin ratio = " << hmin << "\tmax ratio = " << hmax << std::endl;

    cout<<"h->GetMaximum(): "<<h->GetMaximum()<<", h->GetMinimum(): "<<h->GetMinimum()<<endl;

    if ( h->GetMaximum()<1.01 ) //h->SetMaximum(1.01);
        h->SetMaximum(0.99);
    if ( h->GetMinimum()>0.99 ) //h->SetMinimum(0.99);
        h->SetMinimum(1.01);

    return h;
}






int main(int argc, char** argv) {

    // use a default atlas inclusive grid
    //
    //string inputname="atlas2012_top.txt";

    //attempt to use EnvVar to find steeringfile, otherwise use default folder path steering/ in current dir
    string steeringName = "atlas2012_top.txt";
    string steeringPath = "steering";
    string steeringFile = steeringPath+"/"+steeringName;

    string steeringDefaultPath=GetEnv("STEERINGPATH");
    if(steeringDefaultPath.size()>0) {
        steeringFile=steeringDefaultPath+"/"+steeringName;
        cout<<" makegridfromsherpa::main: STEERINGPATH environment varaible found, using path: "<<steeringFile<<endl;
    }
    else {
        cout<<" makegridfromsherpa::main: STEERINGPATH environment varaible not set, using default: "<<steeringFile<<endl;
    }


//attempt to use EnvVar to find PDFsets, otherwise use default folder path PDFsets/ in current dir
    string pdfSetName = "CT10.LHgrid";
    string pdfSetPath = "PDFsets";
    string pdfSetFile = pdfSetPath+"/"+pdfSetName;

    string pdfSetDefaultPath = GetEnv("LHAPATH");

    if(pdfSetDefaultPath.size()>0) {
        pdfSetFile=pdfSetDefaultPath+"/"+pdfSetName;
        cout<<" makegridfromsherpa::main: LHAPATH environment varaible found, using path: "<<pdfSetFile<<endl;
    }
    else {
        cout<<" makegridfromsherpa::main: LHAPATH environment varaible not set, using default: "<<pdfSetFile<<endl;
    }


    /*
    //allow easier user input for different argumentparameters
    if(argc>=1 && argc%2=0) {
        for(int i=0; i<argc; i++) {
            if(argv[i].compare("-s")) {
                steeringFile=argv[i+1];
                std::cout<<" makegridfromsherpa::main: Using steeringfile located: "<<steeringFile<<std::endl;
            }
            if(argv[i].comapre("-e")) {
                nevmax=atoi(argv[i+1]);
                std::cout<<" makegridfromsherpa::main: Using events: "<<nevmax<<std::endl;
            }
            if(argv[i].comapre("-p")) {
                pdfSetFile=argv[i+1];
                std::cout<<" makegridfromsherpa::main: Using PDFset: "<<pdfSetFile<<std::endl;
            }
        }
    }
    else {
        cout<<"Invalid options passed."<<endl;
        exit(0);
    }
    */

    //allow for passing of another steering file name (default=inputname) and a different number of events(default=all_events)
    if ( argc>1 ) steeringFile = string(argv[1]);
    std::cout << " makegridfromsherpa::main: Reading steering file " << steeringFile << std::endl;

    long int nevmax=LONG_MAX; //allow for halt of execution when a user defined maximum number of events is reached
    if ( argc>2 ) {
        nevmax = atoi(argv[2]);
        std::cout << " makegridfromsherpa::main: Reading Number of events " << nevmax << std::endl;
    }


    //starting and ending indexes for histogram loop
    const int startIndex=0;
    const int endIndex=3; //change to choose ntups to go over


    //Create a uniquely named grid for each Type: 0-R, 1-B, 2-RthenB
    MyGrid *mygrid[endIndex];
    for(int i=0; i<endIndex; i++) {
        string version=ntup_names[i];
        cout<<" makegridfromsherpa::main: Creating grid using steeringFile: "<<steeringFile<<", version: "<<version<<endl;
        mygrid[i] = new MyGrid(steeringFile, version);
    }

    const int NGrid=mygrid[0]->GetNGrid(); //NGrid happens to be the same for all Types


    //create histograms to store tests:
    TH1D *htest1[endIndex][NGrid];  //htest1[0][]==R-Type, htest1[1][]==B-Type, htest1[2][]==RB-Type
    TH1D *htestRB[NGrid];           //holds the results of all htest1[i_R] + htest1[i_B]
    TH1D *href[endIndex][NGrid];    //holds histogram references from mygrid[x] to ensure they are equal to htest1[x]

    //maintain a count of the number of events run over for each histogram type we are testing
    int htestEventCount[endIndex];
    for(int i=0; i<endIndex; i++)
        htestEventCount[i]=0;


    //
    //Loop over different histograms to test different *-like* approaches
    //
    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        cout<<" makegridfromsherpa::main: Starting histoIndex loop: "<<histoIndex<<endl;

        //read in each NTuple for each test histogram
        string NtupName=mygrid[histoIndex]->GetInputNtupDir();
        if (histoIndex==i_R || histoIndex==i_B)
            NtupName+="/NTuple"+ntup_names[histoIndex]+"-like";
        else if(histoIndex==i_RB)
            NtupName+="/NTuple_*-like*";
        NtupName+=".root";

        cout<<" makegridfromsherpa::main: Opening "<<NtupName<<endl;


        TChain *fChain= new TChain("t3");

        fChain->Add(TString(NtupName));
        fChain->Print();

        t3 t(fChain);


        // set-up jet algorithm via fastjet wrapper class
        cout<<" makegridfromsherpa::main: Set up jet algorithm:"<<endl;
        fjClustering* jetclus = new fjClustering(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);

        MyEvent* myevent = new MyEvent();


        Long64_t nentries = fChain->GetEntries();
        cout<<" makegridfromsherpa::main: Number of events= "<<nentries<<endl;


        //
        // set-up test histogram
        //
        const int NGrid=mygrid[histoIndex]->GetNGrid();

        for (int igrid=0; igrid<NGrid; igrid++)
        {
            htest1[histoIndex][igrid]=mygrid[histoIndex]->GetHisto(igrid,string("htest1"+ntup_names[histoIndex]));

            htest1[histoIndex][igrid]->SetTitle(string("htest1"+ntup_names[histoIndex]).c_str());
            htest1[histoIndex][igrid]->SetName(string("htest1"+ntup_names[histoIndex]).c_str());
            htest1[histoIndex][igrid]->SetLineColor(histoIndex+1); //1-black(B), 2-red(R), 3-green(RthenB)

            if(debug) cout<<" makegridfromsherpa::main: Got histogram: "<<string("htest1_"+histoIndex)<<", with line color: "<<histoIndex+1<<endl;
        }


        //
        // Event loop
        //
        cout<<" makegridfromsherpa::main: Loop over events for: "<<ntup_names[histoIndex]<<endl;
        Long64_t nbytes = 0, nb = 0;
        
        for (Long64_t jentry=0; jentry<nentries && jentry<nevmax; jentry++)
        {
            t.GetEntry(jentry);

            nbytes += nb;
            //t.Show();

            Int_t np=(Int_t)t.nparticle;


            //inform user every NUMBER of events completed to illustrate progress
            if (jentry%100000==0 && jentry>0) {
                cout<<" Type: "<<ntup_names[histoIndex]<<"(histoIndex:"<<histoIndex<<")"<<", events processed "<<jentry<<" Hund-Thou of ";
                if(argc>2) cout<<nevmax<<". ";
                cout<<"TOTAL EVENTS: "<<nentries<<""<<endl;
            }

            /*
              if (debug) {
              // Documentation see
              // http://sherpa.hepforge.org/doc/SHERPA-MC-2.0.0.html
              //
              cout<<" " <<endl;
              cout<<" Event id = "<<t.id<<endl; //Event ID to identify correlated real sub-events.

              cout<<" Incoming parton id1 = "<<t.id1<<" id2 = "<<t.id2<<endl; // PDG code of incoming parton 1/2
              for (int i=0; i<np; i++)
              cout<<" Outcoming parton kf["<<i<<"]= "<<t.kf[i]<<endl; //

              cout<<" x1 = "<<t.x1<<" x2 = "<<t.x2<<endl;     // Bjorken-x of incoming parton 1/2
              cout<<" x1p= "<<t.x1p<<" x2p= "<<t.x2p<<endl;   // x’ for I-piece of incoming parton 1/2
              // Factorisation and Factorisaton scale.
              cout<<" fac_scale= "<<t.fac_scale<<" ren_scale= "<<t.ren_scale
              <<" alphas= "<<t.alphas<<endl; //" alphasPower= "<<(Char_t)t.alphasPower<<endl;
              //
              // nuwgt Number of additional ME weights for loops and integrated subtraction terms.
              // usr_wgt[nuwgt] Additional ME weights for loops and integrated subtraction terms.
              //8.8.4.1 Computing (differential) cross sections of real correction events with statistical errors
              cout<<" nuwgt= "<<t.nuwgt<<endl;
              // units are GeV and pb
              // weight    Event weight, if sub-event is treated independently.
              // weight2   Event weight, if correlated sub-events are treated as single event.
              // see
              cout<<" weight= "<<t.weight<<" weight2= "<<t.weight2<<endl;
              // me_wgt    ME weight (w/o PDF), corresponds to ’weight’.
              // me_wgt2   ME weight (w/o PDF), corresponds to ’weight2’.
              cout<<" me_wgt= "<<t.me_wgt<<" me_wgt2= "<<t.me_wgt2<<endl;
              cout<<" b_part= "<<t.part<<endl;
              cout<<" b_nparticle= "<<t.nparticle<<endl; //Number of outgoing partons.
              //cout<<" np= "<<np<<endl;
              //Int_t np=t.b_nparticle;
              }
            */


            int iorder=int(t.alphasPower);
            int itype =int(t.part[0]);
            //   int itype1 =int(t.part[1]);
            //cout<<" itype= "<<itype<<" itype1= "<<itype1<<endl;
            //cout<<" iorder= "<<iorder<<" itype= "<<itype<<endl;
            //if (iorder!=1) continue;

            jetclus->ClearJets();
            myevent->ClearEvent();;
            myevent->SetCMS(7000.);
            
            double wgt2_fac = (1.0/ (double)(2.*PI));
            if(iorder==2) {
                wgt2_fac = pow(wgt2_fac,2.0);
            } else { //3
                wgt2_fac = pow(wgt2_fac,3.0);
            }

            
            myevent->SetWeight(t.me_wgt2*wgt2_fac); //dealing entirely with weight2
            myevent->SetXSection(t.weight2);

            myevent->SetOrder(iorder);
            myevent->SetType (itype);
            myevent->SetEventId(t.id==21? 0:t.id); //conversion from sherpa gluon to appl_grid convention
            myevent->SetX1(t.x1);
            myevent->SetX2(t.x2);
            myevent->SetQ2(t.fac_scale*t.fac_scale);


            //
            // fill incoming partons
            //
            double pxin=0.;
            double pyin=0.;
            double pzin=0.;
            double Ein=0.;
            for (int ip=0; ip<np; ip++) {
                pxin+=t.px[ip];
                pyin+=t.py[ip];
                pzin+=t.pz[ip];
                Ein+=t.E[ip];
            }

            double shat=Ein*Ein-(pxin*pxin+pyin*pyin+pzin*pzin);
            shat/=t.x1*t.x2;
            myevent->SetCMS(sqrt(shat));
            //cout<<" pz= "<<pzin<<" E= "<<Ein<<" shat= "<<shat<<" sqrt(shat)= "<<sqrt(shat)<<endl;

            int pid=t.id1;
            if(pid==21) pid=0; //conversion from sherpa gluon to appl_grid convention
            myevent->push_back(pxin,pyin, t.x1*pzin,t.x1*Ein,pid);

            pid=t.id2;
            if(pid==21) pid=0; //conversion from sherpa gluon to appl_grid convention
            myevent->push_back(pxin,pyin,-t.x2*pzin,t.x2*Ein,pid);



            //std::cout<<" TEST: Checking name: "<<steeringFile<<std::endl;
            /*
            TString subProcType = steeringFile.c_str();
            if(subProcType.Contains("-gg")) {
              if(pid!=0) {
                  continue;
              }
            }
            */


            for (int ip=0; ip<np; ip++)
            {
                // Momentum components of the partons  kf  Parton PDG code
                pid=t.kf[ip];
                if(pid==21) pid=0; //conversion from sherpa gluon to appl_grid convention

                if (abs(pid)==6) {
                    myevent->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
                }

                if (abs(pid)==11||abs(pid)==12) {
                    myevent->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
                    continue;
                } else
                    // need to put in Gavins code here
                    // for top that will only work in LO
                    jetclus->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
            }

            //
            // run jet algorithm
            //
            //filling partons not jets
            /*
              jetclus->doClustering();
              if (debug) {
              cout<<" makegridfromsherpa::main: After jet clustering... "<<endl;
              jetclus->PrintJets();
              }

              myevent->push_back(jetclus->GetJets());

              if (debug) {
              cout<<" makegridfromsherpa::main: After pushing back jets, printing my event..."<<endl;
              myevent->Print2();
              }
            */


            //
            // fill the grid with the event for this histogram
            //
            mygrid[histoIndex]->fill(myevent);


            //
            // fill each test histogram
            //
            double obs=myevent->GetInvariantMass12(); // replace by GetObservable from steering
            for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++)
            {
                if (mygrid[histoIndex]->eventcuts(myevent,igrid)==false) continue;

                htest1[histoIndex][igrid]->Fill(obs,t.weight2);
            } //end loop over grid

            htestEventCount[histoIndex]++; //keep count of event for each type, 0-B, 1-R, 2-RthenB
        } //end loop over events




        //
        // get and set up test external histograms
        //
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++)
        {
            href[histoIndex][igrid]=(TH1D*)mygrid[histoIndex]->GetReference(igrid);
            if (!href[histoIndex][igrid]) cout<<" makegridfromsherpa::main: Reference from grid not found ! "<<endl;
            else {
                cout<<" makegridfromsherpa::main: Reference from grid found ! "<<endl;

                //Normalise(href[igrid],evuncorr*yfac,xfac,true);
                TString nameTitle="interal_href"+ntup_names[histoIndex];
                href[histoIndex][igrid]->SetTitle(nameTitle);
                href[histoIndex][igrid]->SetName(nameTitle);
                href[histoIndex][igrid]->SetLineColor(7); //SKY BLUE
                cout<<" makegridfromsherpa::main: Reference grid line color set to: 7"<<endl;
            }
            if (!href[histoIndex][igrid]) cout<<" makegridfromsherpa::main: href not found after norm ! "<<endl;
            else cout<<" makegridfromsherpa::main: Reference found after norm ! "<<endl;

            //cout<<" Print SubprocessRefHistos(igrid)"<<endl;
            //mygrid->PrintSubprocessRefHistos(igrid);
            //mygrid->PrintRefHistos(igrid);
        }


        cout<<" makegridfromsherpa::main: End histIndex loop: "<<histoIndex<<" of "<<(endIndex-1)<<", events this loop: "<<htestEventCount[histoIndex]<<endl;
    } //end of loop over all htest1 histograms



    //write grid to scale the internal reference histos and save them to *.root files
    cout<<" makegridfromsherpa::main: Printing and writing grid "<<endl;
    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        mygrid[histoIndex]->Print();
        
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++)
            mygrid[histoIndex]->ScaleInternalRefHistos(igrid);

        mygrid[histoIndex]->write_grid();
    }
    cout<<" makegridfromsherpa::main: Grid written "<<endl;






    //
    // Scale htest1
    // Add histograms of MyGrid for R-Like(i_R) to B-Like(i_B). Should be equal to both combined(i_RB)
    //
    cout<<" makegridfromsherpa::main: Adding R-Type and B-Type"<<endl;

    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        //scale each histograms by one devided by number of events depending on type
        for(int igrid=0; igrid<NGrid; igrid++) {

            //scale each histograms by one devided by number of events depending on type
            htest1[histoIndex][igrid]->Scale(1.0/htestEventCount[histoIndex]); //hR/nR, hB/nB, hBR/nBR, etc
            href[histoIndex][igrid]->Scale(1.0/htestEventCount[histoIndex]); //hrefR/nR, hrefB/nB,hrefRB/nRB, etc


            //sum the scaled R and B-type together, "should" be the same as htest[i_RB] where R and B are done together
            if(histoIndex==i_RB)
            {
                if(startIndex<=i_R && endIndex>=i_B) //make sure that R and B exist before adding them!
                {
                    htestRB[igrid]=mygrid[0]->GetHisto(igrid,string("htestRplusB_")); //mygrid index=0 is arbitrary, chosen to use GetHisto func
                    htestRB[igrid]->Add(htest1[i_R][igrid]);
                    htestRB[igrid]->Add(htest1[i_B][igrid]);

                    //htestRB[igrid]->Scale(1.0/htestEventCount[i_RB]);

                    htestRB[igrid]->SetLineColor(4); //blue
                }
            }
        }

        cout<<" makegridfromsherpa::main: Added R-Type and B-Type!"<<endl;



        //
        // Normalize after adding R and B-Types across all R and B events
        //
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++) {

            double evtot      =mygrid[histoIndex]->GetTotalEventNumber(igrid);
            double evalltot   =mygrid[histoIndex]->GetTotalEventNumber();
            double evuncorr   =mygrid[histoIndex]->GetUncorrellatedEventNumber(igrid);
            double evalluncorr=mygrid[histoIndex]->GetUncorrellatedEventNumber();
            //double evalltot=mygrid[0]->GetUncorrellatedEventNumber();
            cout<<" makegridfromsherpa::main:**NORM CHECK: evtot: "<<evtot<<", evalltot: "<<evalltot<<", evuncorr: "<<evuncorr<<", evalluncorr: "<<evalluncorr<<endl;
            cout<<" makegridfromsherpa::main:**NORM CHECK: nR: "<<htestEventCount[0]<<", nB: "<<htestEventCount[1]<<", nRB: "<<htestEventCount[2]<<endl;

            //temporarily changed to a default!!
            evtot      =1;
            evalltot   =1;
            evuncorr   =1;
            evalluncorr=1;


            /*          //The evtot from the grid should be equal to what we counted in the event loop for each R, B, and RB Type
            //NOTE: evtot is 0 if mygrid->fill(*) never runs!
            cout<<" makegridfromsherpa::main: Comparing event counts:\n\t evtot: "<<evtot<<" events counted: "<<htestEventCount[histoIndex]+1<<endl; //+1 b/c start count at 0
            if(evtot!=htestEventCount[histoIndex]+1) {
            cout<<" makegridfromsherpa::main: WARNING: Event counts were not the same on histIndex: "<<histoIndex<<endl;
            exit(0);
            }
            */

            cout<<" makegridfromsherpa::main: Norm: evalltot= "<<evalltot<<" evalluncorr= "<<evalluncorr<<endl;
            cout<<" makegridfromsherpa::main: Norm: evtot= "<<evtot<<" evuncorr= "<<evuncorr<<endl;

            MyData *mydata=mygrid[histoIndex]->GetMyData(igrid);
            if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
            double yfac=mydata->GetUnitfbFactor();
            double xfac=mydata->GetUnitGeVFactor();
            cout<<" makegridfromsherpa::main: MyGrid::Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;



            mygrid[histoIndex]->Normalise(htest1[histoIndex][igrid],evtot*yfac,xfac,true);   //normalise 0-B, 1-R, and 2-RthenB Type
            mygrid[histoIndex]->Normalise(href[histoIndex][igrid],evtot*yfac,xfac,true);     //normalise hrefB, hrefR, and hrefRthenB

            if(histoIndex==i_RB)
                mygrid[histoIndex]->Normalise(htestRB[igrid],evtot*yfac,xfac,true);     //normalise histo of scaled R + scaled B together if it exists
        }

    } //end of loop over all htest1 histograms






    //
    // Convolute and print all external histograms
    //
    cout<< "\n makegridfromsherpa::main: Performing convolute"<<endl;
    //NGrid = mygrid[0]->GetNGrid(); //NGrid will be the same for all grids, so grid[0] is arbirary
    TH1D* convGridHistos[endIndex+1][NGrid];
    TH1D* LOconvGridHistos[endIndex+1][NGrid];
    TH1D* subProcConvGridHistos[endIndex+1][NGrid][121];
    
    int nLoops = 1;

    //string pdf_set_name = "PDFsets/CT10.LHgrid"; //hardcoded
    //LHAPDF::initPDFSet(pdf_set_name.c_str(), 0);
    LHAPDF::initPDFSet(pdfSetFile.c_str(), 0);



    std::vector<std::vector<double> > ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));
    //creating an identity matrix needed in convolute function - temp printing to make sure it's correct
    for(int x=0; x<ckm2.size(); x++) {
        for(int y=0; y<ckm2.at(x).size(); y++) {
            if(x==y) ckm2.at(x).at(y)=1;
        }
    }



    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        TFile *fout;

        for(int igrid=0; igrid<NGrid; igrid++)
        {
            mygrid[histoIndex]->GetGrid(igrid)->setckm(ckm2);

            LOconvGridHistos[histoIndex][igrid] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute( evolvepdf_, alphaspdf_, 0 );
            LOconvGridHistos[histoIndex][igrid]->SetName((TString) ("LOconvolute_for" + ntup_names[histoIndex]));
            LOconvGridHistos[histoIndex][igrid]->SetTitle((TString) ("LOconvolute_for" + ntup_names[histoIndex]));
            LOconvGridHistos[histoIndex][igrid]->SetLineColor(kBlue);
            
            for(int isubproc=0; isubproc<mygrid[histoIndex]->GetNSubProcess(igrid);  isubproc++){
                subProcConvGridHistos[histoIndex][igrid][isubproc] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute_subproc(isubproc, evolvepdf_, alphaspdf_, nLoops );
                string sub_proc_hist_name="subProcConvolute_for"+ntup_names[histoIndex]+"_subProc-"+to_string(isubproc);
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetName((TString) (sub_proc_hist_name));
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name));
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetLineColor(kBlue);
            }
            
            convGridHistos[histoIndex][igrid] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute( evolvepdf_, alphaspdf_, nLoops );
            convGridHistos[histoIndex][igrid]->SetName((TString) ("convolute_for" + ntup_names[histoIndex]));
            convGridHistos[histoIndex][igrid]->SetTitle((TString) ("convolute_for" + ntup_names[histoIndex]));
            convGridHistos[histoIndex][igrid]->SetLineColor(kBlue);
            cout<< " makegridfromsherpa::main: Printing after convol for type: "<<ntup_names[histoIndex]<<"(histoIndex:"<<histoIndex<<")"<<endl;


            string filename=mygrid[histoIndex]->GetGridName(igrid);
            filename+=ntup_names[histoIndex]+"_histos.root";
            fout= new TFile(filename.c_str(),"recreate");


            MyData *mydata=mygrid[histoIndex]->GetMyData(igrid);
            if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
            double yfac=mydata->GetUnitfbFactor();
            double xfac=mydata->GetUnitGeVFactor();
            cout<<" makegridfromsherpa::main: Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


            LOconvGridHistos[histoIndex][igrid]->Scale(1.0/htestEventCount[histoIndex]);
            mygrid[histoIndex]->Normalise(LOconvGridHistos[histoIndex][igrid],yfac,xfac,true);
            LOconvGridHistos[histoIndex][igrid]->Write();
            LOconvGridHistos[histoIndex][igrid]->Draw();

            for(int isubproc=0; isubproc<mygrid[histoIndex]->GetNSubProcess(igrid);  isubproc++){
                subProcConvGridHistos[histoIndex][igrid][isubproc]->Scale(1.0/htestEventCount[histoIndex]);
                mygrid[histoIndex]->Normalise(subProcConvGridHistos[histoIndex][igrid][isubproc],yfac,xfac,true);
                subProcConvGridHistos[histoIndex][igrid][isubproc]->Write();
                subProcConvGridHistos[histoIndex][igrid][isubproc]->Draw();
            }
            
            convGridHistos[histoIndex][igrid]->Scale(1.0/htestEventCount[histoIndex]);
            mygrid[histoIndex]->Normalise(convGridHistos[histoIndex][igrid],yfac,xfac,true);
            convGridHistos[histoIndex][igrid]->Write();
            convGridHistos[histoIndex][igrid]->Draw();


            htest1[histoIndex][igrid]->Print("all");
            htest1[histoIndex][igrid]->Write();
            htest1[histoIndex][igrid]->Draw();

            href[histoIndex][igrid]->Print("all");
            href[histoIndex][igrid]->Write();
            href[histoIndex][igrid]->Draw();

            TH1D* ratio1 = divide( convGridHistos[histoIndex][igrid],htest1[histoIndex][igrid] );
            if ( ratio1 ) {
                ratio1->SetTitle(string("ratio1_convolute/htest1"+ntup_names[histoIndex]).c_str());
                ratio1->SetName(string("ratio1_convolute/htest1"+ntup_names[histoIndex]).c_str());
                ratio1->Print("all");
                ratio1->Write();
                ratio1->Draw();
            }
            TH1D* ratio2 = divide( htest1[histoIndex][igrid], convGridHistos[histoIndex][igrid] );
            if ( ratio2 ) {
                ratio2->SetTitle(string("ratio2_htest1/convolute"+ntup_names[histoIndex]).c_str());
                ratio2->SetName(string("ratio2_htest1/convolute"+ntup_names[histoIndex]).c_str());
                ratio2->Print("all");
                ratio2->Write();
                ratio2->Draw();
            }
        }

        fout->Write();
        fout->Close();
    }











    //after performing convolutes for B, R, RthenB types, add grids and histo for R and B to get RplusB and check it's convolute
    cout<<" makegridfromsherpa::main: Adding B-Type mygrid to R-Type mygrid"<<endl;
    mygrid[i_R]->AddGrid(mygrid[i_B]); //STILL NEEDS APPLgrid add implimentation so the internals update as well as MyGrid!!
    mygrid[i_R]->SetGridVersionName(string("_RplusB"));
    mygrid[i_R]->write_grid();



    cout<<" makegridfromsherpa::main: Normalising Internal Reference hIstograms "<<endl;
    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++) {
            mygrid[histoIndex]->NormaliseInternalRefHistos(igrid);
            
            mygrid[histoIndex]->SetGridVersionName(ntup_names[histoIndex]+"_norm");
            mygrid[histoIndex]->write_grid();
        }
    }
    cout<<" makegridfromsherpa::main: Internal Reference hIstograms normalised!"<<endl;



    TH1D* ConvHistoRplusB[mygrid[i_R]->GetNGrid()];
    TH1D* hrefRplusB[mygrid[i_R]->GetNGrid()];

    for(int igrid=0; igrid<mygrid[i_R]->GetNGrid(); igrid++)
    {
        TFile *fout;

        string filename=mygrid[i_R]->GetGridName(igrid);
        filename+="_RplusB_histos.root";
        fout= new TFile(filename.c_str(),"recreate");

        //problem: only performing convolute to R...
        ConvHistoRplusB[igrid] = (TH1D*)mygrid[i_R]->GetGrid(igrid)->convolute( evolvepdf_, alphaspdf_, nLoops );
        ConvHistoRplusB[igrid]->SetName((TString) ("convolute_for_RplusB"));
        ConvHistoRplusB[igrid]->SetTitle((TString) ("convolute_for_RplusB"));
        ConvHistoRplusB[igrid]->SetLineColor(kBlue);

        //NOTE: NEED TO ACCOUNT FOR xfac and yfac being different when adding grid's MyData??
        MyData *mydata=mygrid[i_R]->GetMyData(igrid);
        if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
        double yfac=mydata->GetUnitfbFactor();
        double xfac=mydata->GetUnitGeVFactor();
        cout<<" makegridfromsherpa::main: Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


        ConvHistoRplusB[igrid]->Print("all");
        ConvHistoRplusB[igrid]->Scale(1.0/(htestEventCount[i_R]+htestEventCount[i_B]));
        mygrid[i_R]->Normalise(ConvHistoRplusB[igrid],yfac,xfac,true);

        ConvHistoRplusB[igrid]->Print("all");
        ConvHistoRplusB[igrid]->Write();
        ConvHistoRplusB[igrid]->Draw();

        htestRB[igrid]->Print("all");
        htestRB[igrid]->Write();
        htestRB[igrid]->Draw();


        //NOTE: currently the applgrids are not added, only the MyGrid histos, so to get the reference histo for RplusB 
        // we have to get each individually and add them
        hrefRplusB[igrid] = (TH1D*)mygrid[i_R]->GetReference(igrid); 
        hrefRplusB[igrid]->Add( (TH1D*)mygrid[i_B]->GetReference(igrid) );
        if (!hrefRplusB[igrid]) cout<<" makegridfromsherpa::main: Reference from grid not found ! "<<endl;
        else {
            cout<<" makegridfromsherpa::main: Reference from grid found ! "<<endl;

            TString nameTitle="interal_href_RplusB";
            hrefRplusB[igrid]->SetTitle(nameTitle);
            hrefRplusB[igrid]->SetName(nameTitle);
            hrefRplusB[igrid]->SetLineColor(7); //SKY BLUE
        }

        //note the scaling to total events for R and B
        hrefRplusB[igrid]->Scale(1.0/(htestEventCount[i_R]+htestEventCount[i_B])); //Are these inidividual internal refs scaled??
        mygrid[i_R]->Normalise(hrefRplusB[igrid],yfac,xfac,true);     //normalise hrefB, hrefR, and hrefRB

        hrefRplusB[igrid]->Print("all");
        hrefRplusB[igrid]->Write();
        hrefRplusB[igrid]->Draw();

        TH1D* ratio1 = divide( ConvHistoRplusB[igrid],htestRB[igrid] );
        if ( ratio1 ) {
            //ratio->Scale(1/htestEventCount[i_B]);
            ratio1->SetTitle(string("ratio1_convolute/htestRplusB").c_str());
            ratio1->SetName(string("ratio1_convolute/htestRplusB").c_str());
            ratio1->Print("all");
            ratio1->Write();
            ratio1->Draw();
        }
        TH1D* ratio2 = divide( htestRB[igrid], ConvHistoRplusB[igrid] );
        if ( ratio2 ) {
            //ratio->Scale(1/htestEventCount[i_B]);
            ratio2->SetTitle(string("ratio2_htestRplusB/convolute").c_str());
            ratio2->SetName(string("ratio2_htestRplusB/convolute").c_str());
            ratio2->Print("all");
            ratio2->Write();
            ratio2->Draw();
        }

        fout->Write();
        fout->Close();
    }











    /*
        //
        // Write out external histograms
        // NOTE: Could(should?) be made nicer by looping over all htest1 histograms instead of hardcoding indexes in htest1
        //
        cout<<"\n makegridfromsherpa::main: Writing test histos: "<<endl;
        for (int igrid=0; igrid<mygrid[0]->GetNGrid(); igrid++) {
            string filename=mygrid[0]->GetGridFullFileName(igrid);

            filename.replace( filename.find(".root"), 5, "-histos.root");
            cout<<" makegridfromsherpa::main: Write histos to filename= "<<filename<<endl;


            TFile *f1= new TFile(filename.c_str(),"recreate");

            //TApplication myapp("myapp",0,0);
            //TApplication *theApp = new TApplication("My ROOT Application",0,0);
            //theApp->SetReturnFromRun(true);

            MyData *mydata=mygrid[0]->GetMyData(igrid);
            if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
            else cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] read "<<endl;

            MyFrameData *myframe= new MyFrameData(600,600,mydata);
            if (!myframe) cout<<" makegridfromsherpa::main: myframe not found "<<endl;
            else cout<<" makegridfromsherpa::main: frame created "<<endl;

            mydata->DrawData();
            cout<<" makegridfromsherpa::main: DrawData finished "<<endl;



            //display and print to terminal desired histogram information

            ////terminal data output
            //htest1[i_R][igrid]->Draw("same"); //R-Type
            //htest1[i_B][igrid]->Draw("same"); //B-type
            //htest1[i_RB][igrid]->Draw("same"); //RB-type

            //href[i_R][igrid]->Draw("same"); //ref for R-type
            //href[i_B][igrid]->Draw("same"); //ref for B-type
            //href[i_RB][igrid]->Draw("same"); //ref for RB-Type
            //htestRB[igrid]->Draw("same"); //RplusB


            ////graphical data output
            //htest1[i_R][igrid]->Print("all"); //R-Type
            //htest1[i_B][igrid]->Print("all"); //B-Type
            //htest1[i_RB][igrid]->Print("all"); //RB-Type

            //href[i_R][igrid]->Print("all"); //ref for R-Type
            //href[i_B][igrid]->Print("all"); //ref for B-Type
            //href[i_RB][igrid]->Print("all"); //ref for RB-Type
            //htestRB[igrid]->Print("all"); //RplusB


            gPad->Update();
            //theApp->Run(kTRUE);
        }


        gPad->Print("xsec.pdf");
        system("open xsec.pdf &");

    */
    return 0;
}
