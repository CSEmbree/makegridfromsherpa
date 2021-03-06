/*
 * Title:   makegridfromsherpa
 * Purpose: Given SHERPA ntuples, we test appl_grid LO and NLO convolute
 * Authors: Dr. Carli
 * Helpers: Mark Sutton, Cameron Embree
 */

//general
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

//Dr. Carli's My* classes
#include "MyEvent.h"
#include "MyData.h"
#include "MyFrameData.h"
#include "MyGrid.h"

//root
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>

//extra
#include "LHAPDF/LHAPDF.h"
#include "VariableDefinitions.h"    //added for convolute


#define PI 3.141592653589793238462


/*
 * EXAMPLE Execution:
 *                   ./makegridfromsherpa
 *                   ./makegridfromsherpa <filename> <numevents>
 *                   ./makegridfromsherpa steering/atlas2012_top-config.txt 1000000
 */



//Constructs for more easily running over NTuples.
//Using B-like Ntup first
enum ntup_types {i_B=0,i_R,i_RB};           //indexes for different histograms depending on "NTuple_*-like" in htest1
string ntup_names[]= {"_B","_R","_RthenB"}; //names of each htest1 index type, which comes from the NTuples they run over



long int nevmax=LONG_MAX; //allow for halt of execution when a user defined maximum number of events is reached
double wgt_fac=-1; //factor needed for conversion from SHERPA to ApplGrid when setting the weight 
static const int GLUON = 6; //used for more readable code


bool debug=false;
bool debug2=false;


extern "C" void evolvepdf_(const double& , const double& , double* );
extern "C" double alphaspdf_(const double& Q);


/*
//calculates the NLO weight based on partons id1 and id2
double GetWeightNLO(int id1, int id2, t3 t) {

    //current particle numering convention:
    //0-tbar, 1-bbar, 2-cbar, 3-sbar, 4-ubar, 5-dbar, 6-glu, 7-d, 8-u, 9-s, 10-c, 11-b, 12-t
    //hardcoded temp val in case numbering convention changes
    static const int    GLUON       = 6;
    static const int    BOTTOM_BAR  = 1;
    static const int    BOTTOM      = 11;

    double wgt = 0; //wgt that will be returned by end of NLO weight calc



    //****END -- TEST TO CHECK CORRECT WEIGHT - STEFAN's code
    //LHAPDF::initPDFSet(pdfSetFile.c_str(), 0); //need this? Assuming it is already init???


    //setup xf1 and xf2 to be filled from evolvepdf
    int nWgts=13;
    double xf1[nWgts];
    double xf2[nWgts];
    double xf1p[nWgts];
    double xf2p[nWgts];
    for( int iproc=0 ; iproc<nWgts ; iproc++ ) {
        xf1[iproc] = 0.0; //emptied
        xf2[iproc] = 0.0; //emptied
    }

    //evolve to get weights in xf1, xf2
    getPDF(t.x1,t.fac_scale,xf1);
    getPDF(t.x2,t.fac_scale,xf2);


    //computing f1 from xf1 and f2 from xf2 for later usage (and readability) conveyance
    double f1[nWgts];
    double f2[nWgts];
    double f1p[nWgts];
    double f2p[nWgts];
    for( int iproc=0 ; iproc<nWgts ; iproc++ ) {
        f1[iproc] = xf1[iproc] / t.x1;
        f2[iproc] = xf2[iproc] / t.x2;
    }



    double asf=1;
    //double lr=log(mur2/sqr(p_vars->m_mur));
    //double lf=log(muf2/sqr(p_vars->m_muf));
    double lr=0;
    double lf=0;

    if(debug) {
        for(int i=0; i<nWgts; i++) std::cout<<"xf1["<<i<<"]: "<<xf1[i]<<std::endl; //show all
        for(int i=0; i<nWgts; i++) std::cout<<"xf2["<<i<<"]: "<<xf2[i]<<std::endl; //show all
    }


    //get weights from f1(xf1/x1) of partion id1 and id2 for later usage conveyance
    double fa = f1[id1]; // (xf1/x1)[id1]
    double fb = f2[id2]; // (xf2/x2)[id2]



    //get information needed to compute weight
    double w[9];
    w[0] = t.me_wgt + t.usr_wgts[0] * lr + t.usr_wgts[1] * lr * lr / 2.0;

    bool wnz=false;
    for ( int i=1 ; i<9 ; ++i ) {
        w[i] = t.usr_wgts[i+1] + t.usr_wgts[i+9] * lf;
        if (w[i]==0) wnz=true;
    }

    wgt = w[0] * fa * fb;


    //wgt=t.me_wgt2+t.usr_wgts[0]*lr+t.usr_wgts[1]*lr*lr/2.0;

    if (wnz==true) {
        double faq = 0.;
        double faqx = 0.0;

        double fag = 0.0;
        double fagx = 0.0;

        double fbq = 0.0;
        double fbqx = 0.0;

        double fbg = 0.0;
        double fbgx = 0.0;

        if ( id1 != GLUON ) {
            //QUARK
            faq = fa;
            fag = f1[6];

            getPDF( t.x1/t.x1p, t.fac_scale, xf1p );

            faqx = xf1p[id1] / t.x1;
            fagx = xf1p[6]   / t.x1;
        }
        else {
            //GLU
            fag=fa;
            for ( int i=1 ; i<nWgts-1 ; ++i)
                if( i!=GLUON ) faq += f1[i];

            getPDF( t.x1/t.x1p, t.fac_scale , xf1p );

            fagx = (xf1p[id1] / t.x1);
            for ( int i=1 ; i<nWgts-1 ; ++i )
                if( i != GLUON ) faqx += (xf1p[i]/t.x1);
        }
        if ( id2 != GLUON ) {
            //QUARK
            fbq = fb;
            fbg = f2[GLUON];

            getPDF( t.x2/t.x2p, t.fac_scale, xf2p );

            fbqx = xf2p[id2] / t.x2;
            fbgx = xf2p[GLUON] / t.x2;
        }
        else {
            //GLU
            fbg = fb;
            for ( int i = 1 ; i<nWgts-1 ; ++i)
                if( i != GLUON ) fbq += f2[i];

            getPDF( t.x2/t.x2p , t.fac_scale , xf2p);

            fbgx = (xf2p[id2] / t.x2);
            for ( int i=1 ; i<nWgts-1 ; ++i )
                if( i != GLUON ) fbqx += (xf2p[i] / t.x2);
        }

        //compute weight
        wgt+=(faq*w[1]+faqx*w[2]+fag*w[3]+fagx*w[4])*fb;
        wgt+=(fbq*w[5]+fbqx*w[6]+fbg*w[7]+fbgx*w[8])*fa;

    }
    else {
        std::cout<<" makrgridfromsherpa::GetWeightNLO: ERROR: wnz invalid? wnz:"<<wnz<<std::endl;
        exit(0);
    }


    wgt=wgt*asf;
    std::cout<<" makrgridfromsherpa::GetWeightNLO: RESULT: "<<wgt<<std::endl;




    //reference prints
    std::cout<<"\tid1: "<<id1<<", id2: "<<id2<<std::endl;
    std::cout<<"\txf1["<<id1<<"]: "<<xf1[id1]<<std::endl;
    std::cout<<"\txf2["<<id2<<"]: "<<xf2[id2]<<std::endl;
    std::cout<<"\tt.weight: "<<t.weight<<std::endl;
    std::cout<<"\tt.weight2: "<<t.weight2<<std::endl;
    std::cout<<"\tt.me_wgt: "<<t.me_wgt<<std::endl;
    std::cout<<"\tt.me_wgt2: "<<t.me_wgt2<<std::endl;
    double myweight=(xf1[id1]*xf2[id2]*t.me_wgt)/(t.x1*t.x2);
    std::cout<<"\t(xf1["<<id1<<"]*xf2["<<id2<<"]*t.me_wgt)/(t.x1*t.x2)= "<<myweight<<std::endl;

    //****END -- TEST TO CHECK CORRECT WEIGHT - NEW











    //retreive values for x# and x#prime from tree for later usage (and readability) conveyance
    double x1  = t.x1;
    double x2  = t.x2;
    double x1p = t.x1p;
    double x2p = t.x2p;
    double fac_scale = t.fac_scale;

    //containers for weights from evolvepdf_ calls
    double xf1[nWgts];
    double xf2[nWgts];
    double f1 [nWgts];
    double f2 [nWgts];

    //sum over all f1 and f2
    double sf1;
    double sf2;



    //get information needed to compute weight
    double w[9];
    w[0] = t.me_wgt + t.usr_wgts[0] * lr + t.usr_wgts[1] * lr * lr / 2.0;


    wgt = ( w[0] * wgt_fac ) / npairs;
    htestWgt = ??

    myevent->SetWeight( wgt );
    myevent->SetX1( t.x1 );
    myevent->SetX2( t.x2 );
    mygrid[histoIndex]->fill( myevent );
    for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, htestWgt );
    eventCount[histoIndex]++;


    bool wnz=false;
    for ( int i=1 ; i<9 ; ++i ) {
        w[i] = t.usr_wgts[i+1] + t.usr_wgts[i+9] * lf;
        if (w[i]==0) wnz=true;
    }

    if(wnz)
        //Implimenting Mark's suggestions
        if( id1 != GLUON )
        {
            if( id2 != GLUON ) {
                //QUARK-QUARK
                //needed: fa, fb, fap, fbp, fag, fbg, fagp, fbgp
                //using : xf1, xf2, x1, x2, xf1p, xf2p, x1p, x2p

                getPDF( x1, fac_scale, xf1 );
                getPDF( x2, fac_scale, xf2 );
                fa  = xf1[id1]   / x1;
                fb  = xf2[id2]   / x2;
                fag = xf1[GLUON] / x1;
                fbg = xf2[GLUON] / x2;

                getPDF( x1/x1p, fac_scale, xf1p );
                getPDF( x2/x2p, fac_scale, xf2p );
                fap  = xf1p[id1]   / x1;
                fbp  = xf2p[id2]   / x2;
                fagp = xf1p[GLUON] / x1;
                fbgp = xf2p[GLUON] / x2;


                1  fa   * fb   *w[1];
                2  fap  * fb   *w[2] * (1/x1p);
                3  fag  * fb   *w[3];
                4  fagp * fb   *w[4] * (1/x1p);
                5  fa   * fb   *w[5];
                6  fa   * fbp  *w[6] * (1/x2p);
                7  fa   * fbg  *w[7];
                8  fa   * fbgp *w[8] * (1/x2p);


                myevent->SetX1( t.x1p );
                myevent->SetX2( t.x2 );

                wgt  = w[1] * wgt_fac * (1.0/npairs);
                hwgt = fa*w[1] * fb;
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[2] * (1/x1p) * wgt_fac * (1.0/npairs);
                hwgt = fap*w[2]  * fb;
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[3] * wgt_fac * (1.0/npairs);
                hwgt = fag*w[3]  * fb;
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[4] * (1/x1p) * wgt_fac * (1.0/npairs);
                hwgt = fagp*w[4] * fb;
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );


                myevent->SetX1( t.x1 );
                myevent->SetX2( t.x2p );

                wgt  = w[5] * wgt_fac * (1.0/npairs);
                hwgt = fa * fb*w[5];
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[6] * (1/x2p) * wgt_fac * (1.0/npairs);
                hwgt = fa * fbp*w[6];
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[7] * wgt_fac * (1.0/npairs);
                hwgt = fa * fbg*w[7];
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );

                wgt  = w[8] * (1/x2p) * wgt_fac * (1.0/npairs);
                hwgt = fa * fbgp*w[8];
                myevent->setWeight( wgt );
                mygrid[histoIndex]->fill( myevent );
                htest3->fill( obs, hwgt );
            }
            else {
                //QUARK-GLUON
            }
        }
        else {
            if( id2 != GLUON ) {
                //GLUON-QUARK
                //needed: sf1, sf1p, fg1, fg1p, fb, fbp, fg2, fg2p
                //using : xf1, xf2, x1, x2, xf1p, xf2p, x1p, x2p

                getPDF( x1, fac_scale, xf1 );
                getPDF( x2, fac_scale, xf2 );
                fb  = xf2[id2]   / x2;
                fg1 = xf1[GLUON] / x1;
                fg2 = xf2[GLUON] / x2;
                for( int i=1; i<nWgts-1; i++ ) {
                    sf1 += xf1[i] / x1;
                }

                getPDF( x1/x1p, fac_scale, xf1p );
                getPDF( x2/x2p, fac_scale, xf2p );
                fbp  = xf2p[id2]   / x2;
                fg1p = xf1p[GLUON] / x1;
                fg2p = xf2p[GLUON] / x2;
                for( int i=1; i<nWgts-1; i++ ) {
                    sf1p = xf1p[i] / x1;
                }


                1  sf1  * fb   * w[1];
                2  sf1p * fb   * w[2] * (1/x1p);
                3  fg1  * fb   * w[3];
                4  fg1p * fb   * w[4] * (1/x1p);
                5  f1g  * fb   * w[5];
                6  f1g  * fbp  * w[6] * (1/x2p);
                7  f1g  * f2g  * w[7];
                8  f1g  * f2gp * w[8] * (1/w2p);

            }
            else {
                //GLUON-GLUON
                //needed: fg1, fg2, sf1, sf2, fg1p, fg2p, sf1p, sf2p
                //using : xf1, xf2, x1, x2, xf1p, xf2p, x1p, x2p

                getPDF( x1, fac_scale, xf1 );
                getPDF( x2, fac_scale, xf2 );
                fg1 = xf1[GLUON] / x1;
                fg2 = xf2[GLUON] / x2;
                for( int i=1; i<nWgts-1; i++ ) {
                    sf1 += xf1[i] / x1;
                    sf2 += xf2[i] / x2;
                }


                getPDF( x1/x1p, fac_scale, xf1p );
                getPDF( x2/x2p, fac_scale, xf2p );
                fg1p = xf1p[GLUON] / x1;
                fg2p = xf2p[GLUON] / x2;
                for( int i=1; i<nWgts-1; i++ ) {
                    sf1p += xf1p[i] / x1;
                    sf2p += xf2p[i] / x2;
                }


                1   sf1  * fg2  * w[1];
                2   sf1p * fg2  * w[2] * (1/x1);
                3   fg1  * fg2  * w[3];
                4   fg1p * fg2  * w[3] * (1/x1); //Should this be w[4]???
                5   fg1  * sf2  * w[5];
                6   fg1  * sf2p * w[6] * (1/x2);
                7   fg1  * fg2  * w[7];
                8   fg1  * fg2p * w[8] * (1/x2);
            }
        }
    else {
        std::cout<<" makrgridfromsherpa::GetWeightNLO: ERROR: wnz invalid? wnz:"<<wnz<<std::endl;
        exit(0);
    }


    std::cout<<" makrgridfromsherpa::GetWeightNLO: NLO weight is: "<<wgt<<std::endl;
    return wgt;
}
*/

//skeleton to replace hard coded and repetative NLO fillings
double PrepAndFillEvent(MyEvent *myevent, MyGrid *mygrid, double w, double x1, double x2, double x1p, double x2p, double id1, double id2, double xsec=0) {

    myevent->SetX1( x1 );
    myevent->SetX2( x2 );

    myevent->SetID( 0, id1 );   //setting parton A's ID
    myevent->SetID( 1, id2 );   //setting parton B's ID

    myevent->SetXSection( 0 ); //hack for not filling ref histos


    //determine subprocess that may be needed to divide by
    int currentSubProc = mygrid->GetDecideSubProcess( id1, id2 );
    int npairs         = mygrid->GetNSubProcessPairs( currentSubProc );

    //general weight calculation
    double wgt = ( w * wgt_fac ) / npairs;

    return wgt;
}



//calculates mins & maxs of ren & fac scales to get the correct Q2
void GetRenAndFacMaxAndMins(string NtupName, double *facMin, double *facMax, double *renMin, double *renMax) {

    TChain *fChaintmp= new TChain("t3");

    fChaintmp->Add(TString(NtupName));
    if(debug) fChaintmp->Print();

    t3 ttmp(fChaintmp);


    // determine Q2 and x boundaries from fac and ren scales
    Long64_t nentriestmp = fChaintmp->GetEntries();
    for ( Long64_t jentry=0 ; jentry<nentriestmp && jentry<nevmax ; jentry++ ) {
        ttmp.GetEntry(jentry);
        if (ttmp.fac_scale < *facMin) *facMin = ttmp.fac_scale;
        if (ttmp.fac_scale > *facMax) *facMax = ttmp.fac_scale;
        if (ttmp.ren_scale < *renMin) *renMin = ttmp.ren_scale;
        if (ttmp.ren_scale > *renMax) *renMax = ttmp.ren_scale;
    }

    if(debug) std::cout<<" makegridfromsherpa::GetRenAndFacMaxAndMins: ntuple:"<<NtupName
                           <<", facMin= "<<*facMin
                           <<", facMax= "<<*facMax
                           <<", renMin= "<<*renMin
                           <<", renMax= "<<*renMax<<std::endl;

    delete fChaintmp; //cleanup
}


//makes conversion of id1, id2 and fills f[]
void getPDF(const double& x, const double& Q2, double* f) {

    evolvepdf_(x, Q2, f);

    for(int id=0; id<13; id++) f[id]=x;
}

//allows for getting value of environment variables
string GetEnv( const string & var ) {

    const char* res= getenv( var.c_str() );

    std::string s = res!=NULL? res:"";
    return s;
}

//method from Mark to divide two histograms nicely for ratio comparisons between two histos
TH1D* divide( const TH1D* h1, const TH1D* h2 ) {
    std::cout<<" makegridfromsherpa::divide: dividing two histos..."<<std::endl;

    bool DBG=true;
    if ( h1==NULL || h2==NULL ) return NULL;

    TH1D* h = (TH1D*)h1->Clone();

    if ( DBG ) std::cout << " makegridfromsherpa::divide:histograms h1: " << h1->GetTitle() << ", h2: " << h2->GetTitle() << std::endl;


    for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) {
        double b  = h2->GetBinContent(i);
        double be = h2->GetBinError(i);
        double t  = h1->GetBinContent(i);
        double te = h1->GetBinError(i);

        double r  = ( b!=0 ? t/b : 0 );
        double re = 0;

        h->SetBinContent( i, r );
        h->SetBinError( i, re ) ;
    }

    double hmin = h->GetBinContent(1);
    double hmax = h->GetBinContent(1);

    for ( int i=2 ; i<=h->GetNbinsX() ; i++ ) {
        double d = h->GetBinContent(i);
        if ( hmin>d ) hmin=d;
        if ( hmax<d ) hmax=d;
    }

    if ( DBG ) std::cout << " makegridfromsherpa::divide: \tmin ratio = " << hmin << "\tmax ratio = " << hmax << std::endl;

    cout<<" makegridfromsherpa::divide: h->GetMaximum(): "<<h->GetMaximum()<<", h->GetMinimum(): "<<h->GetMinimum()<<endl;

    if ( h->GetMaximum()<1.01 ) //h->SetMaximum(1.01);
        h->SetMaximum(0.99);
    if ( h->GetMinimum()>0.99 ) //h->SetMinimum(0.99);
        h->SetMinimum(1.01);

    return h;
}








int main(int argc, char** argv) {

    // use a default atlas inclusive grid

    //attempt to use EnvVar to find steeringfile, otherwise use default folder path steering/ in current dir
    string steeringName = "atlas2012_top-config.txt"; //*-config for lumi, without for generic
    string steeringPath = "steering";
    string steeringFile = steeringPath + "/" + steeringName;

    string steeringDefaultPath = GetEnv("STEERINGPATH");

    if( steeringDefaultPath.size() > 0 ) {
        steeringFile=steeringDefaultPath+"/"+steeringName;
        if (debug) cout<<" makegridfromsherpa::main: STEERINGPATH environment variable found, using path: "<<steeringFile<<endl;
    } else {
        if (debug)cout<<" makegridfromsherpa::main: STEERINGPATH environment varaible not set, using default: "<<steeringFile<<endl;
    }



    //attempt to use EnvVar to find PDFsets, otherwise use default folder path PDFsets/ in current dir
    string pdfSetName = "CT10.LHgrid"; //"MSTW2008nlo90cl.LHgrid";
    string pdfSetPath = "PDFsets";
    string pdfSetFile = pdfSetPath + "/" + pdfSetName;

    string pdfSetDefaultPath = GetEnv("LHAPATH");

    if( pdfSetDefaultPath.size() > 0 ) {
        pdfSetFile = pdfSetDefaultPath + "/" + pdfSetName;
        if (debug) cout<<" makegridfromsherpa::main: LHAPATH environment varaible found, using path: "<<pdfSetFile<<endl;
    } else {
        if (debug) cout<<" makegridfromsherpa::main: LHAPATH environment varaible not set, using default: "<<pdfSetFile<<endl;
    }



    //allow for passing of another steering file name (default=inputname) and a different number of events(default uses all_events)
    if ( argc>1 ) {
        steeringFile = string(argv[1]);
        if (debug) std::cout << " makegridfromsherpa::main: Reading steering file " << steeringFile << std::endl;
    }
    if ( argc>2 ) {
        nevmax = atoi(argv[2]);
        if (debug) std::cout << " makegridfromsherpa::main: Reading Number of events " << nevmax << std::endl;
    }


    //**starting and ending indexes for histogram looping over ntuples B and R
    const int startIndex = 0;
    const int endIndex   = 1; //change to choose ntups to go over


    //Create uniquely named grids for each Type: 0-B, 1-R, 2-RthenB
    MyGrid *mygrid[endIndex];
    for(int i=0; i<endIndex; i++) {
        string version=ntup_names[i];
        if (debug) cout<<" makegridfromsherpa::main: Creating grid using steeringFile: "<<steeringFile<<", version: "<<version<<endl;
        mygrid[i] = new MyGrid(steeringFile, version);

        /*
        double facMin =  1.e20;
        double facMax = -1.e20;

        double renMin =  1.e20;
        double renMax = -1.e20;

        string NtupName=mygrid[histoIndex]->GetInputNtupDir();
        if ( histoIndex==i_R || histoIndex==i_B )
            NtupName += "/NTuple" + ntup_names[histoIndex] + "-like";
        else if ( histoIndex==i_RB )
            NtupName += "/NTuple_*-like*";
        NtupName+=".root";


        //get the Max/Min of fac and ren scales across all events to fill grid with correct Q2 Low and Up
        GetRenAndFacMaxAndMins( NtupName, &facMin, &facMax, &renMin, &renMax );


        double q2low = facMin;
        double q2up  = facMax;
        int nQ2bins  = 1;

        if (fabs(facmin-facmax)<1.e-6) {
            q2low   = q2low - q2low / 100.;
            q2up    = q2up + q2up / 100.;
            nQ2bins = 1;
        } else  {
            std::cout<<" makegridfromsherpa::main: NEED to develop an algorithm for this case "<<std::endl;
        }

        mygrid[i]->SetQ2Low (q2low*q2low);
        mygrid[i]->SetQ2Up  (q2up*q2up);
        mygrid[i]->SetQ2bins(nQ2bins);

        if(debug)
            std::cout<<" makegridfromsherpa::main: SetQ2Low: "<<q2low*q2low
                    <<", SetQ: "<<q2up*q2up
                    <<", nQ2bins: "<<nQ2bins<<std::endl;
        */
    }


    const int NGrid=mygrid[0]->GetNGrid(); //NGrid happens to be the same for all Types


    //create histograms to store tests:
    TH1D *htest1[endIndex][NGrid];  //htest1[0][]==R-Type, htest1[1][]==B-Type, htest1[2][]==RB-Type, scaled by uncorr events

    TH1D *htest3[endIndex][NGrid];  //filled in NLO, reference for htest1 to make sure NLO is filled correctly

    TH1D *htestRB[NGrid];           //holds the results of all htest1[i_R] + htest1[i_B]

    TH1D *href[endIndex][NGrid];    //holds histogram references from mygrid[x] to ensure they are equal to htest1[x]

    std::vector<int> nSubCount = std::vector<int>(25,0); //replace 25 by nprocs from lumi

    //maintain a count of the number of events run over for each histogram type we are testing
    int eventCount[endIndex];
    int uncorrEventCount[endIndex];
    for( int i=0 ; i<endIndex ; i++ ) {
        uncorrEventCount[i]=0;
        eventCount[i]=0;
    }


    //
    //Loop over different histograms to test different *-like* approaches
    //
    for( int histoIndex=startIndex ; histoIndex<endIndex ; histoIndex++ )
    {
        if (debug) cout<<" makegridfromsherpa::main: Starting histoIndex loop: "<<histoIndex<<endl;
        TChain *fChain= new TChain("t3");
        
        
        
        //read in each NTuple for each test histogram
        string NtupName=mygrid[histoIndex]->GetInputNtupDir();

        if ( histoIndex==i_R || histoIndex==i_B )
            NtupName += "/NTuple" + ntup_names[histoIndex] + "-like";
        else if(histoIndex==i_RB)
            NtupName += "/NTuple_*-like*";
        NtupName += ".root";

        if (debug) cout<<" makegridfromsherpa::main: Opening "<<NtupName<<endl;


        
        fChain->Add( TString( NtupName ) );
        if (debug) fChain->Print();

        t3 t(fChain);


        
        // set-up jet algorithm via fastjet wrapper class
        if (debug) cout<<" makegridfromsherpa::main: Set up jet algorithm:"<<endl;
        fjClustering* jetclus = new fjClustering( fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best );

        //create a new event for each event. This could be done better no?
        MyEvent* myevent = new MyEvent();


        Long64_t nentries = fChain->GetEntries();
        if (debug) cout<<" makegridfromsherpa::main: Number of events= "<<nentries<<endl;


        //
        // set-up test histogram
        //
        const int NGrid=mygrid[histoIndex]->GetNGrid();

        for ( int igrid=0 ; igrid<NGrid ; igrid++ )
        {
            string htest1name = "htest1" + ntup_names[histoIndex];
            htest1[histoIndex][igrid] = mygrid[histoIndex]->GetHisto( igrid, htest1name );
            htest1[histoIndex][igrid]->SetTitle( htest1name.c_str() );
            htest1[histoIndex][igrid]->SetName ( htest1name.c_str() );
            htest1[histoIndex][igrid]->SetLineColor(histoIndex+1);

            
            string htest3name = "htest3" + ntup_names[histoIndex] + "_NLO";
            htest3[histoIndex][igrid] = mygrid[histoIndex]->GetHisto( igrid, htest3name );
            htest3[histoIndex][igrid]->SetTitle( htest3name.c_str() );
            htest3[histoIndex][igrid]->SetName ( htest3name.c_str() );
            htest3[histoIndex][igrid]->SetLineColor(histoIndex+1);

            if(debug) cout<<" makegridfromsherpa::main: Got histogram: "<<string("htest1_"+histoIndex)<<", with line color: "<<histoIndex+1<<endl;
        }

        if (debug) cout<<" makegridfromsherpa::main: Loop over events for: "<<ntup_names[histoIndex]<<endl;



        //
        // Event loop over all events
        //
        Long64_t nbytes = 0, nb = 0;
        LHAPDF::initPDFSet( pdfSetFile.c_str(), 0 );

        for ( Long64_t jentry=0 ; jentry<nentries && jentry<nevmax; jentry++ )
        {
            t.GetEntry(jentry);
            nbytes += nb;
            //t.Show();

            Int_t np=(Int_t)t.nparticle;
            //inform user every NUMBER of events completed to illustrate progress. May need to be adjusted depending on code exeuction speed
            if (jentry%1000==0 && jentry>0) {
                cout<<" Type: "<<ntup_names[histoIndex]<<"(histoIndex:"<<histoIndex<<")"<<", events processed "<<jentry<<"  of  ";
                if(argc>2) cout<<nevmax<<". ";
                cout<<"TOTAL EVENTS: "<<nentries<<""<<endl;
            }


            if (debug) {
                // Documentation see
                // http://sherpa.hepforge.org/doc/SHERPA-MC-2.0.0.html
                //

                cout<<"\nmakegridfromsherpa::main: SHERPA EVENT INFO >>>>>>>>" <<endl;
                cout<<"  Event id = "<<t.id<<endl; //Event ID to identify correlated real sub-events.

                cout<<"  Incoming parton id1 = "<<t.id1<<", id2 = "<<t.id2<<endl; // PDG code of incoming parton 1/2
                for (int i=0; i<np; i++)
                    cout<<"  Outcoming parton kf["<<i<<"]= "<<t.kf[i]<<endl; //

                cout<<"  x1 = "<<t.x1<<", x2 = "<<t.x2<<endl;     // Bjorken-x of incoming parton 1/2
                cout<<"  x1p= "<<t.x1p<<", x2p= "<<t.x2p<<endl;   // x’ for I-piece of incoming parton 1/2
                // Factorisation and Factorisaton scale.
                cout<<"  fac_scale= "<<t.fac_scale<<", ren_scale= "<<t.ren_scale
                    <<" alphas= "<<t.alphas<<endl; //" alphasPower= "<<(Char_t)t.alphasPower<<endl;
                //
                // nuwgt Number of additional ME weights for loops and integrated subtraction terms.
                // usr_wgt[nuwgt] Additional ME weights for loops and integrated subtraction terms.
                //8.8.4.1 Computing (differential) cross sections of real correction events with statistical errors
                cout<<"  nuwgt= "<<t.nuwgt<<endl;
                // units are GeV and pb
                // weight    Event weight, if sub-event is treated independently.
                // weight2   Event weight, if correlated sub-events are treated as single event.
                // see
                cout<<"  weight= "<<t.weight<<", weight2= "<<t.weight2<<endl;
                // me_wgt    ME weight (w/o PDF), corresponds to ’weight’.
                // me_wgt2   ME weight (w/o PDF), corresponds to ’weight2’.
                cout<<"  me_wgt= "<<t.me_wgt<<", me_wgt2= "<<t.me_wgt2<<endl;
                cout<<"  b_part= "<<t.part<<endl;
                cout<<"  b_nparticle= "<<t.nparticle<<endl; //Number of outgoing partons.
                //cout<<" np= "<<np<<endl;
                //Int_t np=t.b_nparticle;
                cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<endl;
            }





            int iorder=int(t.alphasPower);
            int itype =int(t.part[0]);

            jetclus->ClearJets();
            myevent->ClearEvent();;
            myevent->SetCMS(7000.);



            //pepare coeficients
            int Wsize=13;
            double xf1[Wsize];
            double xf2[Wsize];

            for(int i=0; i<Wsize; i++) {
                xf1[i]=0.0; //reset
                xf2[i]=0.0; //reset
            }


            //use evolvepdf_ to get weights in f1, f2
            getPDF( t.x1, t.fac_scale, xf1 );
            getPDF( t.x2, t.fac_scale, xf2 );


            //particle num convention conversion from SHERPA
            int id1 = t.id1;
            int id2 = t.id2;
            if(t.id1==21) id1=0;
            id1 = id1+6;
            if(t.id2==21) id2=0;
            id2 = id2+6;

            if(debug) std::cout<<" makegridfromsherpa::main: iorder: "<<iorder<<std::endl;


            //******FORCING REJECTION FOR TESTING
            if(iorder<=2) continue; //Accepting only NLO
            //if(iorder!=2) continue; //Accepting only LO

            //if(!(id1==6 && id2==6)) continue; //Accepting only certain particles
            //if(id1==6 || id2==6) continue;
            //******FORCING REJECTION FOR TESTING




            //weights are now set explicitly bellow depending on event order being LO or NLO
            //myevent->SetWeight(wgt);
            myevent->SetXSection( t.weight2 );

            myevent->SetOrder( iorder );
            myevent->SetType( itype );
            myevent->SetEventId( t.id );


            myevent->SetQ2( t.fac_scale * t.fac_scale ); // applgrid takes sqrt(fac_scale) when evaluating pdf and alphas


            //
            // fill incoming partons
            //
            double pxin=0.;
            double pyin=0.;
            double pzin=0.;
            double Ein=0.;

            for ( int ip=0 ; ip<np ; ip++ ) {
                pxin += t.px[ip];
                pyin += t.py[ip];
                pzin += t.pz[ip];
                Ein  += t.E[ip];
            }

            double shat = Ein*Ein - ( pxin*pxin + pyin*pyin + pzin*pzin );
            double s    = shat / ( t.x1*t.x2 );
            myevent->SetCMS( sqrt(s) );



            int pid   = t.id1;
            double ep = ( sqrt(s) / 2.0 );

            if( pid == 21 ) pid = 0; //conversion from sherpa gluon to appl_grid convention
            myevent->push_back( pxin, pyin, t.x1 * ep, t.x1 * ep, pid);
            //cout<<" makegridfromsherpa::main: pid1: "<<pid<<endl;

            pid=t.id2;
            if( pid == 21 ) pid = 0; //conversion from sherpa gluon to appl_grid convention
            myevent->push_back( pxin, pyin, -t.x2 * ep, t.x2 * ep, pid );
            //cout<<" makegridfromsherpa::main: pid2: "<<pid<<endl;




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


            


            //Need subprocess ID to compute the number of pairs as part of a combo. Refer to Lumi_pdf implimentation
            int subProcID = mygrid[histoIndex]->GetDecideSubProcess( t.id1==21? 0:t.id1, t.id2==21? 0:t.id2 );
            double npairs = mygrid[histoIndex]->GetNSubProcessPairs( subProcID );


            //weight multiplication factor needed by ApplGrid when setting Grid weights when filling mygrid
            wgt_fac = pow( ( 2.0 * PI ) / t.alphas, iorder );


            //needed for filling reference histos htest1 and htest2
            double obs = myevent->GetInvariantMass12(); //only needs to be computed once per events with the same id




            int nWgts=13;   //number of flavors tBar(0)...gluon(6)...t(13)
            
            
            //coefficients used in Stefan's code.
            double f1  [nWgts];
            double f2  [nWgts];
            double xf1p[nWgts];
            double xf2p[nWgts];
            for( int iproc=0 ; iproc<nWgts ; iproc++ ) {
                xf1[iproc] = 0.0; //emptied
                xf2[iproc] = 0.0; //emptied
            }


            //evolve to get weights in xf1, xf2
            getPDF( t.x1, t.fac_scale, xf1 );
            getPDF( t.x2, t.fac_scale, xf2 );


            //Filling f1 and f2 now for easy access in Stefan's code
            for( int iproc=0 ; iproc<nWgts ; iproc++ ) {
                f1[iproc] = xf1[iproc] / t.x1;
                f2[iproc] = xf2[iproc] / t.x2;
            }


            //computing f1 from xf1 and f2 from xf2 for later usage (and readability) conveyance
            double fa = f1[id1]; // (xf1/x1)[id1]
            double fb = f2[id2]; // (xf2/x2)[id2]
            double lr = 0.0;     //NOTE: NOT SET CURRENTLY!
            double lf = 0.0;     //NOTE: NOT SET CURRENTLY!


            double wgtCheck = 0;


            //fill weight depending on LO or NLO
            if(iorder==2) {
                //
                // iorder = 2 = LO
                //

                //compute weight to fill my event with when LO
                double wgt = (t.me_wgt2 * wgt_fac) / npairs;

                wgtCheck = t.me_wgt2 * fa * fb;

                std::cout<<" wgtCheck: "<<wgtCheck<<", t.weight2: "<<t.weight<<", iorder: "<<iorder<<std::endl;


                //set event based on LO
                myevent->SetWeight(wgt); //dealing entirely with weight2
                myevent->SetX1(t.x1);
                myevent->SetX2(t.x2);

                mygrid[histoIndex]->fill(myevent);
                eventCount[histoIndex]++;
                //if(debug) std::cout<<" makegridfromsherpa::main: LO Weight set to: "<<wgt<<std::endl;
            }
            else
            {
                //
                // iorder = 3 = NLO
                //
                bool NLOdebug = false;   //for activating NLO debug statements
                int currentSubProc = -1; //holder for result of deciding the sub process of two partons

                
                //w used for setting weights
                double w[9];
                w[0] = t.me_wgt + t.usr_wgts[0] * lr + t.usr_wgts[1] * lr * lr / 2.0; //special case for w[0]??


                bool wnz=false;
                for ( int i=1 ; i<9 ; ++i )
                {
                    w[i] = t.usr_wgts[i+1] + t.usr_wgts[i+9] * lf;
                    if (w[i]==0) wnz=true;
                }


                double  wgt      = 0.0; //weight for event
                double  htestWgt = 0.0; //weight for htest histo
                nSubCount.clear();



                //always fill w[0].
                //Note this will be the first fill for B-like but the only fill for R-like
                htestWgt = ( w[0] * fa * fb );
                wgt      = ( w[0] * wgt_fac ) / npairs;
                double wgtTest  = ( w[0] * fa * fb );


                myevent->SetWeight( wgt );
                myevent->SetX1( t.x1 );
                myevent->SetX2( t.x2 );
                mygrid[histoIndex]->fill( myevent );
                //std::cout<<"TEST: P0: htestWgt: "<<htestWgt<<std::endl;
                for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, htestWgt );

                eventCount[histoIndex]++;



                if(NLOdebug)std::cout<<"TEST: id1: "<<id1<<", id2: "<<id2<<", wnz: "<<wnz<<std::endl;

                if (wnz==true)
                {
                    wgt = 0.0; //reset

                    double  faq  = 0.0;
                    double  faqx = 0.0;

                    double  fag  = 0.0;
                    double  fagx = 0.0;

                    double  fbq  = 0.0;
                    double  fbqx = 0.0;

                    double  fbg  = 0.0;
                    double  fbgx = 0.0;

                    //
                    //calculating coeficients to fill ref histos before filling/setting up mygrid
                    // literally a copy of Stefan's code.
                    //
                    if ( id1 != GLUON ) {
                        //Parton 1 - QUARK
                        faq = fa;
                        fag = f1[6];

                        getPDF( t.x1/t.x1p, t.fac_scale, xf1p );

                        faqx = xf1p[id1] / t.x1;
                        fagx = xf1p[6]   / t.x1;
                    }
                    else {
                        //Parton 1 - GLU
                        fag=fa;
                        for ( int i=1 ; i<nWgts-1 ; ++i)
                            if( i!=GLUON ) faq += f1[i];

                        getPDF( t.x1/t.x1p, t.fac_scale, xf1p );

                        fagx = (xf1p[id1] / t.x1);
                        for ( int i=1 ; i<nWgts-1 ; ++i )
                            if( i != GLUON ) faqx += (xf1p[i]/t.x1);
                    }
                    if ( id2 != GLUON ) {
                        //Parton 2 - QUARK
                        fbq = fb;
                        fbg = f2[GLUON];

                        getPDF( t.x2/t.x2p, t.fac_scale, xf2p );

                        fbqx = xf2p[id2] / t.x2;
                        fbgx = xf2p[GLUON] / t.x2;
                    }
                    else {
                        //Parton 2 - GLU
                        fbg = fb;
                        for ( int i = 1 ; i<nWgts-1 ; ++i)
                            if( i != GLUON ) fbq += f2[i];

                        getPDF( t.x2/t.x2p, t.fac_scale , xf2p);

                        fbgx = (xf2p[id2] / t.x2);
                        for ( int i=1 ; i<nWgts-1 ; ++i )
                            if( i != GLUON ) fbqx += (xf2p[i] / t.x2);
                    }


                    //wgtCheck += ( faq*w[1] + faqx*w[2] + fag*w[3] + fagx*w[4] ) * fb;
                    //wgtCheck += ( fbq*w[5] + fbqx*w[6] + fbg*w[7] + fbgx*w[8] ) * fa;
                    //std::cout<<"wgtCheck: "<<wgtCheck<<" t.weight: "<<t.weight<<std::endl;






                    if(NLOdebug) std::cout<<"TEST: jentry: "<<jentry<<", id1: "<<id1<<" , id2: "<<id2<<std::endl;

                    //
                    //Filling grid and ref histos from coeficients and 'w' just generated
                    //
                    if ( id1 != GLUON ) {
                        //
                        // parton 1 = QUARK
                        //

                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[1] * wgt_fac ) / npairs; //<--weights determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*faq*w[1] );




                        myevent->SetX1( t.x1 / t.x1p );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[2] * ( 1.0 / t.x1p ) * wgt_fac ) / npairs; //<--weights determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*faqx*w[2] );




                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, 0 );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[3]  * wgt_fac ) / npairs; //<--weights determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*fag*w[3] );




                        myevent->SetX1( t.x1 / t.x1p );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, 0 );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[4] * ( 1.0 / t.x1p ) * wgt_fac ) / npairs; //<--weights determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*fagx*w[4] );


                    }
                    else {
                        //
                        // parton 1 = GLU
                        //
                        if(NLOdebug) std::cout<<"TEST: P1 = GLU"<<std::endl;

                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );

                        //count the number of subprocs that make up
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav )
                        {
                            if( iflav != GLUON )
                            {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );

                                if( currentSubProc != -1 )
                                    nSubCount[currentSubProc]++; //count the number of occurences of that particular subproc
                            }
                        }

                        //filling grids for all subprocs
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav )
                        {
                            if( iflav != GLUON )
                            {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );

                                if( currentSubProc != -1 )
                                {
                                    npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );
                                    wgt = w[1] * wgt_fac * ( 1.0 / npairs ) * ( 1.0 / nSubCount[currentSubProc] ); //<--weights determined


                                    myevent->SetID( 0, (iflav-6) );
                                    myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                                    myevent->SetWeight( wgt );
                                    if(NLOdebug) std::cout<<"P1-21: filling ..."<<std::endl;
                                    myevent->SetXSection( 0 ); //hack for not filling ref histos
                                    mygrid[histoIndex]->fill( myevent );
                                }
                            }
                        }

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*faq*w[1]  );

                        for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use


                        //set x1 and x2 for all of the upcoming subprocs event fills
                        myevent->SetX1( t.x1 / t.x1p );
                        myevent->SetX2( t.x2 );


                        //count the number of subprocs that make up
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );

                                if( currentSubProc != -1 )
                                    nSubCount[currentSubProc]++;
                            }
                        }

                        //filling grids for all subprocs
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );
                                if(NLOdebug) std::cout<<"P2-22: flv1: "<<(t.id1==21? 0:t.id1)<<", flv2: "<<(iflav-6)<<", DecideSubProc: "<<currentSubProc<<std::endl;

                                if( currentSubProc != -1 ) {
                                    npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );
                                    wgt = w[2] * wgt_fac * ( 1 / t.x1p ) * ( 1.0 / npairs) * ( 1.0 / nSubCount[currentSubProc] ); //<--weight determined

                                    myevent->SetID( 0, (iflav-6) );
                                    myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                                    myevent->SetWeight( wgt );
                                    if(NLOdebug) std::cout<<"P1-22: filling ..."<<std::endl;
                                    myevent->SetXSection( 0 ); //hack for not filling ref histos
                                    mygrid[histoIndex]->fill( myevent );
                                }
                            }
                        }

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*faqx*w[2] );

                        for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use


                        
                        //glu components 
                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, 0 ); //glu
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = w[3] * wgt_fac * ( 1.0 / npairs);//<--weight determined
                        myevent->SetWeight( wgt );
                        mygrid[histoIndex]->fill( myevent );

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*fag*w[3] );


                        
                        myevent->SetX1( t.x1 / t.x1p );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, 0 ); //glu
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = w[4] * wgt_fac * ( 1.0 / t.x1p ) * ( 1.0 / npairs);//<--weight determined
                        myevent->SetWeight( wgt );
                        mygrid[histoIndex]->fill( myevent );

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fb*fagx*w[4] );

                    }
                    if ( id2 != GLUON ) {
                        //
                        // parton 2 = QUARK
                        //
                        if(NLOdebug) std::cout<<"TEST: P2 = QUARK"<<std::endl;
                        
                        
                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[5] * wgt_fac ) / npairs;//<--weight determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbq*w[5] );



                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 / t.x2p );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[6] * ( 1.0 / t.x2p ) * wgt_fac ) / npairs;//<--weight determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbqx*w[6] );



                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, 0 );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[7]  * wgt_fac ) / npairs;//<--weight determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbg*w[7] );



                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 / t.x2p );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, 0 );
                        myevent->SetXSection( 0 ); //hack for not filling ref histos
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = ( w[8] * ( 1.0 / t.x2p ) * wgt_fac ) / npairs;//<--weight determined
                        myevent->SetWeight( wgt );

                        mygrid[histoIndex]->fill( myevent );
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbgx*w[8] );

                    }
                    else {
                        //
                        // parton 2 = GLU
                        //
                        if(NLOdebug) std::cout<<"TEST: P2 = GLU"<<std::endl;

                        
                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );

                        //count the number of subprocs that make up
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );

                                if( currentSubProc != -1 ) 
                                    nSubCount[currentSubProc]++;
                            }
                        }

                        //filling grids for all subprocs
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );

                                if( currentSubProc != -1 ) {
                                    npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );
                                    wgt = w[5] * wgt_fac * ( 1.0 / npairs ) * ( 1.0 / nSubCount[currentSubProc] );//<--weight determined

                                    myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                                    myevent->SetID( 1, (iflav-6) );
                                    myevent->SetWeight( wgt );
                                    if(NLOdebug) std::cout<<"P2-21: filling ..."<<std::endl;
                                    myevent->SetXSection( 0 ); //hack for not filling ref histos
                                    mygrid[histoIndex]->fill( myevent );
                                }
                            }
                        }
                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbq*w[5] );

                        for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use




                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 / t.x2p );


                        //count the number of subprocs that make up
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {
                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );

                                if( currentSubProc != -1 )
                                    nSubCount[currentSubProc]++;
                            }
                        }

                        //filling grids for all subprocs
                        for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
                            if( iflav != GLUON ) {

                                currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );
                                if(NLOdebug) std::cout<<"P2-22: flv1: "<<(iflav-6)<<", flv2: "<<(t.id2==21? 0:t.id2)<<", DecideSubProc: "<<currentSubProc<<std::endl;

                                if( currentSubProc != -1 ) {
                                    npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );
                                    wgt = w[6] * wgt_fac * ( 1.0 / npairs) * ( 1 / t.x2p ) * ( 1.0 / nSubCount[currentSubProc] );//<--weight determined

                                    myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                                    myevent->SetID( 1, (iflav-6) );
                                    myevent->SetWeight( wgt );
                                    if(NLOdebug) std::cout<<"P2-22: filling ..."<<std::endl;
                                    myevent->SetXSection( 0 ); //hack for not filling ref histos
                                    mygrid[histoIndex]->fill( myevent );
                                }
                            }
                        }

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbqx*w[6] );

                        for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use



                        //glu component
                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, 0 ); //glu
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = w[7] * wgt_fac * ( 1.0 / npairs);//<--weight determined
                        myevent->SetWeight( wgt );
                        mygrid[histoIndex]->fill( myevent );

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbg*w[7] );


                        myevent->SetX1( t.x1 );
                        myevent->SetX2( t.x2 / t.x2p );
                        myevent->SetID( 0, (t.id1==21? 0:t.id1) );
                        myevent->SetID( 1, 0 ); //glu
                        currentSubProc = mygrid[histoIndex]->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
                        npairs = mygrid[histoIndex]->GetNSubProcessPairs( currentSubProc );

                        wgt = w[8] * wgt_fac * ( 1.0 / t.x2p ) * ( 1.0 / npairs);//<--weight determined
                        myevent->SetWeight( wgt );
                        mygrid[histoIndex]->fill( myevent );

                        for( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[histoIndex][igrid]->Fill( obs, fa*fbgx*w[8] );
                    }

                    //PERSONAL REFERENCE FOR PICKING WHAT TO FILL EACH htest at what time (from stefan's code)

                    wgtTest += ( faq*w[1] + faqx*w[2] + fag*w[3] + fagx*w[4] ) * fb;
                    wgtTest += ( fbq*w[5] + fbqx*w[6] + fbg*w[7] + fbgx*w[8] ) * fa;

                    //if(NLOdebug) std::cout<<"TEST: wgtTest: "<<wgtTest<<", t.weight2: "<<t.weight2<<std::endl;

                    if(NLOdebug)
                        std::cout<<"TEST"
                                 <<": faq:"  <<faq
                                 <<", fbq:"  <<fbq
                                 <<", faqx:" <<faqx
                                 <<", fbqx:" <<fbqx
                                 <<", fag:"  <<fag
                                 <<", fbg:"  <<fbg
                                 <<", fagx:" <<fagx
                                 <<", fbgx:" <<fbgx<<std::endl;

                    if(NLOdebug) std::cout<<"test: Current event number finishing "<<jentry<<std::endl;
                }
                else {
                    std::cout<<" makegridfromsherpa::main: ERROR: wnz was false? Running over 'R'? "<<std::endl;
                    exit(0); //TEST temp
                }
            }






            //
            // fill each test histogram
            //
            obs = myevent->GetInvariantMass12(); // to be replaced by GetObservable from steering
            for ( int igrid=0 ; igrid<mygrid[histoIndex]->GetNGrid() ; igrid++ )
            {
                if (mygrid[histoIndex]->eventcuts(myevent,igrid)==false) continue;

                htest1[histoIndex][igrid]->Fill( obs, t.weight2 );
                
                //htest1[histoIndex][igrid]->Fill( obs, wgtCheck ); //testing to make sure wgt is filled correctly
                
            } //end loop over grid

            uncorrEventCount[histoIndex]++; //keep count of event for each type, 0-B, 1-R, 2-RthenB



        } //end loop over events

        if (debug)
            std::cout<<"\nmakegridfromsherpa::main: Finished running over events!\n"<<std::endl;

        //
        // get and set up test external histograms
        //
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++)
        {
            href[histoIndex][igrid]=(TH1D*)mygrid[histoIndex]->GetReference(igrid);
            if (!href[histoIndex][igrid]) {
                cout<<" makegridfromsherpa::main: Reference from grid not found ! "<<endl;
                exit(0);
            }
            else {

                TString nameTitle="internal_href"+ntup_names[histoIndex];
                href[histoIndex][igrid]->SetTitle(nameTitle);
                href[histoIndex][igrid]->SetName(nameTitle);
                href[histoIndex][igrid]->SetLineColor(7); //SKY BLUE

            }

        }


        if (debug) {
            cout<<" makegridfromsherpa::main: End histIndex loop: "
                <<histoIndex<<" of "<<(endIndex-1)
                <<", events this loop: "<<uncorrEventCount[histoIndex]<<endl;
        }

    } //end of loop over all htest1 histograms



    //write grid to scale the internal reference histos and save them to *.root files
    if (debug) cout<<" makegridfromsherpa::main: Scaling and writing mygrid "<<endl;
    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {

        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++)
            mygrid[histoIndex]->ScaleInternalRefHistos(igrid);

        mygrid[histoIndex]->write_grid();

    }
    if (debug) cout<<" makegridfromsherpa::main: Grids written "<<endl;




    //
    // Scale htest histos
    // Add histograms of MyGrid for R-Like(i_R) to B-Like(i_B). Should be equal to both combined(i_RB)
    //
    if (debug) cout<<" makegridfromsherpa::main: Adding R-Type and B-Type"<<endl;

    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        //scale each histograms by one devided by number of events depending on type
        for(int igrid=0; igrid<NGrid; igrid++) {

            //scale each histograms by one devided by number of events depending on type
            htest1[histoIndex][igrid]->Scale( 1.0 / uncorrEventCount[histoIndex] ); //hR/nR, hB/nB, hBR/nBR, etc
            htest3[histoIndex][igrid]->Scale( 1.0 / uncorrEventCount[histoIndex] );

            href[histoIndex][igrid]->Scale(   1.0 / uncorrEventCount[histoIndex] ); //hrefR/nR, hrefB/nB,hrefRB/nRB, etc


            cout<<" makegridfromsherpa::main: htest1 printing..."<<endl;
            htest1[histoIndex][igrid]->Print("all");
            htest3[histoIndex][igrid]->Print("all");


            cout<<" makegridfromsherpa::main: hnorm printing (htest1/binswidth)... "<<endl;
            TH1D *hnorm1=(TH1D*)htest1[histoIndex][igrid]->Clone("hnorm");
            mygrid[histoIndex]->DivideByBinWidth(hnorm1);
            hnorm1->Print("all");


            //sum the scaled R and B-type together, "should" be the same as htest[i_RB] where R and B are done together
            if(histoIndex==i_RB)
            {
                if(startIndex<=i_B && endIndex>=i_R) //make sure that R and B exist before adding them!
                {
                    htestRB[igrid] = mygrid[0]->GetHisto(igrid,string("htestRplusB_")); //mygrid index=0 is arbitrary, chosen to use GetHisto func
                    htestRB[igrid]->Add(htest1[i_R][igrid]);
                    htestRB[igrid]->Add(htest1[i_B][igrid]);
                    htestRB[igrid]->SetLineColor(4); //blue
                }
            }
        }

        if (debug) cout<<" makegridfromsherpa::main: Added R-Type and B-Type!"<<endl;
    }

    /*
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
                cout<<" makegridfromsherpa::main:**NORM CHECK: nR: "<<uncorrEventCount[0]<<", nB: "<<uncorrEventCount[1]<<", nRB: "<<uncorrEventCount[2]<<endl;

                //temporarily changed to a default!!
                evtot      =1;
                evalltot   =1;
                evuncorr   =1;
                evalluncorr=1;

                cout<<" makegridfromsherpa::main: Norm: evalltot= "<<evalltot<<" evalluncorr= "<<evalluncorr<<endl;
                cout<<" makegridfromsherpa::main: Norm: evtot= "<<evtot<<" evuncorr= "<<evuncorr<<endl;

                MyData *mydata=mygrid[histoIndex]->GetMyData(igrid);
                if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
                double yfac=mydata->GetUnitfbFactor();
                double xfac=mydata->GetUnitGeVFactor();
                cout<<" makegridfromsherpa::main: MyGrid::Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


                //update names for normalised ones
                htest1[histoIndex][igrid]->SetTitle(string("htest1"+ntup_names[histoIndex]+"-norm").c_str());
                htest1[histoIndex][igrid]->SetName(string("htest1"+ntup_names[histoIndex]+"-norm").c_str());

                href[histoIndex][igrid]->SetTitle(string("internal_href"+ntup_names[histoIndex]+"-norm").c_str());
                href[histoIndex][igrid]->SetName(string("internal_href"+ntup_names[histoIndex]+"-norm").c_str());


                mygrid[histoIndex]->Normalise(htest1[histoIndex][igrid],evtot*yfac,xfac,true);   //normalise 0-B, 1-R, and 2-RthenB Type
                mygrid[histoIndex]->Normalise(href[histoIndex][igrid],evtot*yfac,xfac,true);     //normalise hrefB, hrefR, and hrefRthenB

                if(histoIndex==i_RB) {
                    htestRB[igrid]->SetTitle(string("htestRplusB_-norm").c_str());
                    htestRB[igrid]->SetName(string("htestRplusB_-norm").c_str());
                    mygrid[histoIndex]->Normalise(htestRB[igrid],evtot*yfac,xfac,true);//normalise histo of scaled R + scaled B together
                }
            }

        } //end of loop over all htest1 histograms
    */



    //
    // Convolute and print all external histograms
    //
    if (debug) cout<< "\n makegridfromsherpa::main: Performing convolute"<<endl;
    //NGrid = mygrid[0]->GetNGrid(); //NGrid will be the same for all grids, so grid[0] is arbirary
    TH1D* convGridHistos[endIndex+1][NGrid];
    TH1D* LOconvGridHistos[endIndex+1][NGrid];
    TH1D* subProcConvGridHistos[endIndex+1][NGrid][121];
    TH1D* LOsubProcConvGridHistos[endIndex+1][NGrid][121];

    int nLoops = 1;


    //creating an identity matrix needed in convolute function - temp printing to make sure it's correct
    std::vector<std::vector<double> > ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));
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
            string filename=mygrid[histoIndex]->GetGridName(igrid);
            filename+=ntup_names[histoIndex]+"_histos.root";
            fout= new TFile(filename.c_str(),"recreate");
            mygrid[histoIndex]->GetGrid(igrid)->setckm(ckm2);



            //NLO convolute for all
            convGridHistos[histoIndex][igrid] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute( getPDF, alphaspdf_, nLoops );
            convGridHistos[histoIndex][igrid]->SetName((TString) ("convolute_for" + ntup_names[histoIndex]));
            convGridHistos[histoIndex][igrid]->SetTitle((TString) ("convolute_for" + ntup_names[histoIndex]));
            convGridHistos[histoIndex][igrid]->SetLineColor(kBlue);
            convGridHistos[histoIndex][igrid]->Write();



            //LO convolute for all
            cout<<" calling convolute: "<<endl;
            LOconvGridHistos[histoIndex][igrid] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute( getPDF, alphaspdf_, 0 );
            LOconvGridHistos[histoIndex][igrid]->SetName((TString) ("LOconvolute_for" + ntup_names[histoIndex]));
            LOconvGridHistos[histoIndex][igrid]->SetTitle((TString) ("LOconvolute_for" + ntup_names[histoIndex]));
            LOconvGridHistos[histoIndex][igrid]->SetLineColor(kBlue);
            LOconvGridHistos[histoIndex][igrid]->Write();



            for(int isubproc=0; isubproc<mygrid[histoIndex]->GetNSubProcess(igrid);  isubproc++) {

                //NLO convolute for subprocs
                subProcConvGridHistos[histoIndex][igrid][isubproc] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute_subproc(isubproc, getPDF, alphaspdf_, nLoops );
                string sub_proc_hist_name="subProc-"+to_string(isubproc)+"-convolute_for"+ntup_names[histoIndex];
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetName((TString) (sub_proc_hist_name));
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name));
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetLineColor(kBlue);



                //LO convolute for subprocs
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc] = (TH1D*)mygrid[histoIndex]->GetGrid(igrid)->convolute_subproc(isubproc, getPDF, alphaspdf_, 0 );
                sub_proc_hist_name="LOsubProc-"+to_string(isubproc)+"-convolute_for"+ntup_names[histoIndex];
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->SetName((TString) (sub_proc_hist_name));
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name));
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->SetLineColor(kBlue);
            }



            //get information needed to perform mygrid normalisation calls
            MyData *mydata=mygrid[histoIndex]->GetMyData(igrid);
            if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
            double yfac=mydata->GetUnitfbFactor();
            double xfac=mydata->GetUnitGeVFactor();
            if (debug) cout<<" makegridfromsherpa::main: Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


            // convolute histograms do not need to be normalised !
            //save a scaled and normalised version
            //LOconvGridHistos[histoIndex][igrid]->Scale(1.0/uncorrEventCount[histoIndex]);
            //LOconvGridHistos[histoIndex][igrid]->Write();

            /*
                  mygrid[histoIndex]->Normalise(LOconvGridHistos[histoIndex][igrid],yfac,xfac,true);
                  LOconvGridHistos[histoIndex][igrid]->SetName((TString) ("LOconvolute_for" + ntup_names[histoIndex]+"-norm"));
                  LOconvGridHistos[histoIndex][igrid]->SetTitle((TString) ("LOconvolute_for" + ntup_names[histoIndex]+"-norm"));
                  LOconvGridHistos[histoIndex][igrid]->Write();
            */

            /*
            // sub process convolute histograms do not need to be normalised !
            for(int isubproc=0; isubproc<mygrid[histoIndex]->GetNSubProcess(igrid);  isubproc++) {

                //NLO convolute for subprocs
                // convolute histograms do not need to be normalised !
                subProcConvGridHistos[histoIndex][igrid][isubproc]->Write();
                string sub_proc_hist_name="subProc-"+to_string(isubproc)+"-convolute_for"+ntup_names[histoIndex];

                mygrid[histoIndex]->DivideByBinWidth(subProcConvGridHistos[histoIndex][igrid][isubproc]);
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetName((TString) (sub_proc_hist_name+"-norm"));
                subProcConvGridHistos[histoIndex][igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name+"-norm"));



                //LO convolute for subprocs
                // convolute histograms do not need to be normalised !
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->Write();
                sub_proc_hist_name="LOsubProc-"+to_string(isubproc)+"-convolute_for"+ntup_names[histoIndex];

                mygrid[histoIndex]->DivideByBinWidth(LOsubProcConvGridHistos[histoIndex][igrid][isubproc]);
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->SetName((TString) (sub_proc_hist_name+"-norm"));
                LOsubProcConvGridHistos[histoIndex][igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name+"-norm"));
            }
            */




            //save a scaled and DivideByBinWidth version
            href[histoIndex][igrid]->Write(); //scaled version

            mygrid[histoIndex]->DivideByBinWidth(href[histoIndex][igrid]);
            string hrefName = "internal_href"+ntup_names[histoIndex]+"-norm";
            href[histoIndex][igrid]->SetTitle( hrefName.c_str() );
            href[histoIndex][igrid]->SetName(  hrefName.c_str() );
            href[histoIndex][igrid]->Write(); //DivideByBinWidth version

            std::cout<<" makegridfromsherpa::main: Printing LO convolute histo..."<<std::endl;
            LOconvGridHistos[histoIndex][igrid]->Print("all");





            //output origional htest1 (uncorr events events) & a DivideByBinWidth version
            htest1[histoIndex][igrid]->Print("all");
            htest1[histoIndex][igrid]->Write();

            string name1 = "hnorm1"+ntup_names[histoIndex]+"_uncorr_Div";
            TH1D *htest1norm=(TH1D*)htest1[histoIndex][igrid]->Clone( name1.c_str() );
            htest1norm->SetTitle( name1.c_str() );
            htest1norm->SetName(  name1.c_str() );

            mygrid[histoIndex]->DivideByBinWidth(htest1norm); //DivideByBinWidth version
            std::cout<<" makegridfromsherpa::main: Printing htest1 normalised..."<<std::endl;
            htest1norm->Print("all");

            TH1D* ratioLO = divide( LOconvGridHistos[histoIndex][igrid], htest1norm ); //LOconvolute/htestnorm
            if ( ratioLO ) {
                string ratioName = "ratio_LO-conv/htest1norm"+ntup_names[histoIndex];
                ratioLO->SetTitle( ratioName.c_str() );
                ratioLO->SetName(  ratioName.c_str() );
                ratioLO->Print("all");

                ratioLO->DrawCopy();
                gPad->Print("ratioLO.pdf");
            }





            //output htest3 & a DivideByBinWidth version
            htest3[histoIndex][igrid]->Print("all");
            htest3[histoIndex][igrid]->Write(); //scaled version

            string name3 = "hnorm3"+ntup_names[histoIndex]+"_NLO_Div";
            TH1D *htest3norm=(TH1D*)htest3[histoIndex][igrid]->Clone( name3.c_str() );
            htest3norm->SetTitle( name3.c_str() );
            htest3norm->SetName(  name3.c_str() );

            mygrid[histoIndex]->DivideByBinWidth(htest3norm); //DivideByBinWidth version
            std::cout<<" makegridfromsherpa::main: Printing htest3 normalised..."<<std::endl;
            htest3norm->Print("all");

            TH1D* ratioNLO = divide( convGridHistos[histoIndex][igrid], htest3norm ); //NLOconvolute/htestnorm
            if ( ratioNLO ) {
                string ratioName = "ratio_NLO-conv/htest3norm"+ntup_names[histoIndex];
                ratioNLO->SetTitle( ratioName.c_str() );
                ratioNLO->SetName(  ratioName.c_str());
                ratioNLO->Print("all");

                ratioNLO->DrawCopy();
                gPad->Print("ratioNLO.pdf");
            }

        }

        fout->Write(); //writing root file
        fout->Close();
    }





    /*
    //write overlays of histograms you are interested in
    string operation = "./overlay";
    string hfile1 = ""; //"atlas2012_5fb_top_mtt_ljet_B_norm-ghistos.root";
    string hfile2 = ""; //"atlas2012_5fb_top_mtt_ljet_R_histos.root";
    string hname1 = ""; //"convolute_for_R";
    string hname2 = ""; //"internal_href_R";

    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        cout<< " makegridfromsherpa::main: Writing overlays for: "<<ntup_names[histoIndex]<<"(histoIndex:"<<histoIndex<<")"<<endl;
        for(int igrid=0; igrid<NGrid; igrid++)
        {
            hfile1="atlas2012_5fb_top_mtt_ljet"+ntup_names[histoIndex]+"_histos.root";
            hfile2="atlas2012_5fb_top_mtt_ljet"+ntup_names[histoIndex]+"_histos.root";
            hname1="convolute_for"+ntup_names[histoIndex];
            hname2="internal_href"+ntup_names[histoIndex];

            string performThis = operation+" "+hfile1+" "+hfile2+" "+hname1+" "+hname2;
            cout<< " makegridfromsherpa::main: performing: "<<performThis<<endl;

            system(performThis.c_str());
        }
    }
    */

    /*
        //Add R and B grids together whenall needed grid spaces exist
        if(startIndex<=i_B && endIndex>i_RB)
        {
            //after performing convolutes for B, R, RthenB types, add grids and histo for R and B to get RplusB and check it's convolute
          if (debug) cout<<" makegridfromsherpa::main: Adding B-Type mygrid to R-Type mygrid"<<endl;
            mygrid[i_R]->AddGrid(mygrid[i_B]); //STILL NEEDS to handle APPLgrid add/combine so the internals update as well as MyGrid!!
            mygrid[i_R]->SetGridVersionName(string("_RplusB"));
            mygrid[i_R]->write_grid();
        }

    */
    if (debug) cout<<" makegridfromsherpa::main: Normalising Internal Reference histograms "<<endl;
    for(int histoIndex=startIndex; histoIndex<endIndex; histoIndex++)
    {
        for (int igrid=0; igrid<mygrid[histoIndex]->GetNGrid(); igrid++) {
            mygrid[histoIndex]->NormaliseInternalRefHistos(igrid);

            mygrid[histoIndex]->SetGridVersionName(ntup_names[histoIndex]+"_norm");
            mygrid[histoIndex]->write_grid();
        }
    }
    if (debug) cout<<" makegridfromsherpa::main: Internal Reference histograms normalised!"<<endl;

    /*
        //Add R and B grids together whenall needed grid spaces exist
        if(startIndex<=i_B && endIndex>i_RB)
        {
            TH1D* ConvHistoRplusB[mygrid[i_R]->GetNGrid()];
            TH1D* hrefRplusB[mygrid[i_R]->GetNGrid()];

            for(int igrid=0; igrid<mygrid[i_R]->GetNGrid(); igrid++)
            {
                TFile *fout;

                string filename=mygrid[i_R]->GetGridName(igrid);
                filename+="_RplusB_histos.root";
                fout= new TFile(filename.c_str(),"recreate");

                //problem: only performing convolute to R...
                ConvHistoRplusB[igrid] = (TH1D*)mygrid[i_R]->GetGrid(igrid)->convolute( getPDF, alphaspdf_, nLoops );
                ConvHistoRplusB[igrid]->SetName((TString) ("convolute_for_RplusB"));
                ConvHistoRplusB[igrid]->SetTitle((TString) ("convolute_for_RplusB"));
                ConvHistoRplusB[igrid]->SetLineColor(kBlue);

                //NOTE: NEED TO ACCOUNT FOR xfac and yfac being different when adding grid's MyData??
                MyData *mydata=mygrid[i_R]->GetMyData(igrid);
                if (!mydata) cout<<" makegridfromsherpa::main: mydata["<<igrid<<"] not found "<<endl;
                double yfac=mydata->GetUnitfbFactor();
                double xfac=mydata->GetUnitGeVFactor();
                if (debug) cout<<" makegridfromsherpa::main: Normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


                //ConvHistoRplusB[igrid]->Print("all");
                ConvHistoRplusB[igrid]->Scale(1.0/(uncorrEventCount[i_R]+uncorrEventCount[i_B]));
                ConvHistoRplusB[igrid]->Write();
                ConvHistoRplusB[igrid]->Draw();
                mygrid[i_R]->Normalise(ConvHistoRplusB[igrid],yfac,xfac,true);
                ConvHistoRplusB[igrid]->SetName((TString) ("convolute_for_RplusB-norm"));
                ConvHistoRplusB[igrid]->SetTitle((TString) ("convolute_for_RplusB-norm"));
                ConvHistoRplusB[igrid]->Write();
                ConvHistoRplusB[igrid]->Draw();

                //htestRB[igrid]->Print("all");
                htestRB[igrid]->Write();
                htestRB[igrid]->Draw();


                //NOTE: currently the applgrids are not added, only the MyGrid histos, so to get the reference histo for RplusB
                // we have to get each individually and add them
                hrefRplusB[igrid] = (TH1D*)mygrid[i_R]->GetReference(igrid);
                hrefRplusB[igrid]->Add( (TH1D*)mygrid[i_B]->GetReference(igrid) );
                if (!hrefRplusB[igrid]) cout<<" makegridfromsherpa::main: Reference from grid not found ! "<<endl;
                else {
    	      if (debug) cout<<" makegridfromsherpa::main: Reference from grid found ! "<<endl;

                    TString nameTitle="internal_href_RplusB";
                    hrefRplusB[igrid]->SetTitle(nameTitle);
                    hrefRplusB[igrid]->SetName(nameTitle);
                    hrefRplusB[igrid]->SetLineColor(7); //SKY BLUE
                }

                //note the scaling to total events for R and B
                hrefRplusB[igrid]->Scale(1.0/(uncorrEventCount[i_R]+uncorrEventCount[i_B])); //Are these inidividual internal refs scaled??
                mygrid[i_R]->Normalise(hrefRplusB[igrid],yfac,xfac,true);     //normalise hrefB, hrefR, and hrefRB

                hrefRplusB[igrid]->Print("all");
                hrefRplusB[igrid]->Write();
                hrefRplusB[igrid]->Draw();

                TH1D* ratio1 = divide( ConvHistoRplusB[igrid],htestRB[igrid] );
                if ( ratio1 ) {
                    //ratio->Scale(1/uncorrEventCount[i_B]);
                    ratio1->SetTitle(string("ratio1_convolute/htestRplusB").c_str());
                    ratio1->SetName(string("ratio1_convolute/htestRplusB").c_str());
                    ratio1->Print("all");
                    ratio1->Write();
                    ratio1->Draw();
                }
                TH1D* ratio2 = divide( htestRB[igrid], ConvHistoRplusB[igrid] );
                if ( ratio2 ) {
                    //ratio->Scale(1/uncorrEventCount[i_B]);
                    ratio2->SetTitle(string("ratio2_htestRplusB/convolute").c_str());
                    ratio2->SetName(string("ratio2_htestRplusB/convolute").c_str());
                    ratio2->Print("all");
                    ratio2->Write();
                    ratio2->Draw();
                }

                fout->Write();
                fout->Close();
            }
        }
    */


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

    for( int histoIndex=0 ; histoIndex<endIndex ; histoIndex++) {
        std::cout<<"uncorr events for "<<ntup_names[histoIndex]<<": "<<uncorrEventCount[histoIndex]<<std::endl;
        std::cout<<"tot events for "   <<ntup_names[histoIndex]<<": "<<eventCount[histoIndex]<<std::endl;
    }

    //system("root -l test.C");

    return 0;
}
