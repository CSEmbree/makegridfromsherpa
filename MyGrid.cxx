//
//   for sherpa
//
#include <iostream>
using namespace std;
#include <string>
//#include "appl_grid/appl_igrid.h"


#include "MyGrid.h"
/******************************************************************
 ** Method Implementations
 ******************************************************************/
MyGrid::MyGrid(string name, string version)
{

    debug=false;
    bookrefsubprocess=true;

    alluncorrevents  =  0; // uncorrelated events passing cuts (subevents are counted as one event)
    allevents      =  0; // all events passing cuts
    uncorrevents.clear();
    events.clear();
    newevent=false;
    eventid=-99;

    // Grid architecture
    nXbins  = 30;
    xLow    = 1.0e-8, xUp = 1.0;
    xorder  = 5;
    nQ2bins = 5; //<--**changed
    q2Low   = 70.*70.;
    q2Up = 7000.*7000.; //squared
    qorder  = 2;

    iorder  = 0;

    // lowest order inalphas
    lowest_order = 2;
    // how many loops
    nloops = 1;

    processnumber=1;

    pdf_function="";
    subprocesssteername="";
    ntupname="";
    steername=name;
    versionname=version; //optional argument

    cout<<" MyGrid::steering "<<steername<<endl;
    gridname.clear();
    //cout<<" before reading "<<endl;
    this->ReadSteering(name);
    //cout<<" after reading "<<endl;

    bool incljets=false;
    cout<<" MyGrid::Initialize fill for inclusive jets"<<endl;

    pdf_function=subprocesssteername; //subprocesssteername was read into from reading steering file

    this->Initialize();


    const int nSub=nsub;

    weight= new double[nSub]; // needs to be put after ReadSteering


    istylemarkerdata.clear();
    icolormarkerdata.clear();
    mydata.clear();
    //isBooked.clear();


    //save important vector of histograms into a new vector for later retreival when adding another to this one
    std::vector<vhref> referenceHistogramVectors;
    referenceHistogramVectors.clear();
    referenceHistogramVectors.push_back(hreference);        //save pointer to hreference
    referenceHistogramVectors.push_back(hreferencefine);    //save pointer to hreferencefine



    /*
    vector< vector< vector < TH1D* > > > histos = this->GetOutputHistograms();
    cout<<"Vec<vec<vec<hist>>> size: "<<histos.size()<<endl;
    exit(0);
    */

    return;
}

void MyGrid::Initialize() {

    //cout<<"start MyGrid::Initialize()"<<endl;

    //MySubProcess *mySub= new MySubProcess(subprocesssteername);
    //this->SetSubProcess(mySub);

    //cout<<" MyGrid:Initialize Number of grids to produce "<<gridname.size()<<endl;
    ///sherpaw_pdf *mypdf=new sherpaw_pdf(pdf_function,mysub);

    cout << " MyGrid::Initialize: subprocesssteername: " <<subprocesssteername <<  endl;


    /*
        mypdf = new  generic_pdf(subprocesssteername); //<--**

        nsub=mypdf->GetSubProcessNumber(); //<--**
        cout<<" TEST: nsub= "<<nsub<<endl;

        mypdf->SetNProc(nsub); //<--**
        cout<<" MyGrid::Initialize nsub= "<<nsub<<endl; //<--**
    */


    // Loop over grid
    for (int  igrid = 0; igrid < gridname.size(); igrid++) {

        cout<<" MyGrid::Initialize Grid Name "<<this->GetGridName(igrid)<<endl;
        //this->SetStyleMarker (20+igrid);
        //this->SetColorMarker (1+igrid);
        this->SetStyleMarker (20);
        this->SetColorMarker (1);


        MyData *mydatatmp= new MyData;
        if (!mydatatmp) cout<<" MyGrid::Initialize MyData could not be created "<<this->GetGridName(igrid)<<endl;
        else cout<<" MyGrid::Initialize MyData created name= "<<this->GetGridName(igrid)<<endl;

        mydatatmp->ReadData(TString(this->GetDataName(igrid)).Data());
        mydatatmp->SetMarker(20,1);
        mydata.push_back(mydatatmp);

        string fname="output/"+this->GetGridFileName(igrid);
        cout<<" MyGrid::Initialize filename= "<<fname<<endl;

        uncorrevents.push_back(0);
        events.push_back(0);
        this->book_grid(igrid);
        isBooked.push_back(true);





        if (bookrefsubprocess) {

            TString hnameref="href";
            hreference.push_back(this->GetHisto(igrid,hnameref.Data()));
            hnameref+="sum2";
            hreferencesum2.push_back(this->GetHisto(igrid,hnameref.Data()));

            TString hnamereftmp="hreftmp";
            hreferencetmp.push_back(this->GetHisto(igrid,hnamereftmp.Data()));
            hnamereftmp+="sum2";
            hreferencesum2tmp.push_back(this->GetHisto(igrid,hnamereftmp.Data()));

            TString hnamereffine="hreffine"+this->GetGridVersionName();
            TH1D * htmp=this->GetHisto(igrid,hnamereffine.Data());
            if (!htmp) cout<<" MyGrid::Initialize Histo not found "<<hnamereffine.Data()<<endl;
            int nbins=htmp->GetNbinsX();
            double xmin=htmp->GetBinLowEdge(1);
            double xmax=htmp->GetBinWidth(nbins)/2.+htmp->GetBinCenter(nbins);
            nbins*=10;
            //cout<<" MyGrid::Initialize nbins= "<<nbins<<" xmin= "<<xmin<<" xmax= "<<xmax<<endl;
            TString htit=htmp->GetTitle()+TString(" fine bins");
            TH1D *hfine= new TH1D(hnamereffine,htit,nbins,xmin,xmax);
            hreferencefine.push_back(hfine);

            hnamereffine+="sum2";
            TH1D * hfinesum2=(TH1D*) hfine->Clone(hnamereffine);
            hreferencefinesum2.push_back(hfinesum2);

            TString hnamereffinetmp="hreffinetmp";
            TH1D * hfinetmp=(TH1D*) hfine->Clone(hnamereffinetmp);
            hreferencefinetmp.push_back(hfinetmp);
            hnamereffinetmp+="sum2";
            TH1D * hfinetmpsum2=(TH1D*) hfinetmp->Clone(hnamereffinetmp);
            hreferencefinesum2tmp.push_back(hfinetmpsum2);

            hrefLOsubprocesseshistos.clear();
            hrefsubprocesseshistos.clear();
            hrefLOsubprocesseshistossum2.clear();
            hrefsubprocesseshistossum2.clear();

            hrefLOsubprocesseshistostmp.clear();
            hrefsubprocesseshistostmp.clear();
            hrefLOsubprocesseshistossum2tmp.clear();
            hrefsubprocesseshistossum2tmp.clear();

            vhref vtmphreflo;
            vhref vtmphreflosum2;
            vhref vtmphref;
            vhref vtmphrefsum2;

            vhref vtmphreflotmp;
            vhref vtmphreflosum2tmp;
            vhref vtmphreftmp;
            vhref vtmphrefsum2tmp;

            for ( int isub=0; isub<nsub; isub++) {
                //cout<<" MyGrid::Initialize: TEST nsub:"<<nsub<<endl; //should be 121
                //exit(0);
                
                
                TString hname="hrefLOsub";
                hname+=isub;
                //cout<<" hname= "<<hname<<endl;
                TH1D *htestlo=this->GetHisto(igrid,hname.Data());
                TString htit=TString(hname)+TString(" LO ")+TString(GetGridVersionName().c_str());//TString(this->GetSubProcessName(isub)); //<--**
                htestlo->SetTitle(htit);
                vtmphreflo.push_back(htestlo);

                TString hnametmp="hreftmpLOsub";
                hnametmp+=isub;
                //cout<<" hname= "<<hname<<endl;
                TH1D *htestlotmp=this->GetHisto(igrid,hnametmp.Data());
                TString htittmp=TString(hnametmp)+TString(" LO ")+TString(GetGridVersionName().c_str());
                htestlotmp->SetTitle(htittmp);
                vtmphreflotmp.push_back(htestlotmp);


                hname+="sum2";
                TH1D* hsum2=(TH1D*)htestlo->Clone(hname);
                vtmphreflosum2.push_back(hsum2);

                hnametmp+="sum2";
                TH1D* hsum2tmp=(TH1D*)htestlotmp->Clone(hnametmp);
                vtmphreflosum2tmp.push_back(hsum2tmp);



                TString hname2="hrefsub";
                hname2+=isub;
                //cout<<" hname2= "<<hname2<<endl;
                TH1D *hnlo=this->GetHisto(igrid,hname2.Data());
                TString htit2=TString(hname2)+TString(" ")+TString(GetGridVersionName().c_str());//TString(this->GetSubProcessName(isub)); //<--**
                hnlo->SetTitle(htit2);
                vtmphref.push_back(hnlo);

                TString hname2tmp="hreftmpsub";
                hname2tmp+=isub;
                //cout<<" hname2= "<<hname2<<endl;
                TH1D *hnlotmp=this->GetHisto(igrid,hname2tmp.Data());
                TString htit2tmp=TString(hname2tmp)+TString(" ")+TString(GetGridVersionName().c_str());
                hnlotmp->SetTitle(htit2tmp);
                vtmphreftmp.push_back(hnlotmp);


                hname2+="sum2";
                TH1D* hnlosum2=(TH1D*)hnlo->Clone(hname2);
                vtmphrefsum2.push_back(hnlosum2);

                hname2tmp+="sum2tmp";
                TH1D* hnlosum2tmp=(TH1D*)hnlotmp->Clone(hname2tmp);
                vtmphrefsum2tmp.push_back(hnlosum2tmp);

            }
            hrefLOsubprocesseshistos.push_back(vtmphreflo);
            hrefLOsubprocesseshistossum2.push_back(vtmphreflosum2);

            hrefsubprocesseshistos.push_back(vtmphref);
            hrefsubprocesseshistossum2.push_back(vtmphrefsum2);

            hrefLOsubprocesseshistostmp.push_back(vtmphreflotmp);
            hrefLOsubprocesseshistossum2tmp.push_back(vtmphreflosum2tmp);

            hrefsubprocesseshistostmp.push_back(vtmphreftmp);
            hrefsubprocesseshistossum2tmp.push_back(vtmphrefsum2tmp);

        }
    }


    //delete mypdf;
    //basic_pdf *mypdfTEST = new basic_pdf(); // TEST

    //mypdf = dynamic_cast<basic_pdf*>( appl::appl_pdf::getpdf(subprocesssteername) ); //TEST
    //mypdf = dynamic_cast<basic_pdf*>( appl::appl_pdf("basic") ); //TEST
    //mypdf->BasicInit(); //TEST



    mypdf = dynamic_cast<generic_pdf*>( appl::appl_pdf::getpdf(subprocesssteername) ); //<--**
    if(!mypdf)
        cout<<" MyGrid::Initialize: Warning: mypdf not found"<<endl;

    //mypdf->PrintSubprocess(); //<--** was used to check that subprocesses had been received


    nsub=mypdf->GetSubProcessNumber(); //<--**
    //nsub=mygrid.back()->subProcesses(0);//<--**appl_grid method. What was does parameter mean? The 0, 1, ...
    if(nsub<1)
        cout<<" MyGrid::Initialize: WARNING: nsub: "<<nsub<<endl;

    //cout<<" finished MyGrid ini "<<endl;


    return;

}


void MyGrid::PrintRefHistos(int igrid) {
    if(!hreference[igrid]) cout<<" MyGrid::PrintRefHistos reference histogram igrid= "<<igrid<<" not found "<<endl;
    hreference[igrid]->Print("all");
}
void MyGrid::PrintSubprocessRefHistos(int igrid) {

    cout<<" MyGrid::PrintSubprocessRefHistos reference histos "
        <<hrefsubprocesseshistos[igrid].size()<<endl;
    for (int i=0; i<nsub; i++) {
        if (!hrefsubprocesseshistos[igrid][i]) cout<<" reference histogram i= "<<i<<" not found "<<endl;
        else {
            cout<<" MyGrid:PrintSubProcessRefHistos: Histo= "<<hrefsubprocesseshistos[igrid][i]->GetTitle()<<endl;
            hrefsubprocesseshistos[igrid][i]->Print("all");
        }
    }
    /*
      cout<<" MyGrid::PrintSubprocessRefHistos LO reference histos "
          <<hrefLOsubprocesseshistos[igrid].size()<<endl;
      for (int i=0; i<hrefLOsubprocesseshistos[igrid].size(); i++){
       if (!hrefLOsubprocesseshistos[igrid][i]) cout<<" LO reference histogram i= "<<i<<" not found "<<endl;
       else {
        cout<<"MyGrid::PrintSubProcessRefHistos: LO Histo= "<<hrefLOsubprocesseshistos[igrid][i]->GetTitle()<<endl;
        hrefLOsubprocesseshistos[igrid][i]->Print("all");
       }
      }
    */
    return;
}

TH1D *MyGrid::TH1NormStatError(TH1D *hsum, TH1D *hsum2, double norm) {
    TH1D* h1new=(TH1D*) hsum->Clone(hsum->GetName());

    Double_t x, y, y2, ey;
    for (Int_t i=0; i<=hsum->GetNbinsX(); i++) {
        // sherpa formula
        //y= hsum ->GetBinContent(i)/norm;
        //y2=hsum2->GetBinContent(i)/norm;
        //ey=sqrt((y2-y*y)/(norm-1));

        y= hsum ->GetBinContent(i);
        y2=hsum2->GetBinContent(i);
        ey=sqrt((y2/norm-(y*y)/(norm*norm))/(norm-1)); //<--**think about


        y/=norm;
        ey/=norm;

        h1new->SetBinContent(i,y);
        h1new->SetBinError(i,ey);

        //if (debug) cout <<"bincenter= "<<hsum->GetBinCenter(i)<<" Binw = " << x << " y= " << y << endl;
    }
    return h1new;
}

void MyGrid::NormaliseInternalRefHistos(int igrid) {

    double yfac=mydata[igrid]->GetUnitfbFactor();
    double xfac=mydata[igrid]->GetUnitGeVFactor();
    cout<<" MyGrid::NormaliseInternalRefHistos: normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


    if(!hreference[igrid]) {
        cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hreference' for igird: "
            <<igrid<<" not found!"<<endl;
        exit(0); //TEST
    }
    else this->Normalise(hreference[igrid],yfac,xfac,true);

    if(!hreferencefine[igrid]) {
        cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hreferencefine' for igird: "
            <<igrid<<" not found!"<<endl;
        exit(0); //TEST
    }
    else this->Normalise(hreferencefine[igrid],yfac,xfac,true);


    for (int iproc=0; iproc<nsub; iproc++) {
        if (!hrefsubprocesseshistos[igrid][iproc]) {
            cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hrefsubprocesseshistos' for igrid: "
                <<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            exit(0); //TEST
        }
        else this->Normalise(hrefsubprocesseshistos[igrid][iproc],yfac,xfac,true);

        //cout<<" MyGrid::write_grid LO igrid= "<<igrid<<" iproc= "<<iproc<<endl;
        if (!hrefLOsubprocesseshistos[igrid][iproc]) {
            cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hrefLOsubprocesseshistos' for igrid: "
                <<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            exit(0); //TEST
        }
        else this->Normalise(hrefLOsubprocesseshistos[igrid][iproc],yfac,xfac,true);
    }
}

void MyGrid::ScaleInternalRefHistos(int igrid) {

        if(!hreference[igrid])
            cout<<" MyGrid::ScaleInternalRefHistos: WARNING: Reference histogram 'hreference' for igird: "<<igrid<<" not found!"<<endl;
        else hreference[igrid]->Scale(1.0/events[igrid]); //href/nhref

        if(!hreferencefine[igrid])
            cout<<" MyGrid::ScaleInternalRefHistos: WARNING: Reference histogram 'hreferencefine' for igird: "<<igrid<<" not found!"<<endl;
        else hreferencefine[igrid]->Scale(1.0/events[igrid]); //hrefFine/nhrefFine


        for (int iproc=0; iproc<nsub; iproc++) {
            if (!hrefsubprocesseshistos[igrid][iproc])
                cout<<" MyGrid::ScaleInternalRefHistos: WARNING: Reference histogram 'hrefsubprocesseshistos' for igrid: "<<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            else
                hrefsubprocesseshistos[igrid][iproc]->Scale(1.0/events[igrid]); //hrefsubprocesseshistos/nhrefsubprocesseshistos

            if (!hrefLOsubprocesseshistos[igrid][iproc])
                cout<<" MyGrid::ScaleInternalRefHistos: WARNING: Reference histogram 'hrefLOsubprocesseshistos' for igrid: "<<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            else
                hrefLOsubprocesseshistos[igrid][iproc]->Scale(1.0/events[igrid]); //hrefLOsubprocesseshistos/nhrefLOsubprocesseshistos
        }
}


void MyGrid::NormRefHistos(int igrid, double norm) {

    //cout<<" MyGrid::NormRefHistos igrid= "<<igrid<<endl;
    double scale=uncorrevents[igrid]*norm;

    TH1D* href=TH1NormStatError(hreference[igrid],hreferencesum2[igrid],scale);
    if (mydata[igrid]->DivideByBinWidth()) {
        this->DivideByBinWidth(href);
    }
    hreference.push_back(href);
    //cout<<" number hrefsubprocess= "<<hrefsubprocesseshistos[igrid].size()<<endl;

    TH1D* hreffine=TH1NormStatError(hreferencefine[igrid],hreferencefinesum2[igrid],scale);
    if (mydata[igrid]->DivideByBinWidth()) {
        this->DivideByBinWidth(hreffine);
    }
    hreferencefine.push_back(hreffine);


    for (int i=0; i<nsub; i++) {
        if (!hrefsubprocesseshistos[igrid][i]) cout<<" MyGrid::NormRefHistos reference histogram i= "<<i<<" not found "<<endl;
        else {
            if (!hrefsubprocesseshistossum2[igrid][i]) cout<<" MyGrid::NormRefHistos reference histogram sum2 i= "<<i<<" not found "<<endl;
            //cout<<" MyGrid::NormRefHistos: Histo= i= "<<i<<" " <<hrefsubprocesseshistos[igrid][i]->GetTitle()<<endl;
            href=TH1NormStatError(hrefsubprocesseshistos[igrid][i],hrefsubprocesseshistossum2[igrid][i],scale);
            if (mydata[igrid]->DivideByBinWidth()) {
                this->DivideByBinWidth(href);
            }
            hrefsubprocesseshistos[igrid].push_back(href);
        }

        //cout<<" MyGrid::NormRefHistos LO igrid= "<<igrid<<" i= "<<i<<endl;
        if (!hrefLOsubprocesseshistos[igrid][i]) cout<<" MyGrid::NormRefHistos LO reference histogram i= "<<i<<" not found "<<endl;
        else {
            //cout<<" MyGrid::NormRefHistos: LO Histo= "<<hrefLOsubprocesseshistos[igrid][i]->GetTitle()<<endl;
            href=TH1NormStatError(hrefLOsubprocesseshistos[igrid][i],hrefLOsubprocesseshistossum2[igrid][i],scale);
            if (mydata[igrid]->DivideByBinWidth()) {
                this->DivideByBinWidth(href);
            }
            hrefLOsubprocesseshistos[igrid].push_back(href);
            //    hrefLOsubprocesseshistos[igrid][i]->Print("all");
        }
    }
    return;
}

//
void MyGrid::ReadSteering(string fname) {

    steername=fname;
    if (debug) cout<<" MyGrid::ReadSteering steering= "<<steername<<endl;

    ifstream infile(steername.c_str(), ios::in);
    if(!infile) { // Check open
        cerr << " MyGrid::ReadSteering: Can't open " << steername <<"\n";
        infile.close();
        exit (1);
    } else {
        cout <<" MyGrid::ReadSteering read data file: " << steername << endl;
    }

// cout<<" filed opened "<<fname<<endl;

    int iline=0;
    int nsyst=1;

    char line[256];
    char lineFirstWord[256]; //same length as 'line' string buffer
    char text[100];
    char name[100];
    float mys;
    string OptionsFileName = "MyGridOptions.txt"; //<--could be made readable from steering file
    OptionHandler *mygridoptions= new OptionHandler(OptionsFileName);
    //mygridoptions->generateResultFile("output/test.txt"); //method currently not implimented, use debug=true instead

    while (infile.good()) {
        if (debug) cout << " good: " << infile.good() << " eof: " << infile.eof() << endl;
        //if (!infile.good()) break;
        infile.getline(line,sizeof(line),'\n');
        if (debug) cout<< "line= "<< line << "\n";

        sscanf(line," %s",lineFirstWord);
        std::string cpp_lineFirstWord(lineFirstWord);

        if(line[0] != '%') {
            /*if(line[0] != 'g' && line[0] != 'n' && line[0] != 's' && line[0] != 'p' &&
               line[0] != 'q'
               ){*/
            if (mygridoptions->isKnownOption(cpp_lineFirstWord)==false) {
                // do something
                //} else if (strstr(line,"nsub")!=0) {
                //char text[100];
                //sscanf(line," %s %d ",text, &nsub);
                //if (debug) printf("**MyGrid Read: Nsub= %d   \n",nsub);
            } else if (strstr(line,"nprocessnumber")!=0) {
                sscanf(line," %s %d ",text, &processnumber);
                if (debug) printf("**MyGrid::Read processnumber= %d   \n",processnumber);
            } else if (strstr(line,"subprocesssteername")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                subprocesssteername=string(name);
                if (debug) printf("**MyGrid::Read subprocesssteername= %s \n",subprocesssteername.c_str());
                //cout<<"     subprocesssteername= "<<subprocesssteername<<endl;
                //   } else if (strstr(line,"pdffunction")!=0) {
                //sscanf(line," %s %[^\n] ",text, name);
                //pdf_function=string(name);
                //if (debug) printf("**MyGrid Read: pdf_function= %s \n",pdf_function.c_str());
                //cout<<"   pdffunction= "<<pdf_function<<endl;

            } else if (strstr(line,"ntupdiroutput")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                ntupdiroutput=name;
            } else if (strstr(line,"ntupdirinput")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                ntupdirinput=name;
            } else if (strstr(line,"ntupname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                cout<<" MyGrid::Read "<<text<<" "<<name<<endl;
                ntupname=name;
            } else if (strstr(line,"nXbins")!=0) {
                sscanf(line," %s %f ",text, &mys);
                cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                nXbins=mys;
            } else if (strstr(line,"xLow")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                xLow=mys;
            } else if (strstr(line,"xUp")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                xUp=mys;
            } else if (strstr(line,"xorder")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<mys<<endl;
                xorder=mys;
            } else if (strstr(line,"nQ2bins")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                nQ2bins=mys;
            } else if (strstr(line,"q2Low")!=0) {
                sscanf(line," %s %f ",text, &mys);
                cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                q2Low=mys;
            } else if (strstr(line,"q2Up")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                q2Up=mys;
            } else if (strstr(line,"qorder")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                qorder=mys;
            } else if (strstr(line,"griddir")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                gridnamedir=name;
            } else if (strstr(line,"gridname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<name<<endl;
                string myname=name;
                gridname.push_back(myname);
            };
        };
        memset(text, '\0', 100); //reset option input buffers for next run
        memset(name, '\0', 100);
    };
    return;
};

void MyGrid::Print() {

    cout<<" MyGrid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    cout<<" MyGrid::steering: "<<steername<<endl;
    cout<<" ntupdirinput= "<<ntupdirinput<<" ntupname=  "<<ntupname<<endl;
    cout<<" ntupdiroutput= "<<ntupdiroutput<<endl;
    cout<<" gridnamedir=  "<<gridnamedir<<endl;
    cout<<" subprocesssteername: "<<subprocesssteername<<endl;
    //cout<<" pdf_function=  "<<pdf_function<<endl;

    cout<<"\n Grid architecture: "<<endl;
    cout<<"\t  nXbins= "<<nXbins<<"  xLow= "<<xLow<<" xUp= "<<xUp<<" xorder= "<<xorder<<endl;
    cout<<"\t nQ2bins= "<<nQ2bins<<" q2Low= "<<q2Low<<" q2Up= "<<q2Up<<" qorder= "<<qorder<<endl;

    cout<<"\t lowest_order= "<<lowest_order<<" nloops= "<<nloops<<endl;

    cout<<"\n Selected Events: "<<endl;
    std::cout<<"all uncorrellated events= "<<alluncorrevents
             <<" all events= "<<allevents<<std::endl;

    cout<<"\n Number of grids N= "<<gridname.size()<<endl;
    for (int  i = 0; i <   gridname.size(); i++)
    {
        cout<<gridname[i]
            <<" uncorrelated events= "<<uncorrevents[i]
            <<" events= "<<events[i]
            <<" style: "<<this->GetStyleMarker(i)
            <<" color: "<<this->GetColorMarker(i)
            <<endl;
        //mydata[i]->Print();
    }


    cout<<"MyGrid <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;

}

bool MyGrid::file_exists(const string& s) {
    if ( FILE* testfile=fopen(s.c_str(),"r") ) {
        fclose(testfile);
        return true;
    }
    else return false;
}

void MyGrid::book_grid(int igrid)  // inital grid booking
{

    time_t _t;
    time(&_t);

    std::cout<<" ***********************************************"<<std::endl;
    std::cout<<" MyGrid::book_grid("<<igrid<<")"<<std::endl;
    std::cout<<" booking the grids " << ctime(&_t) << std::endl;
    //std::cout<<" Grid Name= "<<this->GetGridName(igrid)<<endl;

    double apramval=5.;
    appl::igrid::transformvar(apramval);
    //appl::igrid transformvar(apramval);


    // binning information for the grid constructor
    // forsee to put this in steering at some point
    //
    //cout<<" MyGrid book grid "<<glabel<<endl;

    bool create_new = false;
    // if the file does not exist, create a new grid...
    string filename=this->GetGridFullFileName(igrid);//+this->GetGridVersionName();


    string _filename = filename;
    string newFileName = GetGridVersionName()+".root";
    _filename.replace( _filename.find(".root"), 5, newFileName.c_str());
    filename=_filename;

    std::cout<<" MyGrid::book_grid: TEST: checking if this file exists: "<<filename<<std::endl;


    std::cout << "book_grid() " << igrid << ", grid name " << filename << std::endl;


    if ( !file_exists(filename) ) {
        cout<<"MyGrid::book_grid: file named '"<<filename<<"' does not exist. create_new=true."<<std::endl;
        create_new = true;
    }
    // or if it does exists but root file is a zombie...
    if ( !create_new ) {
        cout<<" MyGrid::book_grid: create test file "<<endl;
        TFile testFile( (TString(filename)));
        if ( testFile.IsZombie() ) {
            cout<<"MyGrid::book_grid: test file is a zombie! Create_new=true"<<std::endl;
            create_new = true;
        }
        testFile.Close();
    }


    std::cout<<" MyGrid::book_grid: TEST: final create_new is: "<<create_new<<" for file: "<<filename<<std::endl;


    //mydata[igrid]->Print();
    int nObsBins=mydata[igrid]->GetNBins();
    double *obsBins=mydata[igrid]->GetBins();

    //  for (int i=0; i<nObsBins+1; i++)
    // cout<<i<<" obsBins= "<<obsBins[i]<<endl;

    if ( create_new ) {
        cout<<" MyGrid::book Creating NEW grid... " << filename << "\tigrid " << igrid << endl;
        cout<<" nObsBins= "<<nObsBins<<endl;
        cout<<" nQ2bins= "<<nQ2bins<<" q2Low= "<<q2Low<<" q2Up= "<<q2Up<<" qorder= "<<qorder<<endl;
        cout<<" nXbins= "<<nXbins<<" xLow= "<<xLow<<" xUp= "<<xUp<<" xorder= "<<xorder<<endl;
        cout<<" pdf_function= "<<pdf_function<<endl;
        cout<<" lowest_order= "<<lowest_order<<" nloops= "<< nloops<<endl;


        appl::grid *tmpgrid = new appl::grid( nObsBins, obsBins,      // obs bins
                                              nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
                                              nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
                                              pdf_function, lowest_order, nloops );

        mygrid.push_back(tmpgrid);

        //cout<<"Mygrid::book_grid: printing appl_grid..."<<endl;
        //tmpgrid->print();

        //double sqrts=7000.; //put here MyData
        //mygrid[igrid]->setCMScale(sqrts ); // done later

        //int nsub=mypdf->GetSubProcessNumber(); //<--**
        ///int nsub=tmpgrid->subProcesses(0); //<--**appl_grid method. What was does parameter mean? The 0, 1, ...
        //int nsub=mygrid.back()->subProcesses(0);//<--**appl_grid method. What was does parameter mean? The 0, 1, ...

        //cout << " MyGrid::book_grid nsub = " << nsub << endl;
        cout << " MyGrid::book_grid reference histo name = " << mygrid[igrid]->getReference()->GetName() << endl;

    } else {
        std::cout <<" MyGrid::book_grid Using existing grid file "<<filename<<std::endl;
        appl::grid *tmpgrid = new appl::grid(filename); //optimise grid x,Q2 bins
        mygrid.push_back(tmpgrid);
        mygrid[igrid]->getReference()->Reset();
        mygrid[igrid]->optimise(nQ2bins, nXbins);
        std::cout<<*(mygrid[igrid])<<std::endl;
    }


    //appl_pdf::printmap();
    //appl_pdf* mypdf = appl_pdf::getpdf("topsubprocesses.dat"); //<--**


    //std::cout << "grid name " << mypdf << std::endl; //<--**
    //std::cout << "grid name " << mypdf->name() << std::endl; //<--**

    nsub=mygrid[igrid]->subProcesses();


    cout<<" MyGrid::book_grid subprocesses= "<<mygrid[igrid]->subProcesses()<<endl;

    std::cout<<" ***********************************************"<<std::endl;

}



void  MyGrid::fill(MyEvent *myevent )
{
    //  if (!isBooked)
    //  {
    //book_grid();
    //    cout<<" Need to book grid first "<<endl;
    //    return;
    //  }

// counter of number of events
    //newevent=false;
    //if (debug) cout<<"MyGrid::fill eventid= "<<eventid<<" GetEvent= "<<myevent->GetEventId()<<endl;

    //if (eventid!=myevent->GetEventId()){
    //eventid=
    //newevent=true;
    //}
    if (debug) cout<<" MyGrid::fill newevent= "<<newevent<<endl;
    this->SetEventId(myevent->GetEventId());
    if (debug&&this->NewEvent()) cout<<" MyGrid::fill new event "<<endl;

    if (newevent) alluncorrevents++;
    //else cout<<" not new event "<<endl;
    allevents++;

    double x1=myevent->GetX1();
    double x2=myevent->GetX2();
    double q2=myevent->GetQ2();
    double xsec=myevent->GetXSection();
    double mewgt=myevent->GetWeight();

    if (debug)
        cout<<" MyGrid::fill x1= "<<x1<<" x2= "<<x2<<" q2= "<<q2
            <<" mewgt= "<<mewgt<<" xsec= "<<xsec<<endl;

    int Ngrids=this->GetNGrid();
    if (debug) cout<<" MyGrid::fill Ngrids= "<<Ngrids<<endl;
    int id1=myevent->GetID(0);
    int id2=myevent->GetID(1);
    if (debug) cout<<" MyGrid::fill incoming partons id1= "<<id1<<" id2= "<<id2<<endl;

//
// Need to prepare weight vector in subprocesses
//

    int nproc=this->GetProcessNumber();
    if (debug) cout<<" MyGrid::fill nproc= "<<nproc<<endl;
    

    if (nproc==1 || nproc==2) { // process 1 and 2 is W
        int icharge=myevent->GetLeptonCharge();
        if (icharge== 1) nproc=1;
        if (icharge==-1) nproc=2;
        if (debug) cout<<" MyGrid::fill nproc= "<<nproc<<" icharge= "<<icharge<<endl;
    }

    //if (!mypdf) cout<<" MyGrid::fill mypdf not initialized "<<endl; //<--**
    //mypdf->SetCurrentProcess(nproc); //<--** Currently does nothing? nothing seen last code checked: 19-Jun-2013

    int myorder=-1;
    if (nproc==1||nproc==2) { // Wminus and Wplus
        myorder=myevent->GetOrder();
        if (myorder==1) myorder=0;
        if (myorder==2) myorder=1;
        if (debug) cout<<" MyGrid::fill myorder= "<<myorder<<endl;
    }

    if (nproc==3 || nproc==4) { // 3: inclusive jet 4: top
        myorder=myevent->GetOrder();
        if (myorder==2) myorder=0;
        if (myorder==3) myorder=1;
        if (debug) cout<<" MyGrid::fill myorder= "<<myorder<<endl;
    }
        

    if (myorder==-1) {
        cout<<" MyGrid::fill do not know what to do myorder= -1"<<endl;
        return;
    }


    //cout<<" MyGrid::fill: Dynamically making generic_pdf for function 'decideSubProcess'"<<endl; //<--** test output
    //generic_pdf* mypdf  = dynamic_cast<generic_pdf*>( appl::appl_pdf::getpdf(pdf_function) ); //<--**

    //std::cout<<"MyGrid::fill: decideSubProcess using: id1: "<<id1<<", id2: "<<id2<<std::endl;
    int iproc=mypdf->decideSubProcess(id1,id2); //<--**



    if (iproc==-1)  {
        cout<<" MyGrid::fill do not know what to do "<<endl;
        //myevent->Print();
        return;
    }

    if (debug) cout<<" MyGrid::fill iproc= "<<iproc<<endl;

    this->ResetWeight();
    this->SetWeight(iproc,mewgt);
    if (debug) this->PrintWeight();

    double *weight=this->GetWeight();

    for(int igrid = 0; igrid < Ngrids; igrid++) {

        //
        // get cuts
        //
        if (debug) cout<<" MyGrid::fill filling igrid= "<<igrid<<endl;
        bool mycut=this->eventcuts(myevent,igrid);
        if (debug) if (mycut==false) cout<<" MyGrid::fill Event rejected !"<<endl;
        if (mycut==false) continue;
        //
        if (debug) cout<<" MyGrid::fill event accepted "<<endl;
        if (newevent) uncorrevents[igrid]++;
        events[igrid]++;

        if (debug) {
            cout<<" MyGrid::fill uncorrelated events= "<<uncorrevents[igrid]<<endl;
            cout<<" MyGrid::fill events= "<<events[igrid]<<" iproc= "<<iproc<<endl;
            this->PrintWeight();
            cout<<" MyGrid::fill x1= "<<x1<<" x2= "<<x2<<" q2= "<<q2
                <<" mewgt= "<<mewgt<<" xsec= "<<xsec<<" myorder= "<<myorder<<endl;
        }

        //double obsmax=this->GetObsMax(igrid);
        //cout<<" obsmax= "<<obsmax<<endl;
        //if (obs>obsmax) {
        //cout<<" cut by obsmax bins= "<<obsmax<<" obs= "<<obs<<endl;
        //return; //speeds up code
        //}

        incljets=mydata[igrid]->isObsInclJets();
        if (debug) {
            if (incljets) cout<<" MyGrid::fill inclusive jet observable "<<endl;
            else          cout<<" MyGrid::fill inclusive jet not set "<<endl;
        }

        double ptmin =mydata[igrid]->GetJetPtCut();
        double rapmin=mydata[igrid]->GetJetMinY();
        double rapmax=mydata[igrid]->GetJetMaxY();
        if (debug) cout<<" MyGrid::fill ptmin= "<<ptmin<<" rapmin= "<<rapmin<<" rapmax= "<<rapmax<<endl;

        std::vector<fastjet::PseudoJet> myseljets=myevent->GetSelectedJets(ptmin,rapmin,rapmax);
        int nseljets=1;
        // if Observable is call inclusive jets loop over jets, otherwise fill grid only once
        double obs=0;
        if (incljets) {
            nseljets=myevent->GetNumberSelJets();
            if (debug) {
                cout<<" MyGrid::fill number of selected jets= "<<nseljets<<endl;
                myevent->PrintSelJets2();
            }
        } else {
            if (debug) {
                cout<<" MyGrid::fill number of selected jets= "<<myevent->GetNumberSelJets()<<endl;
                myevent->PrintSelJets2();
            }
        }

        for (int ijets=0; ijets<nseljets; ijets++) {
            if (incljets) obs=myseljets[ijets].pt();
            else if (nproc!=4)   obs=myevent->GetLeadingJet();
            else { // need to think about better implementation
                if (nproc==4) {
                    if (myevent->GetNumberSelJets()<2)
                        cout<<" MyGrid::fill something is wrong nsel= "
                            <<myevent->GetNumberSelJets()<<endl;
                    else {
                        if (debug) cout<<" MyGrid::fill number of selected jets= "<<nseljets<<endl;
                        obs=myevent->GetInvariantMass12();
                    }
                    if (debug) cout<<" MyGrid::fill obs= "<<obs<<endl;
                }
            }
            if (debug) cout<<" MyGrid::fill "<<" ijets= "<<ijets<<" obs= "<<obs<<endl;
            if (mygrid[igrid]->isOptimised()) {
                if (debug) cout<<" MyGrid::fill filling grid and reference "<<endl;
                if (!mygrid[igrid]) cout<<"MyGrid::fill grid not found igrid= "<<igrid<<endl;
                mygrid[igrid]->fill(x1,x2, q2, obs, weight, myorder );
            } else {
                if (debug) cout<<" MyGrid::fill filling phasespace "<<endl;
                if (!mygrid[igrid]) cout<<" MyGrid::fill grid not found igrid= "<<igrid<<endl;
                mygrid[igrid]->fill_phasespace(x1, x2, q2, obs, weight, myorder );
            }
            mygrid[igrid]->getReference()->Fill( obs,xsec);

            if (bookrefsubprocess) {

                if (myorder==0) {
                    //hrefLOsubprocesseshistostmp[igrid][iproc]->Fill(obs,weight[iproc]);
                    //hrefLOsubprocesseshistossum2tmp[igrid][iproc]->Fill(obs,weight[iproc]*weight[iproc]);
                    hrefLOsubprocesseshistostmp[igrid][iproc]->Fill(obs,xsec*xsec);
                    hrefLOsubprocesseshistossum2tmp[igrid][iproc]->Fill(obs,xsec*xsec);
                }

                //cout<<" MyGrid::fill subProcs for igrid: "<<igrid<<", obs: "<<obs<<", weight[iproc:"<<iproc<<"]: "<<weight[iproc]<<endl;

                //hrefsubprocesseshistostmp[igrid][iproc]->Fill(obs,weight[iproc]);
                //hrefsubprocesseshistossum2tmp[igrid][iproc]->Fill(obs,weight[iproc]*weight[iproc]);
                hrefsubprocesseshistostmp[igrid][iproc]->Fill(obs,xsec*xsec);
                hrefsubprocesseshistossum2tmp[igrid][iproc]->Fill(obs,xsec*xsec);

                //hreferencetmp[igrid]->Fill(obs,weight[iproc]);
                //hreferencesum2tmp[igrid]->Fill(obs,weight[iproc]*weight[iproc]);
                hreferencetmp[igrid]->Fill(obs,xsec*xsec);
                hreferencesum2tmp[igrid]->Fill(obs,xsec*xsec);

                //hreferencefinetmp[igrid]->Fill(obs,weight[iproc]);
                //hreferencefinesum2tmp[igrid]->Fill(obs,weight[iproc]*weight[iproc]);
                hreferencefinetmp[igrid]->Fill(obs,xsec*xsec);
                hreferencefinesum2tmp[igrid]->Fill(obs,xsec*xsec);
            }
        }
        if (this->NewEvent()) {
            //std::cout<<" MyGrid::fill: TEST: NEW EVENT ADD igrid: "<<igrid<<", iproc: "<<iproc<<std::endl;
            //std::cout<<" MyGrid::fill: TEST: Adding hrefLOsubprocesseshistostmp[igrid:"<<igrid<<"][iproc:"<<iproc<<"]: "<<hrefLOsubprocesseshistostmp[igrid][iproc]<<" TO hrefsubprocesseshistos[igrid:"<<igrid<<"][iproc:"<<iproc<<"]: "<<hrefsubprocesseshistos[igrid][iproc]<<std::endl;

            hreference[igrid]->Add(hreferencetmp[igrid]);
            hreferencetmp[igrid]->Reset();
            //hreference[igrid]->Print("all");

            hreferencefine[igrid]->Add(hreferencefinetmp[igrid]);
            hreferencefinetmp[igrid]->Reset();
            //hreference[igrid]->Print("all");

            hreferencesum2[igrid]->Add(hreferencesum2tmp[igrid]);
            hreferencesum2tmp[igrid]->Reset();

            hreferencefinesum2[igrid]->Add(hreferencefinesum2tmp[igrid]);
            hreferencefinesum2tmp[igrid]->Reset();

            hrefLOsubprocesseshistos[igrid][iproc]->Add(hrefLOsubprocesseshistostmp[igrid][iproc]);
            hrefLOsubprocesseshistostmp[igrid][iproc]->Reset();

            hrefLOsubprocesseshistossum2[igrid][iproc]->Add(hrefLOsubprocesseshistossum2tmp[igrid][iproc]);
            hrefLOsubprocesseshistossum2tmp[igrid][iproc]->Reset();

            hrefsubprocesseshistos[igrid][iproc]->Add(hrefsubprocesseshistostmp[igrid][iproc]);
            hrefsubprocesseshistostmp[igrid][iproc]->Reset();

            hrefsubprocesseshistossum2[igrid][iproc]->Add(hrefsubprocesseshistossum2tmp[igrid][iproc]);
            hrefsubprocesseshistossum2tmp[igrid][iproc]->Reset();
        } // new event
    } //loop over grid
    return;
}

void MyGrid::Normalise(TH1D* h1, double yscale, double xscale=1., bool normtot=false)
{

    Double_t x, y, ey;
    Double_t sigtot=0.;

    for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
        y=h1->GetBinContent(i)*yscale;
        x=h1->GetBinWidth(i);
        //x=1.;
        sigtot+=y*x;
    }

//cout<<" MyGrid::Normalise sigtot= "<<sigtot<<endl;

    for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
        x =h1->GetBinWidth(i);
        y =h1->GetBinContent(i)*yscale*x;
        ey=h1->GetBinError(i)  *yscale*x;
        x =h1->GetBinWidth(i)  *xscale;
        if (x!=0) h1->SetBinContent(i,y/x);
        else      h1->SetBinContent(i,0.);

        if (x!=0) h1->SetBinError(i,ey/x);
        else      h1->SetBinError(i,0.);
    }
    if (normtot) {
        if (sigtot!=0.)
            h1->Scale(1. / sigtot);
    }
    return;
}

//
// just normalise to bin width
//
void MyGrid::DivideByBinWidth(TH1D* h)
{
    if (debug)
        cout<<" MyGrid::DivideByBinWidth called "<<endl;

    for ( int ibin=0 ; ibin<=h->GetNbinsX() ; ibin++ )
    {
        double width = h->GetBinWidth(ibin);
        if (debug) cout<<" MyGrid::DivideByBinWidth width= "<<width
                           <<" context= "<<h->GetBinContent(ibin)
                           <<endl;
        h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
        h->SetBinError  ( ibin, h->GetBinError(ibin)/width );
    }
    return;
}

//extern "C" void write_grid_(double& xstotal)   // writes out grid after some events
void MyGrid::write_grid()   // writes out grid after some events
{


    int Ngrids=this->GetNGrid();

    std::cout<<" MyGrid::write_grid Write out the grid ... Ngrids " << Ngrids << std::endl;


    for(int igrid = 0; igrid < Ngrids; igrid++)
    {
        //mygrid[igrid]->getReference()->Print("all");
        cout<<"MyGrid::write_grid saving grid N=" << igrid+1 << "\tof total " << Ngrids << endl;
        //mygrid[igrid]->setNormalised( true );
        mygrid[igrid]->run() = events[igrid];

        cout<<" run stored "<<endl;

        mygrid[igrid]->trim();

        int trim_size = mygrid[igrid]->size();


        //        cout<<" trimmned "<<endl;
        //        mygrid[igrid]->untrim();
        //        int untrim_size = mygrid[igrid]->size();

        //        cout<<"Saved Space ratio: "<<(trim_size/untrim_size*1.)<<endl;

/*
        //cout<<" untrimmed "<<endl;
        double yfac=mydata[igrid]->GetUnitfbFactor();
        double xfac=mydata[igrid]->GetUnitGeVFactor();
        cout<<" MyGrid::write_grid normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;


        //NEW normalization
        if(!hreference[igrid])
            cout<<" MyGrid::write_grid: WARNING: Reference histogram 'hreference' for igird: "<<igrid<<" not found!"<<endl;
        hreference[igrid]->Scale(1.0/events[igrid]); //href/nhref
        //this->Normalise(hreference[igrid],yfac,xfac,true);

        if(!hreferencefine[igrid])
            cout<<" MyGrid::write_grid: WARNING: Reference histogram 'hreferencefine' for igird: "<<igrid<<" not found!"<<endl;
        hreferencefine[igrid]->Scale(1.0/events[igrid]); //hrefFine/nhrefFine
        //this->Normalise(hreferencefine[igrid],yfac,xfac,true);


        for (int iproc=0; iproc<nsub; iproc++) {
            if (!hrefsubprocesseshistos[igrid][iproc])
                cout<<" MyGrid::write_grid: WARNING: Reference histogram 'hrefsubprocesseshistos' for igrid: "<<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            else {
                hrefsubprocesseshistos[igrid][iproc]->Scale(1.0/events[igrid]); //hrefsubprocesseshistos/nhrefsubprocesseshistos
                //this->Normalise(hrefsubprocesseshistos[igrid][iproc],yfac,xfac,true);
            }

            if (!hrefLOsubprocesseshistos[igrid][iproc])
                cout<<" MyGrid::write_grid: WARNING: Reference histogram 'hrefLOsubprocesseshistos' for igrid: "<<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
            else {
                hrefLOsubprocesseshistos[igrid][iproc]->Scale(1.0/events[igrid]); //hrefLOsubprocesseshistos/nhrefLOsubprocesseshistos
                //this->Normalise(hrefLOsubprocesseshistos[igrid][iproc],yfac,xfac,true);
            }
        }
        //NEW normalization
*/



        string filename=this->GetGridFullFileName(igrid);//+this->GetGridVersionName();

        string _filename = filename;
        string newFileName = GetGridVersionName()+".root";
        _filename.replace( _filename.find(".root"), 5, newFileName.c_str());
        filename=_filename;

        cout<<" MyGrid::write_grid: TEST: Writing to filename= " << filename << "\tigrid: " << igrid << endl;
        if (mygrid[igrid])   mygrid[igrid]->Write(filename);
        else                 cout<<" MyGrid::write_grid() mygrid not found ! "<<endl;






        //        TFile *f = new TFile(filename.c_str(),"update");

        //std::cout<<"TEST: GetGridFullFileName: "<<GetGridFullFileName(igrid)<<", filename: "<<filename<<std::endl;

        if(bookrefsubprocess) {
            _filename = filename;
            newFileName = "-ghistos.root";
            //_filename.replace( _filename.find(".root"), 5, "-ghistos.root");
            _filename.replace( _filename.find(".root"), 5, newFileName.c_str());

            TFile *f = new TFile( _filename.c_str(),"recreate");
            f->Print();

            //gDirectory->cd();
            /*
            cout<<" MyGrid::write_grid to filename= "<<filename<<endl;
            if (!mygrid[igrid]) cout<<" MyGrid::write_grid() mygrid not found ! "<<endl;
            else mygrid[igrid]->Write(filename);
            cout<<" TEST mygrid[igrid] is finished writing!"<<endl; //TEST
            */

            if (!hreference[igrid]) {
                cout<<" MyGrid::write_grid() hreference["<<igrid<<"] not found ! "<<endl;
                exit(0); //TEST
            }
            else hreference[igrid]->Write();

            if (!hreferencefine[igrid]) {
                cout<<" MyGrid::write_grid() hreferencefine["<<igrid<<"] not found ! "<<endl;
                exit(0); //TEST
            }
            else hreferencefine[igrid]->Write();

            for (int iproc=0; iproc<nsub; iproc++) {
                if (!(hrefLOsubprocesseshistos[igrid][iproc])) {
                    cout<<" MyGrid::write_grid() hrefLOsubprocesseshistos["<<igrid<<"]["<<iproc<<"] not found ! "<<endl;
                    exit(0); //TEST
                }
                else hrefLOsubprocesseshistos[igrid][iproc]->Write();

                if (!(hrefsubprocesseshistos[igrid][iproc])) {
                    cout<<" MyGrid::write_grid()  hrefsubprocesseshistos["<<igrid<<"]["<<iproc<<"] not found ! "<<endl;
                    exit(0); //TEST
                }
                else hrefsubprocesseshistos[igrid][iproc]->Write();
            }

            cout << " MyGrid::Write() size(trimmed)=" << trim_size << endl;


            f->Close();
        }

        //int nsub = mygrid[igrid]->subProcesses();

        //delete mygrid[igrid]; //<--** delete break the GetReference(...) function?




    }

    time_t _t;
    time(&_t);

    std::cout<<" ***********************************************"<<std::endl;
    //std::cout<<" saved grids " << ctime(&_t) << std::endl;
    std::cout<<" ***************************************************"<<std::endl;
}


//
// ----------------------------------------------
//    analysis
// ----------------------------------------------
//

void MyGrid::getObservable(MyEvent *myevent)
{
    //
    /*
    // calculate observables
    for(int igrid = 0; igrid < Ngrids; igrid++)Observable[igrid] = 0.0; // initialize

    double p3[4] = {evt[3][2],evt[0][2],evt[1][2],evt[2][2]}; // (E,x,y,z)
    double p4[4] = {evt[3][3],evt[0][3],evt[1][3],evt[2][3]};


    double rapidity3 = 0.0;
    rapidity3 = (p3[0] + p3[3])/(p3[0] - p3[3]);
    (rapidity3 < 1e-13) ? rapidity3 = 100.0 : rapidity3 = 0.5*std::log(rapidity3);

    double rapidity4 = 0.0;
    rapidity4 = (p4[0] + p4[3])/(p4[0] - p4[3]);
    (rapidity4 < 1e-13) ? rapidity4 = 100.0 : rapidity4 = 0.5*std::log(rapidity4);

    double rapidity34 = 0.0;                      // rapidity of particle (3+4) in event record
    rapidity34  = (p3[0] + p4[0]) + (p3[3] + p4[3]);
    rapidity34 /= (p3[0] + p4[0]) - (p3[3] + p4[3]);

    (rapidity34 < 1e-13) ? rapidity34 = 100.0 : rapidity34 = 0.5*std::log(rapidity34);

    double pt3 = 0;
    pt3 = std::sqrt( p3[1]*p3[1] + p3[2]*p3[2] );

    double pt4 = 0;
    pt4 = std::sqrt( p4[1]*p4[1] + p4[2]*p4[2] );

    double pt34 = 0;
    pt34 = std::sqrt( std::pow(p3[1] + p4[1],2) + std::pow(p3[2] + p4[2],2) );

    double mass34 = 0;
    mass34 = std::sqrt( std::pow(p3[0] + p4[0],2) - std::pow(p3[1] + p4[1],2) - std::pow(p3[2] + p4[2],2) - std::pow(p3[3] + p4[3],2) );

    Observable[ 0 ] = rapidity3;
    Observable[ 1 ] = pt3;
    Observable[ 2 ] = rapidity4;
    Observable[ 3 ] = pt4;
    //  Observable[ 4 ] = rapidity34;
    //  Observable[ 5 ] = pt34;
    //  Observable[ 6 ] = mass34;

    */
}

//
// copy from MCFM assume that ckm comes from sherpa
//typedef struct {
//  double vsq[ __nf2__ ][ __nf2__ ], vsum[ __nf2__ ];
//} __ckm__;
//  double precision Vsq(-nf:nf,-nf:nf),Vsum(-nf:nf)
//      common/CKM/Vsq,Vsum
// copy from mcfm
// filling in Need/ckmfill.f
//
bool  MyGrid::eventcuts(MyEvent *myevent, int igrid) {
//
// analysis cuts defined in steering
//
    bool mydebug=false;

    bool accept=false;
    bool onejet=false;
    bool oklept=false;
    bool okneut=false;

    int idjetmax=-1;
    double ptold=-999.;
    int idlep=-1;
    int idneu=-1;
    MyData *mydata=this->GetMyData(igrid);
    if (!mydata) cout<<" MyGrid::eventcuts MyData object not found !"<<endl;

    int njets=0;
    for (int  i = 0; i < myevent->GetN(); i++) {
        int id = myevent->GetID(i);
        if (mydebug) std::cout<<" MyGrid::eventcuts";
        double pt = (myevent->Get(i))->pt();
        double rap = fabs((myevent->Get(i))->rap());
        if (mydebug)
            printf(" %d pt= %6.2f rap= %6.2f  id= %d \n",i,pt, rap, id);

//bool isjet=abs(id)<10||abs(id)==21;
        bool isjet = myevent->IsJet(id);

        if (mydebug)
            if( mydata->JetCut()) cout<<" MyGrid::eventcuts JetCut on ! "<<endl;
            else                  cout<<" MyGrid::eventcuts JetCut off !"<<endl;

        if (mydebug&&isjet) cout<<" MyGrid::eventcuts is a jet "<<i<<endl;
        if (mydata->JetCut()&&isjet) {
            if ((rap<mydata->GetJetMaxY()&&rap>mydata->GetJetMinY())
                    && pt>mydata->GetJetPtCut()) {
                onejet=true;
                if (mydebug) cout<<" MyGrid::eventcuts accept jet "<<i<<endl;
                njets++;
                if (pt>ptold) {
                    ptold=pt;
                    idjetmax=i;
                }
            };
        };

        //bool islepton=abs(id)==11||abs(id)==13;
        bool islepton=myevent->IsChargedLepton(id);
        if (mydebug&&islepton) cout<<" MyGrid::eventcuts is a lepton "<<i<<endl;
        if (mydata->LeptonCut()&&islepton) {
            idlep=i;
            if ((rap<mydata->GetLepMaxY()&&rap>mydata->GetLepMinY())
                    &&pt>mydata->GetLepPtCut()) oklept=true;
            if (mydebug)
                if (oklept) cout<<" MyGrid::eventcuts accept lepton cut "<<endl;
                else        cout<<" MyGrid::eventcuts lepton cut not passed "<<endl;
        };

        bool isneutrino=myevent->IsNeutrino(id);
        if (mydebug&&isneutrino) cout<<" MyGrid::eventcuts is a neutrino "<<i<<endl;
        if (mydata->NeutrinoCut()&&isneutrino) {
            idneu=i;
            if (pt>mydata->GetNeuPtCut()) okneut=true;
            if (mydebug) if (okneut) cout<<" MyGrid::eventcuts accept neutrino cut "<<endl;
        };
    };

    if (mydebug)
        if (okneut) cout<<" MyGrid::eventcuts neutrino cut passed "<<endl;
        else if (mydata->NeutrinoCut()) cout<<" MyGrid::eventcuts neutrino cut not passed "<<endl;

    if (mydebug)
        if (oklept) cout<<" MyGrid::eventcuts lepton cut passed "<<endl;
        else if (mydata->LeptonCut()) cout<<" MyGrid::eventcuts lepton cut not passed "<<endl;

    if (mydebug)
        if (onejet) cout<<" MyGrid::eventcuts jet cut passed "<<endl;
        else if (mydata->JetCut()) cout<<" MyGrid::eventcuts jet cut not passed "<<endl;

    if (mydata->JetCut()) accept=onejet;
    else accept=true;

    if (mydebug&&accept) cout<<" MyGrid::eventcuts event passed 1 "<<endl;

    if (accept) {
        if (mydata->NeutrinoCut()&&!okneut) {
            if (debug) cout<<" MyGrid::eventcut cut by neutrino "<<endl;
            return false;
        }
        if (mydata->LeptonCut()  &&!oklept) {
            if (debug) cout<<" MyGrid::eventcut cut by lepton "<<endl;
            return false;
        }
    }

    if (mydebug&&accept) cout<<" MyGrid::eventcuts event passed 2 "<<endl;

    if (accept&&mydata->MtwCut()) {
        if (mydebug) cout<<" MyGrid::eventcuts idlep= "<<idlep<<" idneu= "<<idneu<<endl;
        if (idlep!=-1&&idneu!=-1) {
            fastjet::PseudoJet pw=*myevent->Get(idlep)+*myevent->Get(idneu);;
            double mtw=pw.mt();
            if (mydebug) cout<<" MyGrid::eventcuts mtw= "<<mtw<<endl;
            if (mtw<mydata->GetMtWCut()) return false;
        };
    };

    if (mydebug&&accept) cout<<" MyGrid::eventcuts event passed 3 "<<endl;

    if (accept&&mydata->GetNJetCut()) {
        if (mydebug) cout<<" MyGrid::eventcuts njets= "<<njets<<endl;
        if(njets<mydata->GetNJet()) return false;
        if (mydebug) if (!accept) cout<<" MyGrid::eventcuts event rejected by Njet cut"<<endl;
    }

    if (mydebug&&accept) cout<<" MyGrid::eventcuts event passed 4 "<<endl;

    for (int  i = 0; i < myevent->GetN(); i++) {
// jet are filled after the lepton
        int id=myevent->GetID(i);
        bool isjet=myevent->IsJet(id);
        if (onejet&&isjet&&oklept&&mydata->DeltaLepJetCut()) {
            double dR=999.;
            if (idjetmax<1) cout<<" MyGrid::eventcuts something is wrong idjetmax= "<<idjetmax<<endl;
            if (idlep<1)    cout<<" MyGrid::eventcuts something is wrong idlep= "<<idlep<<endl;
            if (idlep!=-1) dR=(myevent->Get(i))->delta_R(*myevent->Get(idlep));
            if (mydebug) cout<<i<<" MyGrid::eventcuts dR= "<<dR<<" cut= "<<mydata->GetDeltaRLepJet()<<endl;
            if (dR<mydata->GetDeltaRLepJet()) accept=false;
            if (mydebug) if (!accept) cout<<" MyGrid::eventcuts event cut by dR"<<endl;
        };
    };

    if (mydebug&&accept) cout<<" MyGrid::eventcuts event passed 5 "<<endl;

    if (mydebug||debug)
        if (accept) cout<<" MyGrid::eventcuts event accepted !"<<endl;
        else        cout<<" MyGrid::eventcuts event rejected ! "<<endl;

    return accept;
};

void MyGrid::AddGrid(MyGrid *myOtherGrid) {
    cout<<"MyGrid::AddGrid: called..."<<endl;


    //get all of the other grids reference histograms
    std::vector<std::vector<TH1D*> > *otherRefLOsubprocessesHistos;
    otherRefLOsubprocessesHistos = myOtherGrid->GetRefLOsubprocessesHistos();

    std::vector<std::vector<TH1D*> > *otherRefsubprocessesHistos;
    otherRefsubprocessesHistos = myOtherGrid->GetRefsubprocessesHistos();

    std::vector<TH1D*> *otherReferenceHistos;
    otherReferenceHistos = myOtherGrid->GetReferenceHistos();

    std::vector<TH1D*> *otherReferenceFineHistos;
    otherReferenceFineHistos = myOtherGrid->GetReferenceFineHistos();

/*
    //TESTING TO MAKE SURE REFERENCES WERE ADDED CORRECTLY
    //cout<<"TEST: REFERENCE HISTOGRAMS BEFORE..."<<endl;
    int Ngrids=this->GetNGrid();
    for(int igrid = 0; igrid < Ngrids; igrid++)
    {
        cout<<"TEST: MyGrid::AddGrid: hreference print for igrid: "<<igrid<<endl;
        hreference.at(igrid)->Print("all");
        cout<<"TEST: MyGrid::AddGrid: hreferencefine print for igrid: "<<igrid<<endl;
        hreferencefine.at(igrid)->Print("all");
        cout<<"TEST: MyGrid::AddGrid: OTHERhreference print for igrid: "<<igrid<<endl;
        otherReferenceHistos->at(igrid)->Print("all");
        cout<<"TEST: MyGrid::AddGrid: OTHERhreferencefine print for igrid: "<<igrid<<endl;
        otherReferenceFineHistos->at(igrid)->Print("all");
    }
*/

    //make sure all histos can be added together
    if(hreference.size()!=otherReferenceHistos->size())
        std::cout<<" MyGrid::AddGrid: WARNING: hreference are the not the same but will attemp to be added!"<<std::endl;

    if(hreferencefine.size()!=otherReferenceFineHistos->size())
        std::cout<<" MyGrid::AddGrid: WARNING: hreferencefine are the not the same but will attemp to be added!!"<<std::endl;

    if(hrefLOsubprocesseshistos.size()!=otherRefLOsubprocessesHistos->size())
        std::cout<<" MyGrid::AddGrid: WARNING: hrefLOsubprocesseshistos are the not the same but will attemp to be added!"<<std::endl;

    if(hrefsubprocesseshistos.size()!=otherRefsubprocessesHistos->size())
        std::cout<<" MyGrid::AddGrid: WARNING: hrefsubprocesseshistos are the not the same but will attemp to be added!"<<std::endl;


    //add this grid's and the other grid's reference histos
    for(int href_x=0; href_x<hreference.size(); href_x++) {
        hreference.at(href_x)->Add(otherReferenceHistos->at(href_x));
    }

    for(int hreffine_x=0; hreffine_x<hreferencefine.size(); hreffine_x++) {
        hreferencefine.at(hreffine_x)->Add(otherReferenceFineHistos->at(hreffine_x));
    }

    for(int hrefLO_x=0; hrefLO_x<hrefLOsubprocesseshistos.size(); hrefLO_x++) {
        for(int hrefLO_y=0; hrefLO_y<hrefLOsubprocesseshistos.at(hrefLO_x).size(); hrefLO_y++) {
            hrefLOsubprocesseshistos.at(hrefLO_x).at(hrefLO_y)->Add(otherRefLOsubprocessesHistos->at(hrefLO_x).at(hrefLO_y));
        }
    }

    for(int hrefsub_x=0; hrefsub_x<hrefsubprocesseshistos.size(); hrefsub_x++) {
        for(int hrefsub_y=0; hrefsub_y<hrefsubprocesseshistos.at(hrefsub_x).size(); hrefsub_y++) {
            hrefsubprocesseshistos.at(hrefsub_x).at(hrefsub_y)->Add(otherRefsubprocessesHistos->at(hrefsub_x).at(hrefsub_y));
        }
    }

/*
    //TESTING TO MAKE SURE REFERENCES WERE ADDED CORRECTLY
    //cout<<"TEST: REFERENCE HISTOGRAMS AFTER..."<<endl;
    Ngrids=this->GetNGrid();
    for(int igrid = 0; igrid < Ngrids; igrid++)
    {
        cout<<"TEST: MyGrid::AddGrid: hreference print for igrid: "<<igrid<<endl;
        hreference[igrid]->Print("all");
        cout<<"TEST: MyGrid::AddGrid: hreferencefine print for igrid: "<<igrid<<endl;
        hreferencefine[igrid]->Print("all");
    }
*/


    //add event counters
    //cout<<"TEST: THISalluncorrevents: "<<alluncorrevents<<", OTHERalluncorrevents: "<<myOtherGrid->GetUncorrellatedEventNumber()<<std::endl;
    //cout<<"TEST: THISallevents: "<<allevents<<", OTHERmyOtherGrid->GetTotalEventNumber(): "<<myOtherGrid->GetTotalEventNumber()<<std::endl;
    alluncorrevents+=myOtherGrid->GetUncorrellatedEventNumber();
    allevents+=myOtherGrid->GetTotalEventNumber();
    //cout<<"TEST: THISalluncorrevents: "<<alluncorrevents<<std::endl;
    //cout<<"TEST: THISallevents: "<<allevents<<std::endl;


    //check if adding grids is safe
    int numGrids = myOtherGrid->GetNGrid();
    if( numGrids != this->GetNGrid() ) {
        std::cout<<" MyGrid::AddGrid: WARNING: Number of grids are the not the same but will attemp to be added!"<<std::endl;
        exit(0);
    }

    //add grids
    /*
    for(int igrid=0; igrid<numGrids; igrid++) {
        *(GetGrid(igrid)) += *(myOtherGrid->GetGrid(igrid));
    }
    */


    return;
}
