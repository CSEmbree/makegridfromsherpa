//
//   for sherpa
//
using namespace std;
#include "appl_grid/generic_pdf.h"
#include "MyCrossSection.h"


/******************************************************************
 ** Method Implementations
 ******************************************************************/

MyCrossSection::MyCrossSection(char name[100])
{
    debug=false;
    crosssectionname="";

    events.clear();

    mcscalex.clear();
    datascalex.clear();
    scaley.clear();
    frameid.clear();
    divideid.clear();
    framepointer.clear();

    xlabel=0.75;
    ylabel=0.95;

    processnumber=1;

    pdf_function="";
    subprocesssteername="";
    gridnamedir="";
    datanamedir="";
    ntupname="";
    //ntupdiroutput="";
    //ntupdirinput="";

    steername=name;
    if (debug) cout<<" MyCrossSection:: steering file= "<<steername<<endl;
    gridname.clear();
    vardesc.clear();
    dataname.clear();
    corrname.clear();
    refhistlinestyle.clear();
    refhistlinecolor.clear();
    markerstyle.clear();
    markercolor.clear();
    mydata.clear();

    PDFSetCodes_vec.clear();
    do_PDFBand = false;
    do_AlphaS = false;
    do_RenormalizationScale = false;
    do_FactorizationScale = false;
    ErrorSize = -999;
    //PDFBandType_vec.clear();
// PDFErrorSize_vec.clear();

    this->ReadSteering(name);
    //this->ReadSteeringOptions(name);

    this->Initialize();

    //if (debug) this->Print();

}

void MyCrossSection::Initialize() {

    if (debug) cout<<" MyCrossSection:Initialize()"<<endl;

    //generic_pdf *myodf= new generic_pdf(subprocesssteername);
    //this->SetSubProcess(mypdf);

    cout<<" MyCrossSection::Initialize Number of grids to produce "<<gridname.size()<<endl;
    //sherpaw_pdf *mypdf=new sherpaw_pdf(pdf_function,mysub);
    //
    // Loop over grid
    //
    for (int  igrid = 0; igrid < gridname.size(); igrid++) {

        if (debug) cout<<" MyCrossSection:Initialize Data "<<this->GetDataName(igrid)<<endl;
        if (this->GetDataOk(igrid)) {
            MyData *mydatatmp= new MyData;
            mydatatmp->ReadData(TString(this->GetDataName(igrid)).Data());
            //if (strcmp(this->GetCorrelationsName(igrid).c_str()," ")==0){
            //cout<<" name= "<<this->GetCorrelationsName(igrid).size()<<endl;
            if (this->GetCorrelationsName(igrid).size()!=0) {
                mydatatmp->ReadCorrelations(TString(this->GetCorrelationsName(igrid)).Data());
            } else cout<<" MyCrossSection:Initialize no correlation file found ! "<<endl;

            //if (this->GetMarkerStyle(igrid)!=20)
            mydatatmp->SetMarkerStyle(this->GetMarkerStyle(igrid));
            //else markerstyle[igrid]=mydatatmp->GetMarkerStyle();
            //if (this->GetMarkerColor(igrid)!=1)
            mydatatmp->SetMarkerColor(this->GetMarkerColor(igrid));
            //else markercolor[igrid]=mydatatmp->GetMarkerColor();
            mydata.push_back(mydatatmp);
        } else {
            if (debug) cout<<" MyCrossSection::Initialize no data found "<<endl;
        }
        string fname=this->GetGridName(igrid);
        if (debug) cout<<" MyCrossSection::Initialize Grid Name "<<fname<<endl;

        if (this->file_exists(fname)) {
            if (debug) cout<<" MyCrossSection::Initialize grid file "<<fname<<" exist ! "<<endl;
            isBooked.push_back(false);

            //appl::grid g(fname);
            if (debug) cout<<" MyCrossSection::Initialize open "<<fname<<endl;
            appl::grid *tmpgrid = new appl::grid(fname);
            if (debug) std::cout << " MyCrossSection::Initialize grid CMS energy " << tmpgrid->getCMSScale() << std::endl;
            //cout<<" mytest run= "<<tmpgrid->run()<<endl;
            mygrid.push_back(tmpgrid);
            if (debug) std::cout << " MyCrossSection::Initialize tmpgrid pushed" << std::endl;

            events.push_back(int(mygrid[igrid]->run()));
            if (debug) std::cout << " MyCrossSection::Initialize events pushed" << std::endl;
            //mygrid[igrid]->getDocumentation();
            //cout<<" print htest"<<endl;
            //TH1D* htest=this->GetReference(igrid);
            //htest->Print("all");

        } else {
            cout<<" MyCrossSection::Initialize  file not found "<<fname<<endl;
        }
        if (debug) std::cout << " MyCrossSection::Initialize end of igrid loop" << std::endl;

    }

    if (debug) cout<<" MyCrossSection::Initialize finished  "<<endl;
    return;

}

void MyCrossSection::ReadSteering(char fname[100]) {
    //
    // read steering from file fname
    //
    steername=fname;
    if (debug) cout<<" MyCrossSection::ReadSteering steering "<<steername<<endl;

    ifstream infiletmp(steername.c_str(), ios::in);
//if( debug ) std::cout << " MyCrossSection::ReadSteering Why am I not here? " << std::endl;
    if(!infiletmp) { // Check open
        cerr << " MyCrossSection::ReadSteering Can't open " << steername <<"\n";
        infiletmp.close();
        exit (1);
    } else {
        if (debug) cout <<" MyCrossSection::ReadSteering read data file: " << steername << endl;
    }

    int iline=0;
    int nsyst=1;

    char line[1024];
    char text[256];
    char name[256];

    while (infiletmp.good()) {
        //if (debug) cout << " good: " << infile.good() << " eof: " << infile.eof() << endl;
        //if (!infiletmp.good()) break;
        infiletmp.getline(line,sizeof(line),'\n');
        if(line[0] != 'g') continue;
        if (strstr(line,"gridname")!=0) {
            sscanf(line," %s %[^\n] ",text, name);
            if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
            string myname=name;
            gridname.push_back(myname);
        }
    }


    if (debug) cout<<" MyCrossSection::ReadSteering number of grids "<<gridname.size()<<endl;

// set up options vector, one option for each grid
    for (int igrid=0; igrid<gridname.size(); igrid++) {
        if (debug) cout<<" MyCrossSection::ReadSteering gridname["<<igrid<<"]= "<<gridname[igrid]<<endl;
        markerstyle.push_back(20);
        markercolor.push_back(1);
        mcscalex.push_back(1.);
        datascalex.push_back(1.);
        scaley.push_back(1.);
        frameid.push_back(-1);
        divideid.push_back(-1);
        refhistlinestyle.push_back(1.);
        refhistlinecolor.push_back(1.);
    }

    ifstream infile(steername.c_str(), ios::in);

    int igrid=-1;


    //load in all valid cross section options
    string CSOptionsFileName = "CrossSectionOptions.txt"; //<--could be made readable from steering file
    OptionHandler *mycsoptions= new OptionHandler(CSOptionsFileName);
    //mycsoptions->generateResultFile(); //generates *.txt file showing results of OptionHandler::isKnownOption execution
    char lineFirstWord[1024]; //same length as 'line' string buffer

    //read and set all options and data
    while (infile.good()) {
        //if (debug) cout << " good: " << infile.good() << " eof: " << infile.eof() << endl;
        //if (!infile.good()) break;
        infile.getline(line,sizeof(line),'\n');
        std::string cpp_line(line);
        //std::vector<std::string> colon_split_cpp_line;
        //colon_split_cpp_line.clear();
        //split_string(cpp_line, colon_split_cpp_line, ":");
        //if (debug) std::cout << " MyCrossSection:ReadSteering ReadColon-split of line " << cpp_line << " has size " << colon_split_cpp_line.size() << "\n";
//
        sscanf(line," %s",lineFirstWord);
        std::string cpp_lineFirstWord(lineFirstWord);

        if (debug) cout<< " MyCrossSection::ReadSteering line= "<< line << "\n";
        if(line[0] != '%') { //ignore comments
            /*
            if(line[0] != 'g' && line[0] != 'n' && line[0] != 's' && line[0] != 'p' && line[0] != 'm'
                    && line[0] != 'd'   && line[0] != 's' && line[0] != 'v' && line[0] != 'c' && line[0] != 'G'
                    && line[0] != 'r'   && line[0] != 'f'
              ) {
             */
            if(mycsoptions->isKnownOption(cpp_lineFirstWord)==false) {

                // unsupported read, perhaps data?
                // do something
                //} else if (strstr(line,"nsub")!=0) {
                //char text[100];
                //sscanf(line," %s %d ",text, &nsub);
                //if (debug) printf("**MyCrossSection Read: Nsub= %d   \n",nsub);
            } else if (strstr(line,"debug")!=0) {
                debug=true;
            } else if (strstr(line,"nprocessnumber")!=0) {
                sscanf(line," %s %d ",text, &processnumber);
                if (debug) printf(" MyCrossSection::ReadSteering processnumber= %d   \n",processnumber);
            } else if (strstr(line,"subprocesssteername")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                subprocesssteername=string(name);
                if (debug) printf(" MyCrossSection::ReadSteering subprocesssteername= %s \n",subprocesssteername.c_str());
            } else if (strstr(line,"pdffunction")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                pdf_function=string(name);
                if (debug) cout<<" MyCrossSection:ReadSteering pdffunction= "<<pdf_function<<endl;
                //   } else if (strstr(line,"ntupdiroutput")!=0) {
                //sscanf(line," %s %[^\n] ",text, name);
                //if (debug) cout<<"MyCrossSection Read "<<text<<" "<<name<<endl;
                //ntupdiroutput=name;
                //} else if (strstr(line,"ntupdirinput")!=0) {
                //sscanf(line," %s %[^\n] ",text, name);
                //if (debug) cout<<"MyCrossSection Read "<<text<<" "<<name<<endl;
                //ntupdirinput=name;
            } else if (strstr(line,"ntupname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
                ntupname=name;
            } else if (strstr(line,"Gridnamedir")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
                gridnamedir=name;
                if (debug) cout<<" MyCrossSection::ReadSteering gridnamedir= "<<" "<<gridnamedir<<endl;
            } else if (strstr(line,"datanamedir")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
                datanamedir=name;
            } else if (strstr(line,"gridname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
                string myname=name;
                //gridname.push_back(myname);
                igrid++;
                cout<<" MyCrossSection::ReadSteering igrid= "<<igrid<<" "<<name<<endl;
            } else if (strstr(line,"markerstyle")!=0) {
                int mymarker;
                sscanf(line," %s %d ",text, &mymarker);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<mymarker<<endl;
                if (igrid<0) cout<<" MyCrossSection::ReadSteering something wrong ! "<<endl;
                markerstyle[igrid]=mymarker;
            } else if (strstr(line,"markercolor")!=0) {
                int mymarker;
                sscanf(line," %s %d ",text, &mymarker);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<mymarker<<endl;
                markercolor[igrid]=mymarker;
            } else if (strstr(line,"reflinestyle")!=0) {
                int mystyle;
                sscanf(line," %s %d ",text, &mystyle);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<mystyle<<endl;
                refhistlinestyle[igrid]=mystyle;
            } else if (strstr(line,"reflinecolor")!=0) {
                int mycolor;
                sscanf(line," %s %d ",text, &mycolor);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<mycolor<<endl;
                refhistlinecolor[igrid]=mycolor;
            } else if (strstr(line,"datascalex")!=0) {
                double sx;
                sscanf(line," %s %lf ",text, &sx);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<sx<<endl;
                datascalex[igrid]=sx;
            } else if (strstr(line,"mcscalex")!=0) {
                double sx;
                sscanf(line," %s %lf ",text, &sx);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<sx<<endl;
                mcscalex[igrid]=sx;
            } else if (strstr(line,"scaley")!=0) {
                double sy;
                sscanf(line," %s %lf ",text, &sy);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<sy<<endl;
                scaley[igrid]=sy;
            } else if (strstr(line,"frameid")!=0) {
                int frid;
                sscanf(line," %s %d ",text, &frid);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<frid<<" igrid= "<<igrid<<" frameid= "<<frid<<endl;
                frameid[igrid]=frid;
            } else if (strstr(line,"divideid")!=0) {
                int did;
                sscanf(line," %s %d ",text, &did);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<did<<endl;
                divideid[igrid]=did;
            } else if (strstr(line,"vardesc")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection:ReadSteering "<<text<<" "<<name<<endl;
                if(debug ) std::cout << "Filling vardesc" << std::endl;
                string myname=name;
                vardesc.push_back(myname);
            } else if (strstr(line,"dataname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection::ReadSteering "<<text<<" "<<name<<endl;
                string myname=name;
                dataname.push_back(myname);
            } else if (strstr(line,"corrname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyCrossSection:ReadSteering our name: "<<text<<" "<<name<<endl;
                string myname=name;
                corrname.push_back(myname);
            }


            ////// What theory error types will we consider and how will they be displayed?
            for(int pdfi = 0; pdfi < e_n_PDFs; pdfi++) {
                if(cpp_line == PDF_strs[pdfi]) {
                    PDFSetCodes_vec.push_back(pdfi);
                }
            }
            if( cpp_line == "PDFErrorBand" )            do_PDFBand = true;
            if( cpp_line == "AlphaS" )                  do_AlphaS = true;
            if( cpp_line == "RenormalizationScale" )    do_RenormalizationScale = true;
            if( cpp_line == "FactorizationScale" )      do_FactorizationScale = true;
            if( cpp_line == "OneSigma" )                ErrorSize = e_OneSigma;
            if( cpp_line == "90Percent" )               ErrorSize = e_90Percent;

//
//////// Old pdf_calc_class interface
//
//   if (colon_split_cpp_line.size() > 0 ) {
//
//     std::cout << "MyCrossSection Read in PDF choices for line " << line << "?\n";
//     for(int pdfi = 0; pdfi < e_n_PDFs; pdfi++) {
//       std::cout << "compare " << colon_split_cpp_line.at(0) << " with " << PDF_strs[pdfi] << "\n";
//       int comp_val = (colon_split_cpp_line.at(0) == PDF_strs[pdfi] );
//       std::cout << "comparison value: " << comp_val << "\n";
//       if( colon_split_cpp_line.at(0) == PDF_strs[pdfi] ) {
//         if( colon_split_cpp_line.size() != 3 ) {
//           std::cout << "Error, " << line << " seems to communicate PDF info but does not have format SetName:ErrorType:ConfidenceLevel\n";
//           break;
//         }
//         PDFSetCodes_vec.push_back(pdfi);
//         if( colon_split_cpp_line.at(1) == PDFBandType_strs[e_UseAlphaS] ) PDFBandType_vec.push_back(e_UseAlphaS);
//         else if( colon_split_cpp_line.at(1) == PDFBandType_strs[e_UseErrorBand] ) PDFBandType_vec.push_back(e_UseErrorBand);
//         else {
//           std::cout << "Error, PDF error type is " << colon_split_cpp_line.at(1) << ", expected " << PDFBandType_strs[e_UseAlphaS]<< " or " << PDFBandType_strs[e_UseErrorBand] << "\n";
//           PDFSetCodes_vec.erase( PDFSetCodes_vec.end());
//           break;
//         }
//         std::string IsOneSigmaStr = "OneSigma";
//         std::string Is90PerStr = "90Percent";
//         if( colon_split_cpp_line.at(2) == PDFErrorSize_strs[e_OneSigma] ) PDFErrorSize_vec.push_back(e_OneSigma);
//         else if( colon_split_cpp_line.at(2) == PDFErrorSize_strs[e_90Percent] ) PDFErrorSize_vec.push_back(e_90Percent);
//         else {
//           std::cout << "Error, ConfidenceLevel type is " << colon_split_cpp_line.at(2) << ", expected " << PDFErrorSize_strs[e_OneSigma] << " or " << PDFErrorSize_strs[e_90Percent] << "\n";
//           PDFSetCodes_vec.erase( PDFSetCodes_vec.end());
//           PDFBandType_vec.erase( PDFErrorSize_vec.end());
//           break;
//         }
//       }  /// does first string match one of our PDFs?
//     } /// pdfi
//     std::cout << "After reading in line " << line << ", here is a spew of the PDF specifications:" << std::endl;
//     for(int pdfspi = 0; pdfspi < (int) PDFSetCodes_vec.size(); pdfspi++) {
//       std::cout << "pdfspi = " << pdfspi << ", PDFSet = " << PDF_strs[PDFSetCodes_vec.at(pdfspi)] << ", PDF Band Type = " << PDFBandType_strs[PDFBandType_vec.at(pdfspi)] << ", PDF Error Size = " << PDFErrorSize_strs[PDFErrorSize_vec.at(pdfspi)] << std::endl;
//     }
//   } /// colo_split_cpp_line.size()
        };
    };
    return;
};

/*
void MyCrossSection::ReadSteeringOptions(char fname[100]){

 steername=fname;
 if (debug) cout<<" MyCrossSection::ReadSteeringOptions steering "<<steername<<endl;

 ifstream infile(steername.c_str(), ios::in);
 if(!infile){ // Check open
  cerr << "Can't open " << steername <<"\n";
  infile.close();
  exit (1);
 } else {
   //cout <<" MyCrossSection:ReadSteeringOptions read data file: " << steername << endl;
 }

 int iline=0;
 char line[1024];

 //TString myname;
 //TString mynameold=TString("myoldname");

 // set up options vector, one option for each grid
 for (int igrid=0; igrid<gridname.size(); igrid++){
  markerstyle.push_back(20);
  markercolor.push_back(1);
  datascalex.push_back(1.);
  scaley.push_back(1.);
  frameid.push_back(-1);
  divideid.push_back(-1);
  refhistlinestyle.push_back(1.);
  refhistlinecolor.push_back(1.);
 }

 int igrid=-1;
 while (1) {
  if (debug) cout << " MyCrossSection::ReadSteeringOptions good: " << infile.good() << " eof: " << infile.eof() << endl;
  if (!infile.good()) break;
  infile.getline(line,sizeof(line),'\n');
//
  if (debug) cout<< " MyCrossSection::ReadSteeringOptions line= "<< line << "\n";
  char text[100];     char name[100];
  if(line[0] != '%'){
   if (strstr(line,"gridname")!=0) {
    igrid++;
    sscanf(line," %s %[^\n] ",text, name);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<name<<endl;
    if (debug) cout<<" MyCrossSection: ReadSteeringOptions igrid= "<<igrid<<endl;
   } else if (strstr(line,"markerstyle")!=0) {
    int mymarker;
    sscanf(line," %s %d ",text, &mymarker);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<mymarker<<endl;
    markerstyle[igrid]=mymarker;
   } else if (strstr(line,"markercolor")!=0) {
    int mymarker;
    sscanf(line," %s %d ",text, &mymarker);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<mymarker<<endl;
    markercolor[igrid]=mymarker;
   } else if (strstr(line,"reflinestyle")!=0) {
    int mystyle;
    sscanf(line," %s %d ",text, &mystyle);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<mystyle<<endl;
    refhistlinestyle[igrid]=mystyle;
   } else if (strstr(line,"reflinecolor")!=0) {
    int mycolor;
    sscanf(line," %s %d ",text, &mycolor);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<mycolor<<endl;
    refhistlinecolor[igrid]=mycolor;
   } else if (strstr(line,"datascalex")!=0) {
    double sx;
    sscanf(line," %s %lf ",text, &sx);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<sx<<endl;
    datascalex[igrid]=sx;
   } else if (strstr(line,"scaley")!=0) {
    double sy;
    sscanf(line," %s %lf ",text, &sy);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<sy<<endl;
    scaley[igrid]=sy;
   } else if (strstr(line,"frameid")!=0) {
    int frid;
    sscanf(line," %s %d ",text, &frid);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<frid<<endl;
    cout<<" MyCrossSection::ReadSteeringOptions igrid= "<<igrid<<" frameid= "<<frid<<endl;
    frameid[igrid]=frid;
   } else if (strstr(line,"divideid")!=0) {
    int did;
    sscanf(line," %s %d ",text, &did);
    if (debug) cout<<" MyCrossSection::ReadSteeringOptions "<<text<<" "<<did<<endl;
    divideid[igrid]=did;
   }
  }
 }
 return;
}
*/

void MyCrossSection::Print() {

    cout<<" MyCrossSection::Print >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    cout<<" MyCrossSection::Print steering: "<<steername<<endl;
    cout<<" MyCrossSection::Print subprocesssteername: "<<subprocesssteername<<endl;
    cout<<" MyCrossSection::Print pdf_function=  "<<pdf_function<<endl;
    //cout<<" ntupdir= "<<ntupdirinput<<" ntupoutput= "<<ntupdiroutput<<endl;
    cout<<" MyCrossSection::Print directory for gridname=  "<<gridnamedir<<endl;
    cout<<" MyCrossSection::Print directory for datanamedir=  "<<datanamedir<<endl;

    cout<<" MyCrossSection::Print Number of grids N= "<<gridname.size()<<endl;

    for (int  i = 0; i <   gridname.size(); i++)
    {
        if (this->GetDataOk(i)) {
            if (debug) cout<<" print i= "<<i
                               <<" dataname.size()= "<<dataname.size()
                               <<" events.size()=   "<<events.size()
                               <<" gridname.size()= "<<gridname.size()<<endl;

            cout<<i
                <<" grid= "<<gridname[i]
                <<" data= "<<dataname[i]
                <<" events= "<<events[i]
                <<" style: "<<this->GetMarkerStyle(i)
                <<" color: "<<this->GetMarkerColor(i)
                <<" datascalex: "<<this->GetDataScaleX(i)
                <<" mcscalex: "<<this->GetMCScaleX(i)
                <<" scaley: "<<this->GetScaleY(i)
                <<" frameid: "<<this->GetFrameID(i)
                <<" divideid: "<<this->GetDivideID(i)
                <<endl;
            //mydata[i]->Print();
        } else {
            cout<<i
                <<" grid= "<<gridname[i]
                //<<" events= "<<events[i]
                <<" datascalex: "<<this->GetDataScaleX(i)
                <<" mcscalex: "<<this->GetMCScaleX(i)
                <<" scaley: "<<this->GetScaleY(i)
                <<" frameid: "<<this->GetFrameID(i)
                <<" divideid: "<<this->GetDivideID(i)
                <<endl;
        }
    }

    cout<<" MyCrossSection <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;

}

bool MyCrossSection::file_exists(const string& s) {
    if ( FILE* testfile=fopen(s.c_str(),"r") ) {
        fclose(testfile);
        return true;
    }
    else return false;
}

//
// just normalise to bin width
//
void MyCrossSection::Normalise(TH1D* h1, double yscale, double xscale=1., bool normtot=false)
{

//xscale = 1.;
    Double_t x, y, ey;
    Double_t sigtot=0.;

    for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
        y=h1->GetBinContent(i)*yscale;
        x=h1->GetBinWidth(i);
        sigtot+=y*x;
    }
//sigtot *= yscale;

    for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
        x =h1->GetBinWidth(i);
        y =h1->GetBinContent(i)*yscale*x;
        ey=h1->GetBinError(i)  *yscale*x;
        x =h1->GetBinWidth(i)  *xscale;
        if (x!=0) h1->SetBinContent(i,y/x);
        else      h1->SetBinContent(i,0.);
        std::cout << "BinWidth: " << h1->GetBinWidth(i)  << ", xscale: " << xscale << ", x: " << x << ", bin content: " << h1->GetBinContent(i) << "\n";

        if (x!=0) h1->SetBinError(i,ey/x);
        else     h1->SetBinError(i,0.);
    }

    if (normtot) {
        if (debug) cout<<" MyCrossSection::Normalise sigtot= "<<sigtot<<endl;
        if (sigtot!=0.) {
            if (debug) std::cout << " MyCrossSection::Normalise What is the name? " << h1->GetName() << std::endl;
            if (debug) std::cout << " MyCrossSection::Normalise New norm = " << h1->Integral() << std::endl;
        }
        h1->Scale(1. / sigtot);
    }

    /*
    for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ )
      {
        double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);
        h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
      }
    */
    if (debug) std::cout << " MyCrossSection::Normalise return" << std::endl;
    return;
}

void MyCrossSection::DrawErrors(theory_error_calc *my_theory_calcs[], bool use_this_pdf[], TString x_title, float x_min, float x_max, bool first_of_canv, int error_code)
{
    std::cout << "DrawErrors 1" << std::endl;
    TGraphAsymmErrors *Theory_graph_for_draw[e_n_PDFs];
    double x=0.75, y=0.8;
    for(int pdfi = 0; pdfi < e_n_PDFs; pdfi++) {
        std::cout << "pdfi: " << pdfi << ", error code = " << error_code << ", use it? ";
        std::cout << use_this_pdf[pdfi] << std::endl;
        if( !use_this_pdf[pdfi] ) continue;
        if( error_code == e_PDFBandCode ) {
            std::cout << "Do PDFBand" << std::endl;
            std::cout << "PDF_strs = " << PDF_strs[pdfi] << std::endl;
            std::cout << "calc desc = " << my_theory_calcs[pdfi]->calc_desc << std::endl;
            std::cout << "PDFBand name = " << my_theory_calcs[pdfi]->h_PDFBand_results->GetName() << std::endl;
            Theory_graph_for_draw[pdfi] = (TGraphAsymmErrors*) my_theory_calcs[pdfi]->h_PDFBand_results->Clone((TString) (PDF_strs[pdfi] + "_error_bar_graph"));
        } else if( error_code == e_AlphaSCode ) {
            Theory_graph_for_draw[pdfi] = (TGraphAsymmErrors*) my_theory_calcs[pdfi]->h_AlphaS_results->Clone((TString) (PDF_strs[pdfi] + "_error_bar_graph"));
        } else if( error_code == e_RenormalizationScaleCode ) {
            Theory_graph_for_draw[pdfi] = (TGraphAsymmErrors*) my_theory_calcs[pdfi]->h_RenormalizationScale_results->Clone((TString) (PDF_strs[pdfi] + "_error_bar_graph"));
        } else if( error_code == e_FactorizationScaleCode ) {
            Theory_graph_for_draw[pdfi] = (TGraphAsymmErrors*) my_theory_calcs[pdfi]->h_FactorizationScale_results->Clone((TString) (PDF_strs[pdfi] + "_error_bar_graph"));
        } else if( error_code == e_TotErrorCode ) {
            Theory_graph_for_draw[pdfi] = (TGraphAsymmErrors*) my_theory_calcs[pdfi]->h_TotError_results->Clone((TString) (PDF_strs[pdfi] + "_error_bar_graph"));
        }
    }
    std::cout << "DrawErrors 2" << std::endl;
    bool first_time_pdf = true;
    for(int pdfi = 0; pdfi < e_n_PDFs; pdfi++) {
        if( !use_this_pdf[pdfi] ) continue;
        float min_y = 99999;
        float max_y = 0.;
        for(int pi = 0; pi < Theory_graph_for_draw[pdfi]->GetN(); pi++) {
            Double_t x_val;
            Double_t y_val;
            Theory_graph_for_draw[pdfi]->GetPoint(pi, x_val, y_val);
            if( y_val > max_y ) max_y = y_val;
            if( y_val < min_y ) min_y = y_val;
        }
        max_y *= 40.;
        min_y *= 0.2;
        std::cout << "DrawErrors 2.1, pdfi = " << pdfi << std::endl;
        Theory_graph_for_draw[pdfi]->GetYaxis()->SetRangeUser(min_y, max_y);
        if( x_min < x_max) Theory_graph_for_draw[pdfi]->GetXaxis()->SetRangeUser(x_min, x_max);
        Theory_graph_for_draw[pdfi]->GetXaxis()->SetTitle(x_title);
        Theory_graph_for_draw[pdfi]->SetTitle("");
        if( pdfi == e_CT10 ) {
            Theory_graph_for_draw[pdfi]->SetFillColor(kGreen+2);
            Theory_graph_for_draw[pdfi]->SetFillStyle(3005);
        }
        if( pdfi == e_MSTW2008nlo68cl ) {
            Theory_graph_for_draw[pdfi]->SetFillColor(kRed+2);
            Theory_graph_for_draw[pdfi]->SetFillStyle(3004);
        }
        if( pdfi == e_NNPDF23_nlo_as_0118 ) {
            Theory_graph_for_draw[pdfi]->SetFillColor(kBlue+2);
            Theory_graph_for_draw[pdfi]->SetFillStyle(3002);
        }
        if( pdfi == e_HERAPDF15NLO ) {
            Theory_graph_for_draw[pdfi]->SetFillColor(kGreen+10);
            Theory_graph_for_draw[pdfi]->SetFillStyle(3021);
        }
        if( first_time_pdf && first_of_canv ) Theory_graph_for_draw[pdfi]->Draw("A E2");
        else Theory_graph_for_draw[pdfi]->Draw("E2 same");
        first_time_pdf = false;
    }
    std::cout << "DrawErrors 3" << std::endl;
    //mydata[igrid]->DrawLegend(experiment.c_str(),x,y);
    //cout<<" DrawReference"<<endl;
    // this->DrawReference(igrid);
}


TH1D *MyCrossSection::GetReference(int igrid) {
    if (!mygrid[igrid]) cout<<" MyCrossSection::GetReference mygrid["<<igrid<<"] not filled ! "<<endl;
//cout<<" MyCrossSection::GetReference mygrid["<<igrid<<"]  ! "<<mygrid[igrid]<<endl;
    if (debug) std::cout << " MyCrossSection::GetReferenceFunction for " << igrid << " out of " << mygrid.size() << std::endl;
    TH1D *href=(TH1D*)mygrid[igrid]->getReference();
    if (!href) {
        cout<<"MyCrossSection::GetReference: reference histo not found igrid= "<<igrid<<endl;
        cout <<"MyCrossSection::GetReference: grid name=" << this->GetGridName(igrid) << endl;
        exit(0);
    }
    if (debug) std::cout <<" MyCrossSection::GetReferenceFunction Got reference" << std::endl;
    if (this->GetDataOk(igrid)) {
        double yfac=mydata[igrid]->GetUnitfbFactor();
        double xfac=mydata[igrid]->GetUnitGeVFactor();
        if (debug) std::cout <<" MyCrossSection::GetReferenceFunction Got x and y factors" << std::endl;
        if (debug) cout<<" MyCrossSection::GetReference yfac= "<<yfac<<" xfac= "<<xfac<<endl;
        bool normtot=mydata[igrid]->isNormTot();
        if (debug) std::cout <<" MyCrossSection::GetReferenceFunction Normalise the hist normtot= "
                                 <<normtot<< std::endl;
        //href->Print("all");
        if( debug ) std::cout << "For igrid set initial reference normalization factors as xfac = " << xfac << ", yfac = " << yfac << "\n";
        Normalise(href,yfac,xfac,normtot);
        if (debug) std::cout <<" MyCrossSection::GetReferenceFunction Normalised" << std::endl;

        //href->Print("all");
    } else std::cout <<" MyCrossSection::GetReferenceFunction No data found, not normalised " << std::endl;
    if (debug) std::cout <<" MyCrossSection::GetReferenceFunction Return href" << std::endl;
    href->SetMarkerSize(0.01);
    return href;
}

void MyCrossSection::Draw(int igrid) {
    if (debug) std::cout << " MyCrossSection::Draw igrid " << igrid << std::endl;

// mydata[igrid]->DrawExperimentandYear(xlabel,ylabel);
    mydata[igrid]->DrawData();
    if (debug) mydata[igrid]->Print();
    return;
}


void MyCrossSection::DrawinFrame(int iframe) {

    MyFrame *myframe= new MyFrame(600,600);
    framepointer.push_back(myframe);

    std::vector<int> gridid;
    gridid.clear();

    double datascalex=1.2;
    double scaley=1.2;

    const double BIG=1.e30;
    double ymin= BIG;
    double ymax=-BIG;
    double xmin= BIG;
    double xmax=-BIG;

    double x,y;
    TString name="";

    for (int igrid=0; igrid<this->GetNGrid(); igrid++) {
        int jframe=this->GetFrameID(igrid);
        if (debug) cout<<" MyCrossSection::DrawFrame iframe= "<<iframe<<" jframe= "<<jframe<<endl;
        double xscale=this->GetDataScaleX(igrid);
        double yscale=this->GetScaleY(igrid);

        if (jframe==iframe) {
            if (debug) std::cout << " MyCrossSection::DrawinFrame igrid= " << igrid << ", push back this grid id" << std::endl;
            gridid.push_back(igrid);

            if (this->GetDataOk(igrid)) {
                //mydata[igrid]->Print();
                if (xscale!=1||yscale!=1)
                    mydata[igrid]->Scale(xscale,yscale);
                if (debug) mydata[igrid]->Print();
                name+=(TString) mydata[igrid]->GetFileName();

                crosssectionname=string(name.Data());

                std::cout << " MyCrossSection::DrawinFrameFile name = " << crosssectionname << "\n";

                if (mydata[igrid]->GetLogY1()) {
                    myframe->SetLogY1();
                    scaley=2.0;
                }
                if (mydata[igrid]->GetLogX ()) {
                    myframe->SetLogX();
                    datascalex=1.5;
                }

                //cout<<" datascalex= "<<datascalex<<" scaley= "<<scaley<<endl;

                y=mydata[igrid]->GetMinY();
                //cout<<" min y= "<<y<<endl;
                if (ymin>y) ymin=y;

                y=mydata[igrid]->GetMaxY();
                //cout<<" max y= "<<y<<endl;
                if (ymax<y) ymax=y;

                x=mydata[igrid]->GetMinX();
                //cout<<" min x= "<<x<<endl;
                if (xmin>x) xmin=x;

                x=mydata[igrid]->GetMaxX();
                //cout<<" max x= "<<x<<endl;
                if (xmax<x) xmax=x;
            } else {
                std::cout << " MyCrossSection::DrawinFrameFile No data found ! "  << "\n";
                TH1D* href=this->GetReference(igrid);
                if (!href) cout<<" MyCrossSection::DrawinFrameFile reference not found ! "  << "\n";
                int imax=href->GetNbinsX();
                double bw=href->GetBinWidth(0);
                //cout<<" MyCrossSection::DrawinFrameFile imax= "<<imax<<" bw= "<<bw <<endl;
                xmin=href->GetBinCenter(0)   -2*bw;
                xmax=href->GetBinCenter(imax)+2*bw;
                ymin=href->GetMinimum();
                ymax=href->GetMaximum();

                if (debug)
                    cout<<" MyCrossSection::DrawinFrameFile xmin= "<<xmin<<" xmax= "<<xmax
                        <<" ymin= "<<ymin<<" ymax= "<<ymax
                        <<endl;
            }
        }   /// grid frame matches expectations
    }    /// loop over grid frames

    if (debug) cout<<" MyCrossSection::DrawinFrameFile datascalex= "<<datascalex<<" scaley= "<<scaley<<endl;

    if   (xmin>0) xmin/=datascalex;
    else  xmin*=datascalex;

    if   (xmax>0) xmax*=datascalex;
    else  xmax/=datascalex;

    if   (ymin>0) ymin/=scaley;
    else  ymin*=scaley;

    if   (ymax>0) ymax*=scaley;
    else  ymax/=scaley;

    if (debug)
        cout<<" MyCrossSection::DrawFrame iframe= "<<iframe<<" xmin= "
            <<xmin<<" xmax= "<<xmax<< " ymin= "<<ymin<<" ymax= "<<ymax<<endl;

    if (gridid.size()==0) {
        cout<<" MyCrossSection::DrawFrame no grids found ! iframe= "<<iframe<<endl;
        return;
    }

    int igrid=gridid[0];

    if (!myframe) cout<<" MyCrossSection::DrawinFrameFile myframe not found ! "<<endl;
    else if (debug) cout<<" MyCrossSection::DrawinFrameFile SetFrame name= "<<name<<endl;


    y=ylabel;
    TString titx="", tity="";
    if (this->GetDataOk(igrid)) {
        titx=mydata[igrid]->GetTitleX();
        tity=mydata[igrid]->GetTitleY();
    } else {
        titx=this->GetGridName(igrid);
        //tity=this->GetGridName(igrid);
        name=this->GetGridName(igrid);
    }
    myframe->SetFrameName(name);
    myframe->SetSubPad2TitleX(titx);
    myframe->SetSubPad1TitleY(tity);

    if (debug) cout<<" MyCrossSection::DrawinFrameFile MyCrossSection y= "<<y<<endl;

    myframe->SetXmin( xmin);
    myframe->SetXmax( xmax);
    myframe->SetYmin1(ymin);
    myframe->SetYmax1(ymax);
    myframe->SetYmin2(1.1);
    myframe->SetYmax2(0.9);

    if (debug) cout<<" MyCrossSection now DrawFrame  "<<myframe<<endl;

    myframe->DrawFrame();
    if (debug) cout<<" SubPad = "<< myframe->GetSubPad1()<<endl;
    myframe->GetSubPad1()->cd();

    if (debug)
        cout<<" MyCrossSection::DrawFrame number of grids= "<<gridid.size()<<
            " for frame= "<<iframe<<endl;

    for (int i=0; i<gridid.size(); i++) {
        int igrid=gridid[i];
        if (debug)
            cout<<" MyCrossSection::DrawFrame DrawData for i= " << i << "i grid= "<<igrid<<endl;
        //mydata[igrid]->Print();
        TH1D* href=0;
        std::cout << "Is data ok? " << this->GetDataOk(igrid) << std::endl;
        if (this->GetDataOk(igrid)) {
            if (!mydata[igrid]) cout<<" MyCrossSection::DrawinFrame mydata["<<igrid<<"] not found "<<endl;
            else if (debug)     cout<<" MyCrossSection::DrawinFrame mydata["<<igrid<<"]  found "<<endl;

            mydata[igrid]->DrawData();
            //mydata[igrid]->SetMarkerStyle(this->GetMarkerStyle(igrid));
            //mydata[igrid]->SetMarkerColor(this->GetMarkerColor(igrid));
            y-=0.05;
            mydata[igrid]->DrawExperimentandYear(xlabel,y);
            href=this->GetNormalisedReference(igrid);

        } else
            href=this->GetReference(igrid);
        if (debug) cout<<" MyCrossSection::DrawinFrame now draw reference igrid= "<<igrid<<endl;
//this->DrawReference(igrid);

        if (!href) {
            cout<<" MyCrossSection::DrawinFrame hreference not found "<<endl;
            exit;
        }
        href->Draw("same");
        if( debug ) {
            for(int bi = 1; bi <= href->GetNbinsX(); bi++) {
                Double_t x_val;
                Double_t y_val;
                mydata[igrid]->GetTGraphTotErr()->GetPoint(bi-1, x_val, y_val);
                std::cout << "reference center: " << href->GetBinCenter(bi) << ", reference value: " << href->GetBinContent(bi) << ", data center: " << x_val << ", value: " << y_val << "\n";
            }
        }

    }


    myframe->GetSubPad2()->cd();
    myframe-> SetSubPad2TitleOffsetX(0.8);
    TString titlename="Data/";
    titlename+=" NLO QCD"; // better read this in from somewhere
    if (!this->GetDataOk(igrid)) titlename="";
    myframe->GetYAxis2()->SetTitle(titlename.Data());
//myframe->GetYAxis2()->SetTitle("Data/NLO");

    const int NRatioGrids=gridid.size();
    if (debug) cout<< " MyCrossSection::DrawinFrame NRatioGrids= "<<NRatioGrids<<endl;

    TGraphAsymmErrors* ratiotot[NRatioGrids];
    TGraphAsymmErrors* ratiostat[NRatioGrids];

    double Ymin=BIG, Ymax=-BIG, Xmin=BIG, Xmax=-BIG;


    for (int i=0; i<NRatioGrids; i++) {
        if (debug) cout<<" MyCrossSection::DrawinFrame  i= " << i <<endl;
        int igrid=gridid[i];
        if (debug)
            cout<<" MyCrossSection::DrawinFrame  i= " << i << " gridid[i] = " << gridid[i] << ", Draw reference for igrid= "<<igrid<<endl;
        //mydata[igrid]->Print();
        TH1D *href=this->GetNormalisedReference(igrid);

        TGraphAsymmErrors* gref=TH1TOTGraphAsymm(href);
        if (this->GetDataOk(igrid)) {
            ratiotot[i]=myTGraphErrorsDivide(mydata[igrid]->GetTGraphTotErr(),gref);
            ratiostat[i]=myTGraphErrorsDivide(mydata[igrid]->GetTGraphStatOnly(),gref);

            //cout<<" ms= "<<this->GetMarkerStyle(igrid)<<endl;
            //cout<<" mc= "<<this->GetMarkerColor(igrid)<<endl;
            ratiotot[i]->SetName("ratioreference"+name);
            for(int pi = 0; pi < gref->GetN(); pi++) {
                Double_t data_x_val;
                Double_t data_y_val;
                mydata[igrid]->GetTGraphTotErr()->GetPoint(pi, data_x_val, data_y_val);
                Double_t ref_x_val;
                Double_t ref_y_val;
                gref->GetPoint(pi, ref_x_val, ref_y_val);
                Double_t ratiototx_val;
                Double_t ratiototy_val;
                ratiotot[i]->GetPoint(pi, ratiototx_val, ratiototy_val);
                Double_t ratiostatx_val;
                Double_t ratiostaty_val;
                ratiostat[i]->GetPoint(pi, ratiostatx_val, ratiostaty_val);

                if( debug ) std::cout << "calc ratio for center: " << data_x_val << ", data val: " << data_y_val << ", ref val: " << ref_y_val << ", ratiotot val: " << ratiototy_val << ", ratio stat: " << ratiostaty_val << "\n";
            }

            //ratiotot->SetMarkerStyle(20);
            //ratiotot->SetMarkerColor(1);
            ratiotot[i]->SetMarkerStyle(this->GetMarkerStyle(igrid));
            ratiotot[i]->SetMarkerColor(this->GetMarkerColor(igrid));
            ratiostat[i]->SetMarkerStyle(this->GetMarkerStyle(igrid));
            ratiostat[i]->SetMarkerColor(this->GetMarkerColor(igrid));
            //ratiostat->SetMarkerStyle(this->GetMarkerStyle(igrid));
            //ratiostat->SetMarkerColor(this->GetMarkerColor(igrid));

            ratiotot[i]->ComputeRange(xmin,ymin,xmax,ymax);
            if (ymin<Ymin) Ymin=ymin;
            if (ymax>Ymax) Ymax=ymax;
        } else
            cout<<" MyCrossSection::DrawinFrame  data not found ! "<<endl;
    }

    if (debug) cout<<" MyCrossSection::DrawinFrame Ymin= "<<Ymin<<" Ymax= "<<Ymax<<endl;

    myframe->GetYAxis2()->SetRangeUser(Ymin*0.9,Ymax*1.1);

    if (debug) cout<<" MyCrossSection::DrawinFrame gridid.size= "<<gridid.size()<<endl;

    for (int i=0; i<NRatioGrids; i++) {
        if (debug) cout<<" MyCrossSection::DrawinFrame i= "<<i<<endl;
        if (!ratiostat[i])cout<<" MyCrossSection::DrawinFrame ratiostat not found igrid= "<<igrid<<endl;
        if (!ratiotot[i])cout<<" MyCrossSection::DrawinFrame ratiotot not found igrid= "<<igrid<<endl;
        /*
        cout<<" print stat: "<<ratiostat[i]->GetName()<<endl;
        ratiostat[i]->Print();
        cout<<" print tot: "<<ratiostat[i]->GetName()<<endl;
        ratiotot[i]->Print();
        */
        if (this->GetDataOk(igrid)) {
            ratiostat[i]->Draw("p same");
            ratiotot[i]->Draw("p same");
        }
    }
    std::cout << "End of DrawInFrame" << std::endl;
    return;
}

void MyCrossSection::DrawReference(int igrid) {
//
// Draw reference histograms for grid with index igrid
//
    TH1D *href=this->GetReference(igrid);
    if (!href) cout<<" MyCrossSection::DrawReference reference histo not found igrid="<<igrid<<endl;
    //href->Print("all");
    href->SetLineStyle(refhistlinestyle[igrid]);
    href->SetLineColor(refhistlinecolor[igrid]);
    href->Draw("same");
    return;
}

TH1D* MyCrossSection::GetNormalisedReference(int igrid) {
//
// Get reference histograms for grid with index igrid
// and normalise according to datascalex and scaley
//
    TH1D *href=this->GetReference(igrid);
    if (!href) cout<<" MyCrossSection::DrawReference reference histo not found igrid="<<igrid<<endl;
    //href->Print("all");
    href->SetLineStyle(refhistlinestyle[igrid]);
    href->SetLineColor(refhistlinecolor[igrid]);
    double xscale=this->GetMCScaleX(igrid);
    double yscale=this->GetScaleY(igrid);
    if (debug)
        cout<<" MyCrossSection::DrawReference xscale="<<xscale<<" yscale= "<<yscale<<endl;
    if( debug ) {
        std::cout << "Before normalization href contents for igrid " << igrid << " are\n";
        for(int bi = 1; bi <= href->GetNbinsX(); bi++) {
            float center = href->GetBinCenter(bi);
            float content = href->GetBinContent(bi);
            std::cout << "center = " << center << ", content = " << content << "\n";
        }
    }
    this->Normalise(href,yscale,xscale);
    if( debug ) {
        std::cout << "After normalization href contents for igrid " << igrid << " are\n";
        for(int bi = 1; bi <= href->GetNbinsX(); bi++) {
            float center = href->GetBinCenter(bi);
            float content = href->GetBinContent(bi);
            std::cout << "center = " << center << ", content = " << content << "\n";
        }
    }
    return href;
}

TGraphAsymmErrors* MyCrossSection::GetReferenceRatio(int igrid) {
//
// Get reference histograms for grid with index igrid
// and divide data associated to grid
//
    TH1D *href=this->GetReference(igrid);
    TGraphAsymmErrors* gref=TH1TOTGraphAsymm(href);
//gref->Print();
    TGraphAsymmErrors* ratio=myTGraphErrorsDivide(mydata[igrid]->GetTGraphTotErr(),gref);
    return ratio;
};

TGraphAsymmErrors* MyCrossSection::TH1TOTGraphAsymm(TH1 *h1) {
    //
    // convert the histogram h1 into a graph
    //
    const char *name="**TH1TOTGraphAsymm";
    if (!h1) cout <<name<< " histogram not found !" << endl;
//else h1->Print();
//cout<<name<<" n= "<<h1->GetNbinsX()<<endl;
//TGraphErrors* g1;
    TGraphAsymmErrors* g1= new TGraphAsymmErrors();
//if (!g1) cout<<name<<" graph g1 not created "<<endl;
//g1->Print();

    Double_t x, y, ex, ey;
    for (Int_t i=0; i<h1->GetNbinsX(); i++) {
        //cout<<name<<" i= "<<i<<endl;
        y=h1->GetBinContent(i+1);
        ey=h1->GetBinError(i+1);
        x=h1->GetBinCenter(i+1);
        ex=h1->GetBinWidth(i+1)/2.;

        //cout << i<<" x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

        g1->SetPoint(i,x,y);
        g1->SetPointError(i,ex,ex,ey,ey);

    }

//g1->Print();
//cout<<name<<" return "<<endl;
    return g1;
}

TGraphAsymmErrors* MyCrossSection::myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr) {
// Divide two TGraphAsymmErrors: new=g1/g2
//
// noerr=0: put all errors to zero
//       1: add errors from two graph quadrically
//       2: set errors from graph 2 to zero

    const bool debug=false;
    const char *name="**myTGraphErrorsDivide:";

    if (!g1) cout << name<<" g1 does not exist ! " << endl;
    if (!g2) cout << name<<" g2 does not exist ! " << endl;

    Int_t n1=g1->GetN();
    Int_t n2=g2->GetN();

    if (n1!=n2) {
        printf("%s: vector do not have same number of entries !  \n",name);
        cout <<name<< " g1: " << g1->GetName() << " n1= " << n1 << endl;
        cout <<name<< " g2: " << g2->GetName() << " n2= " << n2 << endl;
    }

    TGraphAsymmErrors* g3= new TGraphAsymmErrors();
    if (!g3) cout<<name<<" problem to make new vector ! " << endl;
    g3->SetName       (g1->GetName());
    g3->SetMarkerStyle(g1->GetMarkerStyle());
    g3->SetMarkerColor(g1->GetMarkerColor());
    g3->SetLineColor  (g1->GetLineColor());

    Double_t   x1=0.,   y1=0., x2=0., y2=0.;
    Double_t dx1h=0., dx1l=0.;
    Double_t dy1h=0., dy1l=0.;
    Double_t dy2h=0., dy2l=0.;

// Double_t* X1 = g1->GetX();
//Double_t* Y1 = g1->GetY();
    Double_t* EXhigh1 = g1->GetEXhigh();
    Double_t* EXlow1 =  g1->GetEXlow();
    Double_t* EYhigh1 = g1->GetEYhigh();
    Double_t* EYlow1 =  g1->GetEYlow();

//Double_t* X2 = g2->GetX();
//Double_t* Y2 = g2->GetY();
    Double_t* EXhigh2 = g2->GetEXhigh();
    Double_t* EXlow2 =  g2->GetEXlow();
    Double_t* EYhigh2 = g2->GetEYhigh();
    Double_t* EYlow2 =  g2->GetEYlow();

    Int_t iv=0;
    for (Int_t i1=0; i1<n1; i1++) {  //loop over point of graph1
        Int_t matchcount=0;
        for (Int_t i2=0; i2<n2; i2++) {//loop over point of graph2
            g1->GetPoint(i1,x1,y1);
            g2->GetPoint(i2,x2,y2);
            Double_t emean=(EXhigh1[i1]+EXhigh2[i2]+EXlow1[i1]+EXlow2[i2])/4.;
            if (fabs(x1-x2)>emean) {
                //cout <<name<<" x1 and x2 not the same x1= "<<x1<<" x2= "<<x2<<endl;
            } else { // do something only if x1=x2
                matchcount++;
                //cout <<name<<" x1 and x2 match x1= "<<x1<<" x2= "<<x2<<endl;
                dx1h  = EXhigh1[i1];
                dx1l  = EXlow1[i1];
                if (y1!=0.) dy1h  = EYhigh1[i1]/y1;
                else        dy1h  = 0.;
                if (y2!=0.) dy2h  = EYhigh2[i2]/y2;
                else        dy2h  = 0.;
                if (y1!=0.) dy1l  = EYlow1 [i1]/y1;
                else        dy1l  = 0.;
                if (y2!=0.) dy2l  = EYlow2 [i2]/y2;
                else        dy2l  = 0.;

                if (debug)
                    printf("%s: %d %d dy1=%f %f dy2=%f %f sqrt= %f %f \n",name,i1,i2,dy1l,dy1h,dy2l,dy2h,
                           sqrt(dy1l*dy1l+dy2l*dy2l),sqrt(dy1h*dy1h+dy2h*dy2h));

                if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
                else        g3->SetPoint(iv, x1,y2);
                Double_t el=0.;
                Double_t eh=0.;

                if (noerr==2) {
                    dy2l=0.;
                    dy2h=0.;
                }
                if (y1!=0. && y2!=0.) el=sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
                if (y1!=0. && y2!=0.) eh=sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);


                if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
                if (noerr==0) {
                    g3->SetPointError(iv,dx1l,dx1h,0,0);
                } else {
                    g3->SetPointError(iv,dx1l,dx1h,el,eh);
                }

                iv++;
            }
        }
        //if (matchcount>1) {cout<<name<<" too many x-points matched ! "<<endl; exit (1);}
        if (matchcount>1) {
            cout<<name<<" too many x-points matched ! "<<endl;
        }
    }
    return g3;
}

int MyCrossSection::GetFrameNumber() {
//
// count frames/canvases that are given in the steering files
//
    int nframe=0;
    std::vector<int> frameid2;
    frameid2.clear();
    for (int i=0; i<frameid.size(); i++) {
        int iframe=frameid[i];
        bool found=false;
        for (int i2=0; i2<frameid2.size(); i2++) {
            int iframe2=frameid2[i2];
            if (iframe==iframe2) {
                found=true;
            };
        };
        if (!found) frameid2.push_back(iframe);
        //cout<<" iframe= "<<iframe<<" found= "<<found<<endl;
    };
    if (debug) cout<<" MyCrossSection::GetFrameNumber Number of frames= "<<frameid2.size()<<endl;
    return frameid2.size();
};



void MyCrossSection::split_string(std::string str, std::vector<std::string>& split_results, std::string delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        split_results.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

