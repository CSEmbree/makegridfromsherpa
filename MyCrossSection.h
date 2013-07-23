#ifndef MyCrossSection_H
#define MyCrossSection_H

#include <stdlib.h> // exit()
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <string>
//#include <map>

using std::string;
#include <sys/time.h> 

#include "root.h"
//#include "MySubProcess.h"
#include "MyData.h"
#include "MyFrame.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/mcfmw_pdf.h"
#include "appl_grid/generic_pdf.h"
#include "theory_error_info.h"

#include "OptionHandler.h"

//#include "appl_grid/appl_pdf.h"



class MyCrossSection {
 private:

  bool debug; 

  string crosssectionname;
  string glabel;
  string subprocesssteername;

  //string ntupdirinput;
  //string ntupdiroutput;

  string ntupname;
  string steername;

  string datanamedir;
  string gridnamedir;


  double xlabel;
  double ylabel;

  std::vector<string> gridname;
  std::vector<string> dataname;
  std::vector<string> corrname;
  std::vector<string> vardesc;
  std::vector<int> markerstyle;
  std::vector<int> markercolor;
  std::vector<double> mcscalex; // factor to scale x-values
  std::vector<double> datascalex; // factor to scale x-values
  std::vector<double> scaley; // factor to scale y-values
  std::vector<int> frameid;   // grids/data with same frameid will be overlayed
  std::vector<int> divideid;  // grids/data with same divideid will be divided
  std::vector<int> refhistlinestyle; // line style of reference histogram
  std::vector<int> refhistlinecolor; // line color of reference histogram
  std::vector<long unsigned int> events; // number of events from grid
  std::vector<MyFrame*> framepointer; // pointer to MyFrame


  int processnumber; // process number foresee setting from reaction

  string pdf_function;

  std::vector<appl::grid*> mygrid;   // grid vector  
  //std::vector<appl::grid> mygrid;  // grid vector
  std::vector<bool> isBooked;        // flag grid already booked 
  std::vector<MyData*> mydata;       // information about data from steering file
  generic_pdf *mypdf;

 public:
  std::vector<int> PDFSetCodes_vec;
  bool do_PDFBand;
  bool do_AlphaS;
  bool do_RenormalizationScale;
  bool do_FactorizationScale;
  int ErrorSize;

  MyCrossSection(char name[100]);

  generic_pdf * GetSubProcess() { return mypdf;};
  void SetSubProcess(generic_pdf *subpro) { mypdf=subpro; return;};
  void split_string(std::string str, std::vector<std::string>& split_results, std::string delimiters);

  void Initialize();
  void ReadSteering(char fname[100]);
  //void ReadSteeringOptions(char fname[100]);
  bool file_exists(const string& s);

  int reflinestyle;
  void SetLineColor(int igrid, int rl) {refhistlinecolor[igrid]=rl; return;};
  void SetLineStyle(int igrid, int rl) {refhistlinestyle[igrid]=rl; return;};

  void Normalise(TH1D* h, double yscale, double xscale, bool normtot);
  int GetNGrid(){return gridname.size();};
  //string GetNtupDirInput(){ return ntupdirinput;};
  //string GetNtupDirOutput(){ return ntupdiroutput;};
  string GetNtupName(){ return ntupname;};

  string GetGridName(int igrid){ return gridnamedir+"/"+gridname[igrid];};
  TString GetTStringGridName(int igrid){ return ((TString) (gridnamedir+gridname[igrid]));  };
  TString GetVarDesc(int igrid){ return ((TString) (vardesc.at(igrid))); }
 /*
  string GetGridName(int i){
   TString name=TString(gridname[i]); 
   name.ReplaceAll(".txt","");
   return string(name.Data());
  };

  string GetGridFileName(int i){
    //return string(this->GetGridName(i).Data())+".root";
    return this->GetGridName(i)+".root";
  }


  string GetGridFullFileName(int i){
    //return ntupdirinput+"/"+this->GetGridFileName(i);
    return this->GetGridDir()+"/"+this->GetGridFileName(i);
  }
 */
  MyData *GetMyData(int igrid){ 
   if (igrid>mydata.size()){
    cout<<" MyCrossSection::GetMyData something is wrong mydata too short igrid= "<<igrid<<endl;
    return 0;
   }
   return mydata[igrid];
  };

  int GetProcessNumber() {return processnumber;};
  void SetProcessNumber(int nproc) {processnumber=nproc; return;};

  string GetGridDir(){ return gridnamedir;};
  bool GetDataOk(int igrid){
   bool flag=true;
   //if (debug) cout<<" MyCrossSection::GetDataOk i= "<<i<<" dataname.size()= "<<dataname.size()<<endl;
   if (dataname.size()<=igrid) flag=false;
   //if (debug) cout<<" MyCrossSection::GetDataOk flag= "<<flag<<endl;
   return flag;
  }                
  string GetDataName(int igrid){
   string bad="";
   if (!this->GetDataOk(igrid)) return bad;
   else return datanamedir+"/"+dataname[igrid];
  };
  
  string GetCorrelationsName(int igrid){ 
   if (corrname.size()<=igrid) {
    cout<<" MyCrossSection::GetCorrelationsName corrname vector too short or correlation matrix not found ! igrid= "<<igrid<<endl;
    return "";
   }
   return datanamedir+"/"+corrname[igrid]; 
  };

  int GetMarkerStyle(int igrid) {
   if (markerstyle.size()<=igrid) {
    cout<<" MyCrossSection::GetMarkerStyle something is wrong markserstyle too short igrid= "<<igrid<<endl;
    return 0;
   }
   return markerstyle[igrid];
  };
  
  int GetMarkerColor(int igrid) {
   if (markercolor.size()<=igrid) {
    cout<<" MyCrossSection::GetMarkerColor something is wrong markercolor too short igrid= "<<igrid<<endl;
    return 0;
   }
   return markercolor[igrid];
  };
  
  void SetStyleMarker (int igrid, int ms) {
   if (igrid>markerstyle.size()){
     cout<<" MyCrossSection::SetStyleMarker something is wrong markerstyle too short igrid= "<<igrid<<endl;
   } else markerstyle[igrid]=ms; 
   return;
  };
  
  void SetColorMarker (int igrid, int mc) {
   if (igrid>markerstyle.size()){
     cout<<" MyCrossSection::SetColorMarker something is wrong markercolor too short igrid= "<<igrid<<endl;
   } else markercolor[igrid]=mc; 
   return;
  };

  double GetDataScaleX(int igrid){ 
   if (datascalex.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleX something is wrong markercolor too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return datascalex[igrid];
  };
  
  double GetMCScaleX(int igrid){ 
   if (datascalex.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleX something is wrong markercolor too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return mcscalex[igrid];
  };
  
  double GetScaleY(int igrid){ 
   if (scaley.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleY something is wrong scaley too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return scaley[igrid];
  };
  
  int GetFrameID(int igrid){ 
   if (frameid.size()<=igrid) {
     cout<<" MyCrossSection::GetFrameID something is wrong frameid too short igrid= "<<igrid<<endl;
     return -1;
   }
   return frameid[igrid];
  };
  
  MyFrame *GetMyFrame(int igrid ) {
   if (framepointer.size()<=igrid) {
     cout<<" MyCrossSection::GetMyFrame something is wrong divideid too short igrid= "<<igrid<<endl;
     return 0;
   }
   return framepointer[igrid];
  };
  
  int GetDivideID(int igrid){ 
   if (divideid.size()<=igrid) {
     cout<<" MyCrossSection::GetDivideID something is wrong divideid too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return divideid[igrid];
  };

  int GetFrameNumber();

  appl::grid *GetGrid(int igrid){
    //  if (!mygrid[igrid]) cout<<" MyCrossSection::GetReference mygrid["<<igrid<<"] not filled ! "<<endl;
   return mygrid[igrid];
  };

  TH1D *GetReference(int igrid);
  void Draw(int igrid);
 // void DrawPDFs(theory_error_calc *my_theory_calcs[], bool use_this_pdf[], TString x_title, float x_min, float x_max, bool first_of_canv);
  void DrawErrors(theory_error_calc *my_theory_calcs[], bool use_this_pdf[], TString x_title, float x_min, float x_max, bool first_of_canv, int error_code);
  void DrawData(int igrid){
   mydata[igrid]->DrawData();
   return;
  };
 void SetLabelX(double x) {xlabel=x; return;};
 void SetLabelY(double y) {ylabel=y; return;};

 void DrawReference(int igrid);
 TGraphAsymmErrors* GetReferenceRatio(int igrid);
 TH1D* GetNormalisedReference(int igrid);

 void DrawinFrame(int iframe);

 double GetObsMax(int igrid){
  if (!mydata[igrid])  cout<<" GetObsMax grid not found= mydata["<<igrid<<"]"<<endl;
  double obsmax=mydata[igrid]->GetMaxX(); 
  //if (debug) cout<<" obsmax= "<<obsmax<<endl;
  return obsmax;
 }

  TH1D *GetHisto(int igrid,string name="htest"){
   if (!mydata[igrid]) cout<<" MyCrossSection::GetHisto mydata not filled ! "<<endl;
   int nobs=mydata[igrid]->GetNBins(); 
   if (debug) 
    cout<<" MyCrossSection::GetHisto nobs= "<<nobs<<endl;
   double *xbins=mydata[igrid]->GetBins();

   string hname=name+this->GetGridName(igrid);
   //cout<<" book histogram hname= "<<hname<<endl;
   TH1D *htest = new TH1D(hname.c_str(),hname.c_str(),nobs,xbins);
   if (!htest) cout<<" htest["<<igrid<<"] not found "<<endl;
   //htest->Print("all");
   return htest;
  }; 

  double GetTotalEventNumber(int igrid){
   if (events.size()<=igrid) return 0; 
   return events[igrid];
  };

  TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr=1);
  TGraphAsymmErrors* TH1TOTGraphAsymm(TH1 *h1);

  void Print();
};


#endif
