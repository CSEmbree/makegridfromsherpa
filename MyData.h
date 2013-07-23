#ifndef MyData_H
#define MyData_H

#include "root.h"
//#include "VariableDefinitions.h"
#include "TMatrixT.h"
#include "string.h"
#include <sstream>

#include "OptionHandler.h"

#define TEXT_SIZE 100

class MyData {
private:
    bool debug;

    bool unitsxTeV; // if true units are in TeV, default is GeV
    bool unitsypb;  // if true units are in pb, default is fb
    bool normtot;   // if true y-axis is normalised by sigtot

    bool dividebybinwith; //data are divided by bin width

    //  TGraphAsymmErrors * xtdatavector;
    //  TGraphAsymmErrors * xtdatavectorstat;
    //  TGraphAsymmErrors * xtdatavectorsyst;

    double ymin, ymax, sqrts;
    double scale;

    //std::vector<std::vector<double> > corr_matrix;
    TMatrixT<double>  *cov_matrix;
    void split_string(std::string str, std::vector<std::string>& split_results, std::string delimiters);

    double miny, maxy, minx, maxx;
    double framexmin;
    bool setframexmin;

    TString experiment;
    TString filename;
    TString corrfilename;
    TString jetalgo;
    int jetalgoR;

    int markerstyle;
    int markercolor;
    double markersize;

    double msize;

    int year;

    TString titlex;
    TString titley;

    TString observable;
    TString reaction;

    bool cutMtw;
    double Mtwcut;

    bool neucut;
    double ptneucut;

    bool jetcut;
    double ptjetcut;
    double yjetmin, yjetmax;

    bool lepcut;
    double ptlepcut;
    double ylepmin, ylepmax;

    bool drljcut;
    double deltaRLepJet;

    bool blogy, blogx, bliny, blinx;

    bool njetcut;
    int njets;

    bool incljets;

    //<<<<<<< .mine
    TGraphAsymmErrors * datavector;
    TGraphAsymmErrors * datavectorstat;
    TGraphAsymmErrors * datavectorsyst;
    TGraphAsymmErrors * datavectortoterr;


    //=======

    //>>>>>>> .r147604
    
    std::string GetFileNameFromPath(std::string filepath){
        std::string filename = filepath;
        const size_t last_slash = filename.find_last_of("\\/");
        if (std::string::npos != last_slash)
            filename.erase(0, last_slash + 1);

        const size_t period_ext = filename.rfind('.');
        if (std::string::npos != period_ext)
            filename.erase(period_ext);
            
        return filename;
    };
    

public:
//  std::vector<int> PDFSetCodes_vec;
//  std::vector<int> PDFBandType_vec;
//  std::vector<int> PDFErrorSize_vec;

    MyData();

    bool isNormTot() {
        return normtot;
    };

    double GetUnitGeVFactor() {// GeV is default unit
        double fac=1.;
        if (unitsxTeV) fac=1./1000.;
        return fac;
    }

    double GetUnitfbFactor() {// fb is default unit
        double fac=1.;
        if (unitsypb) fac=1./1000.;
        return fac;
    }

    TMatrixT<double> GetCovarianceMatrix() {
        return (*cov_matrix);
    }
    TGraphAsymmErrors * GetTGraph() {
        return datavector;
    };
    TGraphAsymmErrors * GetTGraphStatOnly() {
        return datavectorstat;
    };
    TGraphAsymmErrors * GetTGraphSyst() {
        return datavectorsyst;
    };
    TGraphAsymmErrors * GetTGraphTotErr() {
        return datavectortoterr;
    };


    //TGraphAsymmErrors * GetTGraphXt() {return xtdatavector;};
    //TGraphAsymmErrors * GetTGraphStatOnlyXt() {return xtdatavectorstat;};
    //TGraphAsymmErrors * GetTGraphSystXt() {return xtdatavectorsyst;};

    void ScaleGraph(TGraphAsymmErrors *g1, double scalex, double scaley);
    void Scale(double scalex, double scaley);

    bool DivideByBinWidth() {
        return dividebybinwith;
    };
    bool GetLogY1() {
        return blogy;
    };
    bool GetLinY1() {
        return bliny;
    };
    bool GetLogX () {
        return blogx;
    };
    bool GetLinX () {
        return blinx;
    };
    bool GetSetFrameXmin () {
        return setframexmin;
    };
    double GetFrameXmin  () {
        return framexmin;
    };

    int GetNJet()  {
        return njets;
    };
    bool GetNJetCut() {
        return njetcut;
    };

    TString  GetFileName() {
        return filename;
    };
    TString  GetObservable() {
        return observable;
    };
    TString  GetReaction()   {
        return reaction;
    };

    TString GetExperiment() {
        return experiment;
    };
    int GetYear() {
        return year;
    };
    void SetExperiment(TString myexp) {
        experiment=myexp;
        return;
    };

    bool isObsInclJets() {
        return incljets;
    };

    // cuts
    bool JetCut() {
        return jetcut;
    };
    bool LeptonCut() {
        return lepcut;
    };
    bool NeutrinoCut() {
        return neucut;
    };
    bool MtwCut() {
        return cutMtw;
    };
    bool DeltaLepJetCut() {
        return drljcut;
    };

    void DrawExperiment(double x=0.9, double y=0.9) {
        Double_t tsize=0.05;
        TMarker *marker = new TMarker(x-(0.6*tsize),y,8);
        marker->SetMarkerColor(markercolor);
        marker->SetNDC();
        marker->SetMarkerStyle(markerstyle);
        marker->Draw();

        TLatex l;
        l.SetTextAlign(12); //l.SetTextSize(tsize);
        l.SetNDC();
        l.SetTextSize(tsize);
        l.DrawLatex(x,y,experiment.Data());

        return;
    };

    void DrawExperimentandYear(double x=0.9, double y=0.9);
    double GetMtWCut()     {
        return Mtwcut;
    };
    double GetJetPtCut()   {
        return ptjetcut;
    };
    double GetLepPtCut()   {
        return ptlepcut;
    };
    double GetNeuPtCut()   {
        return ptneucut;
    };


    double GetJetMinY()    {
        return yjetmin;
    };
    double GetJetMaxY()    {
        return yjetmax;
    };
    double GetLepMinY()    {
        return ylepmin;
    };
    double GetLepMaxY()    {
        return ylepmax;
    };

    double GetDeltaRLepJet() {
        return deltaRLepJet;
    };

    double GetR() {
        return float(jetalgoR)/10.;
    };

    TString GetTitleY() {
        return titley;
    };
    TString GetTitleX() {
        return titlex;
    };

    double GetYMin() {
        return ymin;
    };
    double GetYMax() {
        return ymax;
    };
    double GetMinY() {
        return miny;
    };
    double GetMaxY() {
        return maxy;
    };
    double GetMinX() {
        return minx;
    };
    double GetMaxX() {
        return maxx;
    };
    double GetSQRTS() {
        return sqrts;
    };

    void ReadData(const char fname[200], double scale=1.);
    void ReadCorrelations(const char fname[200], double scale=1.);
    void DrawData();
    void DrawData(char name[100], float x, float y);
    void DrawDataStatOnly();
    //void DrawDataXtStatOnly();
    //void DrawDataXt();
    void DrawLegend(char name[100], float x, float y);

    void PutData(TH1* hdata, double s, double ymin, double ymax );
    //void TransformPtToXt();
    //void DrawXt(char name[100], float x, float y);

    void SetTextSize(double mysize) {
        msize=mysize;
    };
    void SetMarkerStyle(int myMarker) {
        markerstyle=myMarker;
    };
    void SetMarkerColor(int myMarker) {
        markercolor=myMarker;
    };
    void SetMarker(int myMarkerStyle, int  myMarkerColor) {
        markerstyle=myMarkerStyle;
        markercolor=myMarkerColor;
    };
    int GetMarkerStyle() {
        return markerstyle;
    };
    int GetMarkerColor() {
        return markercolor;
    };

    void Print();
    //void PrintXt();


    double *GetBins() {
        // this function is a badly programmed
        // should avoid new in function call...
        //
        const int nbins=1+datavector->GetN();
        double *bins = new double[nbins];
        double *exl=datavector->GetEXlow();
        double *exh=datavector->GetEXhigh();
        double *exp=datavector->GetX();
        for (int i=0; i<datavector->GetN(); i++) {
            bins[i]=exp[i]-exl[i];
            if (i==datavector->GetN()-1) bins[i+1]= exp[i]+exh[i];
        }
        return bins;
    };

    int GetNBins() {
        return datavector->GetN();
    };
};


#endif
