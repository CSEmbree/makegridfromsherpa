
#include "TSystem.h"
//#include "TH1D.h"
//#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
//#include "MyFrame.h"
//#include "MyFrameData.h"
//#include "MyData.h"
//#include "utils.h"
#include "TLegend.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"
#include "appl_grid/appl_timer.h"
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "LHAPDF.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/LHAPDFfw.h"
#include "hoppet_v1.h"
#include "TMath.h"
#include <assert.h>
#include <string>
#include "VariableDefinitions.h"
#include "TMatrixT.h"
//using namespace LHAPDF;



void GetPDF(const double& x, const double& Q, double* xf);
//extern "C" void initpdfsetc_(const char *);
//extern "C" void initpdf_(const int *);
//extern "C" double alphaspdf_(const double *);
//extern "C" void evolvepdf_(const double *, const double *, double *);

//inline void   evolvePDF(const double & x, const double & Q, double * f) {
//  evolvepdf_(&x , &Q, f);
//}

#ifndef VarCodes 
#define VarCodes
//enum enum_PDFs {e_CT10, e_MSTW2008nlo68cl, e_NNPDF23_nlo_as_0118, e_n_PDFs};
//static string PDF_strs[e_n_PDFs] = {"PDFsets/CT10.LHgrid", "PDFsets/MSTW2008nlo68cl.LHgrid", "PDFsets/NNPDF23_nlo_as_0118.LHgrid"};
//static TString PDF_file_name_strs[e_n_PDFs] = {"CT10", "MSTW2008nlo", "NNPDF23nlo"};
//static TString detailed_PDF_strs[e_n_PDFs] = {"CT10", "MSTW2008nlo", "NNPDF23nlo #alpha_{S} 120"};
//static int nPDFMembers[e_n_PDFs] = {53, 41, 101};

enum enum_PDFBandType {e_UseAlphaS, e_UseErrorBand, e_n_PDFBands};
std::string PDFBandType_strs[e_n_PDFBands] = {"UseAlphaS", "UseErrorBand"};
enum enum_PDFErrorSizeType {e_OneSigma, e_90Percent, e_n_PDFErrorSizeTypes};
std::string PDFErrorSize_strs[e_n_PDFErrorSizeTypes] = {"OneSigma", "90Percent"};

#endif

//static int nPDFMembers[e_n_PDFs] = {3, 3, 3};

class theory_error_calc {
  public:
 //   theory_error_calc();
    theory_error_calc(TString grid_name, int pdfi, bool do_i_do_PDFBand, bool do_i_do_AlphaS, bool do_i_do_RS, bool do_i_do_FS, int the_ErrorSizeType, double xscale);
    TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr=1);
    TGraphAsymmErrors* TH1TOTGraphAsymm(TH1 *h1);
    void GetRatioToTH1(TH1D* href);
    //  TGraphAsymmErrors* TGraph_href = TH1TOTGraphAsymm(h_reference);
    bool do_PDFBand;
    bool do_AlphaS;
    bool do_RenormalizationScale;
    bool do_FactorizationScale;
    
    TString calc_desc;

 //   int PDFBandType;
    int ErrorSizeType;
    //void FillPDFErrorSets();
    void CalcSystErrors();
    void CalcPDFBandErrors();
    void CalcAlphaSErrors();
    void CalcRenormalizationScaleErrors();
    void CalcFactorizationScaleErrors();
    void CalcTotErrors();
 //   void CalcThisPDFSyst(std::vector<TH1D*> err_hist, bool do_prenorm);
 //   void CalcThisAlphaSPDFSyst(std::vector<TH1D*> err_hist, bool do_prenorm);
    //void CalcThisPDFSyst(std::vector<TH1D*> err_hist);
    void CalcChi2(TGraphAsymmErrors *g_theory, TGraphAsymmErrors *g_data, TMatrixT<double> data_cov_matrix, double &chi2);
    void InitializeErrorGraphs();
    int pdf_code;

    
    
    TH1D* h_qqbar_prenorm;
    TH1D* h_gg_prenorm;
    TH1D* h_tot_prenorm;
    TH1D* h_qqbar;
    TH1D* h_gg;
    TH1D* h_tot;
    TH1D* h_gg_frac;
    TH1D* h_qqbar_frac;

    TH1 *TH1NormToTot(TH1 *h1, Double_t yscale, Double_t xscale=1.);

    appl::grid *my_grid;
    std::vector<TH1D*> h_errors_RenormalizationScale;
    std::vector<TH1D*> h_errors_FactorizationScale;
    std::vector<TH1D*> h_errors_PDFBand;
    std::vector<TH1D*> h_errors_prenorm;
    std::vector<TH1D*> h_errors_AlphaS;
    std::vector<TH1D*> h_errors_AlphaS_prenorm;
    TGraphAsymmErrors *h_PDFBand_results;
    TGraphAsymmErrors *h_PDFBand_results_ratio_to_ref;
    TGraphAsymmErrors *h_AlphaS_results;
    TGraphAsymmErrors *h_AlphaS_results_ratio_to_ref;
    TGraphAsymmErrors *h_RenormalizationScale_results;
    TGraphAsymmErrors *h_RenormalizationScale_results_ratio_to_ref;
    TGraphAsymmErrors *h_FactorizationScale_results;
    TGraphAsymmErrors *h_FactorizationScale_results_ratio_to_ref;
    TGraphAsymmErrors *h_TotError_results;
    TGraphAsymmErrors *h_TotError_results_ratio_to_ref;
   // TGraphAsymmErrors *h_PDFSystBand_prenorm;
};

