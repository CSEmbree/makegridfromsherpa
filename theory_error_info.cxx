#include "theory_error_info.h"

void GetPDF(const double& x, const double& Q, double* xf) {
  evolvePDF( x, Q, xf);        //// calls LHAPDF
  //hoppeteval_( x, Q, xf);    
  return;
}

//inline void   evolvePDF(const double & x, const double & Q, double * f) {
//  evolvepdf_(&x , &Q, f);
//}

//theory_error_calc::theory_error_calc(appl::grid my_grid)
//theory_error_calc::theory_error_calc()
//theory_error_calc::theory_error_calc(TString grid_name, int pdfi, int the_PDFBandType, int the_ErrorSizeType, double xscale)
theory_error_calc::theory_error_calc(TString grid_name, int pdfi, bool do_i_do_PDFBand, bool do_i_do_AlphaS, bool do_i_do_RS, bool do_i_do_FS, int the_ErrorSizeType, double xscale)
{
  pdf_code = pdfi;
  do_PDFBand = do_i_do_PDFBand;
  do_AlphaS = do_i_do_AlphaS;
  do_RenormalizationScale = do_i_do_RS;
  do_FactorizationScale = do_i_do_FS;
  calc_desc = "theory_errors";
  if( do_PDFBand ) calc_desc += "_PDFBand";
  if( do_AlphaS ) calc_desc += "_AlphaS";
  if( do_RenormalizationScale ) calc_desc += "_RS";
  if( do_FactorizationScale ) calc_desc += "_FS";
  if( !do_PDFBand && !do_AlphaS && !do_RenormalizationScale && !do_FactorizationScale ) {
    std::cout << "All theory uncertainties disabled. Configuration file error? Nothing to do, so exit" << std::endl;
    assert(0);
  }
 // PDFBandType = the_PDFBandType;
  ErrorSizeType = the_ErrorSizeType;

 // const std::string other_grid_name = "test";
 // const string other_str = "str";
  my_grid = new appl::grid(grid_name.Data());
  my_grid->trim();

  static const int nLoops    = 1;
  static const int nFlavours = 5;
  h_errors_PDFBand.clear();
  h_errors_AlphaS.clear();
  h_errors_RenormalizationScale.clear();
  h_errors_FactorizationScale.clear();

  std::cout << "Fill PDF errors for " << PDF_strs[pdf_code] << std::endl;
    //initPDFset(PDF_strs[pdf_code].c_str());
 // PDFsets/cteq66alphas.LHgrid
 // PDFsets/CT10as.LHgrida
 // ////// In this case initPDF(7) => +0.120, and initPDF(3) => +0.116
 // PDFsets/NNPDF23_as_0121_100.LHgrid   === > For up
 // PDFsets/NNPDF23_as_0117_100.LHgrid   === > For down
 // MSTW: something like this:
 // PDFsets/MSTW2008nlo90cl_asmz+90cl.LHgrid
 // PDFsets/MSTW2008nlo90cl_asmz-90cl.LHgrid
  double stev=1000.;
  TH1D* temp_hist;   TH1D* temp_hist_prenorm;
  std::string default_pdf_set_name = (std::string) ("PDFsets/" + PDF_strs[pdf_code] + ".LHgrid");
  if (pdf_code == e_HERAPDF15NLO ) default_pdf_set_name = (std::string) ("PDFsets/" + PDF_strs[pdf_code] + "_EIG.LHgrid");
  LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
  h_qqbar_prenorm = (TH1D*) my_grid->convolute_subproc(6, GetPDF, alphasPDF, nLoops);      ///// maybe also subprocess 5?
  h_qqbar_prenorm->SetName((TString) ("h_qqbar_prenorm_" + calc_desc));
  TH1D* h_qqbar_prenorm2 = (TH1D*) my_grid->convolute_subproc(5, GetPDF, alphasPDF, nLoops);
  h_qqbar_prenorm2->SetName((TString) ("h_qqbar_prenorm_" + calc_desc));
  h_qqbar_prenorm->Add(h_qqbar_prenorm2);
  h_qqbar_prenorm->SetLineColor(FillColorCodes[pdf_code]);  h_qqbar_prenorm->SetMarkerColor(FillColorCodes[pdf_code]);
  h_qqbar = (TH1D*) TH1NormToTot(h_qqbar_prenorm, 1. / 1000., 1000.*xscale/stev);
  h_qqbar->SetName((TString) ("h_qqbar_" + calc_desc));
  h_qqbar->SetLineColor(FillColorCodes[pdf_code]);  h_qqbar->SetMarkerColor(FillColorCodes[pdf_code]);
  h_gg_prenorm = (TH1D*) my_grid->convolute_subproc(0, GetPDF, alphasPDF, nLoops);      
  h_gg_prenorm->SetName((TString) ("h_gg_prenorm_" + calc_desc));
  h_gg_prenorm->SetLineColor(FillColorCodes[pdf_code]);  h_gg_prenorm->SetMarkerColor(FillColorCodes[pdf_code]);
  h_gg = (TH1D*) TH1NormToTot(h_gg_prenorm, 1. / 1000., 1000.*xscale/stev);
  h_gg->SetName((TString) ("h_gg_" + calc_desc));
  h_gg->SetLineColor(FillColorCodes[pdf_code]);  h_gg->SetMarkerColor(FillColorCodes[pdf_code]);
  h_tot_prenorm = (TH1D*) my_grid->convolute(GetPDF, alphasPDF, nLoops);      
  h_tot_prenorm->SetName((TString) ("h_tot_prenorm_" + calc_desc));
  h_tot_prenorm->SetLineColor(FillColorCodes[pdf_code]);  h_tot_prenorm->SetMarkerColor(FillColorCodes[pdf_code]);
  h_tot = (TH1D*) TH1NormToTot(h_tot_prenorm, 1. / 1000., 1000.*xscale/stev);
  h_tot->SetName((TString) ("h_tot_" + calc_desc));
  h_tot->SetLineColor(FillColorCodes[pdf_code]);  h_tot->SetMarkerColor(FillColorCodes[pdf_code]);
  h_gg_frac = (TH1D*) h_gg_prenorm->Clone((TString) ("h_gg_frac_" + calc_desc));
  h_gg_frac->Divide(h_tot_prenorm);
  h_qqbar_frac = (TH1D*) h_qqbar_prenorm->Clone((TString) ("h_qqbar_frac_" + calc_desc));
  h_qqbar_frac->Divide(h_tot_prenorm);
    //TH1F* h_qqbar_prenorm;
    //TH1F* h_gg_prenorm;
    //TH1F* h_qqbar;
    //TH1F* h_gg;

  if( do_RenormalizationScale ) {
    for(int renscali = 0; renscali < e_n_RenScaleVals; renscali++) {
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops, RenScale_vals[renscali], 1.);
      temp_hist_prenorm->SetName((TString) ("h_xsec_rscale_" + RenScale_strs[renscali]));
      temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
      temp_hist->SetName((TString) ("h_xsec_rscale_" + RenScale_strs[renscali] + "_norm"));
      h_errors_RenormalizationScale.push_back(temp_hist);
    }
  }
  if( do_FactorizationScale ) {
    for(int facscali = 0; facscali < e_n_FacScaleVals; facscali++) {
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops, 1., FacScale_vals[facscali]);
      temp_hist_prenorm->SetName((TString) ("h_xsec_rscale_" + FacScale_strs[facscali]));
      temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
      temp_hist->SetName((TString) ("h_xsec_rscale_" + FacScale_strs[facscali] + "_norm"));
      h_errors_FactorizationScale.push_back(temp_hist);
    }
  }

  //// Calculate PDF errors using alpha_S
  if( do_AlphaS ) {
    LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
    temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
    temp_hist_prenorm->SetName((TString) ("h_xsec_default"));
    h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
    if( pdf_code == e_CT10 ) {
      LHAPDF::initPDFSet(((std::string) ("PDFsets/CT10as.LHgrid")).c_str(), 3);  //// alphaS down
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_CT10as116_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
      initPDF(7);    /// alphaS up
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_CT10as120_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);

    } else if( pdf_code == e_MSTW2008nlo68cl ) {     /// MSTW2008nlo
      LHAPDF::initPDFSet(((std::string) ("PDFsets/MSTW2008nlo68cl_asmz-68cl.LHgrid")).c_str(), 0);   /// alphaS down
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_MSTW2008nloAsDown_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
      LHAPDF::initPDFSet(((std::string) ("PDFsets/MSTW2008nlo68cl_asmz+68cl.LHgrid")).c_str(), 0);   /// alphaS up
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_MSTW2008nloAsUp_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
    } else if( pdf_code == e_NNPDF23_nlo_as_0118 ) {   /// NNPDF23
      LHAPDF::initPDFSet(((std::string) ("PDFsets/NNPDF23_nlo_as_0116.LHgrid")).c_str(), 0);  //// alphaS down
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as117_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
      LHAPDF::initPDFSet(((std::string) ("PDFsets/NNPDF23_nlo_as_0120.LHgrid")).c_str(), 0);  //// alphaS up
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as123_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
    } else if( pdf_code == e_HERAPDF15NLO ) {
      LHAPDF::initPDFSet(((std::string) ("PDFsets/HERAPDF15NLO_ALPHAS.LHgrid")).c_str(), 9);  //// alphaS down
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_HERAPDF15NLOas1156_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
      LHAPDF::initPDFSet(((std::string) ("PDFsets/HERAPDF15NLO_ALPHAS.LHgrid")).c_str(), 11);  //// alphaS up
      temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as1196_prenorm"));
      h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
    }
    for(int alphai = 0; alphai < (int) h_errors_AlphaS_prenorm.size(); alphai++) {
      double stev=1000.;
      std::cout << "use an xscale of " << xscale << " #evs assumed: " << stev << "\n";
      temp_hist = (TH1D*) TH1NormToTot(h_errors_AlphaS_prenorm.at(alphai), 1. / 1000., 1000.*xscale/stev);
      TString the_name = h_errors_AlphaS_prenorm.at(alphai)->GetName();
      temp_hist->SetName((TString) (the_name + "_normalized"));
      h_errors_AlphaS.push_back(temp_hist);
    }
  }

  if( do_PDFBand ) {
    LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
    //// Calculate PDF errors using standard PDF error band
    std::cout << "Calc PDF errors\n";
    for(int pdferri = 0; pdferri < nPDFMembers[pdf_code]; pdferri++) {
      //enum enum_PDFBandType {e_UseAlphaS, e_UseErrorBand, e_n_PDFBands};
      std::cout << "pdferri: " << pdferri << " of " << nPDFMembers[pdf_code] << std::endl;
      if( pdf_code != e_HERAPDF15NLO ) initPDF(pdferri);
      else {   //// e_HERAPDF15NLO
        if( pdferri <= 20 ) {
          initPDF(pdferri);
        } else if( pdferri == 21 ) {
          LHAPDF::initPDFSet("PDFsets/HERAPDF15NLO_VAR.LHgrid", 0);
        } else if( pdferri > 21 ) {
          initPDF(pdferri - 21);
        }
      }
      TH1D* temp_hist_prenorm = (TH1D*) my_grid->convolute( GetPDF, alphasPDF, nLoops);
      TH1D* temp_hist;
      TString this_pdf_err_code = "";    this_pdf_err_code += pdferri;
      temp_hist_prenorm->SetName((TString) ("h_xsec_" + PDF_strs[pdf_code] + "_" + this_pdf_err_code + "_prenorm"));
      h_errors_prenorm.push_back(temp_hist_prenorm);
      double stev=1000.;
      TString the_name = temp_hist_prenorm->GetName();
      std::cout << "use an xscale of " << xscale << " #evs assumed: " << stev << "\n";
      temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
      temp_hist->SetName((TString) (the_name + "_normalized"));
      // if( pdf_code == e_CT10 && pdferri == 0 ) {
      //    std::cout << "Dealing with default CT10. Compare grid results with reference before and after normalization.\n";
      //    for(int bi = 1; bi <= hreftheory->GetNbinsX(); bi++) {
      //      float center = hreftheory->GetBinCenter(bi);
      //      cout << "bi: " << bi << ", center: " << center << ", ref: " << hreftheory->GetBinContent(bi) << ", convolution: " << temp_hist_prenorm->GetBinContent(bi) << ", after norm: " << temp_hist->GetBinContent(bi) << "\n";
      //    }
      //  }
      h_errors_PDFBand.push_back(temp_hist);
      //  TH1D *hreftheorynorm= (TH1D*) TH1NormToTot(hreftheory,1./1000.,1./stev); // xscale fb->pb yscale GeV->TeV
    }   /// pdf errors loop
  }  /// do_PDFBand
  std::cout << "Done with PDF errors loop" << std::endl;

}

void theory_error_calc::InitializeErrorGraphs()
{
  std::cout << "Beginning of InitializePDFErrorGraphs" << std::endl;
  TString this_name = (TString) ("syst_band_" + PDF_strs[pdf_code]);
  Double_t x_vals[10];   Double_t x_errs_low[10];    Double_t x_errs_high[10];   Double_t y_vals[10];   Double_t y_errs_low[10];   Double_t y_errs_high[10];

 // Double_t y_vals_prenorm[10];
  std::cout << "PDF: " << PDF_strs[pdf_code] << ", do_AlphaS: " << do_AlphaS << ", do_PDFBand: " << do_PDFBand << std::endl;
  Int_t n_bins = 0;
  std::cout << "Initialize PDFBand_results TGraphs. size of AlphaS = " << h_errors_AlphaS.size() << ", RenScale size = " << h_errors_RenormalizationScale.size() << ", FactScale size = " << h_errors_FactorizationScale.size() << ", PDFBand size = " << h_errors_PDFBand.size()   << std::endl;
  if( do_AlphaS ) n_bins = h_errors_AlphaS.at(0)->GetNbinsX();
  else if( do_RenormalizationScale ) n_bins = h_errors_RenormalizationScale.at(0)->GetNbinsX();
  else if( do_FactorizationScale ) n_bins = h_errors_FactorizationScale.at(0)->GetNbinsX();
  else if( do_PDFBand ) n_bins = h_errors_PDFBand.at(0)->GetNbinsX();
  else assert(0);
  std::cout << "Start initialization" << std::endl;
  for(int pi = 0; pi < n_bins; pi++) {
    x_vals[pi] = 0;    float x_width = 0;
    if( do_AlphaS ) {
      x_vals[pi] = h_errors_AlphaS.at(0)->GetBinCenter(pi+1);
      x_width = h_errors_AlphaS.at(0)->GetBinWidth(pi+1);
      y_vals[pi] =  h_errors_AlphaS.at(0)->GetBinContent(pi+1);
    } else if( do_PDFBand ) {
      x_vals[pi] = h_errors_PDFBand.at(0)->GetBinCenter(pi+1);
      x_width = h_errors_PDFBand.at(0)->GetBinWidth(pi+1);
      y_vals[pi] = h_errors_PDFBand.at(0)->GetBinContent(pi+1);
    } else if( do_RenormalizationScale ) {
      x_vals[pi] = h_errors_RenormalizationScale.at(0)->GetBinCenter(pi+1);
      x_width = h_errors_RenormalizationScale.at(0)->GetBinWidth(pi+1);
      y_vals[pi] = h_errors_RenormalizationScale.at(0)->GetBinContent(pi+1);
    } else if( do_FactorizationScale ) {
      x_vals[pi] = h_errors_FactorizationScale.at(0)->GetBinCenter(pi+1);
      x_width = h_errors_FactorizationScale.at(0)->GetBinWidth(pi+1);
      y_vals[pi] = h_errors_FactorizationScale.at(0)->GetBinContent(pi+1);
    }
    x_errs_low[pi] = x_width / 2.;
    x_errs_high[pi] = x_width / 2.;
    y_errs_low[pi] = 0.;
    y_errs_high[pi] = 0.;
  }
  std::cout << "Make results graphs" << std::endl;
  h_PDFBand_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
  h_PDFBand_results->SetName((TString) (this_name + "_PDFBand_results"));
  h_PDFBand_results->SetFillColor(FillColorCodes[pdf_code]);
  h_PDFBand_results->SetFillStyle(FillStyleCodes[pdf_code]);
  h_AlphaS_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
  h_AlphaS_results->SetName((TString) (this_name + "_AlphaS_results"));
  h_AlphaS_results->SetFillColor(FillColorCodes[pdf_code]);
  h_AlphaS_results->SetFillStyle(FillStyleCodes[pdf_code]);
  h_RenormalizationScale_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
  h_RenormalizationScale_results->SetName((TString) (this_name + "_RenormalizationScale_results"));
  h_RenormalizationScale_results->SetFillColor(FillColorCodes[pdf_code]);
  h_RenormalizationScale_results->SetFillStyle(FillStyleCodes[pdf_code]);
  h_FactorizationScale_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
  h_FactorizationScale_results->SetName((TString) (this_name + "_FactorizationScale_results"));
  h_FactorizationScale_results->SetFillColor(FillColorCodes[pdf_code]);
  h_FactorizationScale_results->SetFillStyle(FillStyleCodes[pdf_code]);
  h_TotError_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
  h_TotError_results->SetName((TString) (this_name + "_TotError_results"));
  h_TotError_results->SetFillColor(FillColorCodes[pdf_code]);
  h_TotError_results->SetFillStyle(FillStyleCodes[pdf_code]);
 // h_PDFBand_results_prenorm = new TGraphAsymmErrors(n_bins, x_vals, y_vals_prenorm, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
//  h_PDFBand_results_prenorm->SetName((TString) (this_name + "_PDFBand_results_prenorm"));
//  h_PDFBand_results_prenorm->SetFillColor(FillColorCodes[pdf_code]);
//  h_PDFBand_results_prenorm->SetFillStyle(FillStyleCodes[pdf_code]);
  std::cout << "End of InitializePDFErrorGraphs" << std::endl; 
}

void theory_error_calc::CalcSystErrors()
{
  std::cout << "Beginning of CalcPDFSystErrors() function for PDF = " << PDF_strs[pdf_code] <<  std::endl;
  if( do_PDFBand ) CalcPDFBandErrors();
  if( do_AlphaS ) CalcAlphaSErrors();
  if( do_RenormalizationScale )  CalcRenormalizationScaleErrors();
  if( do_FactorizationScale )  CalcFactorizationScaleErrors();
  CalcTotErrors();
  std::cout << "End of CalcPDFSystErrors() function for pdf_code = " << pdf_code << std::endl;
}

void theory_error_calc::CalcTotErrors()
{
  for(int pi = 0; pi < h_TotError_results->GetN(); pi++) {
    Double_t PDFBand_err_high = 0.;    Double_t PDFBand_err_low = 0.;
    Double_t AlphaS_err_high = 0.;    Double_t AlphaS_err_low = 0.;
    Double_t RenormalizationScale_err_high = 0.;    Double_t RenormalizationScale_err_low = 0.;
    Double_t FactorizationScale_err_high = 0.;    Double_t FactorizationScale_err_low = 0.;

    Double_t x_val=-999;   Double_t y_val=-999;
    if( do_PDFBand ) {
      PDFBand_err_high = h_PDFBand_results->GetErrorYhigh(pi);
      PDFBand_err_low = h_PDFBand_results->GetErrorYlow(pi);
      h_PDFBand_results->GetPoint(pi, x_val, y_val);
    } 
    if( do_AlphaS ) {
      AlphaS_err_high = h_AlphaS_results->GetErrorYhigh(pi);
      AlphaS_err_low = h_AlphaS_results->GetErrorYlow(pi);
      h_AlphaS_results->GetPoint(pi, x_val, y_val);
    }
    if( do_RenormalizationScale ) {
      RenormalizationScale_err_high = h_RenormalizationScale_results->GetErrorYhigh(pi);
      RenormalizationScale_err_low = h_RenormalizationScale_results->GetErrorYlow(pi);
      h_RenormalizationScale_results->GetPoint(pi, x_val, y_val);
    }
    if( do_FactorizationScale ) {
      FactorizationScale_err_high = h_FactorizationScale_results->GetErrorYhigh(pi);
      FactorizationScale_err_low = h_FactorizationScale_results->GetErrorYlow(pi);
      h_FactorizationScale_results->GetPoint(pi, x_val, y_val);
    }
    Double_t Tot_err_high = TMath::Sqrt(pow(PDFBand_err_high, 2.) + pow(AlphaS_err_high, 2.) + pow(RenormalizationScale_err_high, 2.) + pow(FactorizationScale_err_high, 2.));
    Double_t Tot_err_low = TMath::Sqrt(pow(PDFBand_err_low, 2.) + pow(AlphaS_err_low, 2.) + pow(RenormalizationScale_err_low, 2.) + pow(FactorizationScale_err_low, 2.));

    h_TotError_results->SetPoint(pi, x_val, y_val);
    h_TotError_results->SetPointEYhigh(pi, Tot_err_high);
    h_TotError_results->SetPointEYlow(pi, Tot_err_low);
  }  /// loop over points

}

void theory_error_calc::CalcRenormalizationScaleErrors()
{
  std::cout << "Beginning of CalcPDFRenormalizationScaleErrors for pdf_code = " << pdf_code << std::endl;
  assert(h_errors_RenormalizationScale.size() >= 3);
  for(int bi = 1; bi <= h_errors_RenormalizationScale.at(0)->GetNbinsX(); bi++) {
    double this_default_val = h_errors_RenormalizationScale.at(e_RenScale1p0)->GetBinContent(bi);
    double this_err_down = h_errors_RenormalizationScale.at(e_RenScale0p5)->GetBinContent(bi);
    double this_err_up = h_errors_RenormalizationScale.at(e_RenScale2p0)->GetBinContent(bi);  
    std::cout << "bi = " << bi << ", default val = " << this_default_val << " +" << this_err_up   << " -" << this_err_down << "\n";
 //   double mid_val = 0.5*(this_err_up + this_err_down);
    double error = 0.5*fabs(this_err_up-this_err_down);
    if( ErrorSizeType == e_90Percent ) error *= 1.645;
    Double_t init_x_val;   Double_t init_y_val;
    h_RenormalizationScale_results->GetPoint(bi-1, init_x_val, init_y_val);
    h_RenormalizationScale_results->SetPoint(bi-1, init_x_val, this_default_val);
    h_RenormalizationScale_results->SetPointEYhigh(bi-1, error);
    h_RenormalizationScale_results->SetPointEYlow(bi-1, error);

     
  } /// bi
  std::cout << "End of CalcThisRenormalizationScalePDFSyst for PDF = " << PDF_strs[pdf_code] << std::endl;

}

void theory_error_calc::CalcFactorizationScaleErrors()
{
  std::cout << "Beginning of CalcPDFFactorizationScaleErrors for pdf_code = " << pdf_code << std::endl;
  assert(h_errors_FactorizationScale.size() >= 3);
  for(int bi = 1; bi <= h_errors_FactorizationScale.at(0)->GetNbinsX(); bi++) {
    double this_default_val = h_errors_FactorizationScale.at(e_FacScale1p0)->GetBinContent(bi);
    double this_err_down = h_errors_FactorizationScale.at(e_FacScale0p5)->GetBinContent(bi);
    double this_err_up = h_errors_FactorizationScale.at(e_FacScale2p0)->GetBinContent(bi);  
    std::cout << "bi = " << bi << ", default val = " << this_default_val << " +" << this_err_up   << " -" << this_err_down << "\n";
 //   double mid_val = 0.5*(this_err_up + this_err_down);
    double error = 0.5*fabs(this_err_up-this_err_down);
    if( ErrorSizeType == e_90Percent ) error *= 1.645;
    Double_t init_x_val;   Double_t init_y_val;
    h_FactorizationScale_results->GetPoint(bi-1, init_x_val, init_y_val);
    h_FactorizationScale_results->SetPoint(bi-1, init_x_val, this_default_val);
    h_FactorizationScale_results->SetPointEYhigh(bi-1, error);
    h_FactorizationScale_results->SetPointEYlow(bi-1, error);

     
  } /// bi
  std::cout << "End of CalcThisFactorizationScalePDFSyst for PDF = " << PDF_strs[pdf_code] << std::endl;

}


void theory_error_calc::CalcAlphaSErrors()
{
  std::cout << "Beginning of CalcPDFAlphaSErrors for pdf_code = " << pdf_code << std::endl;
  assert(h_errors_AlphaS.size() == 3);
  for(int bi = 1; bi <= h_errors_AlphaS.at(0)->GetNbinsX(); bi++) {
    double this_default_val = h_errors_AlphaS.at(0)->GetBinContent(bi);
    double this_err_down = h_errors_AlphaS.at(1)->GetBinContent(bi);
    double this_err_up = h_errors_AlphaS.at(2)->GetBinContent(bi);  
    std::cout << "bi = " << bi << ", default val = " << this_default_val << " +" << this_err_up   << " -" << this_err_down << "\n";
 //   double mid_val = 0.5*(this_err_up + this_err_down);
    double error = 0.5*fabs(this_err_up-this_err_down);
    if( ErrorSizeType == e_90Percent ) error *= 1.645;
    Double_t init_x_val;   Double_t init_y_val;
    h_AlphaS_results->GetPoint(bi-1, init_x_val, init_y_val);
    h_AlphaS_results->SetPoint(bi-1, init_x_val, this_default_val);
    h_AlphaS_results->SetPointEYhigh(bi-1, error);
    h_AlphaS_results->SetPointEYlow(bi-1, error);

     
  } /// bi
  std::cout << "End of CalcThisAlphaSPDFSyst for PDF = " << PDF_strs[pdf_code] << std::endl;

}

void theory_error_calc::CalcPDFBandErrors() 
{
  std::cout << "Beginning of CalcPDFBandErrors for pdf_code = " << pdf_code << std::endl;
  for(int bi = 1; bi <= h_errors_PDFBand.at(0)->GetNbinsX(); bi++) {
    double this_err_up = 0.;    double this_err_down = 0.;
    double average = 0.;    //// needed for NNPDF
    if( pdf_code == e_NNPDF23_nlo_as_0118 ) {
      for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri++) {
        average += h_errors_PDFBand.at(pdferri)->GetBinContent(bi);
      }
      average /= h_errors_PDFBand.size()-1;
      for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri++)  {
        this_err_up += pow(h_errors_PDFBand.at(pdferri)->GetBinContent(bi)-average, 2.);
      }
      this_err_up = TMath::Sqrt(this_err_up / (h_errors_PDFBand.size()-1));
      this_err_down = this_err_up;
      if( ErrorSizeType == e_90Percent ) {
        this_err_up *= 1.645;
        this_err_down *= 1.645;
      }

    }  /// NNPDF23
    else if( pdf_code == e_CT10  ) {   /// CT10
      for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size()-1; pdferri += 2) {
        this_err_up += pow( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - h_errors_PDFBand.at(pdferri+1)->GetBinContent(bi), 2.);
      }
      this_err_up = 0.5*TMath::Sqrt(this_err_up);
      if( ErrorSizeType == e_OneSigma ) this_err_up /= 1.645;
      this_err_down = this_err_up;
    } else if( pdf_code == e_MSTW2008nlo68cl ) {     /// MSTW2008nlo
      double central_val = h_errors_PDFBand.at(0)->GetBinContent(bi);
      for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri ++) {
        double mod_val = h_errors_PDFBand.at(pdferri)->GetBinContent(bi);
        if( mod_val > central_val ) this_err_up += pow(mod_val-central_val, 2.);
        else this_err_down += pow(central_val - mod_val, 2.);
      }
      this_err_down = TMath::Sqrt(this_err_down);
      this_err_up = TMath::Sqrt(this_err_up);
      if( ErrorSizeType == e_90Percent ) {
        this_err_up *= 1.645;
        this_err_down *= 1.645;
      }

    } else if( pdf_code == e_HERAPDF15NLO ) {
      double central_val = h_errors_PDFBand.at(0)->GetBinContent(bi);
      for(int pdferri = 1; pdferri < 20; pdferri += 2) {   //// experimental errors
        this_err_up += pow( 0.5*(h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - h_errors_PDFBand.at(pdferri+1)->GetBinContent(bi)), 2.);
      }
      this_err_down = this_err_up;
      for(int pdferri = 21; pdferri < 29; pdferri++) {   /// model errors
        if( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) > central_val ) {
          this_err_up  += pow( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - central_val, 2.);
          this_err_down += 0.;
        }
        else  {
          this_err_up  += 0.;
          this_err_down += pow( central_val - h_errors_PDFBand.at(pdferri)->GetBinContent(bi), 2.);
        }
      }
      double extreme_pos_diff = 0.;     double extreme_neg_diff = 0.;
      for(int pdferri = 29; pdferri < 33; pdferri++) {     //// parameterization errors
        double diff_central = h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - central_val;
        if( diff_central > 0 && diff_central > extreme_pos_diff ) extreme_pos_diff = diff_central;
        if( diff_central < 0 && diff_central < extreme_neg_diff ) extreme_neg_diff = diff_central;
      }
      if( extreme_pos_diff > 0. ) this_err_up += pow(extreme_pos_diff, 2.);
      if( extreme_neg_diff < 0. ) this_err_down += pow(extreme_neg_diff, 2.);
      this_err_up = TMath::Sqrt(this_err_up);
      this_err_down = TMath::Sqrt(this_err_down);
      if( ErrorSizeType == e_90Percent ) {
        this_err_up *= 1.645;
        this_err_down *= 1.645;
      }
  ///// End of HERAPDF15NLO
    } else assert(0);

    ///// To test the plotting
  //  this_err_up *= 10.;
  //  this_err_down *= 10.;
  //  ///////////////////////
  
    h_PDFBand_results->SetPointEYhigh(bi-1, this_err_up);
    h_PDFBand_results->SetPointEYlow(bi-1, this_err_down);
    //     h_PDFBand_results->SetPointEXhigh(bi-1, 0.);
    //     h_PDFBand_results->SetPointEXlow(bi-1, 0.);
    Double_t x_val;   Double_t y_val;
    h_PDFBand_results->GetPoint(bi-1, x_val, y_val);

  }  /// loop over bins
  std::cout << "End of CalcThisPDFSyst for PDF = " << PDF_strs[pdf_code] << std::endl;
}

TH1 *theory_error_calc::TH1NormToTot (TH1 *h1, Double_t yscale, Double_t xscale)
{
  // divide histogram h1 scale by total number and scale bin width x
  // assumes that histogram is divided by bin-width
  std::cout << "renormalize with xscale = " << xscale << ", yscale = " << yscale << std::endl;
  const char *name="**TH1NormToTot";  
  const bool debug=false;
  if (!h1) cout<<" TH1NormTotot h1 not found "<<endl; 
  TH1D* h1new=(TH1D*) h1->Clone(h1->GetName());

  Double_t x, y, ey;
  Double_t sigtot=0.;
  std::cout << "yscale=  " << yscale << std::endl;
  for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
    y=h1->GetBinContent(i)*yscale;
    x=h1->GetBinWidth(i);
    sigtot+=y*x;
  } 
  std::cout << "sigtot: " << sigtot << std::endl;

  if (debug) cout<<name<<" sigtot= "<<sigtot<<endl;

  for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
    x =h1->GetBinWidth(i);
    y =h1->GetBinContent(i)*yscale*x;
    ey=h1->GetBinError(i)  *yscale*x;
    x =h1->GetBinWidth(i)  *xscale;
    if (x!=0) h1new->SetBinContent(i,y/x);
    else      h1new->SetBinContent(i,0.);
    std::cout << "bin " << i << ", center = " << x << ", y-val = " << y << ", fill with y/x = " << h1new->GetBinContent(i) << "\n";

    if (x!=0) h1new->SetBinError(i,ey/x);
    else      h1new->SetBinError(i,0.);
  } 

  if (sigtot!=0.)
    h1new->Scale(1./sigtot);
  std::cout << "Scale by 1 / sigtot to get the final result " << std::endl;
 
  //if (debug) {
    float integral_so_far = 0.;
    for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
      y =h1new->GetBinContent(i);
      x =h1new->GetBinWidth(i);
      integral_so_far += y*x;
      cout << name<<" "<<i<<" bincenter= "<<h1new->GetBinCenter(i)<<" Binw = " << x << " y= " << y << ", integral so far: " << integral_so_far << endl;  
    }
 // }

  // h1new->Print("all");
  return h1new;
}


TGraphAsymmErrors* theory_error_calc::TH1TOTGraphAsymm(TH1 *h1){
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

TGraphAsymmErrors* theory_error_calc::myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr) {
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
    Double_t el=0.; Double_t eh=0.;

    if (noerr==2) {dy2l=0.; dy2h=0.;}
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
  if (matchcount>1) {cout<<name<<" too many x-points matched ! "<<endl; exit (1);}
 }  
 return g3;
}

void theory_error_calc::GetRatioToTH1(TH1D* href)
{
  TGraphAsymmErrors* tgraph_href = TH1TOTGraphAsymm(href);
  TString ratio_to_ref_name = (TString) h_PDFBand_results->GetName(); 
  ratio_to_ref_name += "_ratio_to_ref";
  h_PDFBand_results_ratio_to_ref = myTGraphErrorsDivide(h_PDFBand_results, tgraph_href);
  h_PDFBand_results_ratio_to_ref->SetName(ratio_to_ref_name);
  h_PDFBand_results_ratio_to_ref->SetFillColor(FillColorCodes[pdf_code]);
  h_PDFBand_results_ratio_to_ref->SetFillStyle(FillStyleCodes[pdf_code]);
  ratio_to_ref_name = (TString) h_AlphaS_results->GetName(); 
  ratio_to_ref_name += "_ratio_to_ref";
  h_AlphaS_results_ratio_to_ref = myTGraphErrorsDivide(h_AlphaS_results, tgraph_href);
  h_AlphaS_results_ratio_to_ref->SetName(ratio_to_ref_name);
  h_AlphaS_results_ratio_to_ref->SetFillColor(FillColorCodes[pdf_code]);
  h_AlphaS_results_ratio_to_ref->SetFillStyle(FillStyleCodes[pdf_code]);
  ratio_to_ref_name = (TString) h_RenormalizationScale_results->GetName(); 
  ratio_to_ref_name += "_ratio_to_ref";
  h_RenormalizationScale_results_ratio_to_ref = myTGraphErrorsDivide(h_RenormalizationScale_results, tgraph_href);
  h_RenormalizationScale_results_ratio_to_ref->SetName(ratio_to_ref_name);
  h_RenormalizationScale_results_ratio_to_ref->SetFillColor(FillColorCodes[pdf_code]);
  h_RenormalizationScale_results_ratio_to_ref->SetFillStyle(FillStyleCodes[pdf_code]);
  ratio_to_ref_name = (TString) h_FactorizationScale_results->GetName(); 
  ratio_to_ref_name += "_ratio_to_ref";
  h_FactorizationScale_results_ratio_to_ref = myTGraphErrorsDivide(h_FactorizationScale_results, tgraph_href);
  h_FactorizationScale_results_ratio_to_ref->SetName(ratio_to_ref_name);
  h_FactorizationScale_results_ratio_to_ref->SetFillColor(FillColorCodes[pdf_code]);
  h_FactorizationScale_results_ratio_to_ref->SetFillStyle(FillStyleCodes[pdf_code]);
  ratio_to_ref_name = (TString) h_TotError_results->GetName(); 
  ratio_to_ref_name += "_ratio_to_ref";
  h_TotError_results_ratio_to_ref = myTGraphErrorsDivide(h_TotError_results, tgraph_href);
  h_TotError_results_ratio_to_ref->SetName(ratio_to_ref_name);
  h_TotError_results_ratio_to_ref->SetFillColor(FillColorCodes[pdf_code]);
  h_TotError_results_ratio_to_ref->SetFillStyle(FillStyleCodes[pdf_code]);

}

void theory_error_calc::CalcChi2(TGraphAsymmErrors *g_theory, TGraphAsymmErrors *g_data, TMatrixT<double> data_cov_matrix, double &chi2)
{
  // https://cds.cern.ch/record/1470588/files/ATL-COM-PHYS-2012-1137.pdf
  std::cout << "CalcChi2" << std::endl;
  assert( g_theory->GetN() == g_data->GetN() );
  assert( g_theory->GetN() == data_cov_matrix.GetNcols() );
  assert( g_theory->GetN() == data_cov_matrix.GetNrows() );

  /// Fill in the theory covariance matrix
  TMatrixT<double> theory_cov_matrix(g_theory->GetN(), g_theory->GetN());
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
      if( pi != pi2 ) theory_cov_matrix(pi, pi2) = 0;
      if( pi == pi2 ) {
        double theory_uncertainty = 0.5*(g_theory->GetErrorYhigh(pi) + g_theory->GetErrorYlow(pi));
        theory_cov_matrix(pi, pi2) = theory_uncertainty*theory_uncertainty;
      }
    }
  }
  std::cout << "At beginning, dump contents of data cov matrix: " << std::endl;
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
      std::cout << data_cov_matrix(pi,pi2) << "\t";
    }
    std::cout << "\n";
  }


  TMatrixT<double> tot_cov_matrix = theory_cov_matrix + data_cov_matrix;
  std::cout << "After adding theory, dump contents of cov matrix: " << std::endl;
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
      std::cout << tot_cov_matrix(pi,pi2) << "\t";
    }
    std::cout << "\n";
  }
  TMatrixT<double> invertex_cov_matrix = tot_cov_matrix.Invert();    //// Now it includes the theory errors in the diagonal elements ...
  std::cout << "After inversion, dump contents of cov matrix: " << std::endl;
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
      std::cout << invertex_cov_matrix(pi,pi2) << "\t";
    }
    std::cout << "\n";
  }

  ///// Loop over bins and determine data-theory matrices
  TMatrixT<double> row_data_minus_theory(1, g_theory->GetN());
  TMatrixT<double> col_data_minus_theory(g_theory->GetN(), 1);
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    Double_t data_val;    Double_t theory_val;
    Double_t x_val;
    g_theory->GetPoint(pi, x_val, theory_val);
    g_data->GetPoint(pi, x_val, data_val);
    row_data_minus_theory(0,pi) = data_val - theory_val;
    col_data_minus_theory(pi,0) = data_val - theory_val;
    std::cout << "At " << x_val << ", data = " << data_val << ", theory = " << theory_val << ", content = " << row_data_minus_theory(0,pi) << "\n";
  }  /// pi

  TMatrixT<double> cov_times_col = invertex_cov_matrix*col_data_minus_theory;
  assert( cov_times_col.GetNrows() == g_theory->GetN());
  assert( cov_times_col.GetNcols() == 1);
  std::cout << "After first multiplication matrix is:\n";
  for(int pi = 0; pi < g_theory->GetN(); pi++) {
    std::cout << cov_times_col(pi, 0) << "\t";
  }
  std::cout << "\n";
  TMatrixT<double> result = row_data_minus_theory*cov_times_col;
  assert( result.GetNrows() == 1);
  assert( result.GetNcols() == 1);

  chi2 = result(0,0);
  std::cout << "Final chi2 = " << chi2 << "\n";
  
}
