#ifndef VarCodes 
#define VarCodes
enum enum_PDFs {e_CT10, e_MSTW2008nlo68cl, e_NNPDF23_nlo_as_0118, e_HERAPDF15NLO, e_n_PDFs};
static string PDF_strs[e_n_PDFs] = {"CT10", "MSTW2008nlo68cl", "NNPDF23_nlo_as_0118", "HERAPDF15NLO"};
static TString PDF_file_name_strs[e_n_PDFs] = {"CT10", "MSTW2008nlo", "NNPDF23nlo", "HERAPDF15NLO_EIG"};
static TString detailed_PDF_strs[e_n_PDFs] = {"CT10", "MSTW2008nlo", "NNPDF23nlo #alpha_{S} 118", "HERAPDF 1.5"};
static int nPDFMembers[e_n_PDFs] = {53, 41, 101, 33};
static int FillStyleCodes[e_n_PDFs] = {3005, 3004, 3002, 3021};
static int FillColorCodes[e_n_PDFs] = {kGreen+10, kRed+2, kBlue+2, kGreen+2};

enum enum_PDFBandType {e_UseAlphaS, e_UseErrorBand, e_n_PDFBands};
static std::string PDFBandType_strs[e_n_PDFBands] = {"UseAlphaS", "UseErrorBand"};
enum enum_PDFErrorSizeType {e_OneSigma, e_90Percent, e_n_PDFErrorSizeTypes};
static std::string PDFErrorSize_strs[e_n_PDFErrorSizeTypes] = {"OneSigma", "90Percent"};

enum enum_RenScales {e_RenScale0p5, e_RenScale1p0, e_RenScale2p0, e_n_RenScaleVals};
static double RenScale_vals[e_n_RenScaleVals] = {0.5, 1., 2.};
static TString RenScale_strs[e_n_RenScaleVals] = {"RenScale0p5", "RenScale1p0", "RenScale2p0"};
enum enum_FacScales {e_FacScale0p5, e_FacScale1p0, e_FacScale2p0, e_n_FacScaleVals};
static double FacScale_vals[e_n_FacScaleVals] = {0.5, 1., 2.};
static TString FacScale_strs[e_n_FacScaleVals] = {"FacScale0p5", "FacScale1p0", "FacScale2p0"};

enum enum_ErrorTypes {e_PDFBandCode, e_AlphaSCode, e_RenormalizationScaleCode, e_FactorizationScaleCode, e_TotErrorCode, e_n_ErrorCodes};
static TString ErrorType_strs[e_n_ErrorCodes] = {"PDFBandCode", "AlphaSCode", "RenormalizationScaleCode", "FactorizationScaleCode", "TotErrorCode"};

#endif
