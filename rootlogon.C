{
cout << endl << "Welcome to the ATLAS rootlogon.C" << endl;
//
// based on a style file from BaBar
//

//..BABAR style from RooLogon.C in workdir
TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

// use plain black on white colors
 Int_t icol=0;
 atlasStyle->SetFrameBorderMode(icol);
 atlasStyle->SetCanvasBorderMode(icol);
 atlasStyle->SetPadBorderMode(icol);
 atlasStyle->SetPadColor(icol);
 atlasStyle->SetCanvasColor(icol);
 atlasStyle->SetStatColor(icol);
 //atlasStyle->SetFillColor(icol);

//
// Set the size of the default canvas
// 600x500 looks almost square
//atlasgStyle->SetCanvasDefH(500);
//atlasgStyle->SetCanvasDefW(600);
//atlasgStyle->SetCanvasDefX(10);
//atlasgStyle->SetCanvasDefY(10);

// set the paper & margin sizes
atlasStyle->SetPaperSize(20,26);
//atlasStyle->SetPaperSize(30,39);
atlasStyle->SetPadTopMargin(0.05);
atlasStyle->SetPadRightMargin(0.05);
atlasStyle->SetPadBottomMargin(0.16);
atlasStyle->SetPadLeftMargin(0.12);
//
// shift Canvas up to make space for labels
//atlasStyle->SetPadTopMargin(0.08);
//atlasStyle->SetPadBottomMargin(0.12);
//atlasStyle->SetPadLeftMargin(0.14);
//atlasStyle->SetPadRightMargin(0.1);


// use large fonts
//Int_t font=72;
Int_t font=42; // 42
Double_t tsize=0.05;
atlasStyle->SetTextFont(font);
atlasStyle->SetLabelFont(font,"x");
atlasStyle->SetTitleFont(font,"x");
atlasStyle->SetLabelFont(font,"y");
atlasStyle->SetTitleFont(font,"y");
atlasStyle->SetLabelFont(font,"z");
atlasStyle->SetTitleFont(font,"z");

atlasStyle->SetTextSize(tsize);
atlasStyle->SetLabelSize(tsize,"x");
atlasStyle->SetTitleSize(tsize,"x");
atlasStyle->SetLabelSize(tsize,"y");
atlasStyle->SetTitleSize(tsize,"y");
atlasStyle->SetLabelSize(tsize,"z");
atlasStyle->SetTitleSize(tsize,"z");
//
//atlasStyle->SetLabelOffset(.01, "x");

//use bold lines and markers
atlasStyle->SetMarkerStyle(20);
//y
atlasStyle->SetMarkerSize(1.3);
atlasStyle->SetHistLineWidth(2);

atlasStyle->SetFuncWidth(2);
atlasStyle->SetLineStyleString(7,"[12 12]"); // postscript dashes
//atlasStyle->SetLineScalePS(2);  // 3 is default
//atlasStyle->SetScreenFactor(1); // default is 1
//
//get rid of X error bars and y error bar caps
//atlasStyle->SetErrorX(0.001);
//atlasStyle->SetErrorMsize(5.);// ??
//set the errors in x to be as broad as the binsize
//atlasStyle->SetErrorX(1.);

//do not display any of the standard histogram decorations
atlasStyle->SetOptTitle(0);
//atlasStyle->SetOptStat(1111);
atlasStyle->SetOptStat(0);
//atlasStyle->SetOptFit(1111);
//atlasStyle->SetOptFit(0);
//atlasStyle->SetOptStat(1110);// no histogram title
//atlasStyle->SetStatFormat("6.4f");
//atlasStyle->SetFitFormat("6.4f");
//atlasStyle->SetStatStyle(0); // hollow
//atlasgStyle->SetStatStyle(1001); // filled
//atlasatlasStyle->SetStatBorderSize(0);
//Style->SetStatW(0.25);
//atlasStyle->SetStatH(0.125);
//atlasStyle->SetStatX(0.90);
//atlasStyle->SetStatY(0.90);

// put tick marks on top and RHS of plots
atlasStyle->SetPadTickX(1);
atlasStyle->SetPadTickY(1);
//TGaxis::SetMaxDigits(2); // restrict the number of digits in labels


// no supressed zeroes!
//gStyle->SetHistMinimumZero(kTRUE);



// Legend without shadow
atlasStyle->SetLegendBorderSize(1);

gROOT->SetStyle("Plain");

//atlasgStyle->SetPadTickX(1);
//atlasgStyle->SetPadTickY(1);
//atlasStyle->SetTickLength(0.02,"xyz");

//
gSystem->Load("libPostscript.so");

// gStyle->SetDateX(x);

}
