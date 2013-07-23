void Normalise(TH1D* h1, double yscale, double xscale=1., bool normtot=false) 
{ 
 std::cout<<" HuHu "<<std::endl;

 //xscale = 1.;
 Double_t x, y, ey;
 Double_t sigtot=0.;

 if (!h1) std::cout<<" Normalise histo not found ! "<<std::endl;
 else {
  std::cout<<" Normalise histo to normalise: "<<std::endl;
  h1->Print("all");
 }

 for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
  y=h1->GetBinContent(i)*yscale;
  x=h1->GetBinWidth(i);
  sigtot+=y*x;
 } 
 //sigtot *= yscale;
 std::cout<<" Normalise sigtot= "<<sigtot<<std::endl;

 for (Int_t i=0; i<=h1->GetNbinsX(); i++) {
  x =h1->GetBinWidth(i);
  y =h1->GetBinContent(i)*yscale*x;
  ey=h1->GetBinError(i)  *yscale*x;
  x =h1->GetBinWidth(i)  *xscale;
  if (x!=0) h1->SetBinContent(i,y/x);
  else      h1->SetBinContent(i,0.);
 // std::std::cout << "BinWidth: " << h1->GetBinWidth(i)  << ", xscale: " << xscale << ", x: " << x << ", bin content: " << h1->GetBinContent(i) << "\n";

 if (x!=0) h1->SetBinError(i,ey/x);
  else     h1->SetBinError(i,0.);
 } 
if (normtot){
  h1->Scale(1. / sigtot);
 }
  return;
}
