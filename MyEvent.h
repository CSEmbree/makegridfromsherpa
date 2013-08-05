
#ifndef MyEvent_H
#define MyEvent_H

#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>

#include "root.h"
//#include "fjClustering.h"
#include "fastjet/PseudoJet.hh"

class MyEvent {
 private:

  bool debug; 
  //  Int_t njets;
  double sqrts;
  double xsection; // cross section
  double weight;   // weight without PDF
  int eventid;     // event ID to identify correlated events
  double x1;       // Bjorken x of parton 1
  double x2;       // Bjorken x of parton 2
  double q2;       // factorisation scale
  int iorder;      // order of alphas
  int itype;      // order of alphas
 
  std::vector<fastjet::PseudoJet> event;
  //fastjet::PseudoJet *myjet;
  std::vector<fastjet::PseudoJet> myseljets;

 public:

  MyEvent();
  virtual  ~MyEvent();

  void Print();
  void Print2();
  void ClearEvent();
  void PrintEvent();
  void push_back(double px, double py, double pz, double E);
  void push_back(double px, double py, double pz, double E,int id);
  void push_back(std::vector<fastjet::PseudoJet> myjets);
  //void push_back(std::vector<fastjet::PseudoJet> myjets);
  void push_back(fastjet::PseudoJet myjets);
  //void push_back(fjClustering *myjets);
  int GetN(){ return event.size();};
  void SetCMS(double s) {sqrts=s; return;};
  double GetCMS() {return sqrts;};

  double GetXSection() {return xsection;};
  void   SetXSection(double xsec) {xsection=xsec; return;};
 
  double GetWeight() {return weight;};
  void SetWeight(double myweight) {weight=myweight; return;};

  int GetOrder() {return iorder;};
  void SetOrder(int myorder) {iorder=myorder; return;};

  bool GetBorn   () {bool flag=false; if (itype==73) flag=true; return flag;};
  bool GetVirtual() {bool flag=false; if (itype==66) flag=true; return flag;};
  bool GetReal   () {bool flag=false; if (itype==82) flag=true; return flag;};

  int GetType() {return itype;};
  void SetType(int mytype) {itype=mytype; return;};

  int GetEventId() {return eventid;};
  void SetEventId(int id) {eventid=id; return;};

  double GetX1(){ return x1;};
  double GetX2(){ return x2;};
  void SetX1(double myx){x1=myx; return;};
  void SetX2(double myx){x2=myx; return;};

  double GetQ2(){ return q2;};
  void SetQ2(double myq2){ q2=myq2; return;};
  fastjet::PseudoJet* Get(int id);

  int GetID(int index){return event[index].user_index();};
  int GetLeadingJetID();

  bool IsJet(int id){ return abs(id)<10||abs(id)==0;};
  bool IsChargedLepton(int id){ return abs(id)==11||abs(id)==13;};
  bool IsNeutrino(int id){ return abs(id)==12||abs(id)==14;};

  int GetLeptonCharge(){
  //
  int icharge=0;
  for (int  i = 0; i <   event.size(); i++){
    if (event[i].user_index()== 11||event[i].user_index()== 13) icharge=-1;
    if (event[i].user_index()==-11||event[i].user_index()==-13) icharge= 1;
  };
  if (icharge==0) std::cout<<" MyEvent: lepton not found "<<std::endl;
  return icharge;
 }

 int GetNumberSelJets(){ return myseljets.size();};

 void PrintSelJets(){ 
  std::cout<<" MyGrid::PrintSelJets Number of jets= "<<this->GetNumberSelJets()<<std::endl;
  for (int  i = 0; i <this->GetNumberSelJets(); i++){
    printf("%d px= %6.2f py= %6.2f pz= %6.2f E= %6.2f id= %d \n",i,
         myseljets[i].px(), myseljets[i].py(),myseljets[i].pz(),
         myseljets[i].e(),  myseljets[i].user_index());
  }
  return;
 };

 void PrintSelJets2(){ 
  std::cout<<" MyGrid::PrintSelJets2 Number of Jets= "<<this->GetNumberSelJets()<<std::endl;
  for (int  i = 0; i <this->GetNumberSelJets(); i++){
    printf("%d pt= %6.2f rap= %6.2f id= %d \n",i,
         myseljets[i].pt(), myseljets[i].rap(),myseljets[i].user_index());
  }
  return;
 };

 std::vector<fastjet::PseudoJet> GetSelectedJets(double ptmin, double rapmin, double rapmax);

 double GetInvariantMass12();
 
 double GetLeadingJet(){
  double obs=0.;
  int id=this->GetLeadingJetID();
  if (this->Get(id))
   obs=(this->Get(id))->pt();
  else std::cout<<"GetLeadingJet No leading jet detected "<<std::endl;
  //std::cout<<" MyGrid::GetLeadingJet obs= "<<obs<<std::endl;
  if (obs!=obs) std::cout<<" MyEvent::GetLeadingJet obs is NAN"<<std::endl;
  return obs;
 }
};




#endif
