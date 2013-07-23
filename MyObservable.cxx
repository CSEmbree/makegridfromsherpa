//
//   for sherpa
//
#include <iostream>
using namespace std;
#include <string>

#include "MyObservable.h"


/******************************************************************
 ** Method Implementations
 ******************************************************************/

MyGrid::MyObservable(MyData *data)
{

  debug=true;

  mydata=data;
  //myevent=event;

  incljets=mydata->isObsInclJets(); 
  observable=mydata->GetObservable(); 
  if (debug) cout<<" observable= "<<observable<<endl;

  obsptjet=false;
  obsrapjet=false;
  if (observable.Contains("PTJET") obsptjet=true;
  if (observable.Contains("RAPJET") obsrapjet=true;
  if (obsptjet) cout<<" MyObservable observable is ptjet"<<endl;
  if (obsrapjet)cout<<" MyObservable observable is rapidity jet"<<endl;

  return;
}


//  std::vector<fastjet::PseudoJet*> myseljets=myevent->GetSelectedJets(ptmin,rapmin,rapmax);


double MyObservable::GetJetObservableValue(fastjet::PseudoJet *myjet){
 cout<<" MyObservable::GetJetObservableValue observable= "<<observable<<endl;
 if (!myjet)  cout<<" MyObservable::GetJetObservableValue myjet not found "<<endl;
 double obs=0.;
 if (obsptjet) obs=myjet->pt();
 if (obsrapjet) obs=myjet->rap();
 return obs;
}
double MyObservable::GetEventObservableValue(MyEvent *myevent) {
  double obs=0.;
  cout<<" not implemented "<<endl;
  return obs;
}

bool  MyObservable::eventcuts(MyEvent *myevent){
//
// analysis cuts defined in steering
//
 bool mydebug=false;

 bool accept=false;
 bool oklept=false;
 bool okneut=false;

 int idlep=-1; int idneu=-1;

 if (!mydata) cout<<" MyObservable:cutwithlepton MyData object not found !"<<endl;

 int njets=0;
 for (int  i = 0; i < myevent->GetN(); i++)
 {

  int id=myevent->GetID(i);
  if (mydebug) std::cout<<"MyObservable::cutwithlepton";
  if (mydebug)
    printf("%d pt= %6.2f rap= %6.2f  id= %d \n",i,
         myevent->Get(i)->pt(), myevent->Get(i)->rap(),id);
   //
  double pt =(myevent->Get(i))->pt();
  double rap=fabs((myevent->Get(i))->rap());

  //bool isjet=abs(id)<10||abs(id)==21;
  bool isjet=myevent->IsJet(id);
  if (mydata->JetCut()&&isjet){
    if ((rap<mydata->GetJetMaxY()&&rap>mydata->GetJetMinY())
       && pt>mydata->GetJetPtCut()) accept=true;
   if (mydebug) if (accept) cout<<"MyObservable::cutwithlepton accept jet cut "<<endl;
   if (accept) njets++;
  };

  //bool islepton=abs(id)==11||abs(id)==13;
  bool islepton=myevent->IsChargedLepton(id);
  if (mydata->LeptonCut()&&islepton){
   idlep=i;
   if ((rap<mydata->GetLepMaxY()&&rap>mydata->GetLepMinY())
       &&pt>mydata->GetLepPtCut()) oklept=true;
   if (mydebug) if (oklept) cout<<"MyObservable::cutwithlepton accept lepton cut "<<endl;
  };

  bool isneutrino=myevent->IsNeutrino(id);
  if (mydata->NeutrinoCut()&&isneutrino){
   idneu=i;
   if (pt>mydata->GetNeuPtCut()) okneut=true;
   if (mydebug) if (okneut) cout<<"MyObservable::cutwithlepton accept neutrino cut "<<endl;
  };

// jet are filled after the lepton
  if (accept&&isjet&&oklept&&mydata->DeltaLepJetCut()){ 
   double dR=999.;
   if (idlep!=-1) dR=(myevent->Get(i))->delta_R(*myevent->Get(idlep));
   if (mydebug) cout<<i<<" dR= "<<dR<<" cut= "<<mydata->GetDeltaRLepJet()<<endl;
   if (dR<mydata->GetDeltaRLepJet()) accept=false;
   if (mydebug) if (!accept) cout<<"MyObservable::cut event cut by dR"<<endl;
  };
 };

 if (accept) {
  if (mydata->NeutrinoCut()&&!okneut) accept=false;
  if (mydata->LeptonCut()  &&!oklept) accept=false;
 }

 if (accept&&mydata->MtwCut()){
  if (mydebug) cout<<"MyObservable::cutwithlepton idlep= "<<idlep<<" idneu= "<<idneu<<endl;
  if (idlep!=-1&&idneu!=-1){
   fastjet::PseudoJet pw=*myevent->Get(idlep)+*myevent->Get(idneu);;
   double mtw=pw.mt();
   if (mydebug) cout<<"MyObservable::cutwithlepton mtw= "<<mtw<<endl;
   if (mtw<mydata->GetMtWCut()) accept=false;
  };
 };
 
 if (accept&&mydata->GetNJetCut()){
   if (mydebug) cout<<"MyObservable::cutwithlepton njets= "<<njets<<endl;
  if(njets<mydata->GetNJet()) accept=false;
  if (mydebug) if (!accept) cout<<"MyObservable::cut event rejected by Njet cut"<<endl;
 }

 if (mydebug)
  if (accept) cout<<"MyObservable::cutwithlepton event accepted !"<<endl; 
  else        cout<<"MyObservable::cutwithlepton event rejected ! "<<endl;

 return accept; 
};


