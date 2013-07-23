#include "MyEvent.h"
#include <iostream>
using namespace std;

/******************************************************************
 ** Method Implementations
 ******************************************************************/

MyEvent::MyEvent()
{
  debug=false;
  sqrts=0.; 
}

MyEvent::~MyEvent()
{
  //cout<<" MyEvent:: destructor "<<endl;
}

void MyEvent::ClearEvent()
{
 event.clear();
}

void MyEvent::Print()
{

 std::cout<<"\n MyEvent: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;
 if (debug) std::cout<<" debug flag on "<<std::endl;
 std::cout<<" s= "<<sqrts<<" GeV "<<" eventid= "<<eventid<<std::endl;

 std::cout<<" Xsection= "<<xsection<<" weight= "<<weight<<std::endl; 

 std::cout<<" x1= "<<x1<<" x2= "<<x2<<" q2= "<<q2<<std::endl;

  std::cout<<" Event record N= "<<event.size()<<std::endl;
  for (int  i = 0; i <   event.size(); i++)
  {
    printf("%d px= %6.2f py= %6.2f pz= %6.2f E= %6.2f id= %d \n",i,
	   event[i].px(), event[i].py(),event[i].pz(),event[i].e(),event[i].user_index());
  }
  std::cout<<" "<<std::endl;
  std::cout<<" MyEvent: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<std::endl;

 }

void MyEvent::Print2()
{

  std::cout<<" Event record N= "<<event.size()<<std::endl;
  for (int  i = 0; i <   event.size(); i++)
  {
    printf("%d pt= %6.2f rap= %6.2f id= %d \n",i,
	   event[i].pt(), event[i].rap(),event[i].user_index());
  }
  std::cout<<" "<<std::endl;
  std::cout<<" MyEvent: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<std::endl;

 }

void MyEvent::push_back(double px, double py, double pz, double E, int id)
{
 // 
  fastjet::PseudoJet *myjet2 = new fastjet::PseudoJet(px,py,pz,E);
  fastjet::PseudoJet myjet3 = *myjet2; 
  //myjet->reset(px,py,pz,E);
  //myjet->set_user_index (id);  
   myjet3.set_user_index (id);  
  //cout<<" id= "<<id<<" push_back pt= "<<myjet->pt()<<endl;
  event.push_back(myjet3);
  //event.push_back(fastjet::PseudoJet(px,py,pz,E));
  delete myjet2;
  return;
}

void MyEvent::push_back(double px, double py, double pz, double E)
{
 //
  int id=-1;
  //set_user_index (id);  
  //event.push_back(fastjet::PseudoJet(px,py,pz,E));
  //myjet->reset(px,py,pz,E);
  this->push_back(px,py,pz,E,id);
  return;
}

void MyEvent::push_back(std::vector<fastjet::PseudoJet> myjets)
{
 //
 if (debug) std::cout<<" MyEvent store output jets N= "<<myjets.size()<<std::endl;
 for (int  i = 0; i < myjets.size(); i++){
  //cout<<i<<" pt= "<<myjets[i].pt()<<endl;
  event.push_back(myjets[i]);
 }
 return;
}

void MyEvent::push_back(fastjet::PseudoJet myjet2)
{
 //
 //if (debug) std::cout<<" MyEvent store output jets N= "<<myjets.size()<<std::endl;
  //fastjet::PseudoJet *mypjets=&myjet2;
  //event.push_back(mypjets);
    event.push_back(myjet2);
   //event.push_back(&myjets[i]);
 return;
}


 int MyEvent::GetLeadingJetID(){

   //cout<<" GetLeadingJetID starting "<<endl;
  double ptmax=-99; int id=-1;
  for (int  i = 0; i <   event.size(); i++)
  {
    //printf("%d pt= %6.2f eta= %6.2f E= %6.2f id= %d \n",i,
    //	   event[i]->pt(), event[i]->rapidity(),event[i]->e(),event[i]->user_index());

    if (abs(event[i].user_index())>10 && !event[i].user_index()==21) continue;
    if (ptmax<event[i].pt()){
     ptmax=event[i].pt();
     id=i;
     //cout<<" ptmax= "<<ptmax<<" id= "<<id<<endl;
    }
  };
  //cout<<" leading jet id= "<<id<<endl;
  if (id==-1) {
   cout<<" MyEvent::GetLeadingJetID Leading jet not found ! "<<endl;
   this-> Print();
  }
  //cout<<" GetLeadingJetID finished "<<endl;
  return id;
 };

 fastjet::PseudoJet *MyEvent::Get(int index){
  fastjet::PseudoJet *tmp=0;
  if (index>event.size()) return tmp;
  return &event[index];
 };
 
double MyEvent::GetInvariantMass12(){
 double obs=0.;
 int id1=-1, id2=-1;
 double ptmax1=-999., ptmax2=-999.;
 for (int  i = 0; i < event.size(); i++) {
  if (abs(event[i].user_index())>10 && !event[i].user_index()==21) continue;
  if (event[i].pt()>ptmax2) {ptmax2=ptmax1; ptmax1=event[i].pt(); id2=id1; id1=i;} 

  if (debug)
   cout<<" MyEvent::GetInvariantMass12 i= "<<i
       <<" ptmax1= "<<ptmax1<<" ptmax2= "<<ptmax2
       <<" event.pt= "<<event[i].pt()
       <<" id1= "<<id1<<" id2= "<<id2
       <<endl;
 }

 if (debug){
  cout<<" MyEvent::GetInvariantMass12 id1= "<<id1<<" id2= "<<id2<<endl;
  cout<<" MyEvent::GetInvariantMass12 ptmax1= "<<ptmax1<<" ptmax2= "<<ptmax2<<endl;
 }

 if (id2==-1) {
  cout<<" MyEvent::GetInvariantMass12 something is wrong id2= "<<id2<<endl;
  return obs;
 }
 fastjet::PseudoJet pw=event[id1]+event[id2];
 obs=pw.m();
 if (debug) cout<<" MyEvent::GetInvariantMass12 obs= "<<obs<<endl;

 //delete pw;

 return obs;
}


std::vector<fastjet::PseudoJet> MyEvent::GetSelectedJets(double ptmin, double rapmin, double rapmax){
//
// select jets:
// ptjet>ptmin and etamin<=etajet<etamax
//
 bool mydebug=false;
 myseljets.clear();
 double en, px, py, pz, pt, rap;
 
 for (int  i = 0; i < this->GetN(); i++) {
  int id=this->GetID(i);
  if (mydebug) 
    printf(" MyEvent:GetSelectedJets %d pt= %6.2f rap= %6.2f  id= %d \n",i,
         this->Get(i)->pt(), this->Get(i)->rap(),id);
  if (!this->IsJet(id)) continue;

  en =(this->Get(i))->e();
  px =(this->Get(i))->px();
  py =(this->Get(i))->py();
  pz =(this->Get(i))->pz();
  pt =(this->Get(i))->pt();
  rap=fabs((this->Get(i))->rap());

  if ((rap<=rapmax&&rap>rapmin)&& pt>ptmin){
   //cout<<" need to implement something "<<endl;
   //myjet->reset(px,py,pz,en);
   //myjet->set_user_index (id);
   fastjet::PseudoJet myjet=*this->Get(i);
   if (debug) cout<<" MyEvent::GetSelectedJets selected jets pt= "<<
        myjet.pt()<<" rap= "<<myjet.rap()<<endl;
   myseljets.push_back(myjet);
  }
 }
 return myseljets;
};



