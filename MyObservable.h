#ifndef MyGrid_H
#define MyObservable_H

#include <stdlib.h> // exit()
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <map>

using std::string;
#include <sys/time.h> 

#include "root.h"
#include "MyEvent.h"
#include "MyData.h"

class MyObservable {
 private:

  bool debug; 

  //  std::vector<MyData*> mydata;// information about data from steering file  
  MyData mydata;
  //MyEvent myevent;

  TString observable;
  bool incljets;
 
  bool obsptjet;
  bool obsrapjet;

  double ptmin;
  double rapmin;
  double rapmax;

 public:
  MyObservable(MyData *data);

  double GetEventObservableValue(MyEvent *myevent);
  double GetJetObservableValue(fastjet::PseudoJet myjet);
  bool eventcuts(MyEvent *myevent);

  MyData *GetMyData(int i){ return mydata;};
  void SetMyData(MyData *mydat ){ mydata=mydat; return;};

  MyData *GetMyEvent(int i){ return myevent;};
  void SetMyEvent(MyData *myev ){ myevent=myev; return;};

  string *GetObservable(){ return observable;};
  void SetObservable(string obs){ =myev; return;};

  

 double GetObsMax(){
  if (!mydata)  cout<<" GetObsMax grid not found= mydata"<<endl;
  double obsmax=mydata[igrid]->GetMaxX(); 
  //if (debug) cout<<" obsmax= "<<obsmax<<endl;
  return obsmax;
 }

};


#endif
