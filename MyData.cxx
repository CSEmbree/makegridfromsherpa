
#include "MyData.h"
#include <iostream>
#include <cmath>
using namespace std;

/******************************************************************
 ** Method Implementations
 ******************************************************************/

MyData::MyData()
{
    debug=false;

    normtot=false;
    unitsxTeV=false;
    unitsypb=false;

    blogy=false;
    blogx=false;
    bliny=true;
    blinx=true;
    
    njetcut=false;

    dividebybinwith=false;

    sqrts=0.;
    ymin=0.;
    ymax=0.;
    datavector=0;
    datavectorstat=0;
    datavectorsyst=0;
    datavectortoterr=0;

    //xtdatavector=0;
    //xtdatavectorstat=0;
    //xtdatavectorsyst=0;
    scale=1.;
    msize=0.03;
    markerstyle=0;
    markercolor=0;
    markersize=1.0;
    if (debug)  cout<<" msize= "<<msize<<endl;
    //cout<<" datavector= "<<datavector<<endl;

    jetalgoR=0;

    observable="";
    //TString reaction;
    jetcut=false;
    lepcut=false;
    deltaRLepJet=false;
    neucut=false;
    cutMtw=false;
    drljcut=false;

    Mtwcut=0.;
    ptjetcut=0.;
    ptlepcut=0.;
    ptneucut=0.;
    yjetmin=-999.;
    yjetmax=999.;
    ylepmin=-999.;
    ylepmax=999.;
    deltaRLepJet=0.;
    setframexmin=false;
    framexmin=0.;

    incljets=false;

//  PDFSetCodes_vec.clear();
//  PDFBandType_vec.clear();
//  PDFErrorSize_vec.clear();
}

void MyData::ReadCorrelations(const char fname[200], double myscale) {
    /*
     * Read data correlation file provided by fname
     */
    debug = true;

    if (debug) std::cout << " MyData::ReadCorrelations  Read correlations fname= "<<fname<< std::endl;
    cov_matrix = new TMatrixT<double>(datavectortoterr->GetN(), datavectortoterr->GetN());



    //remove path characters and file extentions from file path to extract JUST file name
    corrfilename.Clear();
    string pathfilename = this->GetFileNameFromPath(string(fname));
    corrfilename=TString(pathfilename.c_str());

    /* //older less generic filepath name extractor
    corrfilename.Clear();
    corrfilename=TString(fname);
    corrfilename.ReplaceAll(TString("data/"),TString(""));
    corrfilename.ReplaceAll(TString(".txt"),TString("")); */


    // std::vector<std::vector<double>> corr_matrix;
    ifstream infile(fname, ios::in);
    if(!infile) { // Check open
        cerr << " MyData::ReadCorrelations Can't open " << fname <<". Make sure this path exists!\n";
        infile.close();
        exit (1);
    } else {
        if (debug) cout <<" MyData::ReadCorrelations: read data correlations file: " << fname << endl;
    }

    //int iline=0;
    //int nsyst=1;
    //int ntab=0;   int nline=0;
    //char ch;
    char line[1024];
    int row_num; //<--never given a value before being used??

    while (infile.good()) {
        //if (!infile.good()) break;
        infile.getline(line,sizeof(line),'\n');
        if( line[0] == '%' ) continue;
        std::string cpp_line(line);
        std::vector<std::string> split_line;
        split_line.clear();
        split_string(cpp_line, split_line, " ");
        if( debug ) std::cout << "Read in correlation matrix, line " << line << std::endl;
        std::vector<double> my_corr_row;
        my_corr_row.clear();
        if( cpp_line.size() > 3 ) {
            if( debug ) std::cout << "Dump numerical contents (should be identical):\n";
            for(int ci = 2; ci < (int) split_line.size(); ci++) {
                int col_num = ci-2;
                double corr_val = atof( split_line.at(ci).c_str() );
                double row_data_err = 0.5*(datavectortoterr->GetErrorYhigh(row_num) + datavectortoterr->GetErrorYlow(row_num));
                double col_data_err = 0.5*(datavectortoterr->GetErrorYhigh(col_num) + datavectortoterr->GetErrorYlow(col_num));
//      datavectortoterr
                double cov_val = corr_val*row_data_err*col_data_err;
                //     my_corr_row.push_back(this_val);
                //     TMatrixT<double>
                (*cov_matrix)(row_num, col_num) = cov_val;
                if( debug )  std::cout << corr_val << " => " << (*cov_matrix)(row_num, col_num)  << "\t||\t";
            }    /// loop over entries
            if( debug ) std::cout << "\n";
            row_num++;
            //   corr_matrix.push_back(my_corr_row);
        }  /// This is a row of the correlation matrix. Read it in.
    }
}

void MyData::ReadData(const char fname[200], double myscale) {

    if (debug) std::cout << " MyData::ReadData" << std::endl;
    bool errorinpercent=false;

    /* //older less generic filepath name extractor
    filename.Clear();
    filename=TString(fname);
    filename.ReplaceAll(TString("data/"),TString(""));
    filename.ReplaceAll(TString(".txt"),TString(""));*/

    filename.Clear();
    string pathfilename = this->GetFileNameFromPath(string(fname));
    filename=TString(pathfilename.c_str());


    datavector = new  TGraphAsymmErrors();
    datavector->SetName(fname);
    if (!datavector) cout<<" MyData::ReadData can not create data vector "<<endl;
//datavector->Print();

    datavectorsyst = new  TGraphAsymmErrors();
    TString fnamesyst(fname);
    fnamesyst+="syst";
    datavectorsyst->SetName(fnamesyst);
    if (!datavectorsyst) cout<<" MyData::ReadData can not create data vector "<<endl;
//datavectorsyst->Print();

    if (debug) std::cout << " MyData::ReadData Creating datavectorstat" << std::endl;
    datavectorstat = new  TGraphAsymmErrors();
    TString fnamestat(fname);
    fnamestat+="stat";
    datavectorstat->SetName(fnamestat);
    if (!datavectorstat) cout<<" MyData::ReadData can not create data vector "<<endl;
    if (debug) std::cout << " MyData::ReadData Created it" << std::endl;

    if (debug) std::cout << " MyData::ReadData Creating datavectortot" << std::endl;
    datavectortoterr = new  TGraphAsymmErrors();
    TString fnametot(fname);
    fnamestat+="tot";
    datavectortoterr->SetName(fnametot);
    if (!datavectortoterr) cout<<" MyData::ReadData can not create data vector "<<endl;
    if (debug) std::cout << " MyData::ReadData Created it" << std::endl;

//datavectorstat->Print();

//filename=new TString(&fname);

    ifstream infile(fname, ios::in);
    if(!infile) { // Check open
        cerr << "Can't open " << fname <<"\n";
        infile.close();
        exit (1);
    } else {
        if (debug) cout <<" MyData::ReadData: read data file: " << fname << endl;
    }
    
    cout <<"TEST: MyData::ReadData: read data file: " << fname << endl;

    int iline=0;
    int nsyst=1;

//int ntab=0;   int nline=0;
//char ch;
    char line[1024];
//void split_string(string str, vector<string>& split_results, string delimiters)
// while (1) {
// //if (debug) cout << " good: " << infile.good() << " eof: " << infile.eof() << endl;
// if (!infile.good()) break;
// infile.get(ch);

// from http://www.gidnetwork.com/b-63.html
//  if ( ch == '\n') {ntab=0; nline++;};
//   float xv; //only count, if not a number
//   int r = scanf("%e", &xv);
////   if (!r) {
//   if (r!=4) {
//    cout<<" ch= "<<ch<<endl;
//    cout<<" ntab= "<<ntab<<endl;
//    if (ch == '\t') ntab++;
//   }
//  }
//  cout<<" ntab= "<<ntab<<" nline= "<<nline<<endl;
//  return;
//
// below take sscanf 1 by 1 to read the systeamtic errors
//
// nt  getCharacter()
// {
//     int ch;         // define the character
//     do              // loop until a good character is read
//     {
//         ch = getchar();             // read a character
//     } while ( (ch == 0x20)  ||      // check for SPACE
//              ((ch >= 0x09)  &&      // check for the
//               (ch <= 0x0D))         //    other whitespace
//             );
//     return ch;
// }
// int  getCharacter()
// {
//     int ch;         // define the character
//     do              // loop until a good character is read
//     {
//         ch = getchar();         // read a character
//     } while ((ch <= 0x20) ||    // check above SPACE
//              (ch >= 0x7F)       // check below the last
//             );                  //    ASCII char
//     return ch;
// }

    miny=9.e20;
    maxy=-9.e20;
    minx=9.e20;
    maxx=-9.e20;
    char text[TEXT_SIZE];
    char name[TEXT_SIZE];
    float myymin, myymax; //used by YRAP, YJET, YLEP
    float mys; //used by SQRTS, PTJET, PTLEP, PTNEU, MTW, DRLJ

    //load in all valid MyData options
    string DataOptionsFileName = "DataOptions.txt"; //<--could be made readable from steering file
    OptionHandler *mydataoptions= new OptionHandler(DataOptionsFileName);
    //mydataoptions->generateResultFile();
    char lineFirstWord[1024]; //same length as 'line' string buffer

    //read and set all options and data
    while (infile.good()) {
        //if (debug) cout << " MyData::ReadData good: " << infile.good() << " eof: " << infile.eof() << endl;
        //if (!infile.good()) break;
        infile.getline(line,sizeof(line),'\n');
        std::string cpp_line(line);
// std::vector<std::string> colon_split_cpp_line;  colon_split_cpp_line.clear();
// split_string(cpp_line, colon_split_cpp_line, ":");
// std::cout << "Colon-split of line " << cpp_line << " has size " << colon_split_cpp_line.size() << "\n";
//

        if (debug) cout<< " MyData::ReadSteering line= "<< line << "\n";

        sscanf(line," %s",lineFirstWord);
        std::string cpp_lineFirstWord(lineFirstWord);

        if (debug) cout<< "line= "<< line << "\n";
        if(line[0] != '%') { //ignore comments
            /*
            if(line[0] != 'Y' && line[0] != 'S' && line[0] != 'E' && line[0] != 'N' && line[0] != 'J' && line[0] != 'Y'
                   && line[0] != 'D' && line[0] != 'P' && line[0] != 'M' && line[0] != 'O' && line[0] != 'R'
                   && line[0] != 't' && line[0] != 'f' && line[0] != 'u' && line[0] != 'n'
             )
             {*/

            if(mydataoptions->isKnownOption(cpp_lineFirstWord)==false) {
                if(strlen(line) != 0) {
                    double xm=0., xl=0., xh=0., y=0., dypstat=0., dymstat=0., dypsyst=0., dymsyst=0., dymtoterr=0., dyptoterr=0.;
                    if (debug) std::cout << "nsyst: " << nsyst << std::endl;
                    if (nsyst!=1) {
                        //here must be an easier way to do that
                        //cout<<" nsyst= "<<nsyst<<endl;
                        char myf[100];
                        char my0[100]="";
                        char my1[100]="%e %e %e %e %e";
                        char my2[100]="%e";
                        for (int i=0; i<nsyst*2; i++) sprintf(my0,"%s %s",my0,my2);
                        sprintf(myf,"%s %s",my1,my0);
                        //cout<<" myf= "<<myf<<endl;
                        float vsyst[6];
                        if (nsyst>100) cout<<" MyData::ReadData too many systematic uncertainties"<<endl;
                        //??for (int i=0; i<nsyst*2; i++) sscanf(line,"%e",&vsyst[i]);
                        // quick hack for the moment
                        if (nsyst==3) {
                            sscanf(line,myf,&xl, &xh, &y, &dypstat, &dymstat, &vsyst[0], &vsyst[1], &vsyst[2], &vsyst[3], &vsyst[4], &vsyst[5]);
                            //float x1,x2,x3,x4,x5,x6;
                            //sscanf(line,myf,&xl, &xh, &y, &dypstat, &dymstat, &x1, &x2, &x3, &x4, &x5, &x6);
//      cout<<" x1= "<<x1<<" x2= "<<x2<<" x3= "<<x3<<" x4= "<<x4<<" x5= "<<x5<<" x6= "<<x6<<endl;
                            //for (int i=0; i<nsyst*2; i++) cout<<" vsyst= "<<vsyst[i]<<endl;
                        }
                    } else {
                        //OLD IMPLIMENTATION
                        //sscanf(line,"%lf %lf %lf %le %le %le %le %le",&xm, &xl, &xh, &y, &dypstat, &dymstat, &dypsyst, &dymsyst);


                        //NEW IMPLIMENTATION
                        //
                        //
                        std::string dataLine = cpp_line;
                        std::stringstream lineStream(dataLine);
                        std::string cell;
                        int inputCount=0;
                        double error;

                        std::vector<double> errors;

                        //std::cout<<"TEST: line: "<<dataLine<<std::endl;

                        while(std::getline(lineStream,cell,' '))
                        {
                            if(cell.compare("")!=0)
                            {
                                //std::cout<<"TEST: cell: "<<cell<<", for index: "<<inputCount<<std::endl;

                                switch(inputCount) {
                                case 0:
                                    sscanf(cell.c_str(), "%lf", &xm);
                                    break;
                                case 1:
                                    sscanf(cell.c_str(), "%lf", &xl);
                                    break;
                                case 2:
                                    sscanf(cell.c_str(), "%lf", &xh);
                                    break;
                                case 3:
                                    sscanf(cell.c_str(), "%lf", &y);
                                    break;
                                default:
                                    //std::cout<<"TEST: default met: "<<cell<<", for index: "<<inputCount<<std::endl;

                                    if(inputCount>3) {
                                        //std::cout<<"TEST: appropreate input count."<<std::endl;

                                        error=0;
                                        sscanf(cell.c_str(), "%lf", &error);
                                        errors.push_back(error);
                                    }
                                    else {
                                        //std::cout<<" MyData::ReadData: ERROR: invalid input encountered."<<std::endl;
                                        exit(0); //TEST
                                    }
                                }
                                inputCount++;

                                //std::cout<<"TEST: errors.size(): "<<errors.size()<<std::endl;
                            }
                        }

                        //HARD coded temporarily to keep the following code the same as with the old implimentation:
                        // simply converting the indexes of the erros in the vector to the appropreate variables
                        //std::cout<<"TEST: errors.size(): "<<errors.size()<<std::endl;

                        for(int i=0; i<errors.size(); i++) {
                            switch(i) {
                            case 0:
                                dypstat = errors.at(i);
                                break;
                            case 1:
                                dymstat = errors.at(i);
                                break;
                            case 2:
                                dypsyst = errors.at(i);
                                break;
                            case 3:
                                dymsyst = errors.at(i);
                                break;
                            default:
                                std::cout<<" MyData::ReadData: More errors were provided then are currently accounted for!"<<std::endl;
                                exit(0); //TEST
                            }
                            //std::cout<<"TEST: Assigned: "<<errors.at(i)<<", for index: "<<i<<std::endl;
                        }
                        /*
                        std::cout<<"TEST: dypstat: "<<dypstat
                                <<", dymstat: "<<dymstat
                                <<", dypsyst: "<<dypsyst
                                <<", dymsyst: "<<dymsyst<<std::endl;
                        */
                        //exit(0); //TEST
                    }




                    if (debug) printf(" MyData::ReadData xm= %f xl= %f xr= %f  y = %e dypstat = %e  dymstat= %e dypsyst= %e dymsyst= %e \n",
                                          xm,xl, xh, y, dypstat, dymstat, dypsyst, dymsyst);
                    dyptoterr = TMath::Sqrt(pow(dypstat, 2.) + pow(dypsyst, 2.));
                    dymtoterr = TMath::Sqrt(pow(dymstat, 2.) + pow(dymsyst, 2.));

                    if (errorinpercent) {
                        //cout<<"before y= "<<y<<" dypstat= "<<dypstat<<" dymstat= "<<dymstat<<" dymsyst= "<<dymsyst<<" dypsyst= "<<dypsyst<<endl;
                        dypstat*=y/100.;
                        dymstat*=y/100.;
                        dypsyst*=y/100.;
                        dymsyst*=y/100.;
                        dyptoterr*=y/100.;
                        dymtoterr*=y/100.;
                        //cout<<"after dypstat= "<<dypstat<<" dymstat= "<<dymstat<<" dymsyst= "<<dymsyst<<" dypstat= "<<dypstat<<endl;
                    }

                    y*=myscale;
                    dypstat*=myscale;
                    dymstat*=myscale;
                    dypsyst*=myscale;
                    dymsyst*=myscale;
                    dyptoterr*=myscale;
                    dymtoterr*=myscale;

                    if (miny>y) miny=y;
                    if (maxy<y) maxy=y;
                    if (minx>xl) minx=xl;
                    if (maxx<xh) maxx=xh;

                    float x= (xh+xl)/2.;
                    float b= (xh-xl)/2.;

                    if (debug) cout<<" MyData::ReadData iline = "<<iline<<" x= "<<x<<" y= "<<y<<", dymsyst: " << dymsyst << ", dypsyst: " << dypsyst << ", dymstat: " << dymstat << ", dypstat: " << dypstat << endl;
                    // y/=(ymax-ymin)*2;
                    //     datavector->SetPoint(iline ,x, y);
                    datavector->SetPoint(iline ,xm, y);
                    datavector->SetPointEXlow(iline ,xm-x+b);
                    datavector->SetPointEXhigh(iline,x+b-xm);

                    datavectorstat->SetPoint(iline ,xm, y);
                    datavectorstat->SetPointEXlow(iline ,xm-x+b);
                    datavectorstat->SetPointEXhigh(iline,x+b-xm);
                    datavectorstat->SetPointEYlow (iline,dymstat);
                    datavectorstat->SetPointEYhigh(iline,dypstat);


                    datavectorsyst->SetPoint(iline ,xm, y);
                    datavectorsyst->SetPointEXlow(iline ,xm-x+b);
                    datavectorsyst->SetPointEXhigh(iline,x+b-xm);
                    datavectorsyst->SetPointEYlow (iline,dymsyst);
                    datavectorsyst->SetPointEYhigh(iline,dypsyst);

                    datavectortoterr->SetPoint(iline ,xm, y);
                    datavectortoterr->SetPointEXlow(iline ,xm-x+b);
                    datavectortoterr->SetPointEXhigh(iline,x+b-xm);
                    datavectortoterr->SetPointEYlow (iline,dymtoterr);
                    datavectortoterr->SetPointEYhigh(iline,dyptoterr);

                    float dym=sqrt(dymstat*dymstat+dymsyst*dymsyst);
                    float dyp=sqrt(dypstat*dypstat+dypsyst*dypsyst);
                    datavector->SetPointEYlow (iline,dym);
                    datavector->SetPointEYhigh(iline,dyp);
                    iline++;
                }
            } else if (strstr(line,"unitxTeV")!=0) {
                unitsxTeV=true;
            } else if (strstr(line,"normtot")!=0) {
                normtot=true;
            } else if (strstr(line,"unitypb")!=0) {
                unitsypb=true;
            } else if (strstr(line,"framelogy")!=0) {
                blogy=true;
                bliny=false;
            } else if (strstr(line,"frameliny")!=0) {
                blogy=false;
                bliny=true;
            } else if (strstr(line,"framelogx")!=0) {
                blogx=true;
                blinx=false;
                //cout<<" set logx "<<endl;
            } else if (strstr(line,"DIVIDEBYBINWIDTH")!=0) {
                dividebybinwith=true;
            } else if (strstr(line,"framelinx")!=0) {
                blogx=false;
                blinx=true;
            } else if (strstr(line,"framexmin")!=0) {
                float myframexmin;
                sscanf(line," %s %f ",text, &myframexmin);
                setframexmin=true;
                framexmin=myframexmin;
            } else if (strstr(line,"YRAP")!=0) {
                sscanf(line," %s %f %f",text, &myymin,&myymax);
                ymin=myymin;
                ymax=myymax;
                if (debug) printf(" MyData::ReadData ymin= %f  ymax = %f  \n",myymin,myymax);
            } else if (strstr(line,"YJET")!=0) {
                jetcut=true;
                sscanf(line," %s %f %f",text, &myymin,&myymax);
                yjetmin=myymin;
                yjetmax=myymax;
                if (debug) printf(" MyData::ReadData ymin= %f  ymax = %f  \n",myymin,myymax);
            } else if (strstr(line,"YLEP")!=0) {
                lepcut=true;
                sscanf(line," %s %f %f",text, &myymin,&myymax);
                ylepmin=myymin;
                ylepmax=myymax;
                if (debug) printf(" MyData::ReadData ymin= %f  ymax = %f  \n",myymin,myymax);
            } else if (strstr(line,"SQRTS")!=0) {
                sscanf(line," %s %f ",text, &mys);
                sqrts=mys;
                if (debug) printf(" MyData::ReadData s= %f   \n",mys);
            } else if (strstr(line,"NSYST")!=0) {
                sscanf(line," %s %d ",text, &nsyst);
                if (debug) printf(" MyData::ReadData Nsyst= %d   \n",nsyst);
            } else if (strstr(line,"ERRORINPERCENT")!=0) {
                errorinpercent=true;
                if (debug) cout<<fname<<" MyData::ReadData errors given in percent "<<endl;
            } else if (strstr(line,"YEAR")!=0) {
                sscanf(line," %s %d ",text, &year);
            } else if (strstr(line,"JETALGO")!=0) {
                char jetname[100];
                sscanf(line," %s %s ",text, jetname);
                jetalgo=jetname;
            } else if (strstr(line,"EXP")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                experiment=TString(name);
            } else if (strstr(line,"REACTION")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                //sscanf(line," %s %s ",text, name);
                reaction=TString(name);
            } else if (strstr(line,"OBS")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                observable=TString(name);
                if (debug) cout<<" MyData::ReadData observable= "<<observable.Data()<<endl;
                if (observable.Contains("INCLJETS")) incljets=true;
                if (incljets) cout<<" MyData::ReadData variable for inclusive jets "<<endl;
                //double *bins = new double[nbins];
            } else if (strstr(line,"titlex")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                titlex=TString(name);
            } else if (strstr(line,"titley")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                titley=TString(name);
            } else if (strstr(line,"JETALGR")!=0) {
                sscanf(line," %s %d ",text, &jetalgoR);
            } else if (strstr(line,"NJETS")!=0) {
                njetcut=true;
                njets=0;
                sscanf(line," %s %d ",text, &njets);
            } else if (strstr(line,"PTJET")!=0) {
                jetcut=true;
                sscanf(line," %s %f ",text, &mys);
                ptjetcut=mys;
            } else if (strstr(line,"PTLEP")!=0) {
                lepcut=true;
                sscanf(line," %s %f ",text, &mys);
                ptlepcut=mys;
            } else if (strstr(line,"PTNEU")!=0) {
                neucut=true;
                sscanf(line," %s %f ",text, &mys);
                ptneucut=mys;
            } else if (strstr(line,"MTW")!=0) {
                cutMtw=true;
                sscanf(line," %s %f ",text, &mys);
                Mtwcut=mys;
            } else if (strstr(line,"DRLJ")!=0) {
                drljcut=true;
                sscanf(line," %s %f ",text, &mys);
                deltaRLepJet=mys;
            }
            //reset variables for next iteration
            memset(text, '\0', TEXT_SIZE);
            memset(name, '\0', TEXT_SIZE);
        }
    }

    return;
}



void MyData::Print() {

    cout<<"\n MyData::Print>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n "<<endl;
    if (debug) cout<<" debug flag on "<<endl;
    cout<<" reaction= "<<reaction.Data()<<endl;
    cout<<" observable= "<<observable.Data()<<endl;
    cout<<" filename= "<<filename.Data()<<endl;
//cout<<" titlex= "<<titlex.Data()<<endl;
//cout<<" titley= "<<titley.Data()<<endl;
    cout<<" s= "<<sqrts<<" GeV "<<endl;
    cout<<" year= "<<year<<endl;

//cout<<" msize= "<<msize<<" markerstyle= "<<markerstyle<<" markercolor= "<<markercolor<<endl;

    if (blogy) cout<<" axis y log-scale"<<endl;
    if (bliny) cout<<" axis y lin-scale"<<endl;
    if (blogx) cout<<" axis x log-scale"<<endl;
    if (blinx) cout<<" axis x lin-scale"<<endl;
    if (ymin!=0 || ymax!=0)
        cout<<" ymin= "<<ymin<<" ymax= "<<ymax<<endl;

    cout<<"\n Cuts: "<<endl;
    cout<<" Jet "<<jetalgo.Data()<<" R = "<< float(jetalgoR)/10.<<endl;
    if (njetcut) cout<<"\t Njets=> "<<njets<<endl;
    if (jetcut) cout<<"\t Ptjet> "<<ptjetcut<<" GeV "<<yjetmin<<" <yjet< "<<yjetmax<<endl;
    if (lepcut) cout<<"\t Ptlep> "<<ptlepcut<<" GeV "<<ylepmin<<" <ylep< "<<ylepmax<<endl;
    if (neucut) cout<<"\t Ptneu> "<<ptneucut<<" GeV "<<" MTw> " <<Mtwcut<<" GeV"<<endl;
    if (cutMtw) cout<<"\t Mtw> "<<Mtwcut<<" GeV "<<endl;
    if (drljcut)cout<<"\t dRlepJet> "<<deltaRLepJet<<" GeV "<<endl;


//
//cout<<" datavector= "<<datavector<<endl;
    if (!datavector) {
        cout<<"\n datavector not found "<<endl;
    } else {
        cout<<" datavector= "<<datavector->GetName()<<endl;
        datavector->Print();
    }

    /*
    if (!datavectorsyst){
     cout<<" datavectorsyst not found "<<endl;
    } else {
     cout<<" datavectorsyst= "<<datavectorsyst->GetName()<<endl;
     datavectorsyst->Print();
    }
    if (!datavectorstat){
     cout<<" datavectorstat not found "<<endl;
    } else {
     cout<<" datavectorstat= "<<datavectorstat->GetName()<<endl;
     datavectorstat->Print();
    }
    */

    cout<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<<endl;

}

/*
void MyData::PrintXt() {

 cout<<" MyData::PrintXt>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<<endl;
 if (debug) cout<<" debug flag on "<<endl;
 cout<<" s= "<<sqrts<<" ymin= "<<ymin<<" ymax= "<<ymax<<endl;
 cout<<jetalgo.Data()<<" R = "<< float(jetalgoR)/10.<<" year= "<<year<<endl;
 cout<<" msize= "<<msize<<" markerstyle= "<<markerstyle<<" markercolor= "<<markercolor<<endl;
 //cout<<" xtdatavector= "<<xtdatavector<<endl;
 if (!xtdatavector){
  cout<<" xtdatavector not found "<<endl;
 } else {
  cout<<" xtdatavector= "<<datavector->GetName()<<endl;
  xtdatavector->Print();
 }
 cout<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<<endl;

}

void MyData::DrawXt(char name[100], float x, float y){

  if (debug) cout<<" DrawXt: name= "<<name<<" x= "<<x<<" y= "<<y<<endl;

  this->TransformPtToXt();

  this->DrawLegend(name,x,y);

  if (debug)
   this->PrintXt();

   this->DrawDataXt();

  return;

}

void MyData::TransformPtToXt(){
//
// transform dsigma/dpt to dsigma/dxt
//
//

 xtdatavector    = new TGraphAsymmErrors();
 xtdatavectorstat= new TGraphAsymmErrors();
 xtdatavectorsyst= new TGraphAsymmErrors();

 TString fname("xt"); fname+=datavector->GetName();
 xtdatavector->SetName(fname);
 TString fnamestat= datavector->GetName(); fnamestat+="stat";
 xtdatavectorstat->SetName(fnamestat);
 TString fnamesyst= datavector->GetName(); fnamesyst+="stat";
 xtdatavectorsyst->SetName(fnamesyst);

 if (!datavector) {cout<<" Transform: datavector not found !"<<endl;}

 double * X1 = datavector->GetX();
 double * Y1 = datavector->GetY();
 double * EXhigh = datavector->GetEXhigh();
 double * EXlow =  datavector->GetEXlow();

 double * EYhigh = datavector->GetEYhigh();
 double * EYlow =  datavector->GetEYlow();

 if (!datavectorsyst) {cout<<" Transform: datavectorsyst not found !"<<endl;}
 double * EYhighsyst = datavectorsyst->GetEYhigh();
 double * EYlowsyst  = datavectorsyst->GetEYlow();

 if (!datavectorstat) {cout<<" Transform: datavectorstat not found !"<<endl;}
 double * EYhighstat = datavectorstat->GetEYhigh();
 double * EYlowstat  = datavectorstat->GetEYlow();


 double ymean=(ymax-ymin)/2.;
 //double fac=2.*exp(+ymean)/sqrts;
 //double fac=2.*exp(-ymin)/sqrts;
 double  fac=2./sqrts;

 if (debug)
  cout<<" Transform: ymean= "<<ymean<<" sqrts= "<<sqrts<<" fac= "<<fac
      <<" number of points= "<<datavector->GetN()<<endl;

 for (int i=0; i<datavector->GetN(); i++) {

  double x=X1[i]*fac;
  double exl=EXlow[i]*fac;
  double exh=EXhigh[i]*fac;

  double facy=X1[i]*X1[i]*X1[i]/(2*3.141592654);
  double y= Y1[i]*facy;
  double eyl= EYlow[i]*facy;
  double eyh= EYhigh[i]*facy;

  double eylsyst= EYlowsyst[i]*facy;
  double eyhsyst= EYhighsyst[i]*facy;

  double eylstat= EYlowstat[i]*facy;
  double eyhstat= EYhighstat[i]*facy;

  //cout<<" i= "<<i<<" pt= "<<X1[i]<<" sigma= "<<Y1[i]<<" fac= "<<fac<<endl;
  //cout<<" i= "<<i<<" x= "<<x<<" y= "<<y<<" facy="<<facy<<endl;

  xtdatavector->SetPoint(i, x, y);
  xtdatavector->SetPointError(i,exl,exh,eyl,eyh);

  xtdatavectorstat->SetPoint(i, x, y);
  xtdatavectorstat->SetPointError(i,exl,exh,eylstat,eyhstat);

  xtdatavectorsyst->SetPoint(i, x, y);
  xtdatavectorsyst->SetPointError(i,exl,exh,eylsyst,eyhsyst);
 }

 if (debug) cout<<" Transform finished ! "<<endl;
 //xtdatavectorsyst->Print();

 return;

}

*/
void MyData::PutData(TH1* hdata, double mys, double myymin, double myymax) {
    //
    // input data in histogram
    //
    if (!hdata) cout<<" MyData:: PutData: hdata not found "<<endl;
    ymin=myymin;
    ymax=myymax;
    sqrts=mys;

    if (debug) cout<<" MyData::PutData: input data s= "<<sqrts
                       <<" ymin= "<<ymin<<" ymax= "<<ymax<<endl;

    if (debug)
        hdata->Print();


    TGraphAsymmErrors* g1= new TGraphAsymmErrors();
    if (!g1) cout<<" MyData:: PutData: graph g1 not created "<<endl;
// else g1->Print();

    double x, y, ex, ey;
    int iv=0;
    for (int i=0; i<hdata->GetNbinsX(); i++) {
        y=hdata->GetBinContent(i+1);
        ey=hdata->GetBinError(i+1);
        x=hdata->GetBinCenter(i+1);
        ex=hdata->GetBinWidth(i+1)/2.;
        if (y!=0 && ey!=0) {
            // cout << i<<" x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;
            g1->SetPoint(iv,x,y);
            g1->SetPointError(iv,ex,ex,ey,ey);
            iv++;
        }
    }

//g1->Print();
    datavector=g1;

// if (debug){
//cout<<"MyData::PutData: after conversion to TGraphAsymmError:"<<endl;
//if (datavector) datavector->Print();
    //}

    //cout<<"MyData::PutData return :"<<endl;

    return;
}

void MyData::DrawData() {

    Double_t x_val_dv;
    Double_t  y_val_dv;
    Double_t x_val_dvs;
    Double_t y_val_dvs;
    Double_t ylowerr;
    Double_t yhigherr;
    for(int pi = 0; pi < datavector->GetN(); pi++) {
        datavector->GetPoint(pi, x_val_dv, y_val_dv);
        datavectorsyst->GetPoint(pi, x_val_dvs, y_val_dvs);
        ylowerr = datavectorsyst->GetErrorYlow(pi);
        yhigherr = datavectorsyst->GetErrorYhigh(pi);
        if (debug) std::cout << "pi: " << pi << ", x_val_dv: " << x_val_dv << ", x_val_vs: " << x_val_dvs << ", y_val_dv: " << y_val_dv << ", dvs: " << y_val_dvs << " -" << ylowerr << " +" << yhigherr << "\n";
        x_val_dv=0;
        y_val_dv=0;
        x_val_dvs=0;
        y_val_dvs=0;
        ylowerr=0;
        yhigherr=0; //reset variables for next iteration
    }
    if (!datavector) cout<<" DrawData: datavector not found "<<endl;
    datavector->SetMarkerStyle(markerstyle);
    datavector->SetMarkerColor(markercolor);
    datavector->SetMarkerSize(markersize);
    datavector->Draw("P same");

    if (!datavectorstat) cout<<" DrawData: datavectorstat not found "<<endl;
    datavectorstat->SetMarkerStyle(markerstyle);
    datavectorstat->SetMarkerColor(markercolor);
    datavectorstat->SetMarkerSize(markersize);
//  datavectorstat->Draw("|| same");
    datavectorstat->Draw("P same");

    if (!datavectortoterr) cout<<" DrawData: datavectortoterr not found "<<endl;
    datavectortoterr->SetMarkerStyle(markerstyle);
    datavectortoterr->SetMarkerColor(markercolor);
    datavectortoterr->SetMarkerSize(markersize);
    datavectortoterr->Draw("P same");

    Double_t x_val;
    Double_t y_val_toterr;
    Double_t y_val_stat;
    Double_t y_val_syst;
    Double_t y_val_def;
    for(int pi = 0; pi < datavectortoterr->GetN(); pi++) {
        datavectortoterr->GetPoint(pi, x_val, y_val_toterr);
        datavectorstat->GetPoint(pi, x_val, y_val_stat);
        datavectorsyst->GetPoint(pi, x_val, y_val_syst);
        datavector->GetPoint(pi, x_val, y_val_def);

        if (debug) std::cout << "pi: " << pi << ", x_val:  "<< x_val << ", tot cont: " << y_val_toterr << ", stat cont: " << y_val_stat << ", syst cont: " << y_val_syst << ", def cont: " << y_val_def << std::endl;
        x_val=0;
        y_val_toterr=0;
        y_val_stat=0;
        y_val_syst=0;
        y_val_def=0; //reset variables for next iteration
    }

    return;
};

void MyData::DrawDataStatOnly() {
    if (!datavectorstat) cout<<" DrawDataStatOnly: datavectorstat not found "<<endl;
    datavectorstat->SetMarkerStyle(markerstyle);
    datavectorstat->SetMarkerColor(markercolor);
    datavectorstat->SetMarkerSize(1.2);
    return;
};
/*
void MyData::DrawDataXtStatOnly(){
   if (!xtdatavectorstat) cout<<" DrawDataStatOnly: datavectorstat not found "<<endl;
   xtdatavectorstat->SetMarkerStyle(markerstyle);
   xtdatavectorstat->SetMarkerColor(markercolor);
   xtdatavectorstat->SetMarkerSize(1.2);
   return;
};

void MyData::DrawDataXt(){

   if (!xtdatavector) cout<<" DrawDataXt: xtdatavector not found "<<endl;
   else if (debug) xtdatavector->Print();

   xtdatavector->SetMarkerStyle(markerstyle);
   xtdatavector->SetMarkerColor(markercolor);
   xtdatavector->SetMarkerSize(1.2);
   xtdatavector->Draw("P");

   if (!xtdatavectorsyst) cout<<" DrawDataXt: xtdatavectorsyst not found "<<endl;
   xtdatavectorsyst->SetMarkerStyle(markerstyle);
   xtdatavectorsyst->SetMarkerColor(markercolor);
   xtdatavectorsyst->SetMarkerSize(1.2);
   xtdatavectorsyst->Draw("||");

   if (!xtdatavectorstat) cout<<" DrawDataXt: xtdatavectorstat not found "<<endl;
   xtdatavectorstat->SetMarkerStyle(markerstyle);
   xtdatavectorstat->SetMarkerColor(markercolor);
   xtdatavectorstat->SetMarkerSize(1.2);

   return;

};
*/

void MyData::DrawData(char name[100], float x, float y) {

    this->DrawLegend(name,x,y);
    this->DrawData();

    return;
}

void MyData::DrawExperimentandYear(double x, double y) {

    Double_t tsize=0.05;
    this->DrawExperiment(x,y);
    x+=0.088;
//y-=0.05;
    TLatex l;
    l.SetTextAlign(12); //l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextSize(tsize);

    char mydate[100];
    sprintf(mydate,"%d",year);
//cout<<" year= "<<year<<" mydate= "<<mydate<<endl;
    l.DrawLatex(x,y,mydate);

    return;
};



void MyData::DrawLegend(char name[100], float x, float y) {

    char text1[100];
    char text2[100];
    char text3[100];
    char text4[100];

    TMarker *marker = new TMarker(x-(0.6*msize),y,8);
    marker->SetMarkerColor(markercolor);
    marker->SetNDC();
    marker->SetMarkerStyle(markerstyle);
    marker->Draw();
    TLatex l1;
    l1.SetTextAlign(12);
    l1.SetNDC();
    l1.SetTextSize(msize);
    TLatex l2;
    l2.SetTextAlign(12);
    l2.SetNDC();
    l2.SetTextSize(msize);
    TLatex l3;
    l3.SetTextAlign(12);
    l3.SetNDC();
    l3.SetTextSize(msize);
    TLatex l4;
    l4.SetTextAlign(12);
    l4.SetNDC();
    l4.SetTextSize(msize);
    sprintf(text1,"%s",name);
    sprintf(text2,"%d",year);
    sprintf(text3,"#sqrt{s} = %4.0f GeV",sqrts);
    double fractpart, intpart;
    fractpart= std::modf(ymax*10., &intpart);

    //cout<<" ymax= "<<ymax<<" fractpart= "<<fractpart<<" intpart= "<<intpart<<endl;

    if (fractpart<0.01) sprintf(text4,"%3.1f< |y| < %3.1f",ymin,ymax);
    else                sprintf(text4,"%3.2f< |y| < %3.2f",ymin,ymax);

    if (std::fabs(ymin)<0.01) {
        if (fractpart<0.01)  sprintf(text4," |y| < %3.1f",ymax);
        else                 sprintf(text4," |y| < %3.2f",ymax);

    }
    //cout<<" test4= "<<text4<<endl;
    //sprintf(text,"%s  %s %3.1f< |y| < %3.1f",name,text1,ymin,ymax);
    l1.DrawLatex(x,y,text1);
    l2.DrawLatex(x+0.09,y,text2);
    l3.DrawLatex(x+0.16,y,text3);
    l4.DrawLatex(x+0.34,y,text4);

    return;
}

//<<<<<<< .mine
void MyData::ScaleGraph(TGraphAsymmErrors *g1, double scalex, double scaley) {

    Double_t* X1 = g1->GetX();
    Double_t* Y1 = g1->GetY();
    Double_t* EXhigh1 = g1->GetEXhigh();
    Double_t* EXlow1 =  g1->GetEXlow();
    Double_t* EYhigh1 = g1->GetEYhigh();
    Double_t* EYlow1 =  g1->GetEYlow();

    for (Int_t i=0; i<g1->GetN(); i++) {
        g1->SetPoint(i,X1[i]*scalex,Y1[i]*scaley);
        g1->SetPointError(i,EXlow1[i]*scalex,EXhigh1[i]*scalex,
                          EYlow1[i]*scaley,EYhigh1[i]*scaley);
    }
    return;
};

void MyData::Scale(double scalex, double scaley) {
    this->ScaleGraph(datavector      ,scalex,scaley);
    this->ScaleGraph(datavectorstat  ,scalex,scaley);
    this->ScaleGraph(datavectorsyst  ,scalex,scaley);
    this->ScaleGraph(datavectortoterr,scalex,scaley);
    miny*=scaley;
    maxy*=scaley;
    minx*=scalex;
    maxx*=scalex;
    return;
};


//=======
//>>>>>>> .r147604

void MyData::split_string(std::string str, std::vector<std::string>& split_results, std::string delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        split_results.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

