/*
 * Name: Cameron Embree
 * contact: CSEmbree@gmail.com
 * Date created: 7-June-2013
 * Date last edited: 9-June-2013
*/
#include "OptionHandler.h"
using namespace std;



//default constructor
OptionHandler::OptionHandler(std::string fName)
{
    debug=false; //allows additional terminal output
    fileName=fName;
    knownOptions.clear();
    resultsFileNameOn=false;
    readOptions(); //reads file name provided by constuctor
}



//open and read the options that are to be supported.
void OptionHandler::readOptions()
{    
    std::ifstream infile(fileName.c_str(), ios::in);
    if(!infile){ // Check open
        cerr<<" OptionHandler::readOptions: can't open "<<fileName<<"\n";
        infile.close();
        exit (1);
    } 
    else {
        if (debug)cout<<" OptionHandler::readOptions: read data correlations file: "<<fileName<<endl;
    }

    char line[1024];
    while(infile.good())
    {
        infile.getline(line,sizeof(line),'\n');
        //don't retreive comments and spaces
        if(isalpha(line[0])) {
            knownOptions.push_back(string(line));
            if (debug) cout<<"\tOptionHandler::readOptions: Retreived Option: '"<<line<<"'"<<endl;
        }
        memset(line, '\0', 1024); //clear buffer for next iteration
    }
}



//public function for open and read the options that are to be supported. 
//Allow for different options file to be read than the one provided to default constructor
void OptionHandler::readOptions(std::string fName)
{    
    fileName=fName; //overwrite any previously set file name
    knownOptions.clear(); //clear out any previously read in options
    
    std::ifstream infile(fileName.c_str(), ios::in);
    if(!infile){ // Check open
        cerr<<" OptionHandler::readOptions: can't open "<<fileName<<"\n";
        infile.close();
        exit (1);
    } 
    else {
        if (debug) cout<<" OptionHandler::readOptions: read data correlations file: " << fileName << endl;
    }

    char line[1024];
    while(infile.good())
    {
        infile.getline(line,sizeof(line),'\n');
        //don't retreive comments and spaces
        if(isalpha(line[0])) {
            knownOptions.push_back(string(line));
            if (debug) cout <<"\tOptionHandler::readOptions: Retreived Option: '" << line << "'" << endl;
        }
        memset(line, '\0', 1024); //clear buffer for next iteration
    }
}



//check if provided argument is a supported option
bool OptionHandler::isKnownOption(std::string optionName)
{            
    //check through each option to check if requested option is supported
    for(int i=0; i<knownOptions.size(); i++){
        if(knownOptions.at(i).compare(optionName)==0) {
            if(debug) std::cout<<"OptionHandler::isKnownOption: option: '"<<optionName<<"' is supported."<<endl;
            if(resultsFileNameOn==true){
                std::string resTrue = optionName.append(" == true\n");
                fputs (resTrue.c_str(),resultFile);
                //resultFile << resTrue.c_str();
                if(debug) std::cout<<"OptionHandler::isKnownOption: results have been printed to file"<<endl;
            }
            return true; //option is supported
        }
    }
    
        
    if(debug) std::cout<<"OptionHandler::isKnownOption: WARNING: The option: '"<<optionName<<"' is unknown. Check support file:'"<<fileName<<"'"<<endl;
    if(resultsFileNameOn==true){
        std::string resFalse = optionName.append(" == false\n");
        fputs (resFalse.c_str(),resultFile);
        //resultFile << resFalse.c_str();
    }
    return false;   //option is unkown
}



//allow user to send results of OptionHandler method 'isKnownOptions' calls to a file
//default file name is "_OH_RESULTS.txt appended to the fileName"
void OptionHandler::generateResultFile(std::string fResultsName)
{
/*
    string fName = fileName;
    
    //allow for altering of default results file name
    if(fResultsName.compare(string(""))!=0){
        fName = fResultsName;
    }
    if(debug) std::cout<<"OptionHandler::generateResultFile: received filename: "<<fName<<endl;    
    
    resultsFileNameOn=true;
    resultsFileName=fName.substr(0,fName.rfind(".txt")).append("_OH_RESULTS.txt");
    resultFile = fopen (resultsFileName.c_str(),"w");
    //resultFile.open(resultsFileName.c_str());
    
    if(debug) std::cout<<"OptionHandler::generateResultFile: outputing results to: "<<resultsFileName<<endl;    
    
    if (resultFile!=NULL){
        //fputs ("test",resultFile);
    }
    else {
        std::cout<<"OptionHandler::generateResultFile: ERROR: issue creating reference file for:"<<fileName<<endl;
    }
    */
    std::cout<<"OptionHandler::generateResultFile: This method is temporarily out of service!"<<endl;
}



void OptionHandler::printOptions()
{
    std::cout<<"OptionHandler::printOptions: Supported options as extracted from: '"<<fileName<<"'"<<endl;
    for(int i=0; i<knownOptions.size(); i++){
        std::cout<<"\t("<<i<<"/"<<knownOptions.size()<<") "<<knownOptions.at(i);
    }
}

