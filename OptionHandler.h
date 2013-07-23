#ifndef OPTION_HANDLER_H
#define OPTION_HANDLER_H
/*
 * Name: Cameron Embree
 * contact: CSEmbree@gmail.com
 * Date created: 7-June-2013
 * Date last edited: 9-June-2013
*/

#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <fstream>
//#include <direct.h>



class OptionHandler {
    public:
        //variables
        //<none>
        
        //methods
        OptionHandler(std::string FileName); //default constructor
        bool isKnownOption(std::string optionName);
        void readOptions(std::string fName);
        void generateResultFile(std::string fResultsName="");
        void printOptions();
        
        std::vector<std::string> getOptions() {
            return knownOptions;
        }
        int getNumOptions() {
            return knownOptions.size();
        }
        std::string getKnownOptionsPathName() {
            return fileName;
        }
        
        
    private:
        //variables
        bool debug;
        std::vector<std::string> knownOptions;
        std::string fileName;
        std::string resultsFileName;
        bool resultsFileNameOn;
        FILE *resultFile;
        //std::ofstream resultFile;
        
        //methods
        void readOptions(); //called by default constructor
};
  
  
#endif
