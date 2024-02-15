
#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
#include <iomanip>
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include <TUUID.h>
#include "TLegend.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "THStack.h"
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <dirent.h>

#define ROIS 4
#define MAXRONLINES 10000000
#define nentriesMAX 200000000

vector<TString> GetFilenames(
    TString target_dir,
){
    DIR *dr_bin = opendir(input); 
    if (dr_bin == NULL)
    { 
        printf("Could not find any .bin files" ); 
    } 
  
    while ((de_bin = readdir(dr_bin)) != NULL){ 
	    string s_bin = de_bin->d_name;
	    //check in which file the letters "b", "i" and "n" occur
	    if(s_bin.find('B') != std::string::npos && s_bin.find('I') != std::string::npos && s_bin.find('N') != std::string::npos){
	      binfiles.push_back(s_bin);
	    }
    }

    int total_number_files = binfiles.size();
    sort(binfiles.begin(), binfiles.end()); //sort the files in ascending order

    cout << "Following .BIN files were found:" << endl;
    
    for(int i=0;i<total_number_files;i++){
      cout << binfiles[i] << endl;
    }
    cout << " " << endl;

    closedir(dr_bin);
}


void BinToRoot(
    TString run,
    bool hist,
    bool drawhist,
    double maxtime = 864000.,
    double mintime = 0.
){
    TString prefixPath = "";
	TString folder = "/home/chris/Projects/TU5TU7/Experiments/Data/Raw/";  //Path of the .bin file which shall be evaluated
	TString suffix = "/RAW/";
	TString target_dir = "/home/chris/Projects/TU5TU7/Experiments/Data/Root/";

    vector<TString> binfiles;

    DIR *dr_bin = opendir(input); 
    if (dr_bin == NULL)
    { 
        printf("Could not find any .bin files" ); 
    } 
  
    while ((de_bin = readdir(dr_bin)) != NULL){ 
	    string s_bin = de_bin->d_name;
	    //check in which file the letters "b", "i" and "n" occur
	    if(s_bin.find('B') != std::string::npos && s_bin.find('I') != std::string::npos && s_bin.find('N') != std::string::npos){
	      binfiles.push_back(s_bin);
	    }
    }

    int total_number_files = binfiles.size();
    sort(binfiles.begin(), binfiles.end()); //sort the files in ascending order

    cout << "Following .BIN files were found:" << endl;
    
    for(int i=0;i<total_number_files;i++){
      cout << binfiles[i] << endl;
    }
    cout << " " << endl;

    closedir(dr_bin);

    TString filename = folder + run + suffix + "DataR" run + ".BIN";
    TString rootfilename = target_dir + run + ".root";
    cout << "filename: " << filename << endl;
    cout << "rootfilename: " << rootfilename << endl;

    // Create new tree

    int16_t Board;
    int16_t Channel;
    int64_t TimeStamp;
    int16_t Energy_ch;
    uint32_t Flags;
    char waveformCode;
    int32_t numSamples;
    int16_t[numSamples] WfSample;

    TFile *f = new TFile(rootfilename, "RECREATE");
    TTree *t = new TTree("t", "t");
    t->Branch("Board", &Board, "Board/S");
    t->Branch("Channel", &Channel, "Channel/S");
    t->Branch("TimeStamp", &TimeStamp, "TimeStamp/L");
    t->Branch("Energy_ch", &Energy_ch, "Energy_ch/S");
    t->Branch("Flags", &Flags, "Flags/i");
    //t->Branch("waveformCode", &waveformCode, "waveformCode/B");
    //t->Branch("numSamples", &numSamples, "numSamples/I");
    t->Branch("WfSample", &WfSample, "WfSample[numSamples]/S");

    for (i=0;i<total_number_files;i++) {
		caenFileName = input + binfiles[i];
		cout << "Open " << caenFileName << endl;
        fin[i] = fopen(caenFileName, "r");
        if (!fin[i]) {
             cout << "Error when opening file " << caenFileName << endl;
             return;
        } else {
            cout << "Opened file " << caenFileName << endl;
        }  // if (!fin[i])
    }
    cout << " " << endl;





    // Create new directory for the root files
    struct stat st = {0};

	if (stat(target_dir + run + "/", &st) == -1) {
		int directory1 = mkdir(target_dir + run, 0777);
		cout << "Created new directory" << endl;
	}
    TFile *rootFile = new TFile(target_dir + run + "/" + run + "_without_timediff.root","RECREATE");
	cout << "created new file " << endl;
    
}