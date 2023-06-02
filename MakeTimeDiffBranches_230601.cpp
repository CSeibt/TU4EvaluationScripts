#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
#include <cstdio>
#include <vector>
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

TFile* File(
    TString run = "run021",
    TString prefixPath = "Test_230526/"
){
	TString folder = "root/";//"";
	TString list = "DataR_";
    TString Suffix = ".bin";
    
    // Variable definitions.
    TString listModeSuffix = list+run+Suffix;
    TString prefix = prefixPath + run + "/" + folder;

    
    // Open ROOT file
    TFile* rootFile = new TFile(prefix+listModeSuffix+".root","UPDATE");
    //if (rootFile->IsZombie()) return;
    
	return rootFile;
	}
	
	
struct Event {
	Short_t adc;
	UShort_t det;
	ULong64_t time;
	Bool_t pileup;
	Bool_t saturation;
	Long64_t TimeDiffBefore0 = 0;
	Short_t EnergyDepBefore0 = 0;
	Long64_t TimeDiffAfter0 = 0;
	Short_t EnergyDepAfter0 = 0;
	Long64_t TimeDiffBefore1 = 0;
	Short_t EnergyDepBefore1 = 0;
	Long64_t TimeDiffAfter1 = 0;
	Short_t EnergyDepAfter1 = 0;
};

vector<Event*> TreeEntries(
	TFile* File
){
	vector<Event*> TreeEvents = {};
	TTree* Tree = dynamic_cast<TTree*>(file->Get
}

