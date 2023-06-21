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

/* 
 * First version of a complex coincidence scrip, which shall return histograms based on complex coincidence cases, 
 * i.e. events in TU4 which have an entry in a significant energy peak in TU5 and vice versa.
 * v1.0.0
*/

typedef vector<Int_t> EThr;

//This Event struct is reduced (piuleup, saturation and Time missing )
struct Event {
	Short_t Adc;
	UShort_t Det;
//	ULong64_t Time;
//	Bool_t Pileup;
//	Bool_t Saturation;
	Long64_t TimeDiffBefore0;
	Short_t EnergyDepBefore0;
	Long64_t TimeDiffAfter0;
	Short_t EnergyDepAfter0;
	Long64_t TimeDiffBefore1;
	Short_t EnergyDepBefore1;
	Long64_t TimeDiffAfter1;
	Short_t EnergyDepAfter1;
};



//Open the root file containing the respective histograms
TFile* GetInputFile(
    TString run = "run001"
){
	TString folder = "../Test_230526/"+run+"/root/";
    // Open ROOT file
    TFile* rootFile = new TFile(folder+"Data_"+run+"_coincidence.root","UPDATE");
    //if (rootFile->IsZombie()) return;
    
	return rootFile;
}

TH1D* Coincidence(
	TTree* Data,
	Int_t Det,
	Int_t Det1,
	Double_t TimeDiffThr,
){
	Event CurrentEvent;
	if (!Data) {cout << "Error: Faild to retrieve TTree from the input file." << endl;
		return NULL;
		}
	
}

void MakeCoincidenceHist(){
	TString Run = "run002";
	
	TFile* File = GetInputFile(Run);
	TTree* Data = (TTree*)File->Get("Data");
	
	Int_t TU5 = 0;
	Int_t TU4 = 1;
	Double_t thrTimeDiff = 5E6
	EThr Xray1 = {};
	EThr Xray2 = {};
	EThr Xray3 = {};
	
	vector<EThr> Ba133Xrays = {Xray1, Xray2, Xray3};
	
	
}
