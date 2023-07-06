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

typedef vector<Short_t> EThr;

//This Event struct is reduced (piuleup, saturation and Time missing )
struct Event {
	Short_t Adc;
	UShort_t Det;
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

/*
 * This function creates an adc histogram of the main detector (Det1) in regard to coincident events in the other detector (Det2)
 * An event is considered for the histogram if:
 * 1. The time difference towards an event in Det2 is within the timediff threshold
 * 2. The energy deposited in Det2 accords to an X-ray peak, i.e. it is within one of the enrgy thresholds
*/
TH1D* Coincidence(
	TTree* Data,
	Int_t Det1,
	Int_t Det2,
	Long64_t TimeDiffThr,
	vector<EThr> EnergyThrs,
	Int_t Channels,
	TString Name
){
	//Histogram and debug output
	TH1D* Histogram = new TH1D(Name, Name+" ADC spectrum; Channel; Counts", Channels, 0, Channels);
	if (!Data) {cout << "Error: Faild to retrieve TTree from the input file." << endl;
		return Histogram;
	}
	//Set Branch adresses to the structs variable pointers
	Event CurrentEvent;
	Data->SetBranchAddress("det", &CurrentEvent.Det);
	Data->SetBranchAddress("adc", &CurrentEvent.Adc);
	Data->SetBranchAddress("TimeDiff_before0", &CurrentEvent.TimeDiffBefore0);
	Data->SetBranchAddress("TimeDiff_after0", &CurrentEvent.TimeDiffAfter0);
	Data->SetBranchAddress("EnergyDep_before0", &CurrentEvent.EnergyDepBefore0);
	Data->SetBranchAddress("EnergyDep_after0", &CurrentEvent.EnergyDepAfter0);
	Data->SetBranchAddress("TimeDiff_before1", &CurrentEvent.TimeDiffBefore1);
	Data->SetBranchAddress("TimeDiff_after1", &CurrentEvent.TimeDiffAfter1);
	Data->SetBranchAddress("EnergyDep_before1", &CurrentEvent.EnergyDepBefore1);
	Data->SetBranchAddress("EnergyDep_after1", &CurrentEvent.EnergyDepAfter1);
	cout << "Branch Addresses set to CurrentEvent members" << endl;
	
	//Console output for user
	ULong64_t EntryNumber = Data->GetEntries();
	ULong64_t PercentEntries = (long)EntryNumber / 100;
	int k = 1;
	
	///Start histogram filling
	/*
	 * Steps and cases for histogram filling:
	 * 1) Det2 check, use correct TimeDiffs and EnergyDeps
	 * 2) Check whether the current event is an Event of Det1
	 * 3) Apply TimeDiff_before and TimeDifff_after cut
	 * 4.1) TimeDiff_before -> EnergyDep_before -> Energy thresholds
	 * 4.2) TimeDiff_after -> EnergyDep_after -> Energy thresholds
	*/
	if (Det2 == 0){		// 1)
		for (ULong64_t i = 0; i < EntryNumber; i++){
			Data->GetEntry(i);
			if (CurrentEvent.Det == Det1){		// 2)
				if (CurrentEvent.TimeDiffBefore0 < TimeDiffThr){		// 3)
					for (int j = 0; j < EnergyThrs.size(); j++){
						if(CurrentEvent.EnergyDepBefore0 > EnergyThrs[j][0] && CurrentEvent.EnergyDepBefore0 < EnergyThrs[j][1]){		// 4.1)
							Histogram->Fill(CurrentEvent.Adc);
						}
					}
				} else if (CurrentEvent.TimeDiffAfter0 < TimeDiffThr){		// 3)
					for (int j = 0; j < EnergyThrs.size(); j++){
						if(CurrentEvent.EnergyDepAfter0 > EnergyThrs[j][0] && CurrentEvent.EnergyDepAfter0 < EnergyThrs[j][1]){			// 4.2)
							Histogram->Fill(CurrentEvent.Adc);
						}
					}
				}
			}
			if (i == k*PercentEntries){		//User information
				cout << k << "% of total events finished" << endl;
				k++;
			}
		}
	}
	if (Det2 == 1){		// 1)
		for (ULong64_t i = 0; i < EntryNumber; i++){
			Data->GetEntry(i);
			if (CurrentEvent.Det == Det1){		// 2)
				if (CurrentEvent.TimeDiffBefore0 < TimeDiffThr){		// 3)
					for (int j = 0; j < EnergyThrs.size(); j++){
						if(CurrentEvent.EnergyDepBefore0 > EnergyThrs[j][0] && CurrentEvent.EnergyDepBefore0 < EnergyThrs[j][1]){		// 4.1)
							Histogram->Fill(CurrentEvent.Adc);
						}
					}
				} else if (CurrentEvent.TimeDiffAfter0 < TimeDiffThr){		// 3)
					for (int j = 0; j < EnergyThrs.size(); j++){
						if(CurrentEvent.EnergyDepAfter0 > EnergyThrs[j][0] && CurrentEvent.EnergyDepAfter0 < EnergyThrs[j][1]){			// 4.2)
							Histogram->Fill(CurrentEvent.Adc);
						}
					}
				}
			}
			if (i == k*PercentEntries){
				cout << k << "% of total events finished" << endl;
				k++;
			}
		}
	}
	return Histogram;
}

/* Main function */
void MakeCoincidenceHist(){
	TString Run = "run002";
	//File and Tree
	TFile* File = GetInputFile(Run);
	TTree* Data = (TTree*)File->Get("Data");
	//Coincidence() parameters
	Int_t nch = 16384;
	Int_t TU5 = 0;
	Int_t TU4 = 1;
	Long64_t thrTimeDiff = 5E6;
	EThr Xray1 = {417, 597};
	EThr Xray2 = {2700, 3173};
	EThr Xray3 = {3489, 3667};
	
	vector<EThr> Ba133Xrays = {Xray1, Xray2, Xray3};
	//Histograms
	TH1D* CoincidenceTest = Coincidence(Data, TU4, TU5, thrTimeDiff, Ba133Xrays, nch, "TU4_adc_coincidence");
	TCanvas* TestCanvas = new TCanvas("TU4Test", "TU4Test", 100, 10, 1200, 700);
	CoincidenceTest->Draw();
	CoincidenceTest->Write();
	
}
