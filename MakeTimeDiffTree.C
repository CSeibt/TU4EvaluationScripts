#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
#include <cstdio>
#include <vector>
#include <algorithm>
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

#define MaxTime ULLONG_MAX;

#define KeepSaturationEvents 0
#define KeepPileUpEvents 1

TF1* calibration_function;

TTree* GetSortedTree(
    TString run = "run001"
){
    // Open ROOT file
    TFile* inputFile = new TFile("../EvaluationData/"+run+".root","READ");
    //if (rootFile->IsZombie()) return;
    
    // Get data tree
    TTree* inputTree = (TTree*)inputFile->Get("Data");
    
	return inputTree;
}
	

struct Event {
	Short_t Adc;
	UShort_t Det;
	ULong64_t Time;
	Bool_t Pileup;
	Bool_t Saturation;
	Double_t Edep;
	Long64_t TimeDiffBefore;
	Double_t EnergyDepBefore;
	Long64_t TimeDiffAfter;
	Double_t EnergyDepAfter;
};

vector<Event*> TreeEntries(
	TTree* inputTree
){
	vector<Event*> TreeEvents = {};
	
	//Get Values from existing Branches
 	Short_t adc;
	UShort_t det;
	ULong64_t time;
	Bool_t pileup;
	Bool_t saturation;
	inputTree->SetBranchAddress("adc", &adc);
	inputTree->SetBranchAddress("det", &det);
	inputTree->SetBranchAddress("Time", &time);
	inputTree->SetBranchAddress("pileup", &pileup);
	inputTree->SetBranchAddress("saturation", &saturation);

	//Write relevant entries in Event vector
	ULong64_t NumberOfEntries = inputTree->GetEntries();
	ULong64_t PercentEntries = (long)NumberOfEntries / 100;
	cout << "Total number of events: " << NumberOfEntries << endl;
	cout << "Start building of event vector ..." << endl;
	for (ULong64_t i = 0; i < NumberOfEntries; i++) {
		inputTree->GetEntry(i);
		if ((KeepPileUpEvents || !pileup) && (KeepSaturationEvents || !saturation)) {
			Event* CurrentEvent = new Event;
			CurrentEvent->Adc = adc;
			CurrentEvent->Det = det;
			CurrentEvent->Time = time;
			CurrentEvent->Pileup = pileup;
			CurrentEvent->Saturation = saturation;
			CurrentEvent->Edep = 0;
			CurrentEvent->TimeDiffBefore = MaxTime;
			CurrentEvent->EnergyDepBefore = 0;
			CurrentEvent->TimeDiffAfter = MaxTime;
			CurrentEvent->EnergyDepAfter = 0;
			TreeEvents.push_back(CurrentEvent);
		}
		if (i % PercentEntries == 0){
			cout << i / PercentEntries << "% of total events finished" << endl;
		}
	}
	return TreeEvents;
}

void AddTimeDifferences(
// Add calibrated energy to the events in the vector.
	vector<Event*> GoodEvents
){
	//Variable initilization
	Event* LastEvent;

	int firstevent = 0;
	int lastevent = GoodEvents.size() - 1;

	// Dummy event for events which may not have a previous event
	LastEvent = new Event;
	LastEvent->Edep = 0;
	LastEvent->Time = GoodEvents[firstevent]->Time;

	// init pseudo-random generator for MC energy calibration
	TRandom3* rnd = new TRandom3(time(NULL));
	if (!calibration_function)
		cout << "No calibration function for detector" << endl;

	//TimeDiffBefore and EnergyDepBefore with definite previous event
	cout << "Energy calibration, TimeDifferenceBefore and EnergyDepBefore..." << endl;
	for (ULong64_t i = firstevent; i <= lastevent; i++){
		// MC energy calibration
		Double_t ch = GoodEvents[i]->Adc + rnd->Uniform(0, 1);
		Double_t E = calibration_function->Eval(ch);
		GoodEvents[i]->Edep = E;

		GoodEvents[i]->TimeDiffBefore = GoodEvents[i]->Time - LastEvent->Time;
		GoodEvents[i]->EnergyDepBefore = LastEvent->Edep;
		LastEvent = GoodEvents[i];
	}
	cout << "Calibrated energy deposition, TimeDifferenceBefore and EnergyDepBefore calculated for every event" << endl;

	LastEvent = new Event;
	LastEvent->Edep = 0;
	LastEvent->Time = GoodEvents[lastevent]->Time;

	//Calculating TimeDiffAfter and EnergyDepAfter by iterating reversed through tree
	cout << "Start determining TimeDiffAfter and EnergyDepAfter" << endl;

	//TimeDiffAfter and EnergyDepAfter with definite next event
	for (int i = lastevent; i >= firstevent; i--){
		GoodEvents[i]->TimeDiffAfter = LastEvent->Time - GoodEvents[i]->Time;
		GoodEvents[i]->EnergyDepAfter = LastEvent->Edep;
		LastEvent = GoodEvents[i];
	}
	cout << "TimeDiffAfter and EnergyDepAfter calculated for every event" << endl;
	return;
}

TFile* CreateNewFile(
	vector<Event*> GoodEvents,
	TString run = "run001"
){
	//Create new root file with Tree DataNew
	TString PathWithFile = "../EvaluationData/Data_" + run + ".root";
	
	cout << "Start creation of new root file ..." << endl; 
	TFile* File = new TFile(PathWithFile, "RECREATE");
	TTree* DataNew = new TTree("Data", "TU5Data");
	cout << "File and Tree created" << endl;
	//Variable & Branch initiation
	Event CurrentEvent;
	
	DataNew->Branch("det",   &CurrentEvent.Det ,   "det/s");
	DataNew->Branch("adc",    &CurrentEvent.Adc,   "ch/S");
	DataNew->Branch("pileup",&CurrentEvent.Pileup,"pileup/O");
	DataNew->Branch("saturation",&CurrentEvent.Saturation,"saturation/O");
	DataNew->Branch("Time", &CurrentEvent.Time, "Time/l");
	DataNew->Branch("Edep", &CurrentEvent.Edep, "Edep/D");
	DataNew->Branch("TimeDiff_before", &CurrentEvent.TimeDiffBefore, "TimeDiff_before/L");
	DataNew->Branch("TimeDiff_after", &CurrentEvent.TimeDiffAfter, "TimeDiff_after/L");
	DataNew->Branch("EnergyDep_before", &CurrentEvent.EnergyDepBefore, "EnergyDep_before/D");
	DataNew->Branch("EnergyDep_after", &CurrentEvent.EnergyDepAfter, "EnergyDep_after/D");
	cout << "Branches created" << endl;
	
	//Iterate over vector and write Events on the Tree
	cout << "Start writing events into the new tree ..." << endl;
	for (ULong64_t i = 0; i < GoodEvents.size(); i++){
		CurrentEvent = *GoodEvents[i];
		//cout << i << " " << CurrentEvent.Adc << " " << CurrentEvent.Edep << endl;
		DataNew->Fill();
	}

	DataNew->Write();
	cout << "All Events are written into tree!" << endl;
	return File;
}


void SetCalibrationHardCoded()
{
	Double_t TU5_p0 = 0.14812;
	Double_t TU5_p1 = 0.00987247;
	calibration_function = new TF1("fCalTU5", "pol1");
	calibration_function->SetParameters(TU5_p0, TU5_p1);
}


void MakeCalibratedTree(	
	TString run = "run038"

){
	TTree* oldTree = GetSortedTree(run);

	//Do all functions implemented above
	vector<Event*> GoodEvents = TreeEntries(oldTree);
	
	SetCalibrationHardCoded();
	AddTimeDifferences(GoodEvents);
	
	TFile* NewFile = CreateNewFile(GoodEvents, run);
	
	NewFile->Write();
	cout << "New root file saved" << endl;
}
