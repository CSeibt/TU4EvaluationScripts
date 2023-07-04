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

#define Ndet 2 // number of detectors for coincidence analysis

#define MaxTime ULLONG_MAX;

TF1* calibration_function[Ndet];

TTree* GetSortedTree(
    TString run = "run001"
){
    // Open ROOT file
    TFile* inputFile = new TFile(run+".root","UPDATE");
    //if (rootFile->IsZombie()) return;
    
    // Get data tree
    TTree* inputTree = (TTree*)inputFile->Get("Data");
    
	return inputTree;
}
	

struct Event {
	Short_t Adc;
	UShort_t Det;
	ULong64_t Time;
	//Bool_t Pileup;
	//Bool_t Saturation;
	Double_t Edep;
	Long64_t TimeDiffBefore[Ndet];
	Double_t EnergyDepBefore[Ndet];
	Long64_t TimeDiffAfter[Ndet];
	Double_t EnergyDepAfter[Ndet];
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
		if (!pileup && !saturation) {
			Event* CurrentEvent = new Event;
			CurrentEvent->Adc = adc;
			CurrentEvent->Det = det;
			CurrentEvent->Time = time;
			//CurrentEvent->Pileup = pileup;
			//CurrentEvent->Saturation = saturation;
			CurrentEvent->Edep = 0;
			for (int det = 0; det < Ndet; det++) {
				CurrentEvent->TimeDiffBefore[det] = MaxTime;
				CurrentEvent->EnergyDepBefore[det] = 0;
				CurrentEvent->TimeDiffAfter[det] = MaxTime;
				CurrentEvent->EnergyDepAfter[det] = 0;
			}
			TreeEvents.push_back(CurrentEvent);
		}
		if (i % PercentEntries == 0){
			cout << i / PercentEntries << "% of total events finished" << endl;
		}
	}
	return TreeEvents;
}

void AddTimeDifferences(
// Add calibrated energy and time differences to the events in the vector.
	vector<Event*> GoodEvents
){
	//Variable initilization
	//Short_t LastAdc[Ndet];
	//ULong64_t LastTime[Ndet];
	Event* LastEventIn[Ndet];
	
	int firstevent = 0;
	int lastevent = GoodEvents.size() - 1;
	
	// Dummy events for events which may not have a previous event
	for (int det = 0; det < Ndet; det++) {
		//LastAdc[det] = 0;
		//LastTime[det] = 0;
		LastEventIn[det] = new Event;
		LastEventIn[det]->Edep = 0;
		LastEventIn[det]->Det = det;
		LastEventIn[det]->Time = GoodEvents[firstevent]->Time;
	}

	// init pseudo-random generator for MC energy calibration
	TRandom3* rnd = new TRandom3(time(NULL));
	for (int det = 0; det < Ndet; det++) {
		if (!calibration_function[det])
			cout << "No calibration function for detector " << det << endl;
	}

	//TimeDiffBefore and EnergyDepBefore with definite previous event
	cout << "Energy calibration, TimeDifferenceBefore and EnergyDepBefore..." << endl;
	for (ULong64_t i = firstevent; i <= lastevent; i++){
		// MC energy calibration
		Double_t ch = GoodEvents[i]->Adc + rnd->Uniform(0, 1);
		Double_t E = calibration_function[GoodEvents[i]->Det]->Eval(ch);
		GoodEvents[i]->Edep = E;

		for (int det = 0; det < Ndet; det++) {
			//GoodEvents[i]->TimeDiffAfter[det] = GoodEvents[i]->Time - LastTime[det];
			//GoodEvents[i]->EnergyDepAfter[det] = LastAdc[det];
			GoodEvents[i]->TimeDiffBefore[det] = GoodEvents[i]->Time - LastEventIn[det]->Time;
			GoodEvents[i]->EnergyDepBefore[det] = LastEventIn[det]->Edep;
		}
		//cout << "event " << i << " det" << GoodEvents[i]->Det << " Adc " << GoodEvents[i]->Adc << " EnergyDepBefore[0] " << GoodEvents[i]->EnergyDepBefore[0] << " EnergyDepBefore[1] " << GoodEvents[i]->EnergyDepBefore[1] << endl;
		//LastTime[GoodEvents[i]->Det] = GoodEvents[i]->Time;
		//LastAdc[GoodEvents[i]->Det] = GoodEvents[i]->Adc;
		LastEventIn[GoodEvents[i]->Det] = GoodEvents[i];
	}
	cout << "TimeDifferenceBefore and EnergyDepBefore calculated for every event" << endl;

	// Dummy events for events which may not have a previous event
	for (int det = 0; det < Ndet; det++) {
		//LastAdc[det] = 0;
		//LastTime[det] = MaxTime;
		LastEventIn[det] = new Event;
		LastEventIn[det]->Edep = 0;
		LastEventIn[det]->Det = det;
		LastEventIn[det]->Time = GoodEvents[lastevent]->Time;
	}
	
	//Calculating TimeDiffAfter and EnergyDepAfter by iterating reversed through tree
	cout << "Start determining TimeDiffAfter and EnergyDepAfter" << endl;

	//TimeDiffAfter and EnergyDepAfter with definite next event
	for (int i = lastevent; i >= firstevent; i--){
		for (int det = 0; det < Ndet; det++) {
			//GoodEvents[i]->TimeDiffAfter[det] = LastTime[det] - GoodEvents[i]->Time;
			//GoodEvents[i]->EnergyDepAfter[det] = LastAdc[det];
			GoodEvents[i]->TimeDiffAfter[det] = LastEventIn[det]->Time - GoodEvents[i]->Time;
			GoodEvents[i]->EnergyDepAfter[det] = LastEventIn[det]->Edep;
		}
		//cout << "event " << i << " det" << GoodEvents[i]->Det << " Adc " << GoodEvents[i]->Adc << " EnergyDepAfter[0] " << GoodEvents[i]->EnergyDepAfter[0] << " EnergyDepAfter[1] " << GoodEvents[i]->EnergyDepAfter[1] << endl;
		//LastTime[GoodEvents[i]->Det] = GoodEvents[i]->Time;
		//LastAdc[GoodEvents[i]->Det] = GoodEvents[i]->Adc;
		LastEventIn[GoodEvents[i]->Det] = GoodEvents[i];
	}
	cout << "TimeDiffAfter and EnergyDepAfter calculated for every event" << endl;

	return;
}

TFile* CreateNewFile(
	vector<Event*> GoodEvents,
	TString run = "run001"
){
	//Create new root file with Tree DataNew
	TString PathWithFile = "Data_" + run + "_coincidence.root";
	
	cout << "Start creation of new root file ..." << endl; 
	TFile* File = new TFile(PathWithFile, "RECREATE");
	TTree* DataNew = new TTree("Data", "TU5TU4Data");
    cout << "File and Tree created" << endl;
	//Variable & Branch initiation
	Event CurrentEvent;
	
	DataNew->Branch("det",   &CurrentEvent.Det ,   "det/s");
    DataNew->Branch("adc",    &CurrentEvent.Adc,   "ch/S");
    //DataNew->Branch("pileup",&CurrentEvent.Pileup,"pileup/O");
    //DataNew->Branch("saturation",&CurrentEvent.Saturation,"saturation/O");
    DataNew->Branch("Time", &CurrentEvent.Time, "Time/l");
	DataNew->Branch("Edep", &CurrentEvent.Edep, "Edep/D");
    for (int det = 0; det < Ndet; det++) {
		TString name = "TimeDiff_before" + to_string(det);
		TString leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.TimeDiffBefore[det], leaflist.Data());
		name = "TimeDiff_after" + to_string(det);
		leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.TimeDiffAfter[det], leaflist.Data());
		name = "EnergyDep_before" + to_string(det);
		leaflist = name + "/D";
		DataNew->Branch(name.Data(), &CurrentEvent.EnergyDepBefore[det], leaflist.Data());
		name = "EnergyDep_after" + to_string(det);
		leaflist = name + "/D";
		DataNew->Branch(name.Data(), &CurrentEvent.EnergyDepAfter[det], leaflist.Data());
	}
	cout << "Branches created" << endl;
	
	//Iterate over vector and write Events on the Tree
	cout << "Start writing events into the new tree ..." << endl;
	for (ULong64_t i = 0; i < GoodEvents.size(); i++){
		CurrentEvent = *GoodEvents[i];
		//cout << i << " " << CurrentEvent.Adc << " " << CurrentEvent.Edep << " " << CurrentEvent.EnergyDepBefore[0] << " " << CurrentEvent.EnergyDepBefore[1] << endl;
		DataNew->Fill();
	}

	DataNew->Write();
	cout << "All Events are written into tree!" << endl;
	return File;
}


void SetCalibrationHardCoded()
{
	Double_t TU5_p0 = 0.125265;
	Double_t TU5_p1 = 0.00988977;
	Double_t TU4_p0 = -0.380873;
	Double_t TU4_p1 = 0.331185;
	calibration_function[0] = new TF1("fCalTU5", "pol1");
	calibration_function[0]->SetParameters(TU5_p0, TU5_p1);
	calibration_function[1] = new TF1("fCalTU4", "pol1");
	calibration_function[1]->SetParameters(TU4_p0, TU4_p1);
}


void MakeTimeDiffTree(	
	TString run = "run001"

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
