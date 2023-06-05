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
	Long64_t TimeDiffBefore[Ndet];
	Short_t EnergyDepBefore[Ndet];
	Long64_t TimeDiffAfter[Ndet];
	Short_t EnergyDepAfter[Ndet];
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
/*
bool OrderFunction(Event* A, Event* B){
	return A->Time < B->Time;
};

vector<Event*> AddTimeDifferences(
	vector<Event*> GoodEvents
){
	//Sorting Process
	cout << "Sorting events with ascending time order ..." << endl;
	sort(GoodEvents.begin(), GoodEvents.end(), OrderFunction);
	cout << "Finished sorting events" << endl;

	//Variable initilization
	Event* LastEventIn0 = NULL;
	Event* LastEventIn1 = NULL;
	ULong64_t i = 0;

	//TimeDiffBefore and EnergyDepBefore for Events which may not have previous Event
	cout << "Start determining TimeDiffBefore and EnergyDepBefore ..." << endl;
	do
	{
		if(LastEventIn0 != NULL) {
			GoodEvents[i]->TimeDiffBefore0 = GoodEvents[i]->Time - LastEventIn0->Time;
			GoodEvents[i]->EnergyDepBefore0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			GoodEvents[i]->TimeDiffBefore1 = GoodEvents[i]->Time - LastEventIn1->Time;
			GoodEvents[i]->EnergyDepBefore1 = LastEventIn1->Adc;
			}
		if(GoodEvents[i]->Det == 0){
			LastEventIn0 = GoodEvents[i];
		}
		if(GoodEvents[i]->Det == 1){
			LastEventIn1 = GoodEvents[i];
		}
		i++;
	} while (i < GoodEvents.size() && (LastEventIn0==NULL || LastEventIn1==NULL));
	//Just a test of whether I did something wrong
	if(LastEventIn0 != NULL && LastEventIn1 != NULL){
		cout << "Two Previous events are found" << endl; 
	}
	else {
		cout << "Error with finding two previous events: INTERRUPT" << endl;
		return NULL;
	}
	//TimeDiffBefore and EnergyDepBefore with definite previous event
	for (i; i < GoodEvents.size(); i++){
		GoodEvents[i]->TimeDiffBefore0 = GoodEvents[i]->Time - LastEventIn0->Time;
		GoodEvents[i]->EnergyDepBefore0 = LastEventIn0->Adc;
		GoodEvents[i]->TimeDiffBefore1 = GoodEvents[i]->Time - LastEventIn1->Time;
		GoodEvents[i]->EnergyDepBefore1 = LastEventIn1->Adc;
		if(GoodEvents[i]->Det == 0){
			LastEventIn0 = GoodEvents[i];
		}
		if(GoodEvents[i]->Det == 1){
			LastEventIn1 = GoodEvents[i];		
		}
	}
	cout << "TimeDifferenceBefore and EnergyDepBefore calculated for every event" << endl;

	//Calculating TimeDiffAfter and EnergyDepAfter by iterating reversed through tree
	LastEventIn0 = NULL;
	LastEventIn1 = NULL;
	i = GoodEvents.size()-1;
	cout << "Start determining TimeDiffAfter and EnergyDepAfter" << endl;

	do
	{
		if(LastEventIn0 != NULL) {
			GoodEvents[i]->TimeDiffAfter0 = LastEventIn0->Time - GoodEvents[i]->Time;
			GoodEvents[i]->EnergyDepAfter0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			GoodEvents[i]->TimeDiffAfter1 = LastEventIn1->Time - GoodEvents[i]->Time;
			GoodEvents[i]->EnergyDepAfter1 = LastEventIn1->Adc;
			}
		if(GoodEvents[i]->Det == 0){
			LastEventIn0 = GoodEvents[i];
		}
		if(GoodEvents[i]->Det == 1){
			LastEventIn1 = GoodEvents[i];
		}
		i--;
	} while (i >= 0 && (LastEventIn0==NULL || LastEventIn1==NULL));
	//Just a test of whether I did something wrong
	if(LastEventIn0 != NULL && LastEventIn1 != NULL){
		cout << "Two following events are found" << endl; 
	}
	else {
		cout << "Error with finding two following events: INTERRUPT" << endl;
		return NULL;
	}
	//TimeDiffBefore and EnergyDepBefore with definite previous event
	for (i; i >= 0; i--){
		GoodEvents[i]->TimeDiffAfter0 = LastEventIn0->Time - GoodEvents[i]->Time;
		GoodEvents[i]->EnergyDepAfter0 = LastEventIn0->Adc;
		GoodEvents[i]->TimeDiffAfter1 = LastEventIn1->Time - GoodEvents[i]->Time;
		GoodEvents[i]->EnergyDepAfter1 = LastEventIn1->Adc;
		if(GoodEvents[i]->Det == 0){
			LastEventIn0 = GoodEvents[i];
		}
		if(GoodEvents[i]->Det == 1){
			LastEventIn1 = GoodEvents[i];		
		}
	}
	cout << "TimeDiffAfter and EnergyDepAfter calculated for every event" << endl;

	return GoodEvents;
}
*/
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
    DataNew->Branch("time",&CurrentEvent.Time,"time/l");
    for (int det = 0; det < Ndet; det++) {
		TString name = "TimeDiff_before" + to_string(det);
		TString leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.TimeDiffBefore[det], leaflist.Data());
		name = "TimeDiff_after" + to_string(det);
		leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.TimeDiffAfter[det], leaflist.Data());
		name = "EnergyDep_before" + to_string(det);
		leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.EnergyDepBefore[det], leaflist.Data());
		name = "EnergyDep_after" + to_string(det);
		leaflist = name + "/L";
		DataNew->Branch(name.Data(), &CurrentEvent.EnergyDepAfter[det], leaflist.Data());
	}
	cout << "Branches created" << endl;
	
	//Iterate over vector and write Events on the Tree
	cout << "Start writing events into the new tree ..." << endl;
	for (ULong64_t i = 0; i < GoodEvents.size(); i++){
		CurrentEvent = *GoodEvents[i];
		DataNew->Fill();
	}

	DataNew->Write();
	cout << "All Events are written into tree!" << endl;
	return File;
}


void MakeTimeDiffTree(	
	TString run = "run001"

){
	TTree* oldTree = GetSortedTree(run);

	//Do all functions implemented above
	vector<Event*> GoodEvents = TreeEntries(oldTree);
	
	//GoodEvents = AddTimeDifferences(GoodEvents);
	
	TFile* NewFile = CreateNewFile(GoodEvents, "run001");
	
	NewFile->Write();
	cout << "New root file saved" << endl;
}
