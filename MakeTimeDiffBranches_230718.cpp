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
#include <sys/resource.h>

/* MakeTimeDiffBranches file for larger root files, test version only */

const ULong_t MaxTime = kMaxInt;

TFile* File(
    TString run = "run021",
    TString prefixPath = "../Test_230526/"
){
	TString folder = "root/";//"";
	TString list = "DataR_";
    TString Suffix = ".bin_without_timediff";
    
    // Variable definitions.
    TString listModeSuffix = list+run+Suffix;
    TString prefix = prefixPath + run + "/" + folder;

    
    // Open ROOT file
    TFile* rootFile = new TFile(prefix+listModeSuffix+".root","UPDATE");
    //if (rootFile->IsZombie()) return;
    
	return rootFile;
	}


struct Event {
	Event(Short_t adc, UChar_t det, ULong_t time) :
		Adc(adc), Det(det), Time(time)
		{}
	Short_t Adc;
	UChar_t Det;
	ULong_t Time;
};


struct NewEvent {
	NewEvent(Short_t adc, UChar_t det, ULong_t time) :
		Adc(adc), Det(det), Time(time),
		TimeDiffBefore0(MaxTime), TimeDiffAfter0(MaxTime),
		TimeDiffBefore1(MaxTime), TimeDiffAfter1(MaxTime),
		EnergyDepBefore0(0), EnergyDepAfter0(0),
		EnergyDepBefore1(0), EnergyDepAfter1(0)
		{}
	Short_t Adc;
	UChar_t Det;
	ULong64_t Time;
	// --- Defined ĺater
	Short_t EnergyDepBefore0;
	Short_t EnergyDepAfter0;
	Short_t EnergyDepBefore1;
	Short_t EnergyDepAfter1;
	Int_t TimeDiffBefore0;
	Int_t TimeDiffAfter0;
	Int_t TimeDiffBefore1;
	Int_t TimeDiffAfter1;
};

/*
struct Event {
	Event(Short_t adc, UChar_t det, ULong_t time, Bool_t pileup, Bool_t saturation) :
		Adc(adc), Det(det), Time(time), Pileup(pileup), Saturation(saturation),
		TimeDiffBefore0(MaxTime), TimeDiffAfter0(MaxTime),
		TimeDiffBefore1(MaxTime), TimeDiffAfter1(MaxTime),
		EnergyDepBefore0(0), EnergyDepAfter0(0),
		EnergyDepBefore1(0), EnergyDepAfter1(0)
		{}
	Short_t Adc;
	UChar_t Det;
	ULong_t Time;
	Bool_t Pileup;
	Bool_t Saturation;
	// --- Defined ĺater
	Short_t EnergyDepBefore0;
	Short_t EnergyDepAfter0;
	Short_t EnergyDepBefore1;
	Short_t EnergyDepAfter1;
	Int_t TimeDiffBefore0;
	Int_t TimeDiffAfter0;
	Int_t TimeDiffBefore1;
	Int_t TimeDiffAfter1;
}; 
*/

vector<Event>* TreeEntries(TFile* File) // Als CPP-Referenz zurück geben (siehe Wikipedia)
{
	//if(!File || File->IsZombie()) { /*throw error*/ }
	cout << "*** Opening TTree from file: " << File->GetName() << " ***" << endl;
	TTree* Tree = (TTree*)File->Get("Data");
	if (!Tree) {
		cout << "  Error: Failed to retrieve the TTree from the input file." << endl;
		File->Close();
		return NULL;
	}

	// Set branches (variables)
 	Short_t adc;
    UShort_t det;
    ULong64_t time;
    Bool_t pileup;
    Bool_t saturation;
    Tree->SetBranchAddress("adc", &adc);
    Tree->SetBranchAddress("det", &det);
    Tree->SetBranchAddress("Time", &time);
    Tree->SetBranchAddress("pileup", &pileup);
    Tree->SetBranchAddress("saturation", &saturation);

	// 
	ULong64_t totalEntries   = Tree->GetEntries();
	ULong64_t percentEntries = static_cast<long>(totalEntries)/100;
	int k = 1;
	cout << "Total number of events: " << totalEntries << endl;
	
	//
	cout << "Start building of event vector ..." << endl;
	vector<Event>* TreeEvents = new vector<Event>();
	//TreeEvents->reserve(totalEntries); // << memory not enough?
	
	for (ULong64_t currentEntryNr = 0; currentEntryNr < totalEntries; currentEntryNr++) {
		Tree->GetEntry(currentEntryNr);
		if (!pileup && !saturation){
			//Event currentEvent(adc, det, time, pileup, saturation);			
			TreeEvents->push_back(Event(adc, static_cast<char>(det), static_cast<unsigned long>(time)));
		//	cout << "Current vector size: " << TreeEvents->size() << endl;
		}
		if (currentEntryNr == k*percentEntries){
			cout << k << "% of total events finished, current capacity of 'TreeEvents': " << TreeEvents->capacity() << endl;
			k++;
		}
	}
	return TreeEvents; 
}

/*
vector<Event*>* TreeEntries(
	TFile* File
){
	vector<Event*>* TreeEvents = new vector<Event*>();
	ULong64_t MaximumSize = TreeEvents->max_size();
	cout << "Maximum size of event vector: " << MaximumSize << endl;
	TTree* Tree = (TTree*)File->Get("Data");
	if (!Tree) {
		cout << "Error: Failed to retrieve the TTree from the input file." << endl;
		File->Close();
		return TreeEvents;
	}

	//Get Values from existing Branches
 	Short_t adc;
    UShort_t det;
    ULong64_t time;
    Bool_t pileup;
    Bool_t saturation;
    Tree->SetBranchAddress("adc", &adc);
    Tree->SetBranchAddress("det", &det);
    Tree->SetBranchAddress("Time", &time);
    Tree->SetBranchAddress("pileup", &pileup);
    Tree->SetBranchAddress("saturation", &saturation);

	//Write relevant entries in Event vector
	ULong64_t EntryNumber = Tree->GetEntries();
	ULong64_t PercentEntries = (long)EntryNumber / 100;
	int k = 1;
	cout << "Total number of events: " << EntryNumber << endl;
	cout << "Start building of event vector ..." << endl;
	
	rlimit RLimits;
	
	for (ULong64_t i = 0; i < EntryNumber; i++) {
		Tree->GetEntry(i);
		if (!pileup && !saturation){
			Event* CurrentEvent = (Event*)malloc(sizeof(StdEvent));
			CurrentEvent->Adc = adc;
			CurrentEvent->Det = det;
			CurrentEvent->Time = time;
			CurrentEvent->Pileup = pileup;
			CurrentEvent->Saturation = saturation;
			TreeEvents->push_back(CurrentEvent);
		}
		if (i == k*PercentEntries){
			int test = getrlimit(RLIMIT_CORE, &RLimits);
//			cout << k << "% of total events finished, vector size used: " << TreeEvents->size() << "/" << MaximumSize << endl;
			cout << k << "% of total events finished, Usable RAM: soft: " << RLimits.rlim_cur << ", hard: " << RLimits.rlim_max << endl;
			k++;
		}
	}
	return TreeEvents; 
}
*/

bool OrderFunction(Event A, Event B){
	return A.Time < B.Time;
};

/*
vector<Event>* AddTimeDifferences(
	vector<Event>* GoodEvents
){
	//Sorting Process
	cout << "Sorting events with ascending time order ..." << endl;
	sort(GoodEvents->begin(), GoodEvents->end(), OrderFunction);
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
			if (GoodEvents->at(i).Time - LastEventIn0->Time < MaxTime) {GoodEvents->at(i).TimeDiffBefore0 = GoodEvents->at(i).Time - LastEventIn0->Time;}
			GoodEvents->at(i).EnergyDepBefore0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			if (GoodEvents->at(i).Time - LastEventIn1->Time < MaxTime) {GoodEvents->at(i).TimeDiffBefore1 = GoodEvents->at(i).Time - LastEventIn1->Time;}
			GoodEvents->at(i).EnergyDepBefore1 = LastEventIn1->Adc;
			}
		if(GoodEvents->at(i).Det == 0){
			LastEventIn0 = &GoodEvents->at(i);
		}
		if(GoodEvents->at(i).Det == 1){
			LastEventIn1 = &GoodEvents->at(i);
		}
		i++;
	} while (i < GoodEvents->size());									// && (LastEventIn0==NULL || LastEventIn1==NULL)
	
	//Just a test of whether I did something wrong
	if(LastEventIn0 != NULL && LastEventIn1 != NULL){
		cout << "Two Previous events are found" << endl; 
	}
	
	else {
		cout << "Error with finding two previous events: INTERRUPT" << endl;
		return GoodEvents;
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
	i = GoodEvents->size();
	cout << "Start determining TimeDiffAfter and EnergyDepAfter" << endl;

	do
	{
		if(LastEventIn0 != NULL) {
			if (LastEventIn0->Time - GoodEvents->at(i-1).Time < MaxTime) {GoodEvents->at(i-1).TimeDiffAfter0 = LastEventIn0->Time - GoodEvents->at(i-1).Time;}
			GoodEvents->at(i-1).EnergyDepAfter0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			if (LastEventIn1->Time - GoodEvents->at(i-1).Time < MaxTime) {GoodEvents->at(i-1).TimeDiffAfter1 = LastEventIn1->Time - GoodEvents->at(i-1).Time;}
			GoodEvents->at(i-1).EnergyDepAfter1 = LastEventIn1->Adc;
			}
		if(GoodEvents->at(i-1).Det == 0){
			LastEventIn0 = &GoodEvents->at(i-1);
		}
		if(GoodEvents->at(i-1).Det == 1){
			LastEventIn1 = &GoodEvents->at(i-1);
		}
		i--;
	} while (i >= 1);													// && (LastEventIn0==NULL || LastEventIn1==NULL)
	
	
	//Just a test of whether I did something wrong
	if(LastEventIn0 != NULL && LastEventIn1 != NULL){
		cout << "Two following events are found" << endl; 
	}
	else {
		cout << "Error with finding two following events: INTERRUPT" << endl;
		return GoodEvents;
	}
	//TimeDiffBefore and EnergyDepBefore with definite previous event
	for (i; i >= 1; i--){
		GoodEvents[i-1]->TimeDiffAfter0 = LastEventIn0->Time - GoodEvents[i-1]->Time;
		GoodEvents[i-1]->EnergyDepAfter0 = LastEventIn0->Adc;
		GoodEvents[i-1]->TimeDiffAfter1 = LastEventIn1->Time - GoodEvents[i-1]->Time;
		GoodEvents[i-1]->EnergyDepAfter1 = LastEventIn1->Adc;
		if(GoodEvents[i-1]->Det == 0){
			LastEventIn0 = GoodEvents[i-1];
		}
		if(GoodEvents[i-1]->Det == 1){
			LastEventIn1 = GoodEvents[i-1];		
		}
	}
		
	
	cout << "TimeDiffAfter and EnergyDepAfter calculated for every event" << endl;

	return GoodEvents;
}
*/

TFile* CreateNewFile(
	vector<Event>* GoodEvents,
	TString FileName,
	TString Run,
	TString PrefixPath = "../Test_230526/" 
){
	//Create new root file with Tree DataNew
	TString Folder = "/root/";
	TString PathWithFile = PrefixPath + Run + Folder + FileName + ".root";
	
	cout << "Start creation of new root file ..." << endl; 
	TFile* File = new TFile(PathWithFile, "RECREATE");
	TTree* DataNew = new TTree("Data", "TU5TU4Data");
    cout << "File and Tree created" << endl;
	//Variable & Branch initiation
	NewEvent CurrentEvent(0, 0, 0);
	
	DataNew->Branch("det",    &CurrentEvent.Det ,   "det/b");
    DataNew->Branch("adc",    &CurrentEvent.Adc,    "ch/S");
    DataNew->Branch("time",	  &CurrentEvent.Time,	"time/l");
	DataNew->Branch("TimeDiff_before0", &CurrentEvent.TimeDiffBefore0, "TimeDiff_before0/I");
	DataNew->Branch("TimeDiff_after0", &CurrentEvent.TimeDiffAfter0, "TimeDiff_after0/I");
	DataNew->Branch("TimeDiff_before1", &CurrentEvent.TimeDiffBefore1, "TimeDiff_before1/I");
	DataNew->Branch("TimeDiff_after1", &CurrentEvent.TimeDiffAfter1, "TimeDiff_after1/I");
	DataNew->Branch("EnergyDep_before0", &CurrentEvent.EnergyDepBefore0, "EnergyDep_before0/S");
	DataNew->Branch("EnergyDep_after0", &CurrentEvent.EnergyDepAfter0, "EnergyDep_after0/S");
	DataNew->Branch("EnergyDep_before1", &CurrentEvent.EnergyDepBefore1, "EnergyDep_before1/S");
	DataNew->Branch("EnergyDep_after1", &CurrentEvent.EnergyDepAfter1, "EnergyDep_after1/S");
	cout << "Branches created" << endl;
	
	//Create Empty Pointers for the Last Events
	Event* LastEventIn0 = NULL;
	Event* LastEventIn1 = NULL;
	
	//Sort The GoodEvents Vector
	cout << "Sorting events with ascending time order ..." << endl;
	sort(GoodEvents->begin(), GoodEvents->end(), OrderFunction);
	cout << "Finished sorting events" << endl;
	
	//Iterate over vector and write Events on the Tree
	ULong64_t totalEntries   = GoodEvents->size();
	ULong64_t percentEntries = static_cast<long>(totalEntries/100);
	int k = 1;
	cout << "Total number of events: " << totalEntries << endl;
	
	cout << "Start writing events into the new tree ..." << endl;
	for (ULong64_t CurrentEventIndex = 0; CurrentEventIndex < GoodEvents->size(); CurrentEventIndex++){
		CurrentEvent.Adc = GoodEvents->at(CurrentEventIndex).Adc;
		CurrentEvent.Det = GoodEvents->at(CurrentEventIndex).Det;
		CurrentEvent.Time = GoodEvents->at(CurrentEventIndex).Time;
			
		// Write TimeDiffBefore0 and EnergyDepBefore0 if necessary
		if (LastEventIn0 != NULL && CurrentEvent.Time - LastEventIn0->Time >= MaxTime) {
			CurrentEvent.EnergyDepBefore0 = LastEventIn0->Adc;
		}
		if (LastEventIn0 != NULL && CurrentEvent.Time - LastEventIn0->Time < MaxTime) {
			CurrentEvent.TimeDiffBefore0 = CurrentEvent.Time - LastEventIn0->Time;
			CurrentEvent.EnergyDepBefore0 = LastEventIn0->Adc;
		}
		
		// Write TimeDiffBefore1 and EnergyDepBefore1 if necessary
		if (LastEventIn1 != NULL && CurrentEvent.Time - LastEventIn1->Time >= MaxTime) {
			CurrentEvent.EnergyDepBefore1 = LastEventIn1->Adc;
		}
		if (LastEventIn1 != NULL && CurrentEvent.Time - LastEventIn1->Time < MaxTime) {
			CurrentEvent.TimeDiffBefore1 = CurrentEvent.Time - LastEventIn1->Time;
			CurrentEvent.EnergyDepBefore1 = LastEventIn1->Adc;
		}
		// Rewrite LastEventInX pointer to current Event
		if (CurrentEvent.Det == 0) {LastEventIn0 = &(GoodEvents->at(CurrentEventIndex)); }
		if (CurrentEvent.Det == 1) {LastEventIn1 = &(GoodEvents->at(CurrentEventIndex)); }
		

		//Loop for NextEventIn0:
		ULong64_t NextEventIn0Index = CurrentEventIndex+1;
		for (NextEventIn0Index; NextEventIn0Index < GoodEvents->size(); NextEventIn0Index++) {
			if (GoodEvents->at(NextEventIn0Index).Det == 0) {
				if (GoodEvents->at(NextEventIn0Index).Time - CurrentEvent.Time < MaxTime){
					CurrentEvent.TimeDiffAfter0 = GoodEvents->at(NextEventIn0Index).Time - CurrentEvent.Time;
				}
				CurrentEvent.EnergyDepAfter0 = GoodEvents->at(NextEventIn0Index).Adc;
				break;
			}
		}
		
		
		//Loop for NextEventIn1:
		ULong64_t NextEventIn1Index = CurrentEventIndex+1;
		for (NextEventIn1Index; NextEventIn1Index < GoodEvents->size(); NextEventIn1Index++) {
			if (GoodEvents->at(NextEventIn1Index).Det == 1) {
				if (GoodEvents->at(NextEventIn1Index).Time - CurrentEvent.Time < MaxTime){
					CurrentEvent.TimeDiffAfter1 = GoodEvents->at(NextEventIn1Index).Time - CurrentEvent.Time;
				}
				CurrentEvent.EnergyDepAfter1 = GoodEvents->at(NextEventIn1Index).Adc;
				break;
			}
		}
		DataNew->Fill();
		if (CurrentEventIndex == k*percentEntries){
			cout << k << "% of GoodEvents finished " << endl;
			k++;
		}
	}

	DataNew->Write();
	cout << "All Events are written into tree!" << endl;
	return File;
}


void MakeTimeDiffBranches_230718(){
	TString Run = "run004";
	TStopwatch* Stopwatch = new TStopwatch();
	Stopwatch->Start(true);
	TFile* OldFile = File(Run);
	if (!OldFile || OldFile->IsZombie()) {return; } 
	//Do all functions implemented above
	vector<Event>* GoodEvents = TreeEntries(OldFile);
	OldFile->Close();
	cout << "Vector filled, old file closed." << endl;
	//GoodEvents = AddTimeDifferences(GoodEvents);
	
	TFile* NewFile = CreateNewFile(GoodEvents, "Data_"+Run+"_coincidence", Run);
	
	NewFile->Write();
	cout << "New root file saved" << endl;
	Stopwatch->Stop();
	cout << "Total time = " << Stopwatch->RealTime() << " *** CPU time = " << Stopwatch->CpuTime() << endl;
}
