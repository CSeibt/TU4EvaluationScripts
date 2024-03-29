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

const ULong_t MaxTime = kMaxULong;

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
	Event(Short_t adc, UChar_t det, ULong_t time, Bool_t pileup, Bool_t saturation) :
		Adc(adc), Det(det), Time(time), Pileup(pileup), Saturation(saturation),
		TimeDiffBefore0(MaxTime), TimeDiffAfter0(MaxTime),
		TimeDiffBefore1(MaxTime), TimeDiffAfter1(MaxTime),
		EnergyDepBefore0(0), EnergyDepAfter0(0),
		EnergyDepBefore1(0), EnergyDepAfter1(0)
		{}
	Short_t Adc;
	UChar_t Det;
	ULong64_t Time;
	Bool_t Pileup;
	Bool_t Saturation;
	// --- Defined ĺater
	Short_t EnergyDepBefore0;
	Short_t EnergyDepAfter0;
	Short_t EnergyDepBefore1;
	Short_t EnergyDepAfter1;
	Long64_t TimeDiffBefore0;
	Long64_t TimeDiffAfter0;
	Long64_t TimeDiffBefore1;
	Long64_t TimeDiffAfter1;
};

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
			TreeEvents->push_back(Event(adc, static_cast<char>(det), static_cast<unsigned long>(time), pileup, saturation));
		//	cout << "Current vector size: " << TreeEvents->size() << endl;
		}
		if (currentEntryNr == k*percentEntries){
			cout << k << "% of total events finished" << endl;
			k++;
		}
	}
	return TreeEvents; 
}



bool OrderFunction(Event A, Event B){
	return A.Time < B.Time;
};

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
			GoodEvents->at(i).TimeDiffBefore0 = GoodEvents->at(i).Time - LastEventIn0->Time;
			GoodEvents->at(i).EnergyDepBefore0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			GoodEvents->at(i).TimeDiffBefore1 = GoodEvents->at(i).Time - LastEventIn1->Time;
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

	cout << "TimeDifferenceBefore and EnergyDepBefore calculated for every event" << endl;

	//Calculating TimeDiffAfter and EnergyDepAfter by iterating reversed through tree
	LastEventIn0 = NULL;
	LastEventIn1 = NULL;
	i = GoodEvents->size();
	cout << "Start determining TimeDiffAfter and EnergyDepAfter" << endl;

	do
	{
		if(LastEventIn0 != NULL) {
			GoodEvents->at(i-1).TimeDiffAfter0 = LastEventIn0->Time - GoodEvents->at(i-1).Time;
			GoodEvents->at(i-1).EnergyDepAfter0 = LastEventIn0->Adc;
			}
		if(LastEventIn1 != NULL) {
			GoodEvents->at(i-1).TimeDiffAfter1 = LastEventIn1->Time - GoodEvents->at(i-1).Time;
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
	
	cout << "TimeDiffAfter and EnergyDepAfter calculated for every event" << endl;

	return GoodEvents;
}

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
	Event CurrentEvent(0, 0, 0, 0, 0);
	
	DataNew->Branch("det",   &CurrentEvent.Det ,   "det/b");
    DataNew->Branch("adc",    &CurrentEvent.Adc,   "ch/S");
    DataNew->Branch("pileup",&CurrentEvent.Pileup,"pileup/O");
    DataNew->Branch("saturation",&CurrentEvent.Saturation,"saturation/O");
    DataNew->Branch("time",&CurrentEvent.Time,"time/l");
	DataNew->Branch("TimeDiff_before0", &CurrentEvent.TimeDiffBefore0, "TimeDiff_before0/L");
	DataNew->Branch("TimeDiff_after0", &CurrentEvent.TimeDiffAfter0, "TimeDiff_after0/L");
	DataNew->Branch("TimeDiff_before1", &CurrentEvent.TimeDiffBefore1, "TimeDiff_before1/L");
	DataNew->Branch("TimeDiff_after1", &CurrentEvent.TimeDiffAfter1, "TimeDiff_after1/L");
	DataNew->Branch("EnergyDep_before0", &CurrentEvent.EnergyDepBefore0, "EnergyDep_before0/S");
	DataNew->Branch("EnergyDep_after0", &CurrentEvent.EnergyDepAfter0, "EnergyDep_after0/S");
	DataNew->Branch("EnergyDep_before1", &CurrentEvent.EnergyDepBefore1, "EnergyDep_before1/S");
	DataNew->Branch("EnergyDep_after1", &CurrentEvent.EnergyDepAfter1, "EnergyDep_after1/S");
	cout << "Branches created" << endl;
	
	//Iterate over vector and write Events on the Tree
	cout << "Start writing events into the new tree ..." << endl;
	for (ULong64_t i = 0; i < GoodEvents->size(); i++){
		CurrentEvent = GoodEvents->at(i);
		DataNew->Fill();
	}

	DataNew->Write();
	cout << "All Events are written into tree!" << endl;
	return File;
}


void MakeTimeDiffBranches(){
	TString Run = "run004";
	TStopwatch* Stopwatch = new TStopwatch();
	Stopwatch->Start(true);
	TFile* OldFile = File(Run);
	if (!OldFile || OldFile->IsZombie()) {return; } 
	//Do all functions implemented above
	vector<Event>* GoodEvents = TreeEntries(OldFile);
	OldFile->Close();
	
	GoodEvents = AddTimeDifferences(GoodEvents);
	
	TFile* NewFile = CreateNewFile(GoodEvents, "Data_"+Run+"_coincidence", Run);
	
	NewFile->Write();
	cout << "New root file saved" << endl;
	Stopwatch->Stop();
	cout << "Total time = " << Stopwatch->RealTime() << " *** CPU time = " << Stopwatch->CpuTime() << endl;
}
