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



//Open the root file containing the respective histograms
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
	
/*/
List of Branches (Histograms)
* det
* adc
* pileup
* saturation
* ADCout
* Time
* TimeDiff_before#
* TimeDiff_after#
* EnergyDep_before#
* EnergyDep_after#
/*/	

TString OptString(
	Int_t det,
	vector<bool> attributes = {false, false, false, false, false},
	//						   Pile-Up, Saturation, Energy, TimeDiff, EnergyDep
	Int_t det2 = 0,														//Detector 2 for time difference/coincidence investigation
	Int_t pileup = 0,													//pileup = 0 -> non-pile-up-events, pileup = 1 -> pile-up-events
	Int_t saturation = 0,												//saturation = 0 -> non-saturation-events, saturation = 1 -> saturation-events
	Double_t energythr = 1,												//energy threshold in channels
	Double_t timediffthr = 1.E4											//time threshold for time difference/coincidence investigation
){
	const int numchannels = 2; 
	vector<TString> detectors = {"TU5","TU4"};    				//Detector names
	
	
	TString parameters = Form("det==%i", det);
	if (attributes[0]){
		parameters = parameters+Form(" && pileup==%i", pileup);
	}
	if (attributes[1]){
		parameters = parameters+Form(" && saturation==%i", saturation);	
	}
	if (attributes[2]){
		parameters = parameters+Form(" && adc==%f", energythr);
	}
	if (attributes[3]){
		if (attributes[4]){
			parameters = parameters+Form(" && ((TimeDiffBefore%i<%f && EnergyDepBefore%i>%f) || (TimeDiffAfter%i<%f && EnergyDepAfter%i>%f))",
			 det2, timediffthr, det2, energythr, det2, timediffthr, det2, energythr);
		}
		else {
			parameters = parameters+Form(" && ((TimeDiffBefore%i<%f) || (TimeDiffAfter%i<%f))",
			 det2, timediffthr, det2, timediffthr);
		}
	}
	return parameters;
}
	
TH1D* Histogram(
	TTree* data,
	TString histtype,
	TString histname,
	TString attributes,
	Int_t det = 0,														//det and det2 only necessary if histtype = TimeDiff
	Int_t det2 = 1
){
	const int numchannels = 2; 
	vector<TString> detectors = {"TU5","TU4"};    				//Detector names
	
	Double_t nch = 16384;
	
	if (histtype=="adc"){
		TH1D* histo = new TH1D(histname, detectors[det]+" ADC spectrum; Channel; Counts", nch, 0, nch);
		data->Draw("adc>>"+histname, attributes, "goff");
		return histo;
	}
	else if (histtype=="Time"){
		TLeaf* leaf = data->GetBranch("Time")->GetLeaf("Time");				//Get Time Leaf
		if (!leaf) cout << "Could not get leaf Time" << endl;
	
		//cout << "Last entry: " << data->GetEntries("Time")-1 << endl;		//For some reason, if this line is missing, leaf->GetValue(0) turns to 0
		data->GetEntries("Time");
		//Ratemax is in thicks 1.0E12 == 1s
		Double_t tmax = (int)(leaf->GetValue(0) * (1.1/1E12));				//Measuring time times 1.1 in seconds = upper limit of histograms
		Double_t dig = pow(10, (int)(log(tmax)/log(10)));					//determine magnitude of tmax (10^x)
		Double_t tborder = dig * (1+(int)(tmax/dig));						//determine histogram border (tborder=a*10^x)
		Int_t tsteps = int(tborder);										//Set bin number = tborder -> 1 bin = 1 s
	
		//If bin number is over 20000, bin number gets reduced to 2000 < bin number <= 20000 (1 bin = 1, 10,... s)
		while(tsteps > 20000){tsteps = tsteps/10;}							
		//cout << tmax << " and " << tborder << " and " << tsteps << endl;
	
		//Fill Histograms with data
		TH1D* histo = new TH1D(histname, detectors[det]+" signal rate; Time / s; Counts", tsteps, 0, tborder);							//Rate Histogram
		data->Draw("Time/1E12>>"+histname, attributes, "goff");		//Time - Events of det
		return histo;
	}
	else if (histtype=="TimeDiff"){
		TH1D* histo = new TH1D(histname, "t("+detectors[det]+")-t("+detectors[det2]+"); dt / ps; Counts", 2001, -4E7-2000, 4E7+2000);
		data->Draw(Form("TimeDiffBefore%i>>"+histname, det2), attributes, "goff");
		data->Draw(Form("-TimeDiffAfter%i>>+"+histname, det2), attributes, "goff");
		return histo;
	}
	else {
		cout << "Wrong histtype!" << endl;
		return 0;
		
	}
}

void TimeDiff(
	TTree* Tree
){
	Long64_t TotalEntries = Tree->GetEntries();
	cout << "Total Entry number = " << TotalEntries << endl;
	TLeaf* Leaf = Tree->GetBranch("Time")->GetLeaf("Time");
	TLeaf* Pileup = Tree->GetBranch("pileup")->GetLeaf("pileup");
	TLeaf* Saturation = Tree->GetBranch("saturation")->GetLeaf("saturation");
	Double_t LastEntry = Leaf->GetValue(0);
	Double_t CurrentEntry = 0;
	for (Long64_t i = 1; i<= 10; i++){
		CurrentEntry = Leaf->GetValue(i);
		Bool_t CurrentPileUp = Pileup->GetValue(i);
		Bool_t CurrentSaturation = Saturation->GetValue(i);
		/*if (CurrentEntry < LastEntry){
			cout << "Time Order Violated!" << endl;
			}*/
		cout << "EntryIndex = " << i << ", CurrentTme = " << CurrentEntry << ", LastTime = " << LastEntry << ", Pileup Flag = " 
		<< CurrentPileUp << ", Saturation Flag = " << CurrentSaturation << endl;
		LastEntry = CurrentEntry;
		}
	cout << "Finished testing Time Ordering" << endl;
}


void AddTimeDiffBranches(
	TFile* file
){	
    // Access the TTree in the input file
    TTree* tree = dynamic_cast<TTree*>(file->Get("Data"));
    if (!tree) {
        std::cout << "Error: Failed to retrieve the TTree from the input file." << std::endl;
        file->Close();
        return;
    }

    // Create new branches for TimeDiffBefore and TimeDiffAfter
    Double_t TimeDiffAfter0 = 0.0;
    Double_t TimeDiffAfter1 = 0.0;    
    Double_t TimeDiffBefore0 = 0.0;
    Double_t TimeDiffBefore1 = 0.0;
    Short_t EnergyDepAfter0 = 0; 
    Short_t EnergyDepAfter1 = 0; 
    Short_t EnergyDepBefore0 = 0; 
    Short_t EnergyDepBefore1 = 0; 

    TBranch* b_TimeDiffAfter0 = tree->Branch("TimeDiffAfter0", &TimeDiffAfter0, "TimeDiffAfter0/D");
    TBranch* b_TimeDiffAfter1 = tree->Branch("TimeDiffAfter1", &TimeDiffAfter1, "TimeDiffAfter1/D");
    TBranch* b_TimeDiffBefore0 = tree->Branch("TimeDiffBefore0", &TimeDiffBefore0, "TimeDiffBefore0/D");
    TBranch* b_TimeDiffBefore1 = tree->Branch("TimeDiffBefore1", &TimeDiffBefore1, "TimeDiffBefore1/D");
    TBranch* b_EnergyDepAfter0 = tree->Branch("EnergyDepAfter0", &EnergyDepAfter0, "EnergyDepAfter0/S");
    TBranch* b_EnergyDepAfter1 = tree->Branch("EnergyDepAfter1", &EnergyDepAfter1, "EnergyDepAfter1/S");
    TBranch* b_EnergyDepBefore0 = tree->Branch("EnergyDepBefore0", &EnergyDepBefore0, "EnergyDepBefore0/S");
    TBranch* b_EnergyDepBefore1 = tree->Branch("EnergyDepBefore1", &EnergyDepBefore1, "EnergyDepBefore1/S");

    // Retrieve the values of the existing branches
    Short_t adc;
    UShort_t det;
    ULong64_t time;
    Bool_t pileup;
    Bool_t saturation;
    tree->SetBranchAddress("adc", &adc);
    tree->SetBranchAddress("det", &det);
    tree->SetBranchAddress("Time", &time);
    tree->SetBranchAddress("pileup", &pileup);
    tree->SetBranchAddress("saturation", &saturation);

    // Initialize variables to store the previous and next valid event times
    Long64_t prevValidTime0 = 0;
    Long64_t nextValidTime0 = 0;
    Long64_t prevValidTime1 = 0;
    Long64_t nextValidTime1 = 0;
    Short_t nextEnergyDep0 = 0;
    Short_t nextEnergyDep1 = 0;
    

    // Loop over the events in the TTree
    ULong64_t numEntries = tree->GetEntries();
    ULong64_t tenthEntries = (long)numEntries / 100;
    Int_t k = 1;
    cout << "Total number of events: " << numEntries << endl;
    cout << "1% = " << tenthEntries << endl;
    cout << "Start Creation of Time Diff and Energy Dep Branches ..." << endl;
    for (ULong64_t i = 0; i < numEntries; ++i) {
        //cout  << i << " out of " << numEntries << "gets calculated" << endl;
        if (i == tenthEntries*k){
			cout << k <<"0% of events done!" << endl;
			k++;}
        tree->GetEntry(i);

        // Find the next event with pile-up and saturation = false and det = 0
        nextValidTime0 = -1;
        nextEnergyDep0 = 0;
        for (ULong64_t j = i + 1; j < numEntries; ++j) {
            tree->GetEntry(j);
            if (!pileup && !saturation && det == 0) {
                nextValidTime0 = time;
                nextEnergyDep0 = adc;
                break;
            }
        }

		// Find the next event with pile-up and saturation = false and det = 1
        nextValidTime1 = -1;
        nextEnergyDep1 = 0;
        for (ULong64_t j = i + 1; j < numEntries; ++j) {
            tree->GetEntry(j);
            if (!pileup && !saturation && det == 1) {
                nextValidTime1 = time;
                nextEnergyDep1 = adc;
                break;
            }
        }
        
        // Find the previous event with pile-up and saturation = false and det = 0
        TimeDiffBefore0 = 0;
        EnergyDepBefore0 = 0;
        for (ULong64_t j = i + 1; j < numEntries; ++j) {
            tree->GetEntry(j);
            if (!pileup && !saturation && det == 0) {
                TimeDiffBefore0 = time - prevValidTime0;
                EnergyDepBefore0 = adc;
                break;
            }
        }
        
        // Find the previous event with pile-up and saturation = false and det = 1
        TimeDiffBefore1 = -1;
        EnergyDepBefore1 = 0;
        for (ULong64_t j = i + 1; j < numEntries; ++j) {
            tree->GetEntry(j);
            if (!pileup && !saturation && det == 1) {
                TimeDiffBefore1 = time - prevValidTime1;
                EnergyDepBefore1 = adc;
                break;
            }
        }
        
        
		
        // Calculate TimeDiffAfter0 and 1
        TimeDiffAfter0 = (nextValidTime0 > 0) ? nextValidTime0 - time : 0;
        TimeDiffAfter1 = (nextValidTime1 > 0) ? nextValidTime1 - time : 0;

		//Calculate the EnergyDepAfter0 and 1
		EnergyDepAfter0 = nextEnergyDep0;
		EnergyDepAfter1 = nextEnergyDep1;
		
        // Fill the branches with the calculated values
        b_TimeDiffAfter0->Fill();
        b_TimeDiffAfter1->Fill();
        b_TimeDiffBefore0->Fill();
        b_TimeDiffBefore1->Fill();
        b_EnergyDepAfter0->Fill();
        b_EnergyDepAfter1->Fill();
        b_EnergyDepBefore0->Fill();
        b_EnergyDepBefore1->Fill();

        // Update the previous valid event time
        if (!pileup && !saturation && det == 0)
            prevValidTime0 = time;
        if (!pileup && !saturation && det == 1)
            prevValidTime1 = time;    
        
    }
/*
    // Write the modified TTree to the output file
    TFile* outFile = TFile::Open(outputFile, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "Error: Failed to create the output file." << std::endl;
        file->Close();
        return;
    }
*/
    tree->Write();
//    outFile->Close();
    file->Close();

//    delete outFile;
    delete file;

    std::cout << "TimeDiffBefore and TimeDiffAfter branches added and saved to the output file." << std::endl;
}


TH1D* HistCanvas(
	TTree* data,
	TString run,
	TString histname,
	TString histtype,
	TString attributes,
	Int_t det0 = 0,
	Int_t det1 = 1
){
	TCanvas* c1 = new TCanvas(histname, histname, 100, 10, 1200, 700);
	gPad->SetLogy(true);
	TH1D* histogram = Histogram(data, histtype, histname, attributes, det0, det1);
	histogram->SetStats(0);
	histogram->Draw();
	
	return histogram;
}

void MakeHist_230526(){
	//General cosmetics
	gStyle->SetOptStat(1001111);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);

    //Variables
    TString run = "run001";			//run
    Int_t TU5 = 0;					//TU5 (X-ray detector)
    Int_t TU4 = 1;					//TU4 (Ge detector)
    Int_t Sz2 = 2;					//Muon Panel 2
    Double_t thrTimediff = 1.E4;	//time difference threshold
    Double_t thrEnergy = 1;			//energy threshold (pile-up)
    Int_t nch = 16384;
    
    
	TFile* file = File(run);											//Open File 
	TTree* data = (TTree*)file->Get("Data");							//Get data from file
	
	vector<bool> attributes = {true, true, false, false, false};

	
	TString options1 = OptString(0, attributes);
	TString options2 = OptString(1, attributes);
	
	
	
	cout << options1 << endl;
	/*
	TH1D* histogram1 = HistCanvas(data, run, "TU5adc", "adc", options1);
	TH1D* histogramtime1 = HistCanvas(data, run, "histo_rate_TU5", "Time", options1);
	
	TH1D* histogram2 = HistCanvas(data, run, "TU4adc", "adc", options2);
	TH1D* histogramtime2 = HistCanvas(data, run, "histo_rate_TU4", "Time", options2);
	
	histogram1->Write();
	histogramtime1->Write();
	
	histogram2->Write();
	histogramtime2->Write();
	*/
	TH1D* TimeDiff01 = HistCanvas(data, run, "TU5TU4TimeDiff", "TimeDiff", options1, TU5, TU4);
	//AddTimeDiffBranches(file);
}
