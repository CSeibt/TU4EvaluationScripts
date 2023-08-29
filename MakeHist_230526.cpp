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

typedef vector<Int_t> EThr; 

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
	
/*/
List of Branches (Histograms)
* det
* adc
* (pileup)      // commented out
* (saturation)  // no pileup or saturation branch created by MakeTimeDiffTree
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
	Double_t timediffthr = 5.0E6											//time threshold for time difference/coincidence investigation
){
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
			parameters = parameters+Form(" && ((TimeDiff_before%i<%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<%f && EnergyDep_after%i>%f))",
			 det2, timediffthr, det2, energythr, det2, timediffthr, det2, energythr);
		}
		else {
			parameters = parameters+Form(" && ((TimeDiff_before%i<%f) || (TimeDiff_after%i<%f))",
			 det2, timediffthr, det2, timediffthr);
		}
	}
	return parameters;
}


TString CoincidenceOptString(
	Int_t Det,															//Main Detector of interest
	Int_t Det2,															//Coincidence detector
	Int_t TimeDiffThr,												//Time Difference threshold. For TU4-TU5-coincidence approx 5000 ns
	vector<EThr> EnergyThresholds								//Vector of 2d arrays according for the energy window the coincidence event must have
){
	TString Parameters = Form("det==%i && pileup==0 %% saturation==0", Det);
	Parameters = Parameters + Form(" && ((TimeDiff_before%i<%i && (", Det2, TimeDiffThr);
	for (int i = 0; i < EnergyThresholds.size(); i++){
		if (i == 0){
			Parameters = Parameters + Form("(EnergyDep_before%i>%i && EnergyDep_before%i<%i)", Det2, EnergyThresholds[i][0], Det2, EnergyThresholds[i][1]);
		}
		else {
			Parameters = Parameters + Form(" || (EnergyDep_before%i>%i && EnergyDep_before%i<%i)", Det2, EnergyThresholds[i][0], Det2, EnergyThresholds[i][1]);
		}
	}
	Parameters = Parameters + Form(")) || (TimeDiff_after%i<%i && (", Det2, TimeDiffThr);
	for (int i = 0; i < EnergyThresholds.size(); i++){
		if (i == 0){
			Parameters = Parameters + Form("(EnergyDep_after%i>%i && EnergyDep_after%i<%i)", Det2, EnergyThresholds[i][0], Det2, EnergyThresholds[i][1]);
		}
		else {
			Parameters = Parameters + Form(" || (EnergyDep_after%i>%i && EnergyDep_after%i<%i)", Det2, EnergyThresholds[i][0], Det2, EnergyThresholds[i][1]);
		}
	}
	Parameters = Parameters + ")))";
	cout << Parameters << endl;
	return Parameters;
}


void GetAxisBorderAndSteps(Double_t xmax, Double_t &tborder, Int_t &tsteps)
{   // Determine axis upper limit and bin number from given maximum entry
    Double_t dig = pow(10, (int)(log(xmax)/log(10)));					//determine magnitude of tmax (10^x)
    tborder = dig * (1+(int)(xmax/dig));						//determine histogram border (tborder=a*10^x)
    tsteps = int(tborder);										//Set bin number = tborder -> 1 bin = 1 s
    // If bin number is over 20000, bin number gets reduced to 2000 < bin number <= 20000 (1 bin = 1, 10,... s)
    while(tsteps > 20000){tsteps = tsteps/10;}
}

TH1D* Histogram(
	TTree* data,
	TString histtype,
	TString histname,
	TString attributes,
	Int_t det = 0,														//det and det2 only necessary if histtype = TimeDiff
	Int_t det2 = 1
){
// Create a histogram of one of the types {"adc", "Time", "TimeDiff"}
	vector<TString> detectors = {"TU5","TU4"};    				//Detector names
	
	Double_t nch = 16384;
	
	if (histtype=="adc"){
		TH1D* histo = new TH1D(histname, detectors[det]+" ADC spectrum; Channel; Counts", nch, 0, nch);
		data->Draw("adc>>"+histname, attributes, "goff");
		return histo;
	}
	else if (histtype=="time"){
        // Determine time axis border and bin number using last time entry
		/*TLeaf* leaf = data->GetBranch("Time")->GetLeaf("Time");				//Get Time Leaf
		if (!leaf) cout << "Could not get leaf Time" << endl;
	    data->GetEntries("Time");
		Double_t tmax = (int)(leaf->GetValue(0) * (1.1/1E12));				//Measuring time times 1.1 in seconds = upper limit of histograms
        Double_t tborder = 1;
        Int_t tsteps = 1;
        GetAxisBorderAndSteps(tmax, tborder, tsteps);
		*/
		//Fill Histograms with data
		TH1D* histo = new TH1D(histname, detectors[det]+" signal rate; time / s; Counts", 90000, 0, 90000);							//Rate Histogram
		data->Draw("time/1E12>>"+histname, attributes, "goff");		//Time - Events of det
		return histo;
	}
	else if (histtype=="TimeDiff"){
		TH1D* histo = new TH1D(histname, "t("+detectors[det]+")-t("+detectors[det2]+"); dt / ps; Counts", 2001, -4E7-2000, 4E7+2000);
		data->Draw(Form("TimeDiff_before%i>>"+histname, det2), attributes, "goff");
		data->Draw(Form("-TimeDiff_after%i>>+"+histname, det2), attributes, "goff");
		return histo;
	}
	else {
		cout << "Wrong histtype!" << endl;
		return 0;
	}
}

TH1D* HistCanvas(
	TTree* data,
	TString histname,
	TString histtype,
	TString attributes,
	Int_t det0 = 0,
	Int_t det1 = 1,
	bool Write = false
){
// Create a histogram of given options and draw it into a canvas.
	TCanvas* c1 = new TCanvas(histname, histname, 100, 10, 1200, 700);
	gPad->SetLogy(true);
	TH1D* histogram = Histogram(data, histtype, histname, attributes, det0, det1);
	histogram->SetStats(0);
	histogram->Draw();
	if (Write) histogram->Write();
	return histogram;
}

void MakeHist_230526(){
	//General cosmetics
	gStyle->SetOptStat(1001111);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);

    //Variables
    TString run = "run004";			//run
    Int_t TU5 = 0;					//TU5 (X-ray detector)
    Int_t TU4 = 1;					//TU4 (Ge detector)
    Int_t thrTimediff = 5.0E6;	//time difference threshold
    Int_t thrEnergy = 1;			//energy threshold (pile-up)
    Int_t nch = 16384;
    
    
	TFile* file = GetInputFile(run);											//Open File 
	TTree* data = (TTree*)file->Get("Data");							//Get data from file
	
	vector<bool> attributes = {false, false, false, false, false};
    
	
	TString options1 = OptString(0, attributes);
	TString options2 = OptString(1, attributes);
	
	vector<bool> CoinAttributes = {false, false, false, true, false};
	TString TU5CoinOptions = OptString(TU5, CoinAttributes, TU4);
	TString TU4CoinOptions = OptString(TU4, CoinAttributes, TU5);
	 

	
	TH1D* histogram_adc1 = HistCanvas(data, "TU5adc", "adc", options1, TU5, TU4, true);
	TH1D* histogram_time1 = HistCanvas(data, "histo_rate_TU5", "time", options1, TU5, TU4, true);
	
	TH1D* histogram_adc2 = HistCanvas(data, "TU4adc", "adc", options2, TU4, TU5, true);
	TH1D* histogram_time2 = HistCanvas(data, "histo_rate_TU4", "time", options2, TU4, TU5, true);
	
	TH1D* TimeDiff01 = HistCanvas(data, "TU5TU4TimeDiff", "TimeDiff", options1, TU5, TU4, true);
	TH1D* TimeDiff02 = HistCanvas(data, "TU4TU5TimeDiff", "TimeDiff", options2, TU4, TU5, true);
	
	
	//TH1D* CoincidenceTU5 = HistCanvas(data, "TU5adcCoincidence", "adc", TU5CoinOptions, TU5, TU4, true);
	//TH1D* CoincidenceTU4 = HistCanvas(data, "TU4adcCoincidence", "adc", TU4CoinOptions, TU4, TU5, true);
	
	
	EThr Xrays1 = {417, 597};
	EThr Xrays2 = {2700, 3173};
	EThr Xrays3 = {3489, 3667};
	
	vector<EThr> Ba133Xrays = {Xrays1};
	
	//TString XrayCoinOptions = CoincidenceOptString(TU4, TU5, thrTimediff, Ba133Xrays);
	//TH1D* XrayCoincidenceTU4 = HistCanvas(data, "TU4_BaXrayCoincidence", "adc", XrayCoinOptions, TU4, TU5);
	//XrayCoincidenceTU4->Write();
	
}
