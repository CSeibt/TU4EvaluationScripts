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
TFile* GetInputFile(
    TString run = "run001"
){
	// Open ROOT file
	TString input = "../EvaluationData/Data_"+run+".root";
	TFile* rootFile = TFile::Open(input.Data(), "READ");
	if (!rootFile) cout << "Could not open " << input << endl;
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
	Double_t timediffthr = 1.E4											//time threshold for time difference/coincidence investigation
){
	TString parameters = Form("det==%i", det);
//	if (attributes[0]){
//		parameters = parameters+Form(" && pileup==%i", pileup);
//	}
//	if (attributes[1]){
//		parameters = parameters+Form(" && saturation==%i", saturation);	
//	}
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

Double_t GetAxisBorder(TTree* data, TString columname)
{   // Determine axis upper limit from maximum tree entry
	Double_t xmax = data->GetMaximum(columname.Data());
	Double_t dig = pow(10, (int)(log(xmax)/log(10)));	// determine magnitude of xmax (10^x)
	Double_t tborder = dig * (1+(int)(xmax/dig));		// determine histogram border (tborder=a*10^x)
	return tborder;
}

TH1D* Histogram(
	TTree* data,
	TString histname,
	TString histtype,
	TString attributes,
	Int_t det = 0,														//det and det2 only necessary if histtype = TimeDiff
	Int_t det2 = 1
){
// Create a histogram of one of the types {"adc", "Edep", "Time", "TimeDiff"}
	vector<TString> detectors = {"TU5","TU4"};    				//Detector names
	Double_t Nch = 16000;
	Double_t Emax[] = {200, 6000};

	if (histtype=="TimeDiff"){
		TH1D* histo = new TH1D(histname, "t("+detectors[det]+")-t("+detectors[det2]+"); dt / ps; Counts", 2001, -4E7-2000, 4E7+2000);
		data->Draw(Form("TimeDiff_before%i>>"+histname, det2), attributes, "goff");
		data->Draw(Form("-TimeDiff_after%i>>+"+histname, det2), attributes, "goff");
		return histo;
	}	
	if (histtype=="adc"){
		// Determine axis border and bin number using maximum entry
		Double_t border = GetAxisBorder(data, histtype);
		Int_t steps = int(border);
		// Fill histogram with data
		TH1D* histo = new TH1D(histname, detectors[det]+" ADC spectrum; Channel; Counts", Nch, 0, Nch);
		data->Draw("adc>>"+histname, attributes, "goff");
		return histo;
	}
	else if (histtype=="Edep"){
		// Determine axis border and bin number using maximum entry
		Double_t border = Emax[det];
		Int_t steps = int(border);
		while (steps < 1000) {steps = steps*10;}
		// Fill histogram with data
		TH1D* histo = new TH1D(histname, detectors[det]+" energy spectrum; E/keV; Counts", steps, 0, border);
		data->Draw("Edep>>"+histname, attributes, "goff");
		return histo;
	}
	else if (histtype=="Time"){
		// Determine axis border and bin number using maximum entry
		Double_t tborder = GetAxisBorder(data, histtype) * 1.E-12;
		Int_t tsteps = int(tborder);										//Set bin number = tborder -> 1 bin = 1 s
		// If bin number is over 20000, bin number gets reduced to 2000 < bin number <= 20000 (1 bin = 1, 10,... s)
		while(tsteps > 20000) {tsteps = tsteps/10;}

		//Fill Histogram with data
		TH1D* histo = new TH1D(histname, detectors[det] + " signal rate; Time / s; Counts", tsteps, 0, tborder);							//Rate Histogram
		data->Draw("Time/1E12>>"+histname, attributes, "goff");		//Time - Events of det
		return histo;
	}
	else {
		cout << "Wrong histtype: " << histtype << endl;
		return 0;
	}
}

TH1D* HistCanvas(
	TTree* data,
	TString histname,
	TString histtype,
	TString attributes,
	Int_t det0 = 0,
	Int_t det1 = 1
){
// Create a histogram of given options and draw it into a canvas.
	TCanvas* c1 = new TCanvas(histname, histname, 100, 10, 1200, 700);
	gPad->SetLogy(true);
	TH1D* histogram = Histogram(data, histname, histtype, attributes, det0, det1);
	histogram->SetStats(0);
	histogram->Draw();
	
	return histogram;
}

// TODO: 
Double_t GetLiveTime(TTree* data)
{ // Calculate the live time of a given dataset. Get deadtime per event by showing to the user the TimeDiff hist, get real time from time branch, use total number of events. 
	return 1;
}

void MakeHist_230526(){
	//General cosmetics
	gStyle->SetOptStat(1001111);
	gStyle->SetOptFit(0);
	gStyle->SetStripDecimals(kFALSE);

	//Variables
	TString run = "run038";			//run
	Int_t TU5 = 0;					//TU5 (X-ray detector)
	Double_t thrTimediff = 1.E4;	//time difference threshold
	Double_t thrEnergy = 1;			//energy threshold (pile-up)
	Int_t nch = 16384;
    
    
	TFile* file = GetInputFile(run);											//Open File 
	TTree* data = (TTree*)file->Get("Data");							//Get data from file
	
	vector<bool> attributes = {false, false, false, false, false};

	
	TString options1 = OptString(0, attributes);


	//TH1D* histogram_adc1 = HistCanvas(data, "TU5adc", "adc", options1, TU5);
	TH1D* histogram_time1 = HistCanvas(data, "histo_rate_TU5", "Time", options1);
	TH1D* histogram_E1 = HistCanvas(data, "TU5Edep", "Edep", options1);

	/*/ Checks with manual selection string
	TString man_sel = "";
	// TU5 spectrum of coincident events 
	man_sel = "det==0 && (TimeDiff_after1<0.5E6 || TimeDiff_before1<4.E6)";
	TH1D* coincident_E1 = HistCanvas(data, "CoincTU5Edep", "Edep", man_sel, TU5, 0);
	//*/
    

}
