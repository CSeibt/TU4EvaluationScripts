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
    TFile* rootFile = new TFile("Data_"+run+"_coincidence.root","UPDATE");
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
	else if (histtype=="Time"){
        // Determine time axis border and bin number using last time entry
		TLeaf* leaf = data->GetBranch("Time")->GetLeaf("Time");				//Get Time Leaf
		if (!leaf) cout << "Could not get leaf Time" << endl;
	    data->GetEntries("Time");
		Double_t tmax = (int)(leaf->GetValue(0) * (1.1/1E12));				//Measuring time times 1.1 in seconds = upper limit of histograms
        Double_t tborder = 1;
        Int_t tsteps = 1;
        GetAxisBorderAndSteps(tmax, tborder, tsteps);

		//Fill Histograms with data
		TH1D* histo = new TH1D(histname, detectors[det]+" signal rate; Time / s; Counts", tsteps, 0, tborder);							//Rate Histogram
		data->Draw("Time/1E12>>"+histname, attributes, "goff");		//Time - Events of det
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
	Int_t det1 = 1
){
// Create a histogram of given options and draw it into a canvas.
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
    TString run = "run002";			//run
    Int_t TU5 = 0;					//TU5 (X-ray detector)
    Int_t TU4 = 1;					//TU4 (Ge detector)
    Double_t thrTimediff = 1.E4;	//time difference threshold
    Double_t thrEnergy = 1;			//energy threshold (pile-up)
    Int_t nch = 16384;
    
    
	TFile* file = GetInputFile(run);											//Open File 
	TTree* data = (TTree*)file->Get("Data");							//Get data from file
	
	vector<bool> attributes = {false, false, false, false, false};

	
	TString options1 = OptString(0, attributes);
	TString options2 = OptString(1, attributes);


	TH1D* histogram_adc1 = HistCanvas(data, "TU5adc", "adc", options1, TU5);
	//TH1D* histogram_time1 = HistCanvas(data, "histo_rate_TU5", "Time", options1);
	
	TH1D* histogram_adc2 = HistCanvas(data, "TU4adc", "adc", options2, TU4);
	//TH1D* histogram_time2 = HistCanvas(data, "histo_rate_TU4", "Time", options2);
	
	/*/ Time difference spectra
	TH1D* TimeDiff11 = Histogram(data, "TU5TU5TimeDiff", "TimeDiff", options1, TU5, TU5);
	TH1D* TimeDiff22 = Histogram(data, "TU4TU4TimeDiff", "TimeDiff", options2, TU4, TU4);
	TH1D* TimeDiff12 = Histogram(data, "TU5TU4TimeDiff", "TimeDiff", options1, TU5, TU4);
	TH1D* TimeDiff21 = Histogram(data, "TU4TU5TimeDiff", "TimeDiff", options2, TU4, TU5);
	TCanvas* cTimeDiff = new TCanvas("TimeDiff", "TimeDiff");
	cTimeDiff->Divide(2,2);
	cTimeDiff->cd(1);
	gPad->SetLogy();
	TimeDiff11->SetStats(0);
	TimeDiff11->Draw();
	cTimeDiff->cd(2);
	gPad->SetLogy();
	TimeDiff12->SetStats(0);
	TimeDiff12->Draw();
	cTimeDiff->cd(3);
	gPad->SetLogy();
	TimeDiff21->SetStats(0);
	TimeDiff21->Draw();
	cTimeDiff->cd(4);
	gPad->SetLogy();
	TimeDiff22->SetStats(0);
	TimeDiff22->Draw();
	//*/

    // Checks with manual selection string
    TString man_sel = "";
    // TU5 spectrum of coincident events 
    man_sel = "det==0 && (TimeDiff_after1<0.5E6 || TimeDiff_before1<4.E6)";
    TH1D* coincident_adc1 = HistCanvas(data, "CoincTU5adc", "adc", man_sel, TU5, 0);
    // TU4 spectrum of coincident events 
    man_sel = "det==1 && ((TimeDiff_after0<4.E6 && (EnergyDep_after0>2700 && EnergyDep_after0<3173)) || (TimeDiff_before0<0.5E6 && (EnergyDep_before0>2700 && EnergyDep_before0<3173)))"; // 2700-3173 & 3489-3667
    TH1D* coincident_adc2 = HistCanvas(data, "CoincTU4adc", "adc", man_sel, TU4, 0);
    
    /*/ Draw coincident spectra together with raw ones
    TCanvas* cAdc1 = new TCanvas("TU5adc", "TU5adc");
    gPad->SetLogy();
    histogram_adc1->Draw();
    histogram_adc1->GetXaxis()->SetRangeUser(0, 5000);
    coincident_adc1->SetLineColor(kRed);
    coincident_adc1->Draw("same");
    TLegend* l1 = new TLegend(0.6, 0.7, 0.85, 0.8);
    l1->AddEntry(histogram_adc1, "all events");
    l1->AddEntry(coincident_adc1, "coincident events");
    l1->Draw();
    
    TCanvas* cAdc2 = new TCanvas("TU4adc", "TU4adc");
    gPad->SetLogy();
    histogram_adc2->Draw();
    histogram_adc2->GetXaxis()->SetRangeUser(0, 10000);
    coincident_adc2->SetLineColor(kRed);
    coincident_adc2->Draw("same");
    TLegend* l2 = new TLegend(0.6, 0.7, 0.85, 0.8);
    l2->AddEntry(histogram_adc2, "all events");
    l2->AddEntry(coincident_adc2, "coincident events");
    l2->Draw();//*/
}
