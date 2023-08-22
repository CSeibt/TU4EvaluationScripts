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

// TODO: TH1D* GetSpectrumCoincidentWithEnergy(TTree* data, Int_t det0, Int_t det1, Double_t minTimeDiff, Double_t maxTimeDiff, Double_t range0low, Double_t range0up, range1low=0, Double_t range1up=0, range2low=0, Double_t range2up=0) {
// Loop through the tree.

TH2D* GetCoincidentEvsE(TTree* data, Int_t det0, Int_t det1, Double_t minTimeDiff, Double_t maxTimeDiff)
// Make a 2D histogram of all det0 events considered coincident with an event in det1. det0 adc channel on x-axis, det1 adc channel on y-axis. 
// Coincidence time interval: minTimeDiff < t(det0)-t(det1) < maxTimeDiff. 
{
	// create histogram
	TString histname = "name";
	TH2D* h = new TH2D(histname.Data(), "; TU4 E/keV; TU5 E/keV", 800, 0, 800, 1000, 0, 100);
	
	UShort_t det = 0;
	Double_t Edep = 0;
	Long64_t TimeDiffBefore = 0;
	Long64_t TimeDiffAfter = 0;
	Double_t EnergyDepBefore = 0;
	Double_t EnergyDepAfter = 0;
	data->SetBranchAddress("det", &det);
	data->SetBranchAddress("Edep", &Edep);
	data->SetBranchAddress(Form("TimeDiff_before%i", det1), &TimeDiffBefore);
	data->SetBranchAddress(Form("TimeDiff_after%i", det1), &TimeDiffAfter);
	data->SetBranchAddress(Form("EnergyDep_before%i", det1), &EnergyDepBefore);
	data->SetBranchAddress(Form("EnergyDep_after%i", det1), &EnergyDepAfter);
	
	// fill histogram
	for (Long64_t event = 0; event < data->GetEntries(); event++) {
		data->GetEvent(event);
		if (det != det0) continue;
		if (-TimeDiffAfter > minTimeDiff) {
			h->Fill(Edep, EnergyDepAfter); 
			//cout << "det" << det1 << " ch " << EnergyDepAfter << " after det" << det0 << " ch " << adc << endl;
		}
		if (TimeDiffBefore < maxTimeDiff) {
			h->Fill(Edep, EnergyDepBefore); 
			//cout << "det" << det1 << " ch " << EnergyDepAfter << " before det" << det0 << " ch " << adc << endl;
		}
		//cout << "det " << det << " adc " << adc << " dtBefore " << TimeDiffBefore << " dtAfter " << TimeDiffAfter << endl;
		if (-TimeDiffAfter > minTimeDiff && TimeDiffBefore < maxTimeDiff) cout << "Event " << event << " counted twice!" << endl;
	}
	//TString varexp = "EnergyDep_after" + to_string(det1) + ":adc";
	//TString selection = "det==" + to_string(det0) + "&&TimeDiff_after"+to_string(det1) + "<" + to_string(maxTimeDiff);
	//data->Draw(varexp.Data(), selection.Data());
	
	cout << h->Integral() << " histogram entries" << endl;
	
	return h;
}

void MakeHist(){
	//General cosmetics
	gStyle->SetOptStat(1001111);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);

    //Variables
    TString run = "run004";			//run
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


	//TH1D* histogram_adc1 = HistCanvas(data, "TU5adc", "adc", options1, TU5);
	TH1D* histogram_E1 = Histogram(data, "TU5Edep", "Edep", options1, TU5);
	//TH1D* histogram_time1 = HistCanvas(data, "histo_rate_TU5", "Time", options1);
	
	//TH1D* histogram_adc2 = HistCanvas(data, "TU4adc", "adc", options2, TU4);
	TH1D* histogram_E2 = Histogram(data, "TU4Edep", "Edep", options2, TU4);
	//TH1D* histogram_time2 = HistCanvas(data, "histo_rate_TU4", "Time", options2);
	
	/*/ Time difference spectra
	TH1D* TimeDiff11 = HistCanvas(data, "TU5TU5TimeDiff", "TimeDiff", options1, TU5, TU5);
	TH1D* TimeDiff22 = HistCanvas(data, "TU4TU4TimeDiff", "TimeDiff", options2, TU4, TU4);
	TH1D* TimeDiff12 = HistCanvas(data, "TU5TU4TimeDiff", "TimeDiff", options1, TU5, TU4); // t(TU5)-t(TU4), all TU5 events filled twice
	TH1D* TimeDiff21 = HistCanvas(data, "TU4TU5TimeDiff", "TimeDiff", options2, TU4, TU5);//*/
	/*TCanvas* cTimeDiff = new TCanvas("TimeDiff", "TimeDiff");
	cTimeDiff->Divide(2,2);
	cTimeDiff->cd(1);
	//cTimeDiff->cd(1)->SetLogy();
	TimeDiff11->SetStats(0);
	TimeDiff11->Draw();
	cTimeDiff->cd(2);
	//cTimeDiff->cd(2)->SetLogy();
	TimeDiff12->SetStats(0);
	TimeDiff12->Draw();
	cTimeDiff->cd(3);
	//cTimeDiff->cd(3)->SetLogy();
	TimeDiff21->SetStats(0);
	TimeDiff21->Draw();
	cTimeDiff->cd(4);
	//cTimeDiff->cd(4)->SetLogy();
	TimeDiff22->SetStats(0);
	TimeDiff22->Draw();
	//*/

    // Checks with manual selection string
    TString man_sel = "";
    // TU5 spectrum of coincident events 
    man_sel = "det==0 && (TimeDiff_after1<0.5E6 || TimeDiff_before1<4.E6)";
    TH1D* coincident_E1 = HistCanvas(data, "CoincTU5Edep", "Edep", man_sel, TU5, 0);
    // TU4 spectrum of coincident events 
    man_sel = "det==1 && ((TimeDiff_after0<4.E6 && (EnergyDep_after0>30 && EnergyDep_after0<32)) || (TimeDiff_before0<0.5E6 && (EnergyDep_before0>30 && EnergyDep_before0<32)))"; // 2700-3173 & 3489-3667
    TH1D* coincident_E2 = HistCanvas(data, "CoincTU4Edep", "Edep", man_sel, TU4, 0);
    //*/
    
    // Draw coincident spectra together with raw ones
    TCanvas* cEdep1 = new TCanvas("TU5Edep", "TU5Edep");
    gPad->SetLogy();
    histogram_E1->Draw();
    histogram_E1->SetStats(0);
    histogram_E1->GetXaxis()->SetRangeUser(0, 100);
    coincident_E1->SetLineColor(kRed);
    coincident_E1->Draw("same");
    TLegend* l1 = new TLegend(0.6, 0.7, 0.85, 0.8);
    l1->AddEntry(histogram_E1, "all events");
    l1->AddEntry(coincident_E1, "coincident events");
    l1->Draw();
    
    TCanvas* cE2 = new TCanvas("TU4Edep", "TU4Edep");
    gPad->SetLogy();
    histogram_E2->Draw();
    histogram_E2->SetStats(0);
    histogram_E2->GetXaxis()->SetRangeUser(0, 800);
    coincident_E2->SetLineColor(kRed);
    coincident_E2->Draw("same");
    TLegend* l2 = new TLegend(0.6, 0.7, 0.85, 0.8);
    l2->AddEntry(histogram_E2, "all events");
    l2->AddEntry(coincident_E2, "coincident events");
    l2->Draw();//*/
    
	// Draw 2D E vs E spectrum
    TH2D* hEvsE = GetCoincidentEvsE(data, TU4, TU5, -4.0E6, 0.5E6);
    TCanvas* cEvsE = new TCanvas("EvsE");
    hEvsE->SetStats(0);
    hEvsE->Draw("colz");//*/
}
