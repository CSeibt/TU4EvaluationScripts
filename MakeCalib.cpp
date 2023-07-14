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
#include "DrawHist.C"

TH1D* RawSpectrum(
	TFile* file,
	TString histname
){
	TH1D* histo = (TH1D*)file->Get(histname);
	return histo;
}

struct Fitresult {
	TF1* function;
	Double_t En;
	Double_t dEn;
};

//Linear background function ignoring peak area
Double_t FitLinBg(
	Double_t* x,
	Double_t* par
){
	Double_t xx = x[0];													//Actual x-value (Double_t* x in input is array)
	
	// Reject Area
	if (xx > par[2] && xx < par[3]){
		TF1::RejectPoint();
		return 0;
	}
	//Fit Function: LinearBG = A + B * x 
	Double_t LinearBG = par[0] + par[1] * xx;
	return LinearBG;
}

/*/
Fitting a function over a variable number of peaks:
* low ... up: channel range for evaluation
* vector peaks: approximated position of each peak in the range: peaks[0]
* is always the respective peak of interest, i.e. peaks[0]~E
* E, DE: Energy of peak of interest and respective uncertainty
* LowBg: Lower border for background; LowPeak: Lower border for peaks
* UpPeak: Upper border for peaks; UpBg: Upper border of background
* KeepData: If true, TFile and TH1D opened here are kept open and TH1D and TF1 are drawn. If false, heap space handled better
/*/
Fitresult FitPeakN(
	TString run, 
	TString HistName, 
	Double_t LowBg, 
	Double_t LowPeak, 
	Double_t UpPeak, 
	Double_t UpBg, 
	vector<Int_t> peaks, 
	Double_t E, 
	Double_t DE, 
	Double_t sig = 7, 
	bool KeepData = true
){
	TFile* file = File(run);											//Open file for respective run
	TH1D* spectrum = RawSpectrum(file, HistName);
	
	Fitresult FitPeak;													//Create Fit-struct for return
	FitPeak.En = E;
	FitPeak.dEn = DE;
	//Fitresult* FitPeakPointer = &FitPeak;								//Create Pointer for Fit-struct
	cout << "Basics set" << endl;
	
	
	stringstream name;
    name << "run" << run << ", " << std::setprecision(4) << E << "keV"; // unique name to prevent replacing objects
    spectrum->SetName(name.str().c_str());
    // replace '_' by ' ' for title
    string title = name.str(); for (int c = 0; c < title.length(); c++) if (title.at(c) == '_') title.at(c) = ' ';
    //cout << "Name:  " << name.str() << endl << "Title: " << title << endl;
    cout << title << endl;
    spectrum->SetTitle(title.c_str());
    
    //Define fit function based on number of peaks (each peak fitted with gausn) ... 		... maybe if-function for empty peaks vector
    TString function = "gausn(0)";
    for (int i = 1; i<peaks.size(); i++) {
		function += Form("+gausn(%i)", 3*i);
	}
	function += Form("+pol1(%i)", (int)peaks.size()*3);
	cout << function << endl;
	
	//Fitting of linear Background w/ FitLinBg-function
	TF1* LinearBg = new TF1("LinearBg", FitLinBg, LowBg, UpBg, 4);
	//LinearBg->
	LinearBg->FixParameter(2, LowPeak);
	LinearBg->FixParameter(3, UpPeak);
	spectrum->Fit("LinearBg", "0LRQ");

	//Fitting the multi gaussian fixing the background variables
	TF1* fit_fcn = new TF1("fit_fcn", function, LowBg, UpBg);
	
    
    Double_t mpos;														//peak position
	Double_t mval;														//peak maximum
	Double_t mint;														//peak integral
	cout << "FitPeakN Parameter Preset:" << endl;
	
    
    for (int j=0; j<peaks.size(); j++) {
		spectrum->GetXaxis()->SetRangeUser(peaks[j]-sig, peaks[j]+sig);
		mpos = spectrum->GetXaxis()->GetBinCenter(spectrum->GetMaximumBin());		//peak position
		mval = spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(mpos));		//peak maximum
		mint = mval * sig * sqrt(2*M_PI);											//peak integral
		fit_fcn->SetParameter(3*j, mint);
		fit_fcn->SetParameter(3*j+1, mpos);
		fit_fcn->SetParameter(3*j+2, sig);
		//Test couts
		cout << Form("p%i = ", 3*j) << mint << endl;
		cout << Form("p%i = ", 3*j+1) << mpos << endl;
		cout << Form("p%i = ", 3*j+2) << sig << endl;
	}
    
    fit_fcn->FixParameter(3*peaks.size(), LinearBg->GetParameter(0));
    fit_fcn->FixParameter(3*peaks.size()+1, LinearBg->GetParameter(1));
    //Test couts
    cout << Form("p%lu = ", 3*peaks.size()) << LinearBg->GetParameter(0) << endl;
    cout << Form("p%lu = ", 3*peaks.size()+1) << LinearBg->GetParameter(1) << endl;
    
    if (!KeepData){ free(LinearBg); }
    
    spectrum->GetXaxis()->SetRangeUser(LowBg, UpBg);
    // Fit
    spectrum->Fit("fit_fcn", "0LRQ");
    FitPeak.function = fit_fcn;
	
	if (KeepData){
		// draw fit
		new TCanvas();
		spectrum->GetXaxis()->SetTitle("Channel");
		spectrum->GetYaxis()->SetTitle("Entries");
		spectrum->SetStats(0);
		spectrum->Draw("hist");
		fit_fcn->Draw("same");
	} else {
		free(spectrum);
		free(file);
	}
	// output fit results
	cout << "Peak position = " 
		 << fit_fcn->GetParameter(1)
		 << " +/- " << fit_fcn->GetParError(1)  << endl;
	Double_t factor = 2*sqrt(2*log(2)); // switch from std dev to FWHM
	cout << "FWHM = " 
		 << factor*fit_fcn->GetParameter(2)
		 << " +/- " << factor*fit_fcn->GetParError(2)  << endl;
	Double_t calib_factor = E / fit_fcn->GetParameter(1); // guess calibration from single point
	cout << "FWHM = "
		 << calib_factor*factor*fit_fcn->GetParameter(2)
		 << " +/- " << calib_factor*factor*fit_fcn->GetParError(2) << " keV" << endl;
	cout << "Peak content = "
		 << fit_fcn->GetParameter(0)
		 << " +/- " << fit_fcn->GetParError(0) << endl;

	// Calibration data point copy format: write to text file
    return FitPeak;
}

// TU5 Calibration Data
//Fit the peaks used for calibration via FitPeak-functions. From these peaks, obtain peak position (plus uncertainty)
//and respective Energy (FitPeak input, plus uncertainty)
//Output: Channel, Energy, Channel uncertainty, Energy uncertainty 		(as respective vectors in their own vector)
vector<vector<Double_t>> CalibDataPointsTU5(
){
	TString run1 = "run003";		//133Ba
	TString run2 = "run004";		//57Co
	TString TU5 = "TU5adc";
	
	//Vectors containing the (approx.) peak positions. NOTE: the first entry shall match the energy in the respective FitPeakN-function
	//133Ba
	vector<Int_t> peaksBa1 = {440, 470, 505, 540, 565};		//4.28		sig=8
	vector<Int_t> peaksBa2 = {470, 440, 505, 540, 565};		//4.61
	vector<Int_t> peaksBa3 = {505, 440, 470, 540, 565};		//4.94
	//vector<Int_t> peaksBa4 = {440, 470, 505, 540, 565};	//5.28
	vector<Int_t> peaksBa5 = {3095, 3135};		//30.63		sig = 15
	vector<Int_t> peaksBa6 = {3135, 3095};		//30.97
	vector<Int_t> peaksBa7 = {3530, 3620};		//34.99
	//57Co
	vector<Int_t> peaksCo1 = {650, 720};		//6.40		//sig = 8
	vector<Int_t> peaksCo2 = {720, 650};		//7.06
	vector<Int_t> peaksCo3 = {1460};			//14.41		//sig = 10
	
	vector<Fitresult> f;
	
	
	//Fit all peaks and obtain the resulting functions in form "struct Fitresult"
	//FitPeakN(TString run, Double_t LowBg, Double_t LowPeak, Double_t UpPeak, Double_t UpBg, vector<Int_t> peaks, Double_t E, Double_t DE, Double_t sig = 7
	//133Ba
	f.push_back(FitPeakN(run1, TU5, 440-132, 440-32, 565+32, 565+132, peaksBa1, 4.2851, 0.0033, 8));
	f.push_back(FitPeakN(run1, 440-132, 440-32, 565+32, 565+132, peaksBa2, 4.6198, 0.0033, 8));
	f.push_back(FitPeakN(run1, 440-132, 440-32, 565+32, 565+132, peaksBa3, 4.936, 0.006, 8));
//	f.push_back(FitPeakN(run1, 440-132, 440-32, 565+32, 565+132, peaksBa4, 4.2851, 0.0033, 8));
	f.push_back(FitPeakN(run1, 2900, 3095-60, 3135+60, 3300, peaksBa5, 30.625, 0.017, 15));
	f.push_back(FitPeakN(run1, 2900, 3096-60, 3135+60, 3300, peaksBa6, 30.973, 0.017, 15));
	f.push_back(FitPeakN(run1, 3300, 3530-60, 3620+60, 3800, peaksBa7, 34.987, 0.017, 15));
	//57Co
	f.push_back(FitPeakN(run2, 650-132, 650-32, 720+32, 720+132, peaksCo1, 6.3995, 0.0016, 8));
	f.push_back(FitPeakN(run2, 650-132, 650-32, 720+32, 720+132, peaksCo2, 7.058, 0.0016, 8));
	f.push_back(FitPeakN(run2, 1340, 1460-40, 1460+40, 1580, peaksCo3, 14.4129, 0.0006, 10));
	
	//Vectors for TGraphErrors of Energy-Channel-Distribution
	vector<Double_t> ch;
	vector<Double_t> E;
	vector<Double_t> dch;
	vector<Double_t> dE;
	
	for (int i=0; i<f.size(); i++){
		Double_t channel = f[i].function->GetParameter(1);
		Double_t channeldev = f[i].function->GetParError(1);
		ch.push_back(channel);
		dch.push_back(channeldev);
		E.push_back(f[i].En);
		dE.push_back(f[i].dEn);
	}
	
	vector<vector<Double_t>> returnedval = {ch, E, dch, dE};
	return returnedval;
}


// TU4 Calibration Data
//Fit the peaks used for calibration via FitPeak-functions. From these peaks, obtain peak position (plus uncertainty)
//and respective Energy (FitPeak input, plus uncertainty)
//Output: Channel, Energy, Channel uncertainty, Energy uncertainty 		(as respective vectors in their own vector)
vector<vector<Double_t>> CalibDataPointsTU4(
){
	TString run1 = "run003";		//133Ba
	TString run2 = "run004";		//57Co
	TString TU4 = "TU4adc";
	
	//Vectors containing the (approx.) peak positions. NOTE: the first entry shall match the energy in the respective FitPeakN-function
	//133Ba
	vector<Int_t> peaksBa1 = {440, 470, 505, 540, 565};		//4.28		sig=8
	vector<Int_t> peaksBa2 = {470, 440, 505, 540, 565};		//4.61
	vector<Int_t> peaksBa3 = {505, 440, 470, 540, 565};		//4.94
	//vector<Int_t> peaksBa4 = {440, 470, 505, 540, 565};	//5.28
	vector<Int_t> peaksBa5 = {3095, 3135};		//30.63		sig = 15
	vector<Int_t> peaksBa6 = {3135, 3095};		//30.97
	vector<Int_t> peaksBa7 = {3530, 3620};		//34.99
	//57Co
	vector<Int_t> peaksCo1 = {650, 720};		//6.40		//sig = 8
	vector<Int_t> peaksCo2 = {720, 650};		//7.06
	vector<Int_t> peaksCo3 = {1460};			//14.41		//sig = 10
	
	vector<Fitresult> f;
	
	
	//Fit all peaks and obtain the resulting functions in form "struct Fitresult"
	//FitPeakN(TString run, Double_t LowBg, Double_t LowPeak, Double_t UpPeak, Double_t UpBg, vector<Int_t> peaks, Double_t E, Double_t DE, Double_t sig = 7
	//133Ba
	f.push_back(FitPeakN(run1, TU4, 440-132, 440-32, 565+32, 565+132, peaksBa1, 4.2851, 0.0033, 8));
	f.push_back(FitPeakN(run1, TU4, 440-132, 440-32, 565+32, 565+132, peaksBa2, 4.6198, 0.0033, 8));
	f.push_back(FitPeakN(run1, TU4, 440-132, 440-32, 565+32, 565+132, peaksBa3, 4.936, 0.006, 8));
//	f.push_back(FitPeakN(run1, 440-132, 440-32, 565+32, 565+132, peaksBa4, 4.2851, 0.0033, 8));
	f.push_back(FitPeakN(run1, TU4, 2900, 3095-60, 3135+60, 3300, peaksBa5, 30.625, 0.017, 15));
	f.push_back(FitPeakN(run1, TU4, 2900, 3096-60, 3135+60, 3300, peaksBa6, 30.973, 0.017, 15));
	f.push_back(FitPeakN(run1, TU4, 3300, 3530-60, 3620+60, 3800, peaksBa7, 34.987, 0.017, 15));
	//57Co
	f.push_back(FitPeakN(run2, TU4, 650-132, 650-32, 720+32, 720+132, peaksCo1, 6.3995, 0.0016, 8));
	f.push_back(FitPeakN(run2, 650-132, 650-32, 720+32, 720+132, peaksCo2, 7.058, 0.0016, 8));
	f.push_back(FitPeakN(run2, 1340, 1460-40, 1460+40, 1580, peaksCo3, 14.4129, 0.0006, 10));
	
	//Vectors for TGraphErrors of Energy-Channel-Distribution
	vector<Double_t> ch;
	vector<Double_t> E;
	vector<Double_t> dch;
	vector<Double_t> dE;
	
	for (int i=0; i<f.size(); i++){
		Double_t channel = f[i].function->GetParameter(1);
		Double_t channeldev = f[i].function->GetParError(1);
		ch.push_back(channel);
		dch.push_back(channeldev);
		E.push_back(f[i].En);
		dE.push_back(f[i].dEn);
	}
	
	vector<vector<Double_t>> returnedval = {ch, E, dch, dE};
	return returnedval;
}



TGraphErrors* CalibGraph(
	vector<Double_t> ch,
	vector<Double_t> E,
	vector<Double_t> dch,
	vector<Double_t> dE
){
	if(ch.size()==E.size() && ch.size()==dch.size() && ch.size()==dE.size()){
		TGraphErrors* graph = new TGraphErrors(ch.size(), &ch[0], &E[0], &dch[0], &dE[0]);	
		graph->SetTitle("Energy Calibration; Channel; E / keV");
		return graph;
	}
	else {cout << "Vector sizes differ!" << endl;
		return 0;
	}
}

TF1* CalibFit(
	TGraphErrors* graph
){
	TF1* Calib_fnc = new TF1("Calibration_function", "[0]+[1]*x"); 
	graph->Fit("Calibration_function", "S");
	return Calib_fnc;
}

TCanvas* CalibCanvas(
	TGraphErrors* graph
){
	if(graph){
		TCanvas* c = new TCanvas("Energy_Calibration", "Energy_Calibration",100, 10, 1200, 700);
		TF1* Calib = CalibFit(graph);
		graph->Draw("ap");
		Double_t chi2 = Calib->GetChisquare();
		Int_t ndf = Calib->GetNDF();
		Double_t dx = graph->GetErrorX(0);
		cout << "dx for point 2: " << dx << endl;
		cout << "chi^2: " << chi2 << endl;
		cout << "reduced chi^2: " << chi2/ndf << endl;
		return c;
		
	}
	else{return 0;}
}

//Table to compare fitted energies with actual peak energies
void Table(
	vector<Double_t> x,
	vector<Double_t> y,
	vector<Double_t> dx,
	vector<Double_t> dy,
	TF1* fnc
){
	int width=20;
	cout << "ch" 
		 << setw(width) << "E"
		 << setw(width) << "Efit"
		 << setw(width) << "dE"
		 << setw(width) << "E-Efit"
		 << setw(width) << "E-Efit/dE" << endl;
	for (int i=0; i<x.size(); i++){
		Double_t Efit = fnc->Eval(x[i]);
		Double_t dEfit = y[i]-Efit;
		Double_t dEfiterror = dEfit/dy[i];
		cout << x[i] 
			 << setw(width) << y[i]
			 << setw(width) << Efit
			 << setw(width) << dy[i]
			 << setw(width) << dEfit
			 << setw(width) << dEfiterror << endl;
	}
}

void Table2(
	vector<Double_t> x,
	vector<Double_t> y,
	vector<Double_t> dx,
	vector<Double_t> dy,
	TF1* fnc
){
	int width=20;
	cout << "E" 
		 << setw(width) << "ch"
		 << setw(width) << "chfit"
		 << setw(width) << "dch"
		 << setw(width) << "ch-chfit"
		 << setw(width) << "ch-chfit/dch" << endl;
	for (int i=0; i<x.size(); i++){
		Double_t Efit = fnc->Eval(x[i]);
		Double_t dEfit = y[i]-Efit;
		Double_t dEfiterror = dEfit/dy[i];
		cout << x[i] 
			 << setw(width) << y[i]
			 << setw(width) << Efit
			 << setw(width) << dy[i]
			 << setw(width) << dEfit
			 << setw(width) << dEfiterror << endl;
	}
}

void MakeCalib(){
	
	
	//Values for energy calibration, currently unused

	
	//Data Points from automatic Fit:
	vector<vector<Double_t>> TU5Data = CalibDataPointsTU5();				//{ch, E, dch, dE}
	vector<vector<Double_t>> TU4Data = CalibDataPointsTU4();				//{ch, E, dch, dE}

	
	
	TGraphErrors* calibgraph = CalibGraph(fitpoints[0], fitpoints[1], fitpoints[2], fitpoints[3]);
	TF1* calibfnc = CalibFit(calibgraph);
	cout << "Print Calibration Function ..." << endl;
	TCanvas* c1 = CalibCanvas(calibgraph);
	
	
//	TFile* CalibFile = new TFile("TU5/Paper/EnergyCalibration_220919_all.root", "RECREATE");		//Root file containing all calibration data 
//	calibgraph->Write();	
//	calibfnc->Write();
//	c1->Write();																			//from 19.09.22, change name if other calib data is used
																					
	//cout << "..." << endl;
	//cout << "Print Table ..." << endl;
	Table(fitpoints[0], fitpoints[1], fitpoints[2], fitpoints[3], calibfnc);
	/*/
	TGraphErrors* calib2 = CalibGraph(E, ch, dE, dch);
	TF1* calibfnc2 = CalibFit(calib2);
	TCanvas* c1 = CalibCanvas(calib2);
	Table2(E, ch, dE, dch, calibfnc2);
	/*/
	//Fitresult fitfunction = FitPeakN("run24", 600, 800, (vector<Int_t>){650, 750}, 7.2, 0.1);
}
