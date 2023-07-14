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


TFile* GetInputFile(
    TString run = "run001"
){
	TString folder = "../Test_230526/"+run+"/root/";
    // Open ROOT file
    TFile* rootFile = new TFile(folder+"Data_"+run+"_coincidence.root","UPDATE");
    
	return rootFile;
}

TH1D* RawSpectrum(
	TFile* file,
	TString histname
){
	TH1D* histo = (TH1D*)file->Get(histname);
	return histo;
}


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

/* 
 Class PeakInit: This Class contains Initial and final data of one peak, which means main and secondary peak positions, background ranges,
 * guessed Sigma value, corresponding energy and peak content and The fitted (multi)gaussian for the peak.
*/

class PeakInit {
		vector<Short_t> Peaks;											//Peaks within the region of interest: first value =  Main Peak, following values = Secondary Peaks
		Short_t PeaksLow;												//Lower range of the peak area = lowest peak - 4 * sigma
		Short_t PeaksHigh;												//Upper range of the peak area = highest peak + 4 * sigma
		Short_t GuessSigma;												//Guess of the Sigma of the peak
		Short_t BgLow;													//Lower range of the bg area below peaks
		Short_t BgHigh;													//Upper range of the bg area above peaks
	public:
		Double_t E;														//Energy of the main peak
		Double_t dE;													//Energy uncertainty
		TF1* PeakFunction;												//Fitted multi gaussian of the peaks
		Double_t NDecays;												//Determined number of decays in the main peak in the measurement time
		Double_t dNDecays;												//Uncertanty of number of decays
		PeakInit (vector<Short_t>, Short_t, Short_t, Double_t, Double_t);
		PeakInit (vector<Short_t>, Short_t, Short_t, Double_t, Double_t, Double_t, Double_t);
		void FitPeaks (TString, TString, bool);
};

PeakInit::PeakInit(vector<Short_t> PeakV, Short_t Sigma, Short_t BgRange, Double_t En, Double_t dEn) {
	Peaks = PeakV;
	GuessSigma = Sigma;
	PeaksLow = *min_element(Peaks.begin(), Peaks.end()) - 4 * Sigma;
	PeaksHigh = *max_element(Peaks.begin(), Peaks.end()) + 4 * Sigma;
	BgLow = PeaksLow - BgRange;
	BgHigh = PeaksHigh + BgRange;
	E = En;
	dE = dEn;
	NDecays = 0;
	dNDecays = 0;
}

PeakInit::PeakInit(vector<Short_t> PeakV, Short_t Sigma, Short_t BgRange, Double_t En, Double_t dEn, Double_t N, Double_t dN) {
	Peaks = PeakV;
	GuessSigma = Sigma;
	PeaksLow = *min_element(Peaks.begin(), Peaks.end()) - 4 * Sigma;
	PeaksHigh = *max_element(Peaks.begin(), Peaks.end()) + 4 * Sigma;
	BgLow = PeaksLow - BgRange;
	BgHigh = PeaksHigh + BgRange;
	E = En;
	dE = dEn;
	NDecays = N;
	dNDecays = dN;
}

void PeakInit::FitPeaks(
	TString Run,
	TString HistName,
	bool KeepData = true

){
	TFile* file = GetInputFile(Run);											//Open file for respective run
	TH1D* spectrum = RawSpectrum(file, HistName);
	
	cout << "Opened TH1D from " << Run << "." <<endl;
	
	
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
    for (int i = 1; i<Peaks.size(); i++) {
		function += Form("+gausn(%i)", 3*i);
	}
	function += Form("+pol1(%i)", (int)Peaks.size()*3);
	cout << function << endl;
	
	//Fitting of linear Background w/ FitLinBg-function
	TF1* LinearBg = new TF1("LinearBg", FitLinBg, BgLow, BgHigh, 4);
	//LinearBg->
	LinearBg->FixParameter(2, PeaksLow);
	LinearBg->FixParameter(3, PeaksHigh);
	spectrum->Fit("LinearBg", "0LRQ");

	//Fitting the multi gaussian fixing the background variables
	PeakFunction = new TF1("fit_fcn", function, BgLow, BgHigh);
	
    
    Double_t mpos;														//peak position
	Double_t mval;														//peak maximum
	Double_t mint;														//peak integral
	cout << "FitPeakN Parameter Preset:" << endl;
	
    
    for (int j=0; j<peaks.size(); j++) {
		spectrum->GetXaxis()->SetRangeUser(Peaks[j]-GuessSigma, Peaks[j]+GuessSigma);
		mpos = spectrum->GetXaxis()->GetBinCenter(spectrum->GetMaximumBin());		//peak position
		mval = spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(mpos));		//peak maximum
		mint = mval * sig * sqrt(2*M_PI);											//peak integral
		PeakFunction->SetParameter(3*j, mint);
		PeakFunction->SetParameter(3*j+1, mpos);
		PeakFunction->SetParameter(3*j+2, sig);
		//Test couts
		cout << Form("p%i = ", 3*j) << mint << endl;
		cout << Form("p%i = ", 3*j+1) << mpos << endl;
		cout << Form("p%i = ", 3*j+2) << sig << endl;
	}
    
    PeakFunction->FixParameter(3*Peaks.size(), LinearBg->GetParameter(0));
    PeakFunction->FixParameter(3*Peaks.size()+1, LinearBg->GetParameter(1));
    //Test couts
    cout << Form("p%lu = ", 3*Peaks.size()) << LinearBg->GetParameter(0) << endl;
    cout << Form("p%lu = ", 3*Peaks.size()+1) << LinearBg->GetParameter(1) << endl;
    
    if (!KeepData){ free(LinearBg); }
    
    spectrum->GetXaxis()->SetRangeUser(LowBg, UpBg);
    // Fit
    spectrum->Fit("fit_fcn", "0LRQ");
	
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
		 << PeakFunction->GetParameter(1)
		 << " +/- " << PeakFunction->GetParError(1)  << endl;
	Double_t factor = 2*sqrt(2*log(2)); // switch from std dev to FWHM
	cout << "FWHM = " 
		 << factor*PeakFunction->GetParameter(2)
		 << " +/- " << factor*PeakFunction->GetParError(2)  << endl;
	Double_t calib_factor = E / PeakFunction->GetParameter(1); // guess calibration from single point
	cout << "FWHM = "
		 << calib_factor*factor*PeakFunction->GetParameter(2)
		 << " +/- " << calib_factor*factor*PeakFunction->GetParError(2) << " keV" << endl;
	cout << "Peak content = "
		 << PeakFunction->GetParameter(0)
		 << " +/- " << PeakFunction->GetParError(0) << endl;

}



class Calibrations {
	private:
		TGraphErrors* ECalGraph;
		TGraphErrors* EResGraph;
		TGraphErrors* EffCalGraph;
		vector<Double_t> FWHM;
		vector<Double_t> dFWHM;
		vector<Double_t> EnergyEff;
		vector<Double_t> dEnergyEff;
		vector<Double_t> Epsilon;
		vector<Double_t> dEpsilon;
	public:
		vector<PeakInit*> PeakVector;
		TF1* EnergyCalibration;
		TF1* EnergyResolution;
		TF1* EfficiencyCalibration;
		Calibrations(vector<PeakInit*>);
		TF1* CalibrateEnergy();
		void DetermineCalibrationValues();
		TF1* CalibrateResolution();
		TF1* CalibrateEfficiency();
};

Calibrations::Calibrations(vector<PeakInit*> Peaks){
	PeakVector = Peaks;
	FWHM = {};
	dFWHM = {};
	EnergyEff = {};
	dEnergyEff = {};
	Epsilon = {};
	dEpsilon = {};
}

TF1* Calibrations::CalibrateEnergy(){
	//Initialize vectors:
	vector<Double_t> Channels = {};
	vector<Double_t> dChannels = {};
	vector<Double_t> Energies = {};
	vector<Double_t> dEnergies = {};
	//Fill Vectors
	for (int i = 0; i < PeakVector.size(); i++){
		Channels.push_back(PeakVector[i]->PeakFunction->GetParameter(1));
		dChannels.push_back(PeakVector[i]->PeakFunction->GetParError(1));
		Energies.push_back(PeakVector[i]->E);
		dEnergies.push_back(PeakVector[i]->dE);
	}
	//Get TGraphErrors
	ECalGraph = new TGraphErrors(Channels.size(), &Channels[0], &Energies[0], &dChannels[0], &dEnergies[0]);
	ECalGraph->SetTitle("Energy Calibration; Channel; E / keV");
	
	//Get EnergyCalibration
	EnergyCalibration = new TF1("Calibration_function", "[0]+[1]*x"); 
	ECalGraph->Fit("Calibration_function", "0S");
	
	//Draw both
	TCanvas* c = new TCanvas("Energy_Calibration", "Energy_Calibration",100, 10, 1200, 700);
	ECalGraph->Draw("ap");
	EnergyCalibration->Draw("same");
	
	//Output of Fit data as crosscheck
	Double_t chi2 = EnergyCalibration->GetChisquare();
	Int_t ndf = EnergyCalibration->GetNDF();
	Double_t dx = ECalGraph->GetErrorX(0);
	cout << "dx for point 2: " << dx << endl;
	cout << "chi^2: " << chi2 << endl;
	cout << "reduced chi^2: " << chi2/ndf << endl;
	return EnergyCalibration;
}

void Calibrations::DetermineCalibrationValues(){
	Double_t Slope = EnergyCalibration->GetParameter(1);
	Double_t dSlope = EnergyCalibration->GetParError(1);
	for (int i = 0; i < PeakVector.size(); i++){
		FWHM.push_back(PeakVector[i]->PeakFunction->GetParameter(2) * Slope * 2.35);
		dFWHM.push_back(sqrt(pow(2.35 * Slope * PeakVector[i]->PeakFunction->GetParError(2), 2) + pow(2.35 * PeakVector[i]->PeakFunction->GetParameter(2) * dSlope, 2)));
		if (PeakVector[i]->NDecays != 0){
			EnergyEff.push_back(PeakVector[i]->E);
			dEnergyEff.push_back(PeakVector[i]->dE);
			Epsilon.push_back(PeakVector[i]->PeakFunction->GetParameter(0) / PeakVector[i]->NDecays);
			dEpsilon.push_back(sqrt(pow(PeakVector[i]->PeakFunction->GetParError(0) / PeakVector[i]->NDecays, 2) 
			+ pow(PeakVector[i]->PeakFunction->GetParameter(0) * PeakVector[i]->dNDecays / pow(PeakVector[i]->NDecays, 2), 2)));
		}
	}
}

TF1* Calibrations::CalibrateResolution(){
	//Initialize vectors:
	vector<Double_t> Energies = {};
	vector<Double_t> dEnergies = {};
	if (FWHM.size() == 0){
		DetermineCalibrationValues();
	}
	//Fill Vectors
	for (int i = 0; i < PeakVector.size(); i++){
		Energies.push_back(PeakVector[i]->E);
		dEnergies.push_back(PeakVector[i]->dE);
	}
	//Get TGraphErrors
	EResGraph = new TGraphErrors(Channels.size(), &Energies[0], &FWHM[0], &dEnergies[0], &dFWHM[0]);
	EResGraph->SetTitle("FWHM(#it{E}); #it{E} / keV; FWHM / keV");
	
	//Get EnergyCalibration
	EnergyResolution = new TF1("FWHM_fnc", "[0]*sqrt([1]+x)"); 
	EResGraph->Fit("FWHM_fnc", "0S");
	
	//Draw both
	TCanvas* c = new TCanvas("Energy_Resolution", "Energy_Resolution",100, 10, 1200, 700);
	EResGraph->Draw("ap");
	EnergyResolution->Draw("same");
	
	//Output of Fit data as crosscheck
	Double_t chi2 = EnergyResolution->GetChisquare();
	Int_t ndf = EnergyResolution->GetNDF();
	Double_t dx = EResGraph->GetErrorX(0);
	cout << "dx for point 2: " << dx << endl;
	cout << "chi^2: " << chi2 << endl;
	cout << "reduced chi^2: " << chi2/ndf << endl;
	return EnergyResolution;
}


TF1* Calibrations::CalibrateEfficiency(){
	//Initialize vectors:
	if (Epsilon.size() == 0){
		DetermineCalibrationValues();
	}
	//Get TGraphErrors
	EffCalGraph = new TGraphErrors(Channels.size(), &EnergyEff[0], &Epsilon[0], &dEnergyEff[0], &dEpsilon[0]);
	EffCalGraph->SetTitle("Efficiency Calibration; #it{E} / keV; #varepsilon / %");
	
	//Get EnergyCalibration
	EfficiencyCalibration = new TF1("Efficiency_fnc", "[0] * exp(-[1]/x) * (1 - exp(-[2]/x))"); 
	EffCalGraph->Fit("Efficiency_fnc", "0S");
	
	//Draw both
	TCanvas* c = new TCanvas("Efficiency_Calibration", "Efficiency_Calibration",100, 10, 1200, 700);
	EffCalGraph->Draw("ap");
	EfficiencyCalibration->Draw("same");
	
	//Output of Fit data as crosscheck
	Double_t chi2 = EfficiencyCalibration->GetChisquare();
	Int_t ndf = EfficiencyCalibration->GetNDF();
	Double_t dx = EffCalGraph->GetErrorX(0);
	cout << "dx for point 2: " << dx << endl;
	cout << "chi^2: " << chi2 << endl;
	cout << "reduced chi^2: " << chi2/ndf << endl;
	return EfficiencyCalibration;

/* Main function */
void MakeCalib_230711(){
	
}
