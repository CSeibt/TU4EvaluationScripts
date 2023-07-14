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
#include "TColor.h"
#include "MakeCalib.cpp"

Int_t ciTUB = TColor::GetFreeColorIndex();								//Color index of TU Blue
TColor* TUBlue = new TColor(ciTUB, 0, 48/255, 93/255);					//TColor of TU Blue by RGB Values

Int_t ciTULB = TColor::GetFreeColorIndex();								//Color index of TU Light Blue
TColor* TULBlue = new TColor(ciTULB, 0, 105/255, 180/255);					//TColor of TU Light Blue by RGB Values

struct RebinPar {
	Double_t EMin;
	Double_t EMax;
	Int_t NofChannels;
};


//Calibrate a Spectrum and Rebin it to get a uniform histogram (e.g. binwidth = 100eV, binwidth=1keV) 
TH1D* CalibAndRebin(
	TH1D* UncalSpectrum,
	TF1* CalibFunction,
	RebinPar RebinParameters,		//Struct of rebin parameters containing Energy limits and Number of channels in rebinned histogram
	TString namesuffix = "",
	bool save = false
){
	//Important!: Uncalibrated hist: channel 1 goes from 0 to 1, maximum channel goes from max_ch-1 to max_ch
	//E = A + B * channel
	cout << "Start Calibration and Rebinning ... \n" << endl;
	Double_t A = CalibFunction->GetParameter(0);
	Double_t B = CalibFunction->GetParameter(1);
	Int_t UncalNofBins = UncalSpectrum->GetXaxis()->GetNbins();
	Double_t UncalMax = UncalSpectrum->GetXaxis()->GetBinUpEdge(UncalNofBins);
	UncalSpectrum->GetXaxis()->SetLimits(A+B*0, A+B*UncalMax);
	
	Int_t Chmin = UncalSpectrum->GetXaxis()->FindBin(RebinParameters.EMin);
	Int_t Chmax = UncalSpectrum->GetXaxis()->FindBin(RebinParameters.EMax);
	if (Chmax-Chmin > RebinParameters.NofChannels){
		TH1D* CalibratedSpectrum = new TH1D("TU5Energy", "TU5 Calibrated Energy Spectrum", 
											RebinParameters.NofChannels, 
											RebinParameters.EMin, RebinParameters.EMax);
		
		cout << "0 keV Bin range from "<< UncalSpectrum->GetXaxis()->GetBinLowEdge(Chmin) << " to " << UncalSpectrum->GetXaxis()->GetBinUpEdge(Chmin) << endl;
		cout << "160 keV Bin range from "<< UncalSpectrum->GetXaxis()->GetBinLowEdge(Chmax) << " to " << UncalSpectrum->GetXaxis()->GetBinUpEdge(Chmax) << endl;
		// Init pseudo random
		TRandom3 *rnd = new TRandom3(time(NULL));	
		
		for (Int_t bin = Chmin; bin <= Chmax; bin++) {
        // get bin position
        Double_t min = UncalSpectrum->GetXaxis()->GetBinLowEdge(bin);
        Double_t max = UncalSpectrum->GetXaxis()->GetBinUpEdge(bin);
        Int_t entries = UncalSpectrum->GetBinContent(bin);
        // loop over uncal bin entries
			for (Int_t i = 0; i < entries; i++) {
				// random position inside bin
				// fill calibrated value
				CalibratedSpectrum->Fill(rnd->Uniform(min, max));
			}
		}
		
		cout << "Integral Uncalibrated Spectrum = " << 	UncalSpectrum->Integral(Chmin, Chmax) << endl;
		cout << "N of Entries in Calibrated Spectrum = " << CalibratedSpectrum->GetEntries() << endl;		
		cout << "Integral Calibrated Spectrum = " << CalibratedSpectrum->Integral(1, RebinParameters.NofChannels) << endl;
		
		TCanvas* c1 = new TCanvas("Uncalibrated_Spectrum"+namesuffix, "Uncalibrated Spectrum"+namesuffix, 100, 10, 1200, 700);
		UncalSpectrum->SetXTitle("E/keV");
		UncalSpectrum->SetYTitle("Counts");
		UncalSpectrum->Draw();
		TCanvas* c2 = new TCanvas("Calibrated_Spectrum"+namesuffix, "Calibrated Spectrum"+namesuffix, 105, 15, 1200, 700);
		CalibratedSpectrum->SetXTitle("E/keV");
		CalibratedSpectrum->SetYTitle("Counts");
		CalibratedSpectrum->Draw();
		cout << "finished Calibration and Rebinning \n" << endl;
		return CalibratedSpectrum;
	}
	else {return 0;}
}

TCanvas* CombinedCanvas(
	TH1D* Hist1,
	Double_t time1,
	TH1D* Hist2,
	Double_t time2
){
	//gPad->SetRightMargin(0.01);
	gStyle->SetTitleSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "X");
	
	gStyle->SetTitleSize(0.05, "Y");
	gStyle->SetLabelSize(0.05, "Y");
	
	
	//gStyle->SetTitleOffset(0.9, "Y");
	//gStyle->SetCanvasDefX(15);

	TCanvas* combined = new TCanvas("Comparison_above_and_below_ground", "Comparison above and below ground", 100, 10, 1550, 675);
	combined->SetCanvasSize(1500, 625);
    Double_t Integral1 = Hist1->Integral(Hist1->FindBin(2), Hist1->FindBin(20));
    Double_t Integral2 = Hist2->Integral(Hist2->FindBin(2), Hist2->FindBin(20));
    
    Double_t Integral3 = Hist1->Integral(Hist1->FindBin(0), Hist1->FindBin(120));
    Double_t Integral4 = Hist2->Integral(Hist2->FindBin(0), Hist2->FindBin(120));
    cout << "2-10 keV" << endl;
    cout << "Integral above ground = " << Integral1*3600*24/(18*time1) << "+/-" << sqrt(Integral1)*3600*24/(18*time1) << "cts/keV*days" << endl;
    cout << "Integral below ground = " << Integral2*3600*24/(18*time2) << "+/-" << sqrt(Integral2)*3600*24/(18*time2) <<"cts/keV*days"<< endl;
    cout << "0-120 keV" << endl;
    cout << "Integral above ground = " << Integral3*3600*24/(time1) << "+/-" << sqrt(Integral3)*3600*24/(time1) << "cts/keV*days" << endl;
    cout << "Integral below ground = " << Integral4*3600*24/(time2) << "+/-" << sqrt(Integral4)*3600*24/(time2) <<"cts/keV*days"<< endl;
    //List of Activity estimations
    int width = 20;
    cout << "E / keV"
		 << setw(width) << "sigma"
		 << setw(width) << "Bgcts(E, sigma)"
		 << setw(width) << "epsilon"
		 << setw(width) << "A estimated" << endl;
	vector<Double_t> E = {4.28508, 6.39948, 7.05798, 14.4129, 30.6251, 30.9728, 35.822, 53.1622, 80.9979};
	vector<Double_t> eps = {0.0038, 0.014, 0.0167, 0.01081, 0.00162, 0.001458, 0.00114, 0.000199, 0.0000288};
	vector<Double_t> sigma = {0.07054, 0.07801, 0.0802, 0.10146, 0.13714, 0.1378, 0.14675, 0.17054, 0.2127};
	vector<Double_t> Aestimated = {};
	for (int i=0; i<E.size(); i++){
		Double_t ctr = Hist2->Integral(Hist2->FindBin(E[i]-3*sigma[i]), Hist2->FindBin(E[i]+3*sigma[i]))/time2;
		Aestimated.push_back(ctr/eps[i]);
		cout << E[i]
			<< setw(width) << sigma[i]
			<< setw(width) << ctr
			<< setw(width) << eps[i]
			<< setw(width) << Aestimated[i] << endl;
		}

	Hist1->Scale(1/time1);
	Hist2->Scale(1/time2);
	Hist1->SetTitle(" ");
	Hist1->SetXTitle("#it{E} / keV");
	Hist1->SetYTitle("Counting rate / s^{-1} #bullet 20 eV");
	Hist1->SetStats(0);
	Hist2->SetTitle("TU5 Background Spectrum");
	Hist2->SetXTitle("#it{E} / keV");
	Hist2->SetYTitle("Counting rate / s^{-1} #cdot 20 eV");
	Hist2->SetStats(0);
	Hist1->SetLineColor(4);
	Hist1->SetMarkerColor(4);
	Hist2->SetLineColor(8);
	Hist2->SetMarkerColor(8);
	
	Hist1->Rebin(5);		//1 bin = 0.1 keV
	Hist2->Rebin(5);
	Hist1->GetXaxis()->SetRangeUser(0, 120);
	Hist1->GetYaxis()->SetRangeUser(0.000001, 0.001);
	Hist2->GetXaxis()->SetRangeUser(0, 120);
	Hist2->GetYaxis()->SetRangeUser(0.000001, 0.001);
	Hist1->Draw();
	Hist2->Draw("same");
	
	TBox* ROI = new TBox(2.0, 1e-6, 20.0, 0.001);
	ROI->SetFillColorAlpha(kRed, 0.3); // Set the fill color and transparency
	ROI->SetLineColor(kRed);
	ROI->SetFillStyle(3002);
	ROI->Draw("same");
	TLegend* leg1 = new TLegend(0.5,0.75,0.9,0.9); //new legend with x_min,y_min,x_max,y_max
    	leg1->SetFillColor(0); //fill color white
    	leg1->SetTextSize(0.04); 
    	gPad->SetLeftMargin(0.15);
    	leg1->AddEntry(Hist1, "background rate above ground","l");
    	leg1->AddEntry(Hist2, "background rate below ground","l");
    	leg1->Draw("same");
    	gPad->Update();
    
	TGraph* Activitygraph = new TGraph(E.size(), &E[0], &Aestimated[0]);
	TCanvas* ActCanvas = new TCanvas("Activity_estimation", "Activity estimation", 0, 0, 1250, 750);
	ActCanvas->SetCanvasSize(1200, 700);

	Activitygraph->SetTitle("Activity estimation; #it{E} / keV; #it{A} / Bq");


	Activitygraph->SetMarkerStyle(33);
	Activitygraph->SetMarkerSize(2.2);
	const char* TUlightblue = "#0069b4";
	Activitygraph->SetMarkerColor(TColor::GetColor(TUlightblue));
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy(1);
	Activitygraph->Draw("ap");
	return combined;
	
	
	}

TCanvas* CombinedCanvasVectors(
	vector<TH1D*> Histograms,
	vector<Double_t> Times
){
	/* Style Parameters */
	vector<Int_t> Colors = {4, 7, 3, 2};
	//gPad->SetRightMargin(0.01);
	/*
	gStyle->SetTitleSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "X");
	
	gStyle->SetTitleSize(0.05, "Y");
	gStyle->SetLabelSize(0.05, "Y");
	*/
	/* Setting up Canvas and Legend */
	cout << "Start creating Canvas with multiple Histograms ..." << endl;
	TCanvas* Combined = new TCanvas("comparison_of_spectra", "comparison_of_spectra", 100, 10, 1250, 750);
	//Combined->SetCanvasSize(1200, 700);
	cout << "Empty Canvas set up" << endl;
	
	TLegend* Legend = new TLegend(0.6, 0.6, 0.9, 0.9);
	cout << "Empty Legend set up" << endl;
	/* Set up and draw all histograms in vector*/
	for (int i=0; i<Histograms.size(); i++){
		cout << " " << endl;
		Double_t Factor = 1 / Times[i];
		cout << "Integral before rescale: " << Histograms[i]->Integral() << endl;
		cout << "Determined scaling factor: Factor = " << Factor << endl;
		Histograms[i]->Scale(Factor);
		cout << "Integral after rescale: " << Histograms[i]->Integral() << endl;
		Histograms[i]->SetTitle(" ");
		Histograms[i]->SetXTitle("#it{E} / keV");
		Histograms[i]->SetYTitle("Counting rate / s^{-1} #bullet 20 eV");
		Histograms[i]->SetStats(0);
		Histograms[i]->SetLineColor(Colors[i]);
		Histograms[i]->SetMarkerColor(Colors[i]);
	
		Histograms[i]->Rebin(5);		//1 bin = 0.1 keV
		Histograms[i]->GetXaxis()->SetRangeUser(0, 120);
		//Histograms[i]->GetYaxis()->SetRangeUser(0.000001, 0.001);
		cout << "Histograms set up" << endl;
		Legend->AddEntry(Histograms[i], Histograms[i]->GetTitle(), "l");
		if (i == 0) {Histograms[i]->Draw();}
		else {Histograms[i]->Draw("same");}
		cout << "Histograms drawn " << endl;
	}
	Legend->Draw("same");
	gPad->Update();
	cout << "Combined Canvas complete \n" << endl;
	return Combined;
}

void CalibrateSpectrum(){
	
	//TString HistoName = "TU5 adc";
	//TString HistoName2 = "TU5_adc_21222829";
	//TString run1 = "run24";	
	
	RebinPar RebinParameters;
	RebinParameters.EMin = 0;
	RebinParameters.EMax = 160;
	RebinParameters.NofChannels = 8000;
	// 1 bin = 0.02 keV
	
	//Calibration data File
	TFile* CalibFile = new TFile("TU5/Paper/Calibration_031122.root", "READ");
	TF1* calibfnc = (TF1*)CalibFile->Get("Calib_fnc");
	
	//TFile* RootFile = File(run1);
	/*
	TFile* combinedfileabove = new TFile("TU5/Paper/Combination/bg_above_ground_021_030_v2.root", "UPDATE");
	TString combinednameabove = "TU5_adc_21222829";
	TH1D* CombinedHistabove = (TH1D*)combinedfileabove->Get(combinednameabove);
	TH1D* CombHistcal1 = CalibAndRebin(CombinedHistabove, calibfnc, RebinParameters, " above ground", false);
	
	TFile* combinedfilebelow = new TFile("TU5/Paper/Combination/bg_above_ground_034_042_v2.root", "UPDATE");
	TString combinednamebelow = "TU5_adc_3437404142";
	TH1D* CombinedHistbelow = (TH1D*)combinedfilebelow->Get(combinednamebelow);
	TH1D* CombHistcal2 = CalibAndRebin(CombinedHistbelow, calibfnc, RebinParameters, " below ground", false);
	
	Double_t time1 = 1789980.0;
	Double_t time2 = 1727820.0;
	
	TCanvas* Combined = CombinedCanvas(CombHistcal1, time1, CombHistcal2, time2);
	
	TFile* BaFile = File("run038");
	TH1D* UncalBa = RawSpectrum(BaFile, "TU5 adc");
	TH1D* BaHist = CalibAndRebin(UncalBa, calibfnc, RebinParameters);
	//TFile* RootFile2 = new TFile("TU5/Paper/Combination/bg_above_ground_021_030_v2.root", "UPDATE");
	//TH1D* UncalHist = (TH1D*)RootFile->Get(HistoName);
	
	//TH1D* CalibratedHist = CalibAndRebin(UncalHist, calibfnc, RebinParameters, false);
	 
	*/
	TFile* FileDiagTop = new TFile("TU5/Paper/Combination/60co_diagonal_top.root", "UPDATE");
	TString HistnameDiagTop = "TU5_adc_diagtop";
	TH1D* HistDiagTop = (TH1D*)FileDiagTop->Get(HistnameDiagTop);
	TH1D* CalHistDiagTop = CalibAndRebin(HistDiagTop, calibfnc, RebinParameters, "_diagonal_top", false);

	TFile* FileDiagBottom = new TFile("TU5/Paper/Combination/60co_diagonal_bottom.root", "UPDATE");
	TString HistnameDiagBottom = "TU5_adc_diagbottom";
	TH1D* HistDiagBottom = (TH1D*)FileDiagBottom->Get(HistnameDiagBottom);
	TH1D* CalHistDiagBottom = CalibAndRebin(HistDiagBottom, calibfnc, RebinParameters, "_diagonal_bottom", false);

	TFile* FileUnderDet = new TFile("TU5/Paper/Combination/60co_under_detector.root", "UPDATE");
	TString HistnameUnderDet = "TU5_adc_underdet";
	TH1D* HistUnderDet = (TH1D*)FileUnderDet->Get(HistnameUnderDet);
	TH1D* CalHistUnderDet = CalibAndRebin(HistUnderDet, calibfnc, RebinParameters, "_under_detector", false);

	TFile* FileBgApril = new TFile("TU5/Paper/Combination/background_april.root", "UPDATE");
	TString HistnameBgApril = "TU5_adc_background";
	TH1D* HistBgApril = (TH1D*)FileBgApril->Get(HistnameBgApril);
	TH1D* CalHistBgApril = CalibAndRebin(HistBgApril, calibfnc, RebinParameters, "_background_april", false);	

	vector<TH1D*> Co60Histograms = {CalHistDiagTop, CalHistDiagBottom, CalHistUnderDet, CalHistBgApril};
	vector<Double_t> Co60Times = {1, 1, 1, 1};
	
	TCanvas* Co60Canvas = CombinedCanvasVectors(Co60Histograms, Co60Times);
}
