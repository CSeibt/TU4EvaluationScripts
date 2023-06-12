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


/*/
DrawHist:
The purpose of this program is to draw the histograms created by 'MakeHist' from the root file.
* 
* The program contains:
The 'File'-function to open the root file
The 'HistList'-function to print a list of all objects contained in the respective file (so you see which histograms are in the file)
The 'DrawHistograms'-function to draw histograms into one canvas. The histograms are obtained by their names, which are entered as a vector
The 'DrawCanvases'-function, which can be used to draw multiple canvases with 'DrawHist'. Input is a vector of string-vectors
The 'DrawHist'-function acting as a main function. Here you can create the string-vectors containing the histogram names for 'DrawHistograms'
/*/

//Open the root File of interest
TFile* File(
    TString run = "run021",
    TString prefixPath = "../"
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


//Show List of Objects in root file
void HistList(
	TFile* rootFile
){
	cout << "The following Objects are in '" << rootFile->GetName() << "':" << endl;
	rootFile->ls();
	}


//Drawing one Canvas containing one or more histograms
TCanvas* DrawHistograms(
	TFile* rootFile,													//File containing the histograms
	vector<TString> histoName,										//List of the names of the histograms to draw
	TString runid
){
	TCanvas* c = new TCanvas(histoName[0]+"_"+runid, histoName[0]+"_"+runid, 100, 10, 1200, 700);		//Histogram canvas
	
	//Histogram colors: Change/Add colors as you pleases
	vector<short> color = {4, 3, 2, 7, 6, 5};
	
	if((TH1D*)rootFile->Get(histoName[0])){								//Only draw histogram if it exists
		cout << "Draw '" << histoName[0] << "' ..." << endl;
		TH1D* hists[histoName.size()];									//Histogram array
		hists[0] = (TH1D*)rootFile->Get(histoName[0]);					//Get first histogram of list
		hists[0]->SetLineColor(color[0]);
		hists[0]->SetStats(1);
		cout << "n of channels = " << hists[0]->GetXaxis()->GetBinUpEdge(16384) << endl;
		//Fit parameter linear fit:
		//double A = 0;
		//double B = 0.00985674;
		
		//cout << "Integral 1-160: " << hists[0]->Integral(100, 16000) << endl;
		//hists[0]->GetXaxis()->SetLimits(A+B*0, A+B*16384.);
		Double_t tborder = pow(2, 59)/1E12;
		Int_t tsteps = (int)tborder/3600;
		//hists[0]->SetTitle("");											//Set Titles if you want
		//hists[0]->GetXaxis()->SetLimits(0, tborder/(3600));
		hists[0]->SetXTitle("Channels");
		hists[0]->SetYTitle("Counts");
		gPad->SetLogy();
		hists[0]->Draw();
		if(histoName.size()>1){											//If more than one histogram name is listed, draw the other histograms in the same canvas
			for(int i=1;i<histoName.size();i++){						//For loop for all histograms except the first
				if((TH1D*)rootFile->Get(histoName[i])){					//Only draw histogram if it exists
					cout << "Draw '" << histoName[i] << "' ..." << endl;
					hists[i] = (TH1D*)rootFile->Get(histoName[i]);		//Get histograms
					hists[i]->SetLineColor(color[i]),
					hists[i]->Draw("same");								//Draw histograms
				}
				else{cout << "File does not contain any '" << histoName[i] << "'" << endl;}
			}
		}
	}
	else{cout << "Empty canvas!" << endl;
		cout << "File does not contain any '" << histoName[0] <<"'" << endl;}
	gPad->Update();
	return c;															//return Canvas
}

/*/	
//Function to draw multiple canvases using 'DrawHist'
vector<TCanvas*> DrawCanvases(
	TFile* rootFile,													//File containing the histograms
	vector<vector<TString>> hists										//Vector containing the vectors of histogram names for 'DrawHist'
){
	vector<TCanvas*> Canvases;											//Canvas vector
	
	for(int i=0; i<hists.size();i++){
		Canvases.push_back(DrawHistograms(rootFile, hists[i]));				//Execute 'DrawHist' for every vector in 'hists'
		}
	return Canvases;													//return Canvas vector
	}
/*/

void DrawHist(){
	//gStyle->SetOptStat(1001111);
	gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);

    //Variables
    TString run1 = "run038";	
    TString run2 = "run039";												
			
    TString path = "TU5/Paper/";									
    
    TFile* file1 = File(run1);
    TFile* file2 = File(run2);



    //HistList(file);
    
    //Histogram Names
    TString histoName = "TU5 adc";
//	TString histoName2 = "histo_veto_TU5_Sz1";
    
    //Vectors of histogram names (for DrawHistograms)
    vector<TString> histograms = {histoName};
//	vector<TString> histograms2 = {histoName2};
    
    //Vector of vectors of histogram names (for DrawCanvases)
//	vector<vector<TString>> hists1 = {histograms};
    
    //DrawHistograms executes
    TCanvas* c1 = DrawHistograms(file1, histograms, run1);
    TCanvas* c2 = DrawHistograms(file2, histograms, run2);
    //TCanvas* c2 = DrawHistograms(file, histograms2);
    
    //Draw Canavses execute
    //vector<TCanvas*> canvases = DrawCanvases(file, hists1);
    }
