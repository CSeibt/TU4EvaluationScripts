#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
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
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <dirent.h>



int test(){

    struct event{
        int32_t myint;
        vector<int32_t>* myvect; 
    } MyEvent;

    TTree* tree = new TTree("myTree", "My Tree");
    tree->Branch("myint", &MyEvent.myint, "myint/I");
    tree->Branch("myvect", &MyEvent.myvect);

    for (int i = 0; i < 5; i++){
        MyEvent.myint = i;
        MyEvent.myvect->push_back(i);
        tree->Fill();
    }

    // Read events from tree and draw vectors into a TGraph
    TCanvas* canvas = new TCanvas("canvas", "My Canvas");
    TGraph* graph = new TGraph();
    int nPoints = 0;

    tree->SetBranchAddress("myint", &MyEvent.myint);
    tree->SetBranchAddress("myvect", &MyEvent.myvect);

    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        for (int j = 0; j < MyEvent.myvect->size(); j++) {
            graph->SetPoint(nPoints++, MyEvent.myint, MyEvent.myvect->at(j));
        }
    }

    graph->Draw("AC");
    canvas->Update();
    
    //canvas->SaveAs("output.png"); 


    

    return 0;
}