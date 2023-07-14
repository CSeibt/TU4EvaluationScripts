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

#define ROIS 4
#define MAXRONLINES 10000000
#define nentriesMAX 2000000000

void BinToRoot(
	//TString folder = "/cephfs/projekte/astro/Felsenkeller/3He+4He/2ndCampaign/Analysis/TUBunker/Christoph/Runs/Run004/",
	TString run = "run015",
	bool hist = true,				//Fill histograms with data?
	bool drawh = false,				//Draw histograms and save canvases?
	double maxtime = 864000.,
	double mintime = 0.    
	){
		
	//Folders separated into directory + run + folder in run (RAW)
	TString prefixPath = "";
	TString folder = "/ZIH.fast/users/felsdaq/TUBunker/TU5/TU5_TU4_230524_coincidence/DAQ/";  //Path of the .bin file which shall be evaluated
	TString suffix = "/RAW/";

    // Channel  0 = TU5
    // Channel  1 =  
    // Channel  2 = Sz1
    // Channel  3 = 
    // Channel  4 = Sz2 
    // Channel  5 = 
    // Channel  6 = 
    // Channel  7 = 

    //int Deadtimecalc = 1; // Set "=0" if no deadtimecalc is needed. Set "=1" if it is needed, but then also adjust "runID" and "rundate" in the void
    //double digitizer = 179; // Put in 179 or 397 in order to know the arrayposition for the dead time extraction
    int total_number_events = 0;

    const int numchannels = 8;
    //const int numfiles = 1;
    //const int CentralChannel = 0; //Kapsel G = 6
    const int max_num_jumps = 10;

    // Variable definitions.
    
    TString detName[8]= {"TU5","TU4","","","","","",""};
    int det_ch[8] = {0,3,1,4,2,5,6,7}; // Reordering of detector channels: TU5 = det0, TU4=det1, Sz2=det2 etc.
    
    TString input = prefixPath + folder + run + suffix;
    TStopwatch *rootTime = new TStopwatch();
    int i,j,k,m,p,line=0,first,second;
    rootTime->Start();
    double time_ratio = 0;
    int counter[numchannels];


	//Definition of branch variables
    Short_t adc=0;
    UShort_t det=0;
    UInt_t extra;
    Bool_t saturation = false, pileup = false, ADCout=false;
    ULong64_t LostTrig[numchannels], TotTrig[numchannels];
    ULong64_t timeAdopt=0;
    ULong64_t timeLast[numchannels];
    
    TString caenFileName;

    //int mirror=1;
    int maxentries = 0;
    double channelRealTime[numchannels];
    long int linesRead[numchannels];

    TObjString *info;
    bool jump = false;
    int num_time_jumps[numchannels];
    long int time_jumps_at[numchannels][max_num_jumps];
    long int time_shifts[numchannels][max_num_jumps];

	// Variable initialization.
    for (i=0;i<numchannels;i++) {
        timeLast[i] = mintime*2.5e8;
		channelRealTime[i] = 0.;
		linesRead[i] = 0;
        num_time_jumps[i] = 0;
        LostTrig[i] = 0;
        TotTrig[i] = 0;
		//cout << "mintime: " << mintime << endl;
        
        for (j=0;j<max_num_jumps;j++) {
			time_jumps_at[i][j] = 0;
			time_shifts[i][j] = 0;
		}
    }
	cout << "Variable initialization done." << endl;
	
    int nentries=0;
    int 	*index 		= new int[nentriesMAX];
    uint64_t 	*timeAdoptW 	= new uint64_t[nentriesMAX];
    uint16_t 	*detW		= new uint16_t[nentriesMAX];
    int16_t 	*adcW		= new int16_t[nentriesMAX];
    uint16_t 	*saturationW	= new uint16_t[nentriesMAX];
    uint16_t 	*pileupW	= new uint16_t[nentriesMAX];
    uint32_t 	*extrasW	= new uint32_t[nentriesMAX];

	cout << "Temporary variablez defined" << endl;

    // Count the .bin files
    // First search for all files in the directory "folder" with the extension "bin"
    TString binfiles[100];
    int total_number_files=0;
    int b=0;
    struct dirent *de_bin;
    DIR *dr_bin = opendir(input); 
    if (dr_bin == NULL)
    { 
        printf("Could not find any .bin files" ); 
    } 
  
    while ((de_bin = readdir(dr_bin)) != NULL){ 
	    string s_bin = de_bin->d_name;
	    //check in which file the letters "b", "i" and "n" occur
	    if(s_bin.find('b') != std::string::npos && s_bin.find('i') != std::string::npos && s_bin.find('n') != std::string::npos){
	      binfiles[b] = de_bin->d_name;
	      b++;
	      total_number_files++;
	    }
    }
    
    
    cout << "Sorting the files..." << endl;
    sort(binfiles,binfiles+total_number_files);
    
    closedir(dr_bin);
    cout << " " << endl;
    cout << "Following .bin files were found:" << endl;
    
    for(int i=0;i<total_number_files;i++){
      cout << binfiles[i] << endl;
    }
    cout << " " << endl;
    if(total_number_files>30){cout << "ATTENTION!!! Large amount of files were fount. Not all files were read. Please increase the array for total_number_files" << endl;}

    FILE *fin[total_number_files];  
      

    //Create run folder if it does not exist and create a ROOT file.
    struct stat st = {0};

	if (stat("Test_230526/"+run+"/root", &st) == -1) {
		int directory1 = mkdir("../Test_230526/"+run, 0777);
		int directory2 = mkdir("../Test_230526/"+run+"/root", 0777);
		cout << "Created new directory" << endl;
	}
    TFile *rootFile = new TFile("../Test_230526/" + run + "/root/"+binfiles[0]+"_without_timediff.root","RECREATE");
	cout << "created new file " << endl;
    // Plant the data tree.
    TTree *Data = new TTree("Data","TU5 data");
    Data->Branch("det",   &det ,   "det/s");
    Data->Branch("adc",    &adc,   "ch/S");
    Data->Branch("pileup",&pileup,"pileup/O");
    Data->Branch("saturation",&saturation,"saturation/O");
    Data->Branch("ADCout",&ADCout,"ADCout/O");
    Data->Branch("Time",&timeAdopt,"Time/l");
    
    
    Data->Branch("Extra",&extra,"Extra/i");


	//
    int16_t boardN;
    int16_t channelN;
    
    int64_t readTimeTag[1];
    uint64_t readTimeTagLast[numchannels]; 
    uint64_t readTimeTagLastWrap[numchannels];
    
    int16_t readEnergy[1];
    Int_t  wrap[numchannels];
    int headersize=0, wrapped=0;
    
    int32_t readExtras[1];
    UInt_t header[1],left,right;
    
    int32_t numSamples;
    int16_t samples;
   
    int NrOfFilesEnded = 0;
    
    
    // General cosmetics.
    gStyle->SetOptStat(1001111);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);
    
    // Miscellaneous variables.
    double Ratemin=0., Ratemax=(int)maxtime+100, Ratebinsize=1.;    //Ratemax is in thicks 1.0E12 == 1s
    int RateChannels = (Ratemax-Ratemin)/Ratebinsize;

	//Histogram arrays
	
    TH1I *histo_rate[numchannels], *histo_rate_pileup[numchannels], *histo_rate_saturation[numchannels];
    TH1D *histo_raw[numchannels];
    
    // Create histograms (rate, timing, raw_energy)
	for (i=0;i<numchannels;i++) {
		histo_rate[i] = new TH1I(Form("histo_rate[%d]",i), Form("Rate %s; Time [s];s^(-1)",detName[i].Data()),
						RateChannels, Ratemin, Ratemax);
		if(!hist) {histo_rate[i]->SetDirectory(NULL);}
		histo_rate_pileup[i] = new TH1I(Form("histo_rate_pileup[%d]",i), Form("Rate %s; Time [s];s^(-1)",detName[i].Data()),
						RateChannels, Ratemin, Ratemax);
		if(!hist) {histo_rate_pileup[i]->SetDirectory(NULL);}
		histo_rate_saturation[i] = new TH1I(Form("histo_rate_saturation[%d]",i), Form("Rate %s; Time [s];s^(-1)",detName[i].Data()),
						RateChannels, Ratemin, Ratemax);
		if(!hist) {histo_rate_saturation[i]->SetDirectory(NULL);}
                
		}
		histo_raw[i] = new TH1D(Form("histo_raw[%d]",i), Form("Raw %s; Channels;Counts",detName[i].Data()),
						16384, 0, 16384);
		if(!hist) {histo_raw[i]->SetDirectory(NULL);}

		
    
	
	
    // Reading the files
    // Attention!!! The input_ch will be converted
    for (i=0;i<total_number_files;i++) {
		caenFileName = input + binfiles[i];
        fin[i] = fopen(caenFileName, "r");
        if (!fin[i]) {
             cout << "Error when opening file " << caenFileName << endl;
             return;
        } else {
            cout << "Opened file " << caenFileName << endl;
        }  // if (!fin[i])
    }
    cout << " " << endl;


    do{
		nentries = 0;
		//fill the last read value in the array if any
		for (int filen=0;filen<total_number_files;filen++) {
			if(timeAdoptW[nentriesMAX-1-filen]!=0){
				cout << "Fill the last read value in the array... " <<endl;
				timeAdoptW[nentries]	= timeAdoptW[nentriesMAX-1-filen];
				detW[nentries] 		= detW[nentriesMAX-1-filen];
				adcW[nentries] 		= adcW[nentriesMAX-1-filen];
				pileupW[nentries] 	= pileupW[nentriesMAX-1-filen];
                saturationW[nentries] 	= saturationW[nentriesMAX-1-filen];
				extrasW[nentries]	= extrasW[nentriesMAX-1-filen];
				linesRead[filen]++;
				nentries ++ ;
			}
		}
		//cout << "N: " << nentries << endl;
		for (int filen=0;filen<total_number_files;filen++) {
			//cout << "test " << endl;
			//zero the temp variables
			timeAdoptW[nentriesMAX-1-filen] = 0;
			readTimeTag[0] 		= 0;
			
			cout << "Reading data from file... " << endl;
			// Read in listmode file, line by line.
			while (fread(&boardN, 1, sizeof(boardN), fin[filen]) 
					    &&fread(&channelN, 1, sizeof(channelN), fin[filen]) 
					    && (sizeof(readTimeTag) == fread(&readTimeTag, 1, sizeof(readTimeTag), fin[filen]))) 
					    {	
				
				fread(&readEnergy, 1, sizeof(readEnergy), fin[filen]);
				fread(&readExtras, 1, sizeof(readExtras), fin[filen]);
				fread(&numSamples, 1, sizeof(numSamples), fin[filen]);
				fread(&samples ,numSamples, sizeof(samples), fin[filen]);
				//cout << numSamples<< endl;
				
				
				if(channelN<8 && channelN>=0){
				    timeAdoptW[nentries]	= readTimeTag[0];
				    detW[nentries] 		= det_ch[channelN];
				    adcW[nentries] 		= readEnergy[0] & 32767;
				    saturationW[nentries]   = (readExtras[0] & 1024) / 1024;
                    pileupW[nentries]   = (readExtras[0] & 32768) / 32768;
				    extrasW[nentries]	= readExtras[0];
				    linesRead[filen]++;
				    //if(nentries<100){ cout<< channelN << " " << readExtras[0] << " " <<  (readExtras[0] & 32768) / 32768 << endl;}
				    //if(total_number_events<1000000 && (channelN == 0)){cout << boardN << " " << channelN  << " " << timeAdoptW[nentries]  <<" " << adcW[nentries]  <<endl;}
				    
				    if (nentries < nentriesMAX-1-numchannels){
					    nentries ++ ;
				    } else {
					    cout << "!!!! nentries= " << nentries << " > " << nentriesMAX-numchannels << "  TEMP array is too small!!!!" << endl;
				    }
				}
				total_number_events ++;
			} // while ()
			
		} // for
		
		//Check if we reached the end of all files, if so than stop (otherwise the last warp would also be sorted again also when no entries are in)
		for (i=0;i<total_number_files;i++) {
			NrOfFilesEnded += feof(fin[i]);
		}
		cout << " "<<endl;
		cout << "Finished reading data from file. Entries: " <<nentries << endl;
		cout << "Sorting starts..."<<endl;
		cout << " "<<endl;
		
		// Sort index of time stamps by time, upwards
		TMath::Sort(nentries,timeAdoptW,index,kFALSE);
				
		cout << "Sorting finished." << endl;
		cout << "Events are now filled into preliminary data tree..."<<endl;
		cout << " "<<endl;
		
		// Fill events into prelim data tree and rate, timing, raw_spetrum hists.
		for (p=0;p<nentries;p++){
			det = detW[index[p]];

			timeAdopt = timeAdoptW[index[p]];									


			adc	= adcW[index[p]];

			//if(p<1000){cout << det << " " << adc << " " << EnergyDep_before[det] << endl;}
			pileup = pileupW[index[p]];
            saturation = saturationW[index[p]];
			extra = extrasW[index[p]];

			ADCout =   extra & 16;
			if(extra & 32) LostTrig[det]++;
			if(extra & 64) TotTrig[det]++;
			if(TotTrig[det]%64 == 0) TotTrig[det]++;                //Number od events stored in the tree is always higher than caluated from the TotTrig.
																	// it seems, that every 64th TotTrig does not elevate the totTrig bit...
																	

			if(timeAdopt > mintime*1.0E12){
				if((extrasW[index[p]] & 1)){
					jump = true;
					num_time_jumps[det]++;
					time_jumps_at[det][num_time_jumps[det]-1] = timeLast[det];
					time_shifts[det][num_time_jumps[det]-1] = timeAdopt-timeLast[det];
					cout << "There was blind period in " << detName[det].Data() << " at " << timeLast[det]/1.0E12 << "s with length of " << (timeAdopt-timeLast[det])/1.0E12 << "s" << endl;
				}
				
				// Fill the rate hist, only if histograms shall be filled
				if(hist){
					if(!pileup && !saturation){histo_rate[det]-> Fill(1.0*timeAdopt/1.0E12);}  
					else if (!saturation) {histo_rate_pileup[det]-> Fill(1.0*timeAdopt/1.0E12);}
					else {histo_rate_saturation[det]-> Fill(1.0*timeAdopt/1.0E12);}
				}
				
				//if(!pileup){if(EnergyDep_after[det]>20 && EnergyDep_before[det]>20) {histo_rate[det]-> Fill(1.0*timeAdopt/1.0E12);}}   
				//else {histo_rate_pileup[det]-> Fill(1.0*timeAdopt/1.0E12);}
				

				//Fill raw hist, only if histograms shall be filled
				if(hist){
					histo_raw[det]-> Fill(adc);
					}
				
				
				//save the previous timestamp for a given detector, and fill the tree
				//timeLast[det] = timeAdopt;
				//EnergyDepLast[det] = adc;
				//if(p<100){cout << EnergyDep_before2 << endl;}
				Data->Fill();
			}
			
			}
		
		// Debugging output
		if (p%500000== 0) {
			cout << Form("%.1f",100.*p/nentries) <<"% done. " << "Time= " << timeAdopt/1.0E12 << " sec, now altogether " << Data->GetEntries()-1 << " preliminary entries."  << endl;
		}
		
        
		
	
		  
    }while (NrOfFilesEnded<total_number_files && timeAdopt/1.0E12<=maxtime); 



   	for (i=0;i<numchannels;i++) {
		channelRealTime[i] =  timeLast[i]/1.0E12;
	}
	
    // Output for safety.
    cout << "--------------------------------------------------------------|" << endl;
    cout << "File Nr.	  | Lines       | Realtime[s] " << endl;
    cout << "--------------------------------------------------------------|" << endl;
    
    for (i=0;i<numchannels;i++) {
        if (linesRead[i]>0) {
            cout << Form("%s		| %011ld | %.5f   |",detName[i].Data(),linesRead[i],channelRealTime[i]) << endl;
            fclose(fin[i]); // no longer needed
        }
    }
    
    maxentries = Data->GetEntries();
    cout << "--------------------------------------------------------------|" << endl;
    cout << "Number of events:	" << Form("%011d|",maxentries) << endl;
    cout << "Total number of events:	" << Form("%011d|",total_number_events) << endl;
    cout << "--------------------------------------------------------------|" << endl;

    //Put realtime, event number and lost eventnumber information in UserInfo of the tree
    TString info_string;
    for (m=0;m<numchannels;m++) {
        info_string = Form("%s_Real_time: %f s,_Tot_events: %lld ,_Lost_events: %lld", detName[m].Data(), channelRealTime[m], TotTrig[m]*1024, LostTrig[m]*1024);
        info = new TObjString(info_string);
        Data->GetUserInfo()->Add(info);
    }
    
    //Put information in UserInfo of the tree if there was a jump
    if(jump){
        for (m=0;m<numchannels;m++) {
            info_string = Form("Number_of_jump(s)_in_%s: %d", detName[m].Data(),num_time_jumps[m]);
            info = new TObjString(info_string);
            Data->GetUserInfo()->Add(info);
            for (k=0;k<num_time_jumps[m];k++) {
                info_string = Form("at %ld with_length_of %ld", time_jumps_at[m][k], time_shifts[m][k]);
                info = new TObjString(info_string);
                Data->GetUserInfo()->Add(info);
            }
        }
    }
    
    rootFile->Write();
    cout << "Time needed for converting the data to tree: Real time (in seconds)= " << rootTime->RealTime() << " CPU time (in seconds)= " << rootTime->CpuTime() << endl;

   
   //Draw the respective plots if they shall be drawn     

    if (drawh == true){
    
		//Create the rate plots
		TCanvas *c2 = new TCanvas("Rate","TU5 Rates, "+binfiles[0],150,60,1200,700);
		c2->Divide(1,3);
    
		
		//Plot rates
		for (i=0;i<1;i++) {
			c2->cd(i+1);
			histo_rate_saturation[i]->SetAxisRange(mintime,channelRealTime[i],"X");
			histo_rate_saturation[i]->SetLineColor(kGreen);
			histo_rate_saturation[i]->Draw();
			histo_rate[i]->Draw("same");
			histo_rate_pileup[i]->SetLineColor(kRed);
			histo_rate_pileup[i]->Draw("same");
		}

		//Plot rates
		for (i=1;i<3;i++) {
			c2->cd(i+1);
			histo_rate[i]->SetAxisRange(mintime,channelRealTime[i],"X");
			histo_rate[i]->Draw();
			histo_rate_pileup[i]->SetLineColor(kRed);
			histo_rate_pileup[i]->Draw("same");
			histo_rate_saturation[i]->SetLineColor(kGreen);
			histo_rate_saturation[i]->Draw("same");
		}
    

		//if there was jump, its parameters is written on the rate plot
		TLatex notice;
		c2->cd(7);
		if(jump){
			notice.DrawLatexNDC(.0,.9, "Time jump(s): ");
			for (m=0;m<7;m++) {
				if(num_time_jumps[m]>0){
					notice.DrawLatexNDC(.0+0.25*m,.85, Form("%d in %s", num_time_jumps[m], detName[m].Data()));
					for (k=0;k<num_time_jumps[m];k++) {
						notice.DrawLatexNDC(.0+0.25*m,.8-(k*0.05), Form("%0.2fs", time_jumps_at[m][k]/1.0E12));
					}
				}
			}
		}
 
    
		gPad->Update();
		//c2->Print(binfiles[0] + "_rates.pdf");
		rootFile->Add(c2);
    
   
    
    
		//Create and plot the raw spectra
		TCanvas *c3 = new TCanvas("Raw_spectra","TU5 setup - raw spectra, "+binfiles[0],200,110,1200,700);
		c3->Divide(1,3);

		for (i=0;i<3;i++) {
			c3->cd(i+1);
			gPad->SetLogy();
			histo_raw[i]->Draw();
		}

		for (m=0;m<7;m++) {
			if(num_time_jumps[m]>0){
				cout << Form("There is %d time jump(s) in %s", num_time_jumps[m], detName[m].Data()) << endl;
				for (k=0;k<num_time_jumps[m];k++){
					cout << Form("at %0.2fs with length of %0.2fs", time_jumps_at[m][k]/1.0E12, time_shifts[m][k]/1.0E12) << endl;
				}
			}
		}

		rootFile->Add(c3);
		//rootFile->Add(c5);
	
	/*/
	if(!hist){cout << "hist = false" << endl;
		for(i=0;i<numchannels;i++){
			rootFile->Delete(Form("histo_raw[%i]", i));
		}
	}
	* /*/
    rootFile->Write();
	}
    //Save raw spectra to file
    /*/
    ofstream outfile(binfiles[0]+"_raw_spectra.pdf");
    outfile << "Channel_Realtimes[s]:" << "\t";
    for(i = 0;i<numchannels;i++){
		outfile << Form("%.12e",channelRealTime[i]) << "\t";
    }
    outfile << endl;
    outfile << "Ch/Det name" << "\t";
    
    for(i = 0;i<numchannels;i++){
		outfile << detName[i].Data() << "\t";
    }
    
    outfile << endl;
    
    for(k = 1;k<=16384;k++){
		outfile << (k-1) << "\t";
		for(i = 0;i<numchannels;i++){
			outfile << histo_raw[i]->GetBinContent(k) << "\t";
		}
	 
	 outfile << endl;
     
    }
    
    outfile.close();
    * /*/


    cout << "ROOT script effort: Real time (in seconds)= " << rootTime->RealTime() << " CPU time (in seconds)= " << rootTime->CpuTime() << endl;

}


void MakeTree_230526(){
	//Main function
	TString run = "run003";				//run id
	BinToRoot(run, false, false);
	//		  run, hist,  drawh, more information within 'BinToRoot'
}


