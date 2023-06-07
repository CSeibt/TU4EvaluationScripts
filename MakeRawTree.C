
#include <dirent.h>

#define nentriesMAX 20000000

vector<FILE*> GetInputFiles(TString input)
{
    TString caenFileName;
    // Count the .bin files
    // First search for all files in the directory "folder" with the extension "bin"
    TString binfiles[100];
    int total_number_files = 0;
    int b = 0;
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
    
    for(int i = 0; i < total_number_files; i++){
      cout << binfiles[i] << endl;
    }
    cout << " " << endl;

    vector<FILE*> fin;
    // Reading the files
    // Attention!!! The input_ch will be converted
    for (int i=0;i<total_number_files;i++) {
		caenFileName = input + binfiles[i];
        FILE *binFile = fopen(caenFileName, "r");
        if (!binFile) {
             cout << "Error when opening file " << caenFileName << endl;
             return fin;
        } else {
            fin.push_back(binFile);
            cout << "Opened file " << caenFileName << endl;
        }
    }
    cout << " " << endl;
    return fin;
}

TTree* MakePrelimDataTree(TString run, vector<FILE*> fin)
{

	double maxtime = 864000.;
	double mintime = 0.;

    int total_number_events = 0;

    const int numchannels = 8;
    const int max_num_jumps = 10;

    // Variable definitions.
    
    TString detName[8]= {"TU5","TU4","","","","","",""};
    int det_ch[8] = {0,3,1,4,2,5,6,7}; // Reordering of detector channels: TU5 = det0, TU4=det1, Sz2=det2 etc.
    
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
    ULong64_t timeLast[numchannels];
    Long64_t  EnergyDepLast[numchannels];
    ULong64_t timeAdopt=0;
    
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
		EnergyDepLast[i] = 0.;
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
	
    int nentries = 0;
    int 	*index 		= new int[nentriesMAX];
    uint64_t 	*timeAdoptW 	= new uint64_t[nentriesMAX];
    uint16_t 	*detW		= new uint16_t[nentriesMAX];
    int16_t 	*adcW		= new int16_t[nentriesMAX];
    uint16_t 	*saturationW	= new uint16_t[nentriesMAX];
    uint16_t 	*pileupW	= new uint16_t[nentriesMAX];
    uint32_t 	*extrasW	= new uint32_t[nentriesMAX];

	cout << "Temporary variables defined" << endl;

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

    do{
		nentries = 0;
		//fill the last read value in the array if any
		for (int filen=0; filen<fin.size(); filen++) {
			if(timeAdoptW[nentriesMAX-1-filen]!=0){
				cout << "Fill the last read value in the array... " << endl;
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
		for (int filen=0;filen<fin.size();filen++) {
			//cout << "test " << endl;
			//zero the temp variables
			timeAdoptW[nentriesMAX-1-filen] = 0;
			readTimeTag[0] 		= 0;
			
			cout << "Reading data from file... " << endl;
			// Read in listmode file, line by line.
			while (fread(&boardN, 1, sizeof(boardN), fin.at(filen)) 
					    &&fread(&channelN, 1, sizeof(channelN), fin.at(filen)) 
					    && (sizeof(readTimeTag) == fread(&readTimeTag, 1, sizeof(readTimeTag), fin.at(filen)))) 
					    {	
				
				fread(&readEnergy, 1, sizeof(readEnergy), fin.at(filen));
				fread(&readExtras, 1, sizeof(readExtras), fin.at(filen));
				fread(&numSamples, 1, sizeof(numSamples), fin.at(filen));
				fread(&samples ,numSamples, sizeof(samples), fin.at(filen));
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
		for (i=0;i<fin.size();i++) {
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
				
				Data->Fill();
			}
			
			}
		
		// Debugging output
		if (p%500000== 0) {
			cout << Form("%.1f",100.*p/nentries) <<"% done. " << "Time= " << timeAdopt/1.0E12 << " sec, now altogether " << Data->GetEntries()-1 << " preliminary entries."  << endl;
		}
		
    }while (NrOfFilesEnded<fin.size() && timeAdopt/1.0E12<=maxtime); 

    return Data;
}

void MakeRawTree(
    TString run = "run001",
    //TString folder = "/ZIH.fast/users/felsdaq/TUBunker/TU5/TU5_TU4_230524_coincidence/DAQ/",  //Path of the .bin file which shall be evaluated
    TString folder = "/home/hans/Uni/EC/TU5_TU4_coincidence/DAQ/"  //Path of the .bin file which shall be evaluated
)
{
	TString suffix = "/RAW/";
    TString input = folder + run + suffix;
    cout << "input: " << input << endl;

    vector<FILE*> fin = GetInputFiles(input);
    
    TFile *rootFile = new TFile(run+".root","RECREATE");
	cout << "created new file " << endl;
    TTree* prelimDataTree = MakePrelimDataTree(run, fin);
    rootFile->cd("/");
    prelimDataTree->Write("Data", TObject::kOverwrite);
    rootFile->Save();
    rootFile->Close();
}