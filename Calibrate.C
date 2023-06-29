#define V 0
#define MAXNROI 50

Double_t ROI[MAXNROI][5];
int nROI = 0;
Double_t g_guess_p1;

void PrintROIs()
{
    cout << nROI << " ROI: " << endl;
    cout << "E \tE0 \tE1 \tE2 \tE3" << endl;
    for (int i = 0; i < nROI; i++) {
        for (int j = 0; j < 5; j++)
            cout << ROI[i][j] << "\t";
        cout << endl;
    }
}

void ClearROIs()
{
    nROI = 0;
}

void AddROI(Double_t physE, Double_t E0, Double_t E1, Double_t E2, Double_t E3)
{ // Add a peak region of interest. All peak-evaluating functions will use it.
  // (E0, E1) left background interval
  // (E1, E2) peak interval
  // (E2, E3) right background interval
    ROI[nROI][0] = physE;
    ROI[nROI][1] = E0;
    ROI[nROI][2] = E1;
    ROI[nROI][3] = E2;
    ROI[nROI][4] = E3;
    nROI++;
    //PrintROIs();
}

Double_t GetPhysE(int iROI)
{ // Returns a peak's energy position literature value
    return ROI[iROI][0];
}

Double_t GetLimE(int iROI, int index)
{ // Returns a ROI's limits. index 0: left background lower limit, ...
    return ROI[iROI][1+index];
}

TF1* FitPeak(TH1D* spectrum, Double_t guess_p1, Double_t E, Double_t E0, Double_t E3, ofstream* output = 0)
{ // fit peak of uncalibrated spectrum in gven interval (low, up) of energies. 
  // Assume linear calibration with guessed factor p1=E/ch to convert energies into channels

    // convert limits from channels into energies
    Double_t low = E0 / guess_p1;
    Double_t up  = E3 / guess_p1;
    if (V) cout << "Fitting in channel interval " << low << " to " << up << endl;

    // Fit function: gaus + pol1
    TF1* fit_fcn = new TF1("fit_fcn", "gaus(0)+pol1(3)", low, up);

    // start values
    Double_t lbg = spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(low));
    Double_t rbg = spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(up));
    Double_t sbg = (rbg-lbg)/(up-low); // linear slope of background
    spectrum->GetXaxis()->SetRangeUser(low, up);
    Double_t mpos = spectrum->GetXaxis()->GetBinCenter(spectrum->GetMaximumBin());
    Double_t mval = spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(mpos));
    fit_fcn->SetParameters( mval - 0.5*(lbg+rbg),   // gauss constant
                            mpos,                   // gauss mean
                            1,                     // gauss sigma
                            lbg-sbg*low,            // pol1 constant
                            sbg);                   // linear
    
    // check start value quality
    //for (int par = 0; par < 5; par++) cout << "par " << par << " = " << fit_fcn->GetParameter(par) << endl;
    //new TCanvas();
    //spectrum->Draw("hist");
    //fit_fcn->Draw("same");

    // Fit
    spectrum->Fit("fit_fcn", "0LRQ");

    // draw fit
    string drawname = "fit_" + to_string(E);
    new TCanvas();
    gPad->SetLogy(0);
    TH1D* hdraw = (TH1D*)spectrum->Clone(drawname.c_str());
    hdraw->SetStats(0);
    hdraw->Draw("hist");
    fit_fcn->Draw("same");

    // output fit results
    if (V) cout << "Peak position = " 
         << fit_fcn->GetParameter(1)
         << " +/- " << fit_fcn->GetParError(1) << " (simplified uncertainty!)" << endl
         << "Peak width (sigma) = " << fit_fcn->GetParameter(2) << endl;
    Double_t factor = 2*sqrt(2*log(2)); // switch from std dev to FWHM
    
    if (V) cout << "FWHM = " 
         << factor * guess_p1*fit_fcn->GetParameter(2)
         << " +/- " << factor * guess_p1*fit_fcn->GetParError(2) << " (simplified uncertainty!)" << endl;
    // get peak content from fit
    //Double_t cnts_fit = sqrt(2*M_PI)*fit_fcn->GetParameter(0);
    //Double_t Dcnts_fit = sqrt(2*M_PI)*fit_fcn->GetParError(0);
    //if (V) cout << "Peak content (fit) = "
    //     << cnts_fit
    //     << " +/- " << Dcnts_fit << " (simplified uncertainty!)" << endl;

    // output results data line for calibration and FWHM graphs
    // format: channel energy FWHM Dchannel Denergy DFWHM
    if (output)
        (*output) << fit_fcn->GetParameter(1) << "\t" << E << "\t" << factor * guess_p1*fit_fcn->GetParameter(2) << "\t"
                  << fit_fcn->GetParError(1)  << "\t" << 0.01 << "\t" << factor * guess_p1*fit_fcn->GetParError(2) << endl;
    
    return fit_fcn;
}

void FillCalibratedHist(TH1D* hCh, TF1* fCal, TH1D* cal_hist)
{ // Apply calibration. Fill calibrated hist with Monte-Carlo rebinned entries of the uncalibrated spectrum 
    // init pseudo-random generator
    TRandom3* rnd = new TRandom3(time(NULL));
    for (Int_t bin = 1; bin <= hCh->GetNbinsX(); bin++) { // loop over adc channels
        Double_t bin_entries = hCh->GetBinContent(bin);
        Double_t low = hCh->GetXaxis()->GetBinLowEdge(bin);
        Double_t up = hCh->GetXaxis()->GetBinUpEdge(bin);
        for (Int_t entry = 0; entry < bin_entries; entry++) { // loop over channel entries
            // Define a pseudo-random rational channel number inside the bin boundaries
            Double_t rat_ch = rnd->Uniform(low, up);
            // Convert rational channel number into energy using calibration function
            Double_t E = fCal->Eval(rat_ch);
            // Fill energy histogram
            cal_hist->Fill(E);
        }
    }
}

TF1* Calibration(TH1D* hCh, Double_t guess_p1, string calib_name, TH1D* cal_hist = 0, Bool_t draw = true, ofstream* output = 0)
{ // Fit the calibration function for a spectrum using the global ROIs
    // Prepare calibration peaks data collection
    string calib_data_filename = "data_" + calib_name + ".txt";
    std::ofstream ofs(calib_data_filename);
    ofs << "# channel energy FWHM Dchannel Denergy DFWHM" << endl;

    // Prepare output string builder in case we need to write calibration fit results to a file
    std::stringstream width_string_builder;

    // Fit the peaks and write results into text file
    for (int peak = 0; peak < nROI; peak++) {
        Double_t physE = GetPhysE(peak);
        Double_t E0 = GetLimE(peak, 0);
        Double_t E3 = GetLimE(peak, 3);

        // Fit the peak /////////////////////////////////////////
        if (V) cout << "Calibration invoking peak fitter " << hCh->Integral() << "\t" << guess_p1 << "\t" << physE << "\t" << E0 << "\t" << E3 << endl;        
        TF1* fPeak = FitPeak(hCh, guess_p1, physE, E0, E3, &ofs);
        /////////////////////////////////////////////////////////

        if (output) { // build up output string
            Double_t factor = 2*sqrt(2*log(2)); // switch from std dev to FWHM
            width_string_builder << factor * guess_p1 * fPeak->GetParameter(2) << "\t"
                                 << factor * guess_p1 * fPeak->GetParError(2) << "\t";
        }
    }

    // Open the peaks' properties as graph
    TGraphErrors* gCal = new TGraphErrors(calib_data_filename.c_str(), "%lg %lg %*lg %lg %lg");
    string graphtitle = calib_name + " calibration; channel; E / keV";
    gCal->SetTitle(graphtitle.c_str());

    // Do the calibration fit
    TF1* fCal = new TF1("fCal", "pol1", 0, 16000);
    gCal->Fit("fCal", "0Q");

    if (draw) {
        // Draw
        string canvas_name = "Calibration_"+calib_name;
        new TCanvas(canvas_name.c_str());
        gCal->SetMarkerStyle(33);
        gCal->SetMarkerSize(2);
        gCal->Draw("ap");
        fCal->Draw("same");
        gCal->GetXaxis()->SetRangeUser(0, 12000);
        gCal->GetYaxis()->SetRangeUser(0, 3000);

        TGraphErrors* gWth = new TGraphErrors(calib_data_filename.c_str(), "%*lg %lg %lg %*lg %lg %lg");
        //Double_t m = fCal->GetParameter(1);
        //ScaleXY(gWth, 1.0, m);
        gWth->SetTitle("Width; E / keV; FWHM / keV");

        canvas_name = "Width_"+calib_name;
        new TCanvas(canvas_name.c_str());
        gWth->SetMarkerStyle(32);
        gWth->SetMarkerSize(2);
        gWth->Draw("ap");
        gWth->GetXaxis()->SetRangeUser(0, 3000);
        gWth->GetYaxis()->SetRangeUser(0, 4);
    }

    if (output) {
        (*output) << fCal->GetParameter(0) << "\t" << fCal->GetParError(0) << "\t"
                  << fCal->GetParameter(1) << "\t" << fCal->GetParError(1) << "\t"
                  << width_string_builder.str() << endl;
    } else {
        cout << calib_name << ": E/keV = " << fCal->GetParameter(1) << " * ch + " << fCal->GetParameter(0) << endl;
    }

    if (cal_hist) {
        FillCalibratedHist(hCh, fCal, cal_hist);
    }

    return fCal;
}

void CalibrateTU5_133Ba(TH1D* hCh)
// Do calibration for the TU5 Detector. 
// Includes all hard-coded preferences: ROIs, guessed calibration factor
{
    ClearROIs();
    AddROI(4.47, 3.5, 4.0, 4.5, 4.5);
    //AddROI(37.255, 36.8, 36.8, 38.3, 40);
    AddROI(80.9979, 60, 80, 82, 100);
    PrintROIs();

    // guess calibration
    Double_t guess_p1 = 0.0099;

    Calibration(hCh, guess_p1, "TU5");
}

void CalibrateTU4_133Ba(TH1D* hCh)
// Do calibration for the TU5 Detector. 
// Includes all hard-coded preferences: ROIs, guessed calibration factor, detector name
{
    ClearROIs();
    AddROI(80.9979, 70, 70, 90, 90);
    AddROI(276.3989, 260, 260, 290, 290);
    AddROI(302.8508, 290, 290, 320, 320);
    AddROI(356.0129, 340, 340, 370, 370);
    PrintROIs();

    // guess calibration
    Double_t guess_p1 = 0.33;

    Calibration(hCh, guess_p1, "TU4");
}

void Calibrate(
//    string filename = "/home/hans/Uni/EC/TU5_TU4_coincidence/TU4EvaluationScripts/run002.root"
    string filename = "run001.root"
)
{
    // Open root file
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file) cout << "Could not open " << filename << endl;
    // Get data tree
    TTree* data = (TTree*)file->Get("Data");

    // TU5 calibration
    // Get uncalibrated adc histograms
    Int_t Nch = 16000;
    TString histname = "hChTU5";
    TH1D* hChTU5 = new TH1D(histname.Data(), ";TU5 ADC channel; counts", Nch, 0, Nch);
    TString varexp = "adc>>"+histname;
    TString selection = "det==0 && pileup==0 && saturation==0";
    data->Draw(varexp.Data(), selection.Data(), "goff");
    // calibration
    CalibrateTU5_133Ba(hChTU5);

    // TU4 calibration
    // Get uncalibrated adc histograms
    histname = "hChTU4";
    TH1D* hChTU4 = new TH1D(histname.Data(), ";TU4 ADC channel; counts", Nch, 0, Nch);
    varexp = "adc>>"+histname;
    selection = "det==1 && pileup==0 && saturation==0";
    data->Draw(varexp.Data(), selection.Data(), "goff");
    // calibration
    CalibrateTU4_133Ba(hChTU4);
}