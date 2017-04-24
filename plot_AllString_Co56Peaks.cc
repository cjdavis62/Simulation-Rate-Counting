/* 
This script will take in the values from the efficiency file created by Fitting.py and then generate the root file with the rates for each peak

This script requires the input of the total time from the calibration. This is seen in the Time histogram from the g4cuore processed file.

Run the script as
>root
>.L plot_AllString_Co56Peaks.cc
>plot_AllString_Co56Peaks(Time)

Written by: Christopher Davis
christopher.davis@yale.edu
*/


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <utility>
#include "THStack.h"
#include "TMath.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TAxis.h"


using std::cout;
using std::cin;
using std::endl;
using namespace RooFit;


void plot_AllString_Co56Peaks(int Time) {
  
  ifstream EfficiencyFile;
  EfficiencyFile.open("Co56_PeakFits/Output/Efficiency.dat");
  string line;
 
  std::vector<Double_t> Peak;
  std::vector<Double_t> Efficiency;

  while (getline(EfficiencyFile, line))
    {
      istringstream ss(line);
      
      Double_t var1;
      Double_t var2;
      
      ss >> var1 >> var2;
      //    cout << var1 << " " << var2 << "\n";
      
      Peak.push_back(var1);
      Efficiency.push_back(var2);
    }
  cout << "Efficiencies:" << endl;
  for (std::vector<Double_t>::const_iterator i = Efficiency.begin(); i != Efficiency.end(); ++i)
    {
      std::cout << *i << " " << endl;;
    }
  
  TCut multiplicity = "Multiplicity == 1";

  // open the file and the tree
  TFile* f1 = new TFile("/data-mgm/cuore/scratch/simulation_scratch/Co56/Co56_w_r126_g4cuore.root");
  TTree* t1 = (TTree*)f1->Get("outTree");
  
  int nbins = 988;
  int energy_bins = 200;
  double time_scaling = 3600.0/double(Time); // scales from events to events per hour
  cout << "Time Scaling: " << time_scaling << endl;
  
  double eventsToCalibrate = 50; // How many events to require for a channel to be calibrated

  // Empty histograms to be filled by tree
  TH1F* Peak2599 = new TH1F("Peak2599", "Peak2599", nbins, 0, 988);
  TH1F* Peak847 = new TH1F("Peak847", "Peak847", nbins, 0, 988);
  TH1F* Peak1238 = new TH1F("Peak1238", "Peak1238", nbins, 0, 988);
  TH1F* Peak511 = new TH1F("Peak511", "Peak511", nbins, 0, 988);
  TH1F* Peak1771 = new TH1F("Peak1771", "Peak1771", nbins, 0, 988);
  TH1F* Peak1037 = new TH1F("Peak1037", "Peak1037", nbins, 0, 988);
  TH1F* Peak3254 = new TH1F("Peak3254", "Peak3254", nbins, 0, 988);
  TH1F* Peak2035 = new TH1F("Peak2035", "Peak2035", nbins, 0, 988);
  TH1F* Peak1360 = new TH1F("Peak1360", "Peak1360", nbins, 0, 988);
  TH1F* Peak3202 = new TH1F("Peak3202", "Peak3202", nbins, 0, 988);
  TH1F* Peak3451 = new TH1F("Peak3451", "Peak3451", nbins, 0, 988);

  // Stacked histogram to contain the number of events at each peak
  THStack *hs = new THStack("hs", "All Peaks");

  // Make cuts at the energies of each peak
  TCut cut3 = "Ener1 > 2589";
  TCut cut4 = "Ener1 < 2609";
  TCut cut2599 = cut3 && cut4 && multiplicity;

  TCut cut5 = "Ener1 > 837";
  TCut cut6 = "Ener1 < 857";
  TCut cut847 = cut5 && cut6 && multiplicity;

  TCut cut7 = "Ener1 > 1228";
  TCut cut8 = "Ener1 < 1248";
  TCut cut1238 = cut7 && cut8 && multiplicity;

  TCut cut9 = "Ener1 > 501";
  TCut cut10 = "Ener1 < 523";
  TCut cut511 = cut9 && cut10 && multiplicity;

  TCut cut11 = "Ener1 > 1761";
  TCut cut12 = "Ener1 < 1781";
  TCut cut1771 = cut11 && cut12 && multiplicity;

  TCut cut13 = "Ener1 > 1027";
  TCut cut14 = "Ener1 < 1047";
  TCut cut1037 = cut13 && cut14 && multiplicity;

  TCut cut15 = "Ener1 > 3244";
  TCut cut16 = "Ener1 < 3264";
  TCut cut3254 = cut15 && cut16 && multiplicity;
  
  TCut cut17 = "Ener1 > 2000";
  TCut cut18 = "Ener1 < 2045";
  TCut cut2035 = cut17 && cut18 && multiplicity;
  
  TCut cut19 = "Ener1 > 1350";
  TCut cut20 = "Ener1 < 1370";
  TCut cut1360 = cut19 && cut20 && multiplicity;

  TCut cut21 = "Ener1 > 3192";
  TCut cut22 = "Ener1 < 3212";
  TCut cut3202 = cut21 && cut22 && multiplicity;
  
  TCut cut23 = "Ener1 > 3441";
  TCut cut24 = "Ener1 < 3461";
  TCut cut3451 = cut23 && cut24 && multiplicity;
 
  TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);
  c4->cd();
  cout << "Drawing 2599:" << endl;
  t1->Draw("Channel >> Peak2599", cut2599);
  cout << "Drawing 847:" << endl;
  t1->Draw("Channel >> Peak847", cut847);
  cout << "Drawing 1238:" << endl;
  t1->Draw("Channel >> Peak1238", cut1238);
  cout << "Drawing 511:" << endl;
  t1->Draw("Channel >> Peak511", cut511);
  cout << "Drawing 1771:" << endl;
  t1->Draw("Channel >> Peak1771", cut1771);
  cout << "Drawing 1037:" << endl;
  t1->Draw("Channel >> Peak1037", cut1037);
  cout << "Drawing 3254:" << endl;
  t1->Draw("Channel >> Peak3254", cut3254);
  cout << "Drawing 2035:" << endl;
  t1->Draw("Channel >> Peak2035", cut2035);
  cout << "Drawing 1360:" << endl;
  t1->Draw("Channel >> Peak1360", cut1360);
  cout << "Drawing 3202:" << endl;
  t1->Draw("Channel >> Peak3202", cut3202);
  cout << "Drawing 3451:" << endl;
  t1->Draw("Channel >> Peak3451", cut3451);

   // Reduce each peak by their efficiency
  Peak2599->Scale(Efficiency[2]);
  Peak847->Scale(Efficiency[3]);
  Peak1238->Scale(Efficiency[4]);
  Peak511->Scale(Efficiency[5]);
  Peak1771->Scale(Efficiency[6]);
  Peak1037->Scale(Efficiency[7]);
  Peak3254->Scale(Efficiency[0]);
  Peak2035->Scale(Efficiency[1]);
  Peak1360->Scale(Efficiency[8]);
  Peak3202->Scale(Efficiency[9]);
  Peak3451->Scale(Efficiency[10]);
  
  Peak1771->SetLineColor(kAzure);
  Peak2599->SetLineColor(kSpring);
  Peak847->SetLineColor(kMagenta);
  Peak1238->SetLineColor(kOrange);
  Peak511->SetLineColor(kRed);
  Peak1037->SetLineColor(kCyan);
  Peak3254->SetLineColor(kPink);
  Peak2035->SetLineColor(kViolet);
  Peak1360->SetLineColor(kGreen);
  Peak3202->SetLineColor(kBlack);
  Peak3451->SetLineColor(kYellow);

  hs->Add(Peak2599);
  hs->Add(Peak847);
  hs->Add(Peak1238);
  hs->Add(Peak511);
  hs->Add(Peak1771);
  hs->Add(Peak1037);
  hs->Add(Peak3254);
  hs->Add(Peak2035);
  hs->Add(Peak1360);
  hs->Add(Peak3202);
  hs->Add(Peak3451);

  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("Channel");
  hs->GetYaxis()->SetTitle("Peak Events/Channel");
  hs->GetXaxis()->SetNdivisions(13,19, kFALSE);


  //create tree of events

  TTree *tree = new TTree("tree","An example of a ROOT tree");

  Int_t Channel;

  Double_t Rate_2599;
  Double_t Rate_847;
  Double_t Rate_1238;
  Double_t Rate_511;
  Double_t Rate_1771;
  Double_t Rate_1037;
  Double_t Rate_3254;
  Double_t Rate_2035;
  Double_t Rate_1360;
  Double_t Rate_3202;
  Double_t Rate_3451;

  Double_t Rate_Min; // rate for 1 peak calibrated
  Double_t Rate_Two; // rate for 2 peaks calibrated
  Double_t Rate_Three; // rate for 3 peaks calibrated
  Double_t Rate_Four; // rate for 4 peaks calibrated
  Double_t Rate_Five; // rate for 5 peaks calibrated
  Double_t Rate_Six; // rate for 6 peaks calibrated
  Double_t Rate_Seven; // rate for 7 peaks calibrated
  Double_t Rate_Eight; // rate for 8 peaks calibrated
  Double_t Rate_Nine; // rate for 9 peaks calibrated
  Double_t Rate_Ten; // rate for 10 peaks calibrated
  Double_t Rate_Max; //rate for 11 peaks calibrated

  Double_t Time_2599;
  Double_t Time_847;
  Double_t Time_1238;
  Double_t Time_511;
  Double_t Time_1771;
  Double_t Time_1037;
  Double_t Time_3254;
  Double_t Time_2035;
  Double_t Time_1360;
  Double_t Time_3202;
  Double_t Time_3451;

  Double_t Time_Min; //time for 1 peak calibrated
  Double_t Time_Two;
  Double_t Time_Three;
  Double_t Time_Four;
  Double_t Time_Five;
  Double_t Time_Six;
  Double_t Time_Seven;
  Double_t Time_Eight;
  Double_t Time_Nine;
  Double_t Time_Ten;
  Double_t Time_Max;

  Double_t Events_2599;
  Double_t Events_847;
  Double_t Events_1238;
  Double_t Events_511;
  Double_t Events_1771;
  Double_t Events_1037;
  Double_t Events_3254;
  Double_t Events_2035;
  Double_t Events_1360;
  Double_t Events_3202;
  Double_t Events_3451;

  Double_t Events_Min;
  Double_t Events_Two;
  Double_t Events_Three;
  Double_t Events_Four;
  Double_t Events_Five;
  Double_t Events_Six;
  Double_t Events_Seven;
  Double_t Events_Eight;
  Double_t Events_Nine;
  Double_t Events_Ten;
  Double_t Events_Max;


  // Make the Branches for the tree (to be filled later)
  tree->Branch("Rate_2599", &Rate_2599, "Rate_2599/D");
  tree->Branch("CalibrationTime_2599", &Time_2599, "CalibrationTime_2599/D");
  tree->Branch("Rate_847", &Rate_847, "Rate_847/D");
  tree->Branch("CalibrationTime_847", &Time_847, "CalibrationTime_847/D");
  tree->Branch("Rate_1238", &Rate_1238, "Rate_1238/D");
  tree->Branch("CalibrationTime_1238", &Time_1238, "CalibrationTime_1238/D");
  tree->Branch("Rate_511", &Rate_511, "Rate_511/D");
  tree->Branch("CalibrationTime_511", &Time_511, "CalibrationTime_511/D");
  tree->Branch("Rate_1771", &Rate_1771, "Rate_1771/D");
  tree->Branch("CalibrationTime_1771", &Time_1771, "CalibrationTime_1771/D");
  tree->Branch("Rate_1037", &Rate_1037, "Rate_1037/D");
  tree->Branch("CalibrationTime_1037", &Time_1037, "CalibrationTime_1037/D");
   tree->Branch("Rate_3254", &Rate_3254, "Rate_3254/D");
  tree->Branch("CalibrationTime_3254", &Time_3254, "CalibrationTime_3254/D");
  tree->Branch("Rate_2035", &Rate_2035, "Rate_2035/D");
  tree->Branch("CalibrationTime_2035", &Time_2035, "CalibrationTime_2035/D");
  tree->Branch("Rate_1360", &Rate_1360, "Rate_1360/D");
  tree->Branch("CalibrationTime_1360", &Time_1360, "CalibrationTime_1360/D");
  tree->Branch("Rate_3202", &Rate_3202, "Rate_3202/D");
  tree->Branch("CalibrationTime_3202", &Time_3202, "CalibrationTime_3202/D");
  tree->Branch("Rate_3451", &Rate_3451, "Rate_3451/D");
  tree->Branch("CalibrationTime_3451", &Time_3451, "CalibrationTime_3451/D");
  tree->Branch("Rate_Min", &Rate_Min, "Rate_Min/D");
  tree->Branch("CalibrationTime_Min", &Time_Min, "CalibrationTime_Min/D");
  tree->Branch("Rate_Two", &Rate_Two, "Rate_Two/D");
  tree->Branch("CalibrationTime_Two", &Time_Two, "CalibrationTime_Two/D");
  tree->Branch("Rate_Three", &Rate_Three, "Rate_Three/D");
  tree->Branch("CalibrationTime_Three", &Time_Three, "CalibrationTime_Three/D");
  tree->Branch("Rate_Four", &Rate_Four, "Rate_Four/D");
  tree->Branch("CalibrationTime_Four", &Time_Four, "CalibrationTime_Four/D");
  tree->Branch("Rate_Five", &Rate_Five, "Rate_Five/D");
  tree->Branch("CalibrationTime_Five", &Time_Five, "CalibrationTime_Five/D");
  tree->Branch("Rate_Six", &Rate_Six, "Rate_Six/D");
  tree->Branch("CalibrationTime_Six", &Time_Six, "CalibrationTime_Six/D");
  tree->Branch("Rate_Seven", &Rate_Seven, "Rate_Seven/D");
  tree->Branch("CalibrationTime_Seven", &Time_Seven, "CalibrationTime_Seven/D");
  tree->Branch("Rate_Eight", &Rate_Eight, "Rate_Eight/D");
  tree->Branch("CalibrationTime_Eight", &Time_Eight, "CalibrationTime_Eight/D");
  tree->Branch("Rate_Nine", &Rate_Nine, "Rate_Nine/D");
  tree->Branch("CalibrationTime_Nine", &Time_Nine, "CalibrationTime_Nine/D");
  tree->Branch("Rate_Ten", &Rate_Ten, "Rate_Ten/D");
  tree->Branch("CalibrationTime_Ten", &Time_Ten, "CalibrationTime_Ten/D");
  tree->Branch("Rate_Max", &Rate_Max, "Rate_Max/D");
  tree->Branch("CalibrationTime_Max", &Time_Max, "CalibrationTime_Max/D");

  tree->Branch("Channel", &Channel, "Channel/I");

  // Loop over the values in each bin to fill the tree
  for (int n = 0; n < nbins; n++) {
    
    if (n % 247 == 0) {
      cout << (n / 247) * 25 << "% through loop" << endl;
    }
    
    //Get channel
    Channel = TMath::Ceil(Peak2599->GetBinCenter(n+1));

    //get # of events for each peak and find max
    Events_2599 = Peak2599->GetBinContent(n+1);
    Events_847 = Peak847->GetBinContent(n+1);
    Events_1238 = Peak1238->GetBinContent(n+1);
    Events_511 = Peak511->GetBinContent(n+1);
    Events_1771 = Peak1771->GetBinContent(n+1);
    Events_1037 = Peak1037->GetBinContent(n+1);
    Events_3254 = Peak3254->GetBinContent(n+1);
    Events_2035 = Peak2035->GetBinContent(n+1);
    Events_1360 = Peak1360->GetBinContent(n+1);
    Events_3202 = Peak3202->GetBinContent(n+1);
    Events_3451 = Peak3451->GetBinContent(n+1);
    

    Double_t Events[11] = {Events_2599, Events_847, Events_1238, Events_1771, Events_1037, Events_3254, Events_2035, Events_1360, Events_3202, Events_3451}; 

    //cout << Events[11] << endl;

    // sort algorithm. There are definitely better ways to do this. But the ones I see are for c++0x/11, which is sad. So here we are with the dumb sort :`(
    if (Events_511 >= Events_2599) {
      Events_Max = Events_511;
      Events_Two = Events_2599;
    }
    else {
      Events_Max = Events_2599;
      Events_Two = Events_511;
    }
    if (Events_847 >= Events_Max) {
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_847;
    }
    else if (Events_847 >= Events_Two) {
      Events_Three = Events_Two;
      Events_Two = Events_847;
    }
    else {
      Events_Three = Events_847;
    }
    if (Events_1238 >= Events_Max) {
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_1238;
    }
    else if (Events_1238 >= Events_Two) {
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_1238;
    }
    else if (Events_1238 >= Events_Three) {
      Events_Four = Events_Three;
      Events_Three = Events_1238;
    }
    else {
      Events_Four = Events_1238;
    }
    if (Events_1771 >= Events_Max) {
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_1771;
    }
    else if (Events_1771 >= Events_Two) {
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_1771;
    }
    else if (Events_1771 >= Events_Three) {
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_1771;
    }
    else if (Events_1771 >= Events_Four) {
      Events_Five = Events_Four;
      Events_Four = Events_1771;
    }
    else {
      Events_Five = Events_1771;
    }
    if (Events_1037 >= Events_Max) {
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_1037;
    }
    else if (Events_1037 >= Events_Two) {
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_1037;
    }
    else if(Events_1037 >= Events_Three) {
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_1037;
    }
    else if(Events_1037 >= Events_Four) {
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_1037;
    }
    else if(Events_1037 >= Events_Five) {
      Events_Six = Events_Five;
      Events_Five = Events_1037;
    }
    else {
      Events_Six = Events_1037;
    }
    if(Events_3254 >= Events_Max) {
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_3254;
    }
    else if (Events_3254 >= Events_Two) {
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_3254;
    }
    else if (Events_3254 >= Events_Three) {
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Events_3254;
    }
    else if (Events_3254 >= Events_Four) {
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_3254;
    }
    else if (Events_3254 >= Events_Five) {
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_3254;
    }
    else if (Events_3254 >= Events_Six) {
      Events_Seven = Events_Six;
      Events_Six = Events_3254;
    }
    else {
      Events_Seven = Events_3254;
    }
    if (Events_2035 >= Events_Max) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_2035;
    }
    else if (Events_2035 >= Events_Two) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_2035;
    }
    else if (Events_2035 >= Events_Three) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_2035;
    }
    else if (Events_2035 >= Events_Four) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_2035;
    }
    else if (Events_2035 >= Events_Five) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_2035;
    }
    else if (Events_2035 >= Events_Six) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_2035;
    }
    else if (Events_2035 >= Events_Seven) {
      Events_Eight = Events_Seven;
      Events_Seven = Events_2035;
    }
    else {
      Events_Eight = Events_2035;
    }
    if (Events_1360 >= Events_Max) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_1360;
    }
    else if (Events_1360 >= Events_Two) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_1360;
    }
    else if (Events_1360 >= Events_Three) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_1360;
    }
    else if (Events_1360 >= Events_Four) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_1360;
    }
    else if (Events_1360 >= Events_Five) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_1360;
    }
    else if (Events_1360 >= Events_Six) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_1360;
    }
    else if (Events_1360 >= Events_Seven) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_1360;
    }
    else if (Events_1360 >= Events_Eight) {
      Events_Nine = Events_Eight;
      Events_Eight = Events_1360;
    }
    else {
      Events_Nine = Events_1360;
    }
    if (Events_3202 >= Events_Max) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_3202;
    }
    else if (Events_3202 >= Events_Two) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_3202;
    }
    else if (Events_3202 >= Events_Three) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_3202;
    }
    else if (Events_3202 >= Events_Four) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_3202;
    }
    else if (Events_3202 >= Events_Five) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_3202;
    }
    else if (Events_3202 >= Events_Six) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_3202;
    }
    else if (Events_3202 >= Events_Seven) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_3202;
    }
    else if (Events_3202 >= Events_Eight) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_3202;
    }
    else if (Events_3202 >= Events_Nine) {
      Events_Ten = Events_Nine;
      Events_Nine = Events_3202;
    }
    else {
      Events_Ten = Events_3202;
    }

    if (Events_3451 >= Events_Max) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_Max;
      Events_Max = Events_3451;
    }
    else if (Events_3451 >= Events_Two) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_Two;
      Events_Two = Events_3451;
    }
    else if (Events_3451 >= Events_Three) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_Three;
      Events_Three = Events_3451;
    }
    else if (Events_3451 >= Events_Four) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_Four;
      Events_Four = Events_3451;
    }
    else if (Events_3451 >= Events_Five) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_Five;
      Events_Five = Events_3451;
    }
    else if (Events_3451 >= Events_Six) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_Six;
      Events_Six = Events_3451;
    }
    else if (Events_3451 >= Events_Seven) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_Seven;
      Events_Seven = Events_3451;
    }
    else if (Events_3451 >= Events_Eight) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_Eight;
      Events_Eight = Events_3451;
    }
    else if (Events_3451 >= Events_Nine) {
      Events_Min = Events_Ten;
      Events_Ten = Events_Nine;
      Events_Nine = Events_3451;
    }
    else if (Events_3451 >= Events_Ten) {
      Events_Min = Events_Ten;
      Events_Ten = Events_3451;
    }
    else {
      Events_Min = Events_3451;
    }


    Rate_2599 = Events_2599 * time_scaling * 0.2778;
    Rate_847 = Events_847 * time_scaling * 0.2778;
    Rate_1238 = Events_1238 * time_scaling * 0.2778;
    Rate_511 = Events_511 * time_scaling * 0.2778;
    Rate_1771 = Events_1771 * time_scaling * 0.2778;
    Rate_1037 = Events_1037 * time_scaling * 0.2778;      
    Rate_3254 = Events_3254 * time_scaling * 0.2778;
    Rate_2035 = Events_2035 * time_scaling * 0.2778;
    Rate_1360 = Events_1360 * time_scaling * 0.2778;
    Rate_3202 = Events_3202 * time_scaling * 0.2778;
    Rate_3451 = Events_3451 * time_scaling * 0.2778;


    //cout << "Rate 2599: " << Rate_2599 << " Events_2599: " << Events_2599 << endl;

    Rate_Max = Events_Max * time_scaling * 0.2778;
    Rate_Two = Events_Two * time_scaling * 0.2778;
    Rate_Three = Events_Three * time_scaling * 0.2778;
    Rate_Four = Events_Four * time_scaling * 0.2778;
    Rate_Five = Events_Five * time_scaling * 0.2778;
    Rate_Six = Events_Six * time_scaling * 0.2778;
    Rate_Seven = Events_Seven * time_scaling * 0.2778;
    Rate_Eight = Events_Eight * time_scaling * 0.2778;
    Rate_Nine = Events_Nine * time_scaling * 0.2778;
    Rate_Ten = Events_Ten * time_scaling * 0.2778;   
    Rate_Min = Events_Min * time_scaling * 0.2778;


    Time_2599 = eventsToCalibrate / (86.0 * Rate_2599);
    Time_847 = eventsToCalibrate / (86.0 * Rate_847);
    Time_1238 = eventsToCalibrate / (86.0 * Rate_1238);
    Time_511 = eventsToCalibrate / (86.0 * Rate_511);
    Time_1771 = eventsToCalibrate / (86.0 * Rate_1771);
    Time_1037 = eventsToCalibrate / (86.0 * Rate_1037);
    Time_3254 = eventsToCalibrate / (86.0 * Rate_3254);
    Time_2035 = eventsToCalibrate / (86.0 * Rate_2035);
    Time_1360 = eventsToCalibrate / (86.0 * Rate_1360);
    Time_3202 = eventsToCalibrate / (86.0 * Rate_3202);
    Time_3451 = eventsToCalibrate / (86.0 * Rate_3451);

    Time_Max = eventsToCalibrate / (86.0 * Rate_Min);
    Time_Two = eventsToCalibrate / (86.0 * Rate_Two);
    Time_Three = eventsToCalibrate / (86.0 * Rate_Three);
    Time_Four = eventsToCalibrate / (86.0 * Rate_Four);
    Time_Five = eventsToCalibrate / (86.0 * Rate_Five);
    Time_Six = eventsToCalibrate / (86.0 * Rate_Six);
    Time_Seven = eventsToCalibrate / (86.0 * Rate_Seven);
    Time_Eight = eventsToCalibrate / (86.0 * Rate_Eight);
    Time_Nine = eventsToCalibrate / (86.0 * Rate_Nine);
    Time_Ten = eventsToCalibrate / (86.0 * Rate_Ten);
    Time_Min = eventsToCalibrate / (86.0 * Rate_Max);


    tree->Fill();
  }

  TFile* output = new TFile("AllPeaks_rate.root","Recreate");
  tree->Write();
  output->Close();
}
