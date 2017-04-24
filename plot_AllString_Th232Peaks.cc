/* This script takes in a g4cuore processed ROOT File and performs fits to the peaks
These peaks are then analyzed for the ratio of signal/background events
With the number of signal events, it is then calculated how long it would take to get 50 events per channel per peak

This script in particular is for Th232 spectra in a file called AllString_g4cuore.root

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

// Calculate the ratio of signal events to the total number of events
Double_t Acceptance_lineargaus(Double_t offset, Double_t linear, Double_t amplitude, Double_t mean, Double_t sigma, Double_t peak_window)
{
  // Integrate the Background
  Double_t background = offset * peak_window + (0.5) * linear * ((peak_window/2.0 + mean)**2 - ((-peak_window/2.0 + mean)**2));

  // Integrate the Peak
  Double_t signal = 2 * amplitude * (sigma * TMath::Sqrt(TMath::Pi()/2)) * TMath::Erf((peak_window/2.0)/(sigma / TMath::Sqrt(2.0)));

  cout << signal << "\t" << background << endl;
  
  return (signal / (background+signal));
}

// Calculate the raio of signal events to the total number of events
Double_t Acceptance_lineardoublegaus(Double_t offset, Double_t linear, Double_t amplitude1, Double_t mean1, Double_t sigma1, Double_t amplitude2, Double_t mean2, Double_t sigma2, Double_t peak_window)
{
  // Integrate the Background

  Double_t background = offset * peak_window + (0.5) * linear * ((peak_window/2.0 + mean1)**2 - ((-peak_window/2.0 + mean1)**2));

  // Integrate the main peak
  Double_t signal1 = 2 * amplitude1 * (sigma1 * TMath::Sqrt(TMath::Pi()/2)) * TMath::Erf((peak_window/2.0)/(sigma1 / TMath::Sqrt(2.0)));
  // Integrate the subpeak
  Double_t signal2 = 0.5 * amplitude2 * ((sigma2 * TMath::Sqrt(TMath::Pi()/2)) * TMath::Erf((peak_window/2.0 - mean2 + mean1)/(sigma2 / TMath::Sqrt(2.0))) - (sigma2 * TMath::Sqrt(TMath::Pi()/2)) * TMath::Erf((-peak_window/2.0 -mean2 +mean1)/(sigma2 / TMath::Sqrt(2.0))));
 
  Double_t signal = signal1 + signal2;
  cout << signal1 << "\t" << signal2 << "\t" << signal << "\t" << background << endl;
  return (signal / (background+signal));
}

void plot_AllString_Th232Peaks(double Time) {

  int nbins = 988;
  int energy_bins = 200;
  double time_scaling = 3600.0/Time; // scales from events to events per hour

  double eventsToCalibrate = 50; // How many events to require for a channel to be calibrated

  Double_t peak_window = 20; // 20 keV window
  Double_t peak_window_338 = 30; // 30 keV window for 338 double peak
  int energy_bins_338 = int(double(peak_window_338 / peak_window) * energy_bins); // 338 keV is a bit special here
  
  TCut multiplicity = "Multiplicity == 1";

  TFile* f1 = new TFile("AllString_g4cuore.root");
  TTree* t1 = (TTree*)f1->Get("outTree");

  TH1F* Peak2615 = new TH1F("Peak2615", "Peak2615", nbins, 0, 988);
  TH1F* Peak969 = new TH1F("Peak969", "Peak969", nbins, 0, 988);
  TH1F* Peak911 = new TH1F("Peak911", "Peak911", nbins, 0, 988);
  TH1F* Peak583 = new TH1F("Peak583", "Peak583", nbins, 0, 988);
  TH1F* Peak338 = new TH1F("Peak338", "Peak338", nbins, 0, 988);
  TH1F* Peak239 = new TH1F("Peak239", "Peak239", nbins, 0, 988);
  
  TH1F* Energy2615 = new TH1F("Energy2615", "Energy2615", energy_bins, (2615 - peak_window/2), (2615 + peak_window/2));
  TH1F* Energy969 = new TH1F("Energy969", "Energy969", energy_bins, (969 - peak_window/2), (968 + peak_window/2));
  TH1F* Energy911 = new TH1F("Energy911", "Energy911", energy_bins, (911 - peak_window/2), (911 + peak_window/2));
  TH1F* Energy583 = new TH1F("Energy583", "Energy583", energy_bins, (583 - peak_window/2), (583 + peak_window/2));
  TH1F* Energy338 = new TH1F("Energy338", "Energy338", energy_bins_338, (338 - peak_window_338/2), (338 + peak_window_338/2));
  TH1F* Energy239 = new TH1F("Energy239", "Energy239", energy_bins, (239 - peak_window/2), (239 + peak_window/2));
    
  
  // Fit Types
  TF1 * lineargaus = new TF1("linear+gaus", "[0] + [4] * x + [1] * exp(-0.5*((x-[2])/[3])**2)", 0, 5);
  TF2 * lineardoublegaus = new TF1("linear+doublegaus", "[0] + [7] * x + [1] * exp(-0.5*((x-[2])/[3])**2) + [4]*[1] * exp(-0.5*((x-[5]*[2])/[6])**2)", 0, 5);

  // begin RooFits

  // Single Gaussian
  RooRealVar x("x", "x", 0, 4000);
  RooRealVar mean("mean", "mean of gaussian", 10, 0, 4000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0, 10);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  // Second Gaussian
  RooRealVar mean2("mean2", "mean of secondary gaussian", 0, 0, 4000);
  RooRealVar sigma2("sigma2", "sigma of secondary gaussian", 1, 0, 10);
  
  RooGaussian gauss2("gauss2", "secondary gaussian PDF", x, mean2, sigma2);
  
  // Linear background
  RooRealVar a0("a0", "a0", 1, -10000, 10000);
  RooRealVar a1("a1", "a1", 1, -10000, 10000);
  RooPolynomial p2("p2", "p2", x, RooArgList(a0, a1), 0);
 
  //Combined gauss + linear
  RooRealVar gaussfrac("gfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf gausslin("gausslin", "gauss+p2", RooArgList(gauss, p2), RooArgList(gaussfrac));

  //Combined doublegaus
  RooRealVar gaus1frac("g1frac", "fraction of main gaussian", 0.8, 0, 1);
  RooAddPdf doublegaus("doublegaus", "gauss+gauss2", RooArgList(gauss, gauss2), RooArgList(gaus1frac));

  //Combined doublegaus + linear
  RooRealVar doublegausfrac("doublegausfrac", "fraction of gaussians", 0.8, 0, 1);
  RooAddPdf doublegausslin("doublegausslin", "doublegaus+p2", RooArgList(doublegaus, p2), RooArgList(doublegausfrac));

  TCanvas * c5 = new TCanvas("c5", "Roofit", 1200, 1000);
  c5->cd();
  c5->Divide(3,2);
  c5->cd(1);
  t1->Draw("Ener1 >> Energy2615", multiplicity, "goff");
  x.setRange(2615 - peak_window/2, 2615 + peak_window/2);
  RooDataHist data2615("data2615", "2615 peak", x, Energy2615);

  RooPlot* frame2615 = x.frame(Title("RooPlot of x"));
  
  data2615.plotOn(frame2615);
  
  mean.setVal(2615);
  mean.setRange(2614, 2616);
  
  gausslin.fitTo(data2615);
  gausslin.plotOn(frame2615);
  gausslin.plotOn(frame2615, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();
  
  Double_t efficiency_2615 = gaussfrac.getVal();

  frame2615->Draw();

  TCanvas * c4 = new TCanvas("c4", "ROOT fit", 1200, 1000);
  c4->cd();
  c4->Divide(3,2);
  c4->cd(1);
  //  t1->Draw("Ener1 >> Energy2615", multiplicity);
  lineargaus->SetParameter(2, 2615);
  lineargaus->SetParameter(3, 5);
  Energy2615->Fit("linear+gaus");
  Double_t offset_2615 = lineargaus->GetParameter(0);
  Double_t linear_2615 = lineargaus->GetParameter(4);
  Double_t amplitude_2615 = lineargaus->GetParameter(1);
  Double_t mean_2615 = lineargaus->GetParameter(2);
  Double_t sigma_2615 = lineargaus->GetParameter(3);

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  //Double_t efficiency_2615 = Acceptance_lineargaus(offset_2615, linear_2615, amplitude_2615, mean_2615, sigma_2615, peak_window);
  //cout << efficiency_2615 << endl;

  c4->cd(2);
  t1->Draw("Ener1 >> Energy969", multiplicity);
  lineardoublegaus->SetParameter(0, 9876.4);
  lineardoublegaus->SetParameter(1, 6471.7);
  lineardoublegaus->SetParameter(2, 969);
  lineardoublegaus->SetParameter(3, 2.1);
  lineardoublegaus->FixParameter(4, (4.99 / 15.8));
  lineardoublegaus->FixParameter(5, 964.766 / 968.971);
  lineardoublegaus->SetParameter(6, 1);
  lineardoublegaus->SetParameter(7, -9.8625);
  Energy969->Fit("linear+doublegaus");
  Double_t offset_969 = lineardoublegaus->GetParameter(0);
  Double_t linear_969 = lineardoublegaus->GetParameter(7);
  Double_t amplitude1_969 = lineardoublegaus->GetParameter(1);
  Double_t mean1_969 = lineardoublegaus->GetParameter(2);
  Double_t sigma1_969 = lineardoublegaus->GetParameter(3);
  Double_t amplitude2_969 = lineardoublegaus->GetParameter(1) * lineardoublegaus->GetParameter(4);
  Double_t mean2_969 = lineardoublegaus->GetParameter(5) * lineardoublegaus->GetParameter(2);
  Double_t sigma2_969 = lineardoublegaus->GetParameter(6);


  for (int i = 0; i < 8; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  //  Double_t efficiency_969 = Acceptance_lineardoublegaus(offset_969, linear_969,  amplitude1_969, mean1_969, sigma1_969, amplitude2_969, mean2_969, sigma2_969, peak_window);
  // cout << efficiency_969;

  c5->cd(2);

  x.setRange(969 - peak_window/2, 969 + peak_window/2);
  RooDataHist data969("data969", "969 peak", x, Energy969);

  RooPlot* frame969 = x.frame(Title("RooPlot of x"));
  
  data969.plotOn(frame969);
  
  mean.setVal(969);
  mean.setRange(967, 972);

  mean2.setVal(965);
  mean2.setRange(964,967);
  
  doublegausslin.fitTo(data969);
  doublegausslin.plotOn(frame969);
  doublegausslin.plotOn(frame969, Components(p2), LineStyle(kDashed));

  gaus1frac.Print();
  doublegausfrac.Print();
  mean.Print();
  sigma.Print();
  mean2.Print();
  sigma2.Print();
  a0.Print();
  a1.Print();
  
  Double_t efficiency_969 = doublegausfrac.getVal();

  frame969->Draw();

  c4->cd(3);
  t1->Draw("Ener1 >> Energy911", multiplicity);
  lineardoublegaus->SetParameter(2, 911);
  lineardoublegaus->SetParameter(3, 5);
  lineardoublegaus->FixParameter(4, (0.77 / 25.8));
  lineardoublegaus->FixParameter(5, (904.2 / 911.204));
  lineardoublegaus->SetParameter(6, 1);
  Energy911->Fit("linear+doublegaus");
  Double_t offset_911 = lineardoublegaus->GetParameter(0);
  Double_t linear_911 = lineardoublegaus->GetParameter(7);
  Double_t amplitude1_911 = lineardoublegaus->GetParameter(1);
  Double_t mean1_911 = lineardoublegaus->GetParameter(2);
  Double_t sigma1_911 = lineardoublegaus->GetParameter(3);
  Double_t amplitude2_911 = lineardoublegaus->GetParameter(1) * lineardoublegaus->GetParameter(4);
  Double_t mean2_911 = lineardoublegaus->GetParameter(5) * lineardoublegaus->GetParameter(2);
  Double_t sigma2_911 = lineardoublegaus->GetParameter(6);

  //  Double_t efficiency_911 = Acceptance_lineardoublegaus(offset_911, linear_911,  amplitude1_911, mean1_911, sigma1_911, amplitude2_911, mean2_911, sigma2_911, peak_window);

  for (int i = 0; i < 8; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  c5->cd(3);

  x.setRange(911 - peak_window/2, 911 + peak_window/2);
  RooDataHist data911("data911", "911 peak", x, Energy911);

  RooPlot* frame911 = x.frame(Title("RooPlot of x"));
  
  data911.plotOn(frame911);
  
  mean.setVal(911);
  mean.setRange(911-2, 911+2);

  mean2.setVal(905);
  mean2.setRange(903,910);
  
  doublegausslin.fitTo(data911);
  doublegausslin.plotOn(frame911);
  doublegausslin.plotOn(frame911, Components(p2), LineStyle(kDashed));

  gaus1frac.Print();
  doublegausfrac.Print();
  mean.Print();
  sigma.Print();
  mean2.Print();
  sigma2.Print();
  a0.Print();
  a1.Print();
  
  Double_t efficiency_911 = doublegausfrac.getVal();

  frame911->Draw();
  
  c4->cd(4);

  t1->Draw("Ener1 >> Energy583", multiplicity, "goff");

  lineargaus->SetParameter(2,583);
  lineargaus->SetParameter(3,5);
  Energy583->Fit("linear+gaus");
  Double_t offset_583 = lineargaus->GetParameter(0);
  Double_t linear_583 = lineargaus->GetParameter(4);
  Double_t amplitude_583 = lineargaus->GetParameter(1);
  Double_t mean_583 = lineargaus->GetParameter(2);
  Double_t sigma_583 = lineargaus->GetParameter(3);

  // Double_t efficiency_583 = Acceptance_lineargaus(offset_583, linear_583, amplitude_583, mean_583, sigma_583, peak_window);

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c5->cd(4);
  x.setRange(583 - peak_window/2, 583 + peak_window/2);
  RooDataHist data583("data583", "583 peak", x, Energy583);

  RooPlot* frame583 = x.frame(Title("RooPlot of x"));
  
  data583.plotOn(frame583);
  
  mean.setVal(583);
  mean.setRange(583-2, 583+2);
  
  gausslin.fitTo(data583);
  gausslin.plotOn(frame583);
  gausslin.plotOn(frame583, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_583 = gaussfrac.getVal();
  
  frame583->Draw();
  
  c5->cd(5);

  t1->Draw("Ener1 >> Energy338", multiplicity);

  x.setRange(338 - peak_window_338/2, 338 + peak_window_338/2);
  RooDataHist data338("data338", "338 peak", x, Energy338);

  RooPlot* frame338 = x.frame(Title("RooPlot of x"));
  
  data338.plotOn(frame338);
  
  mean.setVal(338);
  mean.setRange(338-2, 338+2);

  mean2.setVal(329);
  mean2.setRange(328,330);
  
  
  doublegausslin.fitTo(data338);
  doublegausslin.plotOn(frame338);
  doublegausslin.plotOn(frame338, Components(p2), LineStyle(kDashed));

  gaus1frac.Print();
  doublegausfrac.Print();
  mean.Print();
  sigma.Print();
  mean2.Print();
  sigma2.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_338 = doublegausfrac.getVal();
  
  frame338->Draw();
  

  c4->cd(5);
  //t1->Draw("Ener1 >> Energy338", multiplicity);
  lineardoublegaus->SetParameter(2, 338);
  lineardoublegaus->SetParameter(3, 5);
  lineardoublegaus->FixParameter(4, (2.95 /  11.27));
  lineardoublegaus->FixParameter(5, (328.0/338.320));
  lineardoublegaus->SetParameter(6, 1);
  Energy338->Fit("linear+doublegaus");
  Double_t offset_338 = lineardoublegaus->GetParameter(0);
  Double_t linear_338 = lineardoublegaus->GetParameter(7);
  Double_t amplitude1_338 = lineardoublegaus->GetParameter(1);
  Double_t mean1_338 = lineardoublegaus->GetParameter(2);
  Double_t sigma1_338 = lineardoublegaus->GetParameter(3);
  Double_t amplitude2_338 = lineardoublegaus->GetParameter(1) * lineardoublegaus->GetParameter(4);
  Double_t mean2_338 = lineardoublegaus->GetParameter(5) * lineardoublegaus->GetParameter(2);
  Double_t sigma2_338 = lineardoublegaus->GetParameter(6);

  //Double_t efficiency_338 = Acceptance_lineardoublegaus(offset_338, linear_338,  amplitude1_338, mean1_338, sigma1_338, amplitude2_338, mean2_338, sigma2_338, peak_window);

  for (int i = 0; i < 8; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  c4->cd(6);
  t1->Draw("Ener1 >> Energy239", multiplicity);
  lineargaus->SetParameter(2,239);
  lineargaus->SetParameter(3,5);
  Energy239->Fit("linear+gaus");
  Double_t offset_239 = lineargaus->GetParameter(0);
  Double_t linear_239 = lineargaus->GetParameter(4);
  Double_t amplitude_239 = lineargaus->GetParameter(1);
  Double_t mean_239 = lineargaus->GetParameter(2);
  Double_t sigma_239 = lineargaus->GetParameter(3);
  
  //Double_t efficiency_239 = Acceptance_lineargaus(offset_239, linear_239, amplitude_239, mean_239, sigma_239, peak_window);

  c5->cd(6);
  x.setRange(239 - peak_window/2, 239 + peak_window/2);
  RooDataHist data239("data239", "239 peak", x, Energy239);

  RooPlot* frame239 = x.frame(Title("RooPlot of x"));
  
  data239.plotOn(frame239);
  
  mean.setVal(239);
  mean.setRange(239-2, 239+2);
  
  gausslin.fitTo(data239);
  gausslin.plotOn(frame239);
  gausslin.plotOn(frame239, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_239 = gaussfrac.getVal();
  cout << efficiency_239 << endl;
  cout << efficiency_338 << endl;
  cout << efficiency_583 << endl;
  cout << efficiency_911 << endl;
  cout << efficiency_969 << endl;
  cout << efficiency_2615 << endl;
  
  
  frame239->Draw();
  
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }
  
  THStack *hs = new THStack("hs", "All Peaks");
  
  TCut cut3 = "Ener1 > 2605";
  TCut cut4 = "Ener1 < 2625";
  TCut cut2615 = cut3 && cut4 && multiplicity;

  TCut cut5 = "Ener1 > 959";
  TCut cut6 = "Ener1 < 979";
  TCut cut969 = cut5 && cut6 && multiplicity;

  TCut cut7 = "Ener1 > 901";
  TCut cut8 = "Ener1 < 921";
  TCut cut911 = cut7 && cut8 && multiplicity;

  TCut cut9 = "Ener1 > 573";
  TCut cut10 = "Ener1 < 593";
  TCut cut583 = cut9 && cut10 && multiplicity;

  TCut cut11 = "Ener1 > 323";
  TCut cut12 = "Ener1 < 353";
  TCut cut338 = cut11 && cut12 && multiplicity;

  TCut cut13 = "Ener1 > 229";
  TCut cut14 = "Ener1 < 249";
  TCut cut239 = cut13 && cut14 && multiplicity;

  TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
  c2->cd();

  t1->Draw("Channel >> Peak2615", cut2615);
  t1->Draw("Channel >> Peak969", cut969);
  t1->Draw("Channel >> Peak911", cut911);
  t1->Draw("Channel >> Peak583", cut583);
  t1->Draw("Channel >> Peak338", cut338);
  t1->Draw("Channel >> Peak239", cut239);

  TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
  c3->cd();
  Peak239->Draw();
  c2->cd();

  Peak2615->Scale(efficiency_2615);
  Peak969->Scale(efficiency_969);
  Peak911->Scale(efficiency_911);
  Peak583->Scale(efficiency_583);
  Peak338->Scale(efficiency_338);
  Peak239->Scale(efficiency_239);
  
  //FullString->Scale(0.06923); // reduce to counts per hour
  //FullString->Scale(0.2778); // reduce to mHz
  /*
  FullString->SetTitle("Channel | Ener1>501 && Ener1<521");
  FullString->GetXaxis()->SetTitle("Channel");
  FullString->GetYaxis()->SetTitle("Counts per hour per Channel");
  FullString->SetFillColor(kBlue);
  */


  Peak338->SetLineColor(kAzure);
  Peak2615->SetLineColor(kSpring);
  Peak969->SetLineColor(kMagenta);
  Peak911->SetLineColor(kOrange);
  Peak583->SetLineColor(kRed);
  Peak239->SetLineColor(kCyan);

  hs->Add(Peak2615);
  hs->Add(Peak969);
  hs->Add(Peak911);
  hs->Add(Peak583);
  hs->Add(Peak338);
  hs->Add(Peak239);

  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("Channel");
  hs->GetYaxis()->SetTitle("Peak Events/Channel");
  hs->GetXaxis()->SetNdivisions(13,19, kFALSE);


  //create tree of events

  TTree *tree = new TTree("tree","An example of a ROOT tree");

  Int_t Channel;

  Double_t Rate_2615;
  Double_t Rate_969;
  Double_t Rate_911;
  Double_t Rate_583;
  Double_t Rate_338;
  Double_t Rate_239;

  Double_t Rate_Max; //rate for 4 peaks calibrated
  Double_t Rate_Min; //rate for 1 peak calibrated
  Double_t Rate_Two; //rate for 2 peaks calibrated
  Double_t Rate_Three; //rate for 3 peaks calibrated
  Double_t Rate_Four; //rate for 4 peaks calibrated
  Double_t Rate_Five; //rate for 5 peaks calibrated

  Double_t Time_2615;
  Double_t Time_969;
  Double_t Time_911;
  Double_t Time_583;
  Double_t Time_338;
  Double_t Time_239;

  Double_t Time_Min; //time for 1 peak calibrated
  Double_t Time_Max;
  Double_t Time_Two;
  Double_t Time_Three;
  Double_t Time_Four;
  Double_t Time_Five;

  Double_t Events_2615;
  Double_t Events_969;
  Double_t Events_911;
  Double_t Events_583;
  Double_t Events_338;
  Double_t Events_239;

  Double_t Events_Max;
  Double_t Events_Min;
  Double_t Events_Two;
  Double_t Events_Three;
  Double_t Events_Four;
  Double_t Events_Five;

  Int_t Calib_Peak583 = 0;
  Int_t Calib_Peak2615 = 0;
  Int_t Calib_Peak969 = 0;
  Int_t Calib_Peak911 = 0;
  Int_t Calib_Peak583and2615 = 0;
  Int_t Calib_Peak583and969 = 0;
  Int_t Calib_Peak583and911 = 0;
  Int_t Calib_Peak2615and969 = 0;
  Int_t Calib_Peak2615and911 = 0;
  Int_t Calib_Peak969and911 = 0;
  Int_t Calib_Peak583and2615and969 = 0;
  Int_t Calib_Peak583and2615and911 = 0;
  Int_t Calib_Peak583and969and911 = 0;
  Int_t Calib_Peak2615and969and911 = 0;

  tree->Branch("Rate_2615", &Rate_2615, "Rate_2615/D");
  tree->Branch("CalibrationTime_2615", &Time_2615, "CalibrationTime_2615/D");
  tree->Branch("Rate_969", &Rate_969, "Rate_969/D");
  tree->Branch("CalibrationTime_969", &Time_969, "CalibrationTime_969/D");
  tree->Branch("Rate_911", &Rate_911, "Rate_911/D");
  tree->Branch("CalibrationTime_911", &Time_911, "CalibrationTime_591/D");
  tree->Branch("Rate_583", &Rate_583, "Rate_583/D");
  tree->Branch("CalibrationTime_583", &Time_583, "CalibrationTime_591/D");
  tree->Branch("Rate_338", &Rate_338, "Rate_338/D");
  tree->Branch("CalibrationTime_338", &Time_338, "CalibrationTime_338/D");
  tree->Branch("Rate_239", &Rate_239, "Rate_239/D");
  tree->Branch("CalibrationTime_239", &Time_239, "CalibrationTime_239/D");

  tree->Branch("Rate_Max", &Rate_Max, "Rate_Max/D");
  tree->Branch("CalibrationTime_Min", &Time_Min, "CalibrationTime_Min/D");
  tree->Branch("Rate_Min", &Rate_Min, "Rate_Min/D");
  tree->Branch("CalibrationTime_Max", &Time_Max, "CalibrationTime_Max/D");
  tree->Branch("Rate_Two", &Rate_Two, "Rate_Two/D");
  tree->Branch("CalibrationTime_Two", &Time_Two, "CalibrationTime_Two/D");
  tree->Branch("Rate_Three", &Rate_Three, "Rate_Three/D");
  tree->Branch("CalibrationTime_Three", &Time_Three, "CalibrationTime_Three/D");
  tree->Branch("Rate_Four", &Rate_Four, "Rate_Four/D");
  tree->Branch("CalibrationTime_Four", &Time_Four, "CalibrationTime_Four/D");
  tree->Branch("Rate_Five", &Rate_Five, "Rate_Five/D");
  tree->Branch("CalibrationTime_Five", &Time_Five, "CalibrationTime_Five/D");

  tree->Branch("Channel", &Channel, "Channel/I");

    for (int n = 0; n < nbins; n++) {

      //Get channel
      Channel = TMath::CeilNint(Peak2615->GetBinCenter(n+1));
      //get # of events for each peak and find max
      Events_2615 = Peak2615->GetBinContent(n+1);
      Events_969 = Peak969->GetBinContent(n+1);
      Events_911 = Peak911->GetBinContent(n+1);
      Events_583 = Peak583->GetBinContent(n+1);
      Events_338 = Peak338->GetBinContent(n+1);
      Events_239 = Peak239->GetBinContent(n+1);

      Double_t Events[6] = {Events_2615, Events_969, Events_911, Events_338, Events_239};

      // sort algorithm
      if (Events_583 >= Events_2615) {
	Events_Max = Events_583;
	Events_Two = Events_2615;
      }
      else {
	Events_Max = Events_2615;
	Events_Two = Events_583;
      }
      if (Events_969 >= Events_Max) {
	Events_Three = Events_Two;
	Events_Two = Events_Max;
	Events_Max = Events_969;
      }
      else if (Events_969 >= Events_Two) {
	Events_Three = Events_Two;
	Events_Two = Events_969;
      }
      else {
	Events_Three = Events_969;
      }
      if (Events_911 >= Events_Max) {
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_Max;
	Events_Max = Events_911;
      }
      else if (Events_911 >= Events_Two) {
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_911;
      }
      else if (Events_911 >= Events_Three) {
	Events_Four = Events_Three;
	Events_Three = Events_911;
      }
      else {
	Events_Four = Events_911;
      }
      if (Events_338 >= Events_Max) {
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_Max;
	Events_Max = Events_338;
      }
      else if (Events_338 >= Events_Two) {
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_338;
      }
      else if (Events_338 >= Events_Three) {
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_338;
      }
      else if (Events_338 >= Events_Four) {
	Events_Five = Events_Four;
	Events_Four = Events_338;
      }
      else {
	Events_Five = Events_338;
      }
      if (Events_239 >= Events_Max) {
	Events_Min = Events_Five;
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_Max;
	Events_Max = Events_239;
      }
      else if (Events_239 >= Events_Two) {
	Events_Min = Events_Five;
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_Two;
	Events_Two = Events_239;
      }
      else if(Events_239 >= Events_Three) {
	Events_Min = Events_Five;
	Events_Five = Events_Four;
	Events_Four = Events_Three;
	Events_Three = Events_239;
      }
      else if(Events_239 >= Events_Four) {
	Events_Min = Events_Five;
	Events_Five = Events_Four;
	Events_Four = Events_239;
      }
      else if(Events_239 >= Events_Five) {
	Events_Min = Events_Five;
	Events_Five = Events_239;
      }
      else {
	Events_Min = Events_239;
      }
      /*
      cout << "*************" << endl;
      cout << "Events_2615: " << Events_2615 << endl;
      cout << "Events_969: " << Events_969 << endl;
      cout << "Events_911: " << Events_911 << endl;
      cout << "Events_583: " << Events_583 << endl;
      cout << "Events_338: " << Events_338 << endl;
      cout << "Events_239: " << Events_239 << endl;

      cout << "Events_Max: " << Events_Max << endl;
      cout << "Events_Two: " << Events_Two << endl;
      cout << "Events_Three: " << Events_Three << endl;
      cout << "Events_Four: " << Events_Four << endl;
      cout << "Events_Five: " << Events_Five << endl;
      cout << "Events_Min: " << Events_Min << endl;
      */
      Rate_2615 = Events_2615 * time_scaling * 0.2778;
      Rate_969 = Events_969 * time_scaling * 0.2778;
      Rate_911 = Events_911 * time_scaling * 0.2778;
      Rate_583 = Events_583 * time_scaling * 0.2778;
      Rate_338 = Events_338 * time_scaling * 0.2778;
      Rate_239 = Events_239 * time_scaling * 0.2778;

      Rate_Max = Events_Max * time_scaling * 0.2778;
      Rate_Min = Events_Min * time_scaling * 0.2778;
      Rate_Two = Events_Two * time_scaling * 0.2778;
      Rate_Three = Events_Three * time_scaling * 0.2778;
      Rate_Four = Events_Four * time_scaling * 0.2778;
      Rate_Five = Events_Five * time_scaling * 0.2778;

      Time_2615 = eventsToCalibrate / (86.0 * Rate_2615);
      Time_969 = eventsToCalibrate / (86.0 * Rate_969);
      Time_911 = eventsToCalibrate / (86.0 * Rate_911);
      Time_583 = eventsToCalibrate / (86.0 * Rate_583);
      Time_338 = eventsToCalibrate / (86.0 * Rate_338);
      Time_239 = eventsToCalibrate / (86.0 * Rate_239);
 
      Time_Min = eventsToCalibrate / (86.0 * Rate_Max);
      Time_Max = eventsToCalibrate / (86.0 * Rate_Min);
      Time_Two = eventsToCalibrate / (86.0 * Rate_Two);
      Time_Three = eventsToCalibrate / (86.0 * Rate_Three);
      Time_Four = eventsToCalibrate / (86 * Rate_Four);
      Time_Five = eventsToCalibrate / (86 * Rate_Five);


      if (Rate_583 <= 0.1 || Rate_2615 <= 0.1 || Rate_969 <= 0.1 || Rate_911 <= 0.1 || Rate_338 <= 0.1 || Rate_239 <= 0.1) {
	cout << "Rate Error " << endl;
	cout << "Rate 583: " << Rate_583 << " Rate 2615: " << Rate_2615 << " Rate 969: " << Rate_969 << " Rate 911: " << Rate_911 << " Rate 338: " << Rate_338 << " Rate 239: " << Rate_239 <<endl;
      }

      tree->Fill();
    }
    /*
      cout << "# of channels calibrated by peak" << endl;
      cout << "583: " << Calib_Peak583 << endl;
      cout << "2615: " << Calib_Peak2615 << endl;
      cout << "969: " << Calib_Peak969 << endl;
      cout << "911: " << Calib_Peak911 << endl;
      cout << "583 and 2615: " << Calib_Peak583and2615 << endl;
      cout << "583 and 969: " << Calib_Peak583and969 << endl;
      cout << "583 and 911: " << Calib_Peak583and911 << endl;
      cout << "2615 and 969: " << Calib_Peak2615and969 << endl;
      cout << "2615 and 911: " << Calib_Peak2615and911 << endl;
      cout << "969 and 911: " << Calib_Peak969and911 << endl;
      cout << "583 and 2615 and 969: " << Calib_Peak583and2615and969 << endl;
      cout << "583 and 2615 and 911: " << Calib_Peak583and2615and911 << endl;
      cout << "583 and 969 and 911: " << Calib_Peak583and969and911 << endl;
      cout << "2615 and 969 and 911: " << Calib_Peak2615and969and911 << endl;
    */


    TFile* output = new TFile("AllPeaks_rate.root","Recreate");
    tree->Write();
    output->Close();

}
