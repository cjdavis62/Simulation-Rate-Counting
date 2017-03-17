/* This script takes in a g4cuore processed ROOT File and performs fits to the peaks
These peaks are then analyzed for the ratio of signal/background events
With the number of signal events, it is then calculated how long it would take to get 50 events per channel per peak

This script in particular is for Co56 spectra in a file called AllString_g4cuore.root

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

void peak__peak_number_() {

  gStyle->SetOptFit(1);
  
  int nbins = 988;
  int energy_bins = 200;
  double time_scaling = 0.01115; // scale from events to events per hour

  double eventsToCalibrate = 50; // How many events to require for a channel to be calibrated

  Double_t peak_window_left = _left_window_; 
  Double_t peak_window_right = _right_window_;

  // looking at the M1 Spectrum
  TCut multiplicity = "Multiplicity == 1";

  // open the file and the tree
  //TFile* f1 = new TFile("/data-mgm/cuore/simulation/Byron/Cobalt56/StringsCobalt_g4cuore.root");
  TFile* f1 = new TFile("_filename_");
  TTree* t1 = (TTree*)f1->Get("outTree");

  // Per channel histograms for each peak
  TH1F* Peak_peaknumber_ = new TH1F("Peak_peak_number_", "Peak_peak_number_", nbins, 0, 988);
   
  // Energy histograms around each peak
  TH1F* Energy_peak_number_ = new TH1F("Energy_peak_number_", "Energy_peak_number_", energy_bins, (_peak_number_ -peak_window_left), (_peak_number_ + peak_window_right));
   
  // Fit Types
  TF1 * lineargaus = new TF1("linear+gaus", "[0] + [4] * x + [1] * exp(-0.5*((x-[2])/[3])**2)", 0, 5);
  TF2 * lineardoublegaus = new TF1("linear+doublegaus", "[0] + [7] * x + [1] * exp(-0.5*((x-[2])/[3])**2) + [4]*[1] * exp(-0.5*((x-[5]*[2])/[6])**2)", 0, 5);

  // begin RooFits

  // Single Gaussian
  RooRealVar x("x", "x", 0, 4000);
  RooRealVar mean("mean", "mean of gaussian", 10, 0, 4000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0, 5);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  // Second Gaussian
  RooRealVar mean2("mean2", "mean of secondary gaussian", 0, 0, 4000);
  RooRealVar sigma2("sigma2", "sigma of secondary gaussian", 1, 0, 10);
  
  RooGaussian gauss2("gauss2", "secondary gaussian PDF", x, mean2, sigma2);
  
  // Linear background
  RooRealVar a0("a0", "a0", 1, 0, 50000);
  RooRealVar a1("a1", "a1", 1, -1000, 1000);
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

  // Make another histogram with 'standard' root fitting. This is to help cross-check fits
  TCanvas * c2 = new TCanvas("c2", "Fits", 1200, 1000);
  c2->Divide(1,2);
  // Create a canvas for Roofit plots and fits
  c2->cd(1);
  
  t1->Draw("Ener1 >> Energy_peak_number_", multiplicity);
  lineargaus->SetParameter(2,_peak_number_);
  lineargaus->SetParameter(3,5);
  Energy_peak_number_->Fit("linear+gaus");
  Double_t offset__peak_number_ = lineargaus->GetParameter(0);
  Double_t linear__peak_number_ = lineargaus->GetParameter(4);
  Double_t amplitude__peak_number_ = lineargaus->GetParameter(1);
  Double_t mean__peak_number_ = lineargaus->GetParameter(2);
  Double_t sigma__peak_number_ = lineargaus->GetParameter(3);

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c2->cd(2);
  x.setRange(_peak_number_ - peak_window_left, _peak_number_ + peak_window_right);
  RooDataHist data_peak_number_("data_peak_number_", "_peak_number_ peak", x, Energy_peak_number_);

  RooPlot* frame_peak_number_ = x.frame(Title("_peak_number_"));
  
  data_peak_number_.plotOn(frame_peak_number_);

  mean.setRange(mean__peak_number_-1, mean__peak_number_+1);
  mean.setVal(mean__peak_number_);
  a0.setVal(offset__peak_number_);
  a1.setVal(linear__peak_number_);
  
  gausslin.fitTo(data_peak_number_);
  gausslin.plotOn(frame_peak_number_);
  gausslin.plotOn(frame_peak_number_, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  frame_peak_number_->Draw();

  Double_t efficiency__peak_number = gaussfrac.getVal();

  c2->SaveAs("Output/Plot__peak_number_.pdf");

  ofstream OutFile;
  OutFile.open("Output/Ratio__peak_number_.dat");
  OutFile << "_peak_number_\t " << efficiency__peak_number << "\n";
}
