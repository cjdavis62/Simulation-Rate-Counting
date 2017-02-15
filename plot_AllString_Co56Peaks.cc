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

void plot_AllString_Co56Peaks() {

  int nbins = 988;
  int energy_bins = 200;
  double time_scaling = 0.01115; // scale from events to events per hour

  double eventsToCalibrate = 50; // How many events to require for a channel to be calibrated

  Double_t peak_window = 20; // 20 keV window
  Double_t peak_window_2035_left = 35; // Give more space to the left to include the subpeak
  Double_t peak_window_2035_right = peak_window/2; // Give the normal space to the right

  // looking at the M1 Spectrum
  TCut multiplicity = "Multiplicity == 1";

  // open the file and the tree
  TFile* f1 = new TFile("AllString_g4cuore.root");
  TTree* t1 = (TTree*)f1->Get("outTree");


  // Per channel histograms for each peak
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
  
  // Energy histograms around each peak
  TH1F* Energy2599 = new TH1F("Energy2599", "Energy2599", energy_bins, (2599 - peak_window/2), (2599 + peak_window/2));
  TH1F* Energy847 = new TH1F("Energy847", "Energy847", energy_bins, (847 - peak_window/2), (847 + peak_window/2));
  TH1F* Energy1238 = new TH1F("Energy1238", "Energy1238", energy_bins, (1238 - peak_window/2), (1238 + peak_window/2));
  TH1F* Energy511 = new TH1F("Energy511", "Energy511", energy_bins, (511 - peak_window/2), (511 + peak_window/2));
  TH1F* Energy1771 = new TH1F("Energy1771", "Energy1771", energy_bins, (1771 - peak_window/2), (1771 + peak_window/2));
  TH1F* Energy1037 = new TH1F("Energy1037", "Energy1037", energy_bins, (1037 - peak_window/2), (1037 + peak_window/2));
  TH1F* Energy3254 = new TH1F("Energy3254", "Energy3254", energy_bins, (3254 - peak_window/2), (3254 + peak_window/2));
  TH1F* Energy2035 = new TH1F("Energy2035", "Energy2035", energy_bins, (2035 - peak_window_2035_left), (2035 + peak_window_2035_right));
  TH1F* Energy1360 = new TH1F("Energy1360", "Energy1360", energy_bins, (1360 - peak_window/2), (1360 + peak_window/2));
  TH1F* Energy3202 = new TH1F("Energy3202", "Energy3202", energy_bins, (3202 - peak_window/2), (3202 + peak_window/2));
  TH1F* Energy3451 = new TH1F("Energy3451", "Energy3451", energy_bins, (3451 - peak_window/2), (3451 + peak_window/2));
    
  
  // Fit Types
  TF1 * lineargaus = new TF1("linear+gaus", "[0] + [4] * x + [1] * exp(-0.5*((x-[2])/[3])**2)", 0, 5);
  TF2 * lineardoublegaus = new TF1("linear+doublegaus", "[0] + [7] * x + [1] * exp(-0.5*((x-[2])/[3])**2) + [4]*[1] * exp(-0.5*((x-[5]*[2])/[6])**2)", 0, 5);

  // begin RooFits

  // Single Gaussian
  RooRealVar x("x", "x", 0, 3000);
  RooRealVar mean("mean", "mean of gaussian", 10, 0, 3000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0, 10);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  // Second Gaussian
  RooRealVar mean2("mean2", "mean of secondary gaussian", 0, 0, 3000);
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

  // Create a canvas for Roofit plots and fits
  TCanvas * c1 = new TCanvas("c1", "Roofit", 1200, 1000);
  c1->cd();
  // Make the canvas into 12 plot spaces (will use 11)
  c1->Divide(4,3);
  c1->cd(1);
  // Fill the histograms with the data
  t1->Draw("Ener1 >> Energy2599", multiplicity, "goff");

  // Set the x range
  x.setRange(2599 - peak_window/2, 2599 + peak_window/2);
  RooDataHist data2599("data2599", "2599 peak", x, Energy2599);

  RooPlot* frame2599 = x.frame(Title("RooPlot of x"));
  
  data2599.plotOn(frame2599);
  
  // Help the fit find the mean because these things are dumb
  mean.setVal(2599);
  mean.setRange(2599-1, 2599+1);
  
  // Do the fit
  gausslin.fitTo(data2599);
  gausslin.plotOn(frame2599);
  gausslin.plotOn(frame2599, Components(p2), LineStyle(kDashed));


  // Print the params
  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();
  
  // Save the efficiency
  Double_t efficiency_2599 = gaussfrac.getVal();

  // Draw the histogram
  frame2599->Draw();

  // Make another histogram with 'standard' root fitting. This is to help cross-check fits
  TCanvas * c2 = new TCanvas("c2", "ROOT fit", 1200, 1000);
  c2->Divide(4,3);
  c2->cd(1);
  // Draw the histogram
  Energy2599->Draw();
  // help the fits because fits algorithms are dumb
  lineargaus->SetParameter(2, 2599);
  lineargaus->SetParameter(3, 5);
  Energy2599->Fit("linear+gaus");

  // Save the parameters
  Double_t offset_2599 = lineargaus->GetParameter(0);
  Double_t linear_2599 = lineargaus->GetParameter(4);
  Double_t amplitude_2599 = lineargaus->GetParameter(1);
  Double_t mean_2599 = lineargaus->GetParameter(2);
  Double_t sigma_2599 = lineargaus->GetParameter(3);

  // Need to release parameters (so next fit can start with a tabula rasa)
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  // Move to the next part of the canvas
  c2->cd(2);
  t1->Draw("Ener1 >> Energy847", multiplicity, "goff");
  Energy847->Draw();
  lineargaus->SetParameter(2, 847);
  lineargaus->SetParameter(3, 5);
  Energy847->Fit("linear+gaus");
  Double_t offset_847 = lineargaus->GetParameter(0);
  Double_t linear_847 = lineargaus->GetParameter(4);
  Double_t amplitude_847 = lineargaus->GetParameter(1);
  Double_t mean_847 = lineargaus->GetParameter(2);
  Double_t sigma_847 = lineargaus->GetParameter(3);

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c1->cd(2);

  x.setRange(847 - peak_window/2, 847 + peak_window/2);
  RooDataHist data847("data847", "847 peak", x, Energy847);

  RooPlot* frame847 = x.frame(Title("RooPlot of x"));
  
  data847.plotOn(frame847);
  
  mean.setVal(847);
  mean.setRange(847-3, 847+3);

  gausslin.fitTo(data847);
  gausslin.plotOn(frame847);
  gausslin.plotOn(frame847, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();
  
  Double_t efficiency_847 = gaussfrac.getVal();

  frame847->Draw();

  c2->cd(3);
  t1->Draw("Ener1 >> Energy1238", multiplicity);
  lineargaus->SetParameter(2, 1238);
  lineargaus->SetParameter(3, 5);
  Energy1238->Fit("linear+gaus");
  Double_t offset_1238 = lineargaus->GetParameter(0);
  Double_t linear_1238 = lineargaus->GetParameter(4);
  Double_t amplitude_1238 = lineargaus->GetParameter(1);
  Double_t mean_1238 = lineargaus->GetParameter(2);
  Double_t sigma_1238 = lineargaus->GetParameter(3);

  for (int i = 0; i < 5; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  c1->cd(3);

  x.setRange(1238 - peak_window/2, 1238 + peak_window/2);
  RooDataHist data1238("data1238", "1238 peak", x, Energy1238);

  RooPlot* frame1238 = x.frame(Title("RooPlot of x"));
  
  data1238.plotOn(frame1238);
  
  mean.setVal(1238);
  mean.setRange(1238-2, 1238+2);

  gausslin.fitTo(data1238);
  gausslin.plotOn(frame1238);
  gausslin.plotOn(frame1238, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();
  
  Double_t efficiency_1238 = gaussfrac.getVal();

  frame1238->Draw();
  
  c2->cd(4);

  t1->Draw("Ener1 >> Energy511", multiplicity, "goff");


  lineargaus->SetParameter(2,511);
  lineargaus->SetParameter(3,5);
  Energy511->Fit("linear+gaus");
  Double_t offset_511 = lineargaus->GetParameter(0);
  Double_t linear_511 = lineargaus->GetParameter(4);
  Double_t amplitude_511 = lineargaus->GetParameter(1);
  Double_t mean_511 = lineargaus->GetParameter(2);
  Double_t sigma_511 = lineargaus->GetParameter(3);

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c1->cd(4);
  x.setRange(511 - peak_window/2, 511 + peak_window/2);
  RooDataHist data511("data511", "511 peak", x, Energy511);

  RooPlot* frame511 = x.frame(Title("RooPlot of x"));
  
  data511.plotOn(frame511);
  
  mean.setVal(511);
  mean.setRange(511-2, 511+2);
  
  gausslin.fitTo(data511);
  gausslin.plotOn(frame511);
  gausslin.plotOn(frame511, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_511 = gaussfrac.getVal();
  
  frame511->Draw();
  
  c1->cd(5);

  t1->Draw("Ener1 >> Energy1771", multiplicity);

  x.setRange(1771 - peak_window/2, 1771 + peak_window/2);
  RooDataHist data1771("data1771", "1771 peak", x, Energy1771);

  RooPlot* frame1771 = x.frame(Title("RooPlot of x"));
  
  data1771.plotOn(frame1771);
  
  mean.setVal(1771);
  mean.setRange(1771-2, 1771+2);

  gausslin.fitTo(data1771);
  gausslin.plotOn(frame1771);
  gausslin.plotOn(frame1771, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_1771 = doublegausfrac.getVal();
  
  frame1771->Draw();

  c2->cd(5);
  lineargaus->SetParameter(2, 1771);
  lineargaus->SetParameter(3, 5);
  Energy1771->Fit("linear+gaus");
  Double_t offset_1771 = lineargaus->GetParameter(0);
  Double_t linear_1771 = lineargaus->GetParameter(4);
  Double_t amplitude_1771 = lineargaus->GetParameter(1);
  Double_t mean_1771 = lineargaus->GetParameter(2);
  Double_t sigma_1771 = lineargaus->GetParameter(3);
  
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c2->cd(6);
  t1->Draw("Ener1 >> Energy1037", multiplicity);
  lineargaus->SetParameter(2,1037);
  lineargaus->SetParameter(3,5);
  Energy1037->Fit("linear+gaus");
  Double_t offset_1037 = lineargaus->GetParameter(0);
  Double_t linear_1037 = lineargaus->GetParameter(4);
  Double_t amplitude_1037 = lineargaus->GetParameter(1);
  Double_t mean_1037 = lineargaus->GetParameter(2);
  Double_t sigma_1037 = lineargaus->GetParameter(3);
  
  c1->cd(6);
  x.setRange(1037 - peak_window/2, 1037 + peak_window/2);
  RooDataHist data1037("data1037", "1037 peak", x, Energy1037);

  RooPlot* frame1037 = x.frame(Title("RooPlot of x"));
  
  data1037.plotOn(frame1037);
  
  mean.setVal(1037);
  mean.setRange(1037-2, 1037+2);
  
  gausslin.fitTo(data1037);
  gausslin.plotOn(frame1037);
  gausslin.plotOn(frame1037, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_1037 = gaussfrac.getVal();
    
  frame1037->Draw();
  
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }

  c2->cd(7);
  t1->Draw("Ener1 >> Energy3254", multiplicity);
  //  lineardoublegaus->SetParameter(0, 9876.4);
  //  lineardoublegaus->SetParameter(1, 6471.7);
  lineardoublegaus->SetParameter(2, 3254);
  lineardoublegaus->SetParameter(3, 2.1);
  lineardoublegaus->FixParameter(4, (1.876 / 7.923));
  lineardoublegaus->FixParameter(5, 3273.079 / 3253.503);
  lineardoublegaus->SetParameter(6, 1);
  //  lineardoublegaus->SetParameter(7, -9.8625);
  Energy3254->Fit("linear+doublegaus");
  Double_t offset_3254 = lineardoublegaus->GetParameter(0);
  Double_t linear_3254 = lineardoublegaus->GetParameter(7);
  Double_t amplitude1_3254 = lineardoublegaus->GetParameter(1);
  Double_t mean1_3254 = lineardoublegaus->GetParameter(2);
  Double_t sigma1_3254 = lineardoublegaus->GetParameter(3);
  Double_t amplitude2_3254 = lineardoublegaus->GetParameter(1) * lineardoublegaus->GetParameter(4);
  Double_t mean2_3254 = lineardoublegaus->GetParameter(5) * lineardoublegaus->GetParameter(2);
  Double_t sigma2_3254 = lineardoublegaus->GetParameter(6);


  for (int i = 0; i < 8; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  c1->cd(7);

  x.setRange(3254 - peak_window/2, 3254 + peak_window/2);
  RooDataHist data3254("data3254", "3254 peak", x, Energy3254);

  RooPlot* frame3254 = x.frame(Title("RooPlot of x"));
  
  data3254.plotOn(frame3254);
  
  mean.setVal(3254);
  mean.setRange(3254-2, 3254+2);

  mean2.setVal(3273);
  mean2.setRange(3273-2,3273+2);
  
  doublegausslin.fitTo(data3254);
  doublegausslin.plotOn(frame3254);
  doublegausslin.plotOn(frame3254, Components(p2), LineStyle(kDashed));

  gaus1frac.Print();
  doublegausfrac.Print();
  mean.Print();
  sigma.Print();
  mean2.Print();
  sigma2.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_3254 = doublegausfrac.getVal();

  frame3254->Draw();

  c2->cd(8);
  t1->Draw("Ener1 >> Energy2035", multiplicity);
  //lineardoublegaus->SetParameter(0, 9876.4);
  //lineardoublegaus->SetParameter(1, 6471.7);
  lineardoublegaus->SetParameter(2, 2035);
  lineardoublegaus->SetParameter(3, 2.1);
  lineardoublegaus->FixParameter(4, (3.016 / 7.77));
  lineardoublegaus->FixParameter(5, 2015.215 / 2034.791);
  lineardoublegaus->SetParameter(6, 1);
  //lineardoublegaus->SetParameter(7, -9.8625);
  Energy2035->Fit("linear+doublegaus");
  Double_t offset_2035 = lineardoublegaus->GetParameter(0);
  Double_t linear_2035 = lineardoublegaus->GetParameter(7);
  Double_t amplitude1_2035 = lineardoublegaus->GetParameter(1);
  Double_t mean1_2035 = lineardoublegaus->GetParameter(2);
  Double_t sigma1_2035 = lineardoublegaus->GetParameter(3);
  Double_t amplitude2_2035 = lineardoublegaus->GetParameter(1) * lineardoublegaus->GetParameter(4);
  Double_t mean2_2035 = lineardoublegaus->GetParameter(5) * lineardoublegaus->GetParameter(2);
  Double_t sigma2_2035 = lineardoublegaus->GetParameter(6);


  for (int i = 0; i < 8; i++) {
    lineardoublegaus->ReleaseParameter(i);
  }

  c1->cd(8);

  x.setRange(2035 - peak_window_2035_left, 2035 + peak_window_2035_right);
  RooDataHist data2035("data2035", "2035 peak", x, Energy2035);

  RooPlot* frame2035 = x.frame(Title("RooPlot of x"));
  
  data2035.plotOn(frame2035);
  
  mean.setVal(2035);
  mean.setRange(967, 972);

  mean2.setVal(2015);
  mean2.setRange(2015-2,2015+2);
  
  doublegausslin.fitTo(data2035);
  doublegausslin.plotOn(frame2035);
  doublegausslin.plotOn(frame2035, Components(p2), LineStyle(kDashed));

  gaus1frac.Print();
  doublegausfrac.Print();
  mean.Print();
  sigma.Print();
  mean2.Print();
  sigma2.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_2035 = doublegausfrac.getVal();

  frame2035->Draw();

  c2->cd(9);
  t1->Draw("Ener1 >> Energy1360", multiplicity);
  lineargaus->SetParameter(2,1360);
  lineargaus->SetParameter(3,5);
  Energy1360->Fit("linear+gaus");
  Double_t offset_1360 = lineargaus->GetParameter(0);
  Double_t linear_1360 = lineargaus->GetParameter(4);
  Double_t amplitude_1360 = lineargaus->GetParameter(1);
  Double_t mean_1360 = lineargaus->GetParameter(2);
  Double_t sigma_1360 = lineargaus->GetParameter(3);
  
  c1->cd(9);
  x.setRange(1360 - peak_window/2, 1360 + peak_window/2);
  RooDataHist data1360("data1360", "1360 peak", x, Energy1360);

  RooPlot* frame1360 = x.frame(Title("RooPlot of x"));
  
  data1360.plotOn(frame1360);
  
  mean.setVal(1360);
  mean.setRange(1360-2, 1360+2);
  
  gausslin.fitTo(data1360);
  gausslin.plotOn(frame1360);
  gausslin.plotOn(frame1360, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_1360 = gaussfrac.getVal();
  
  frame1360->Draw();

  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }
  
  c2->cd(10);
  t1->Draw("Ener1 >> Energy3202", multiplicity);
  lineargaus->SetParameter(2,3202);
  lineargaus->SetParameter(3,5);
  Energy3202->Fit("linear+gaus");
  Double_t offset_3202 = lineargaus->GetParameter(0);
  Double_t linear_3202 = lineargaus->GetParameter(4);
  Double_t amplitude_3202 = lineargaus->GetParameter(1);
  Double_t mean_3202 = lineargaus->GetParameter(2);
  Double_t sigma_3202 = lineargaus->GetParameter(3);
  
  c1->cd(10);
  x.setRange(3202 - peak_window/2, 3202 + peak_window/2);
  RooDataHist data3202("data3202", "3202 peak", x, Energy3202);

  RooPlot* frame3202 = x.frame(Title("RooPlot of x"));
  
  data3202.plotOn(frame3202);
  
  mean.setVal(3202);
  mean.setRange(3202-2, 3202+2);
  
  gausslin.fitTo(data3202);
  gausslin.plotOn(frame3202);
  gausslin.plotOn(frame3202, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_3202 = gaussfrac.getVal();
  
  frame3202->Draw();
  
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }
  
  // begin new stuffs
  c2->cd(11);
  t1->Draw("Ener1 >> Energy3451", multiplicity);
  lineargaus->SetParameter(2,3451);
  lineargaus->SetParameter(3,5);
  Energy3451->Fit("linear+gaus");
  Double_t offset_3451 = lineargaus->GetParameter(0);
  Double_t linear_3451 = lineargaus->GetParameter(4);
  Double_t amplitude_3451 = lineargaus->GetParameter(1);
  Double_t mean_3451 = lineargaus->GetParameter(2);
  Double_t sigma_3451 = lineargaus->GetParameter(3);
  
  c1->cd(11);
  x.setRange(3451 - peak_window/2, 3451 + peak_window/2);
  RooDataHist data3451("data3451", "3451 peak", x, Energy3451);

  RooPlot* frame3451 = x.frame(Title("RooPlot of x"));
  
  data3451.plotOn(frame3451);
  
  mean.setVal(3451);
  mean.setRange(3451-2, 3451+2);
  
  gausslin.fitTo(data3451);
  gausslin.plotOn(frame3451);
  gausslin.plotOn(frame3451, Components(p2), LineStyle(kDashed));

  gaussfrac.Print();
  mean.Print();
  sigma.Print();
  a0.Print();
  a1.Print();

  Double_t efficiency_3451 = gaussfrac.getVal();
    
  frame3451->Draw();
  
  for (int i = 0; i < 5; i++) {
    lineargaus->ReleaseParameter(i);
  }
  

  // Print all the efficiencies
  cout << "Efficiency 2599: " << efficiency_2599 << endl;
  cout << "Efficiency 847: "  << efficiency_847  << endl;
  cout << "Efficiency 1238: " << efficiency_1238 << endl;
  cout << "Efficiency 511: "  << efficiency_511  << endl;
  cout << "Efficiency 1771: " << efficiency_1771 << endl;
  cout << "Efficiency 1037: " << efficiency_1037 << endl;
  cout << "Efficiency 3254: " << efficiency_3254 << endl;
  cout << "Efficiency 2035: " << efficiency_2035 << endl;
  cout << "Efficiency 1360: " << efficiency_1360 << endl;
  cout << "Efficiency 3202: " << efficiency_3202 << endl;
  cout << "Efficiency 3451: " << efficiency_3451 << endl;


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

  t1->Draw("Channel >> Peak2599", cut2599, "goff");
  t1->Draw("Channel >> Peak847", cut847, "goff");
  t1->Draw("Channel >> Peak1238", cut1238, "goff");
  t1->Draw("Channel >> Peak511", cut511, "goff");
  t1->Draw("Channel >> Peak1771", cut1771, "goff");
  t1->Draw("Channel >> Peak1037", cut1037, "goff");
  t1->Draw("Channel >> Peak3254", cut3254, "goff");
  t1->Draw("Channel >> Peak2035", cut2035, "goff");
  t1->Draw("Channel >> Peak1360", cut1360, "goff");
  t1->Draw("Channel >> Peak3202", cut3202, "goff");
  t1->Draw("Channel >> Peak3451", cut3451, "goff");

   // Reduce each peak by their efficiency
  Peak2599->Scale(efficiency_2599);
  Peak847->Scale(efficiency_847);
  Peak1238->Scale(efficiency_1238);
  Peak511->Scale(efficiency_511);
  Peak1771->Scale(efficiency_1771);
  Peak1037->Scale(efficiency_1037);
  Peak3254->Scale(efficiency_3254);
  Peak2035->Scale(efficiency_2035);
  Peak1360->Scale(efficiency_1360);
  Peak3202->Scale(efficiency_3202);
  Peak3451->Scale(efficiency_3451);
  
  //FullString->Scale(0.06923); // reduce to counts per hour
  //FullString->Scale(0.2778); // reduce to mHz
  /*
  FullString->SetTitle("Channel | Ener1>501 && Ener1<521");
  FullString->GetXaxis()->SetTitle("Channel");
  FullString->GetYaxis()->SetTitle("Counts per hour per Channel");
  FullString->SetFillColor(kBlue);
  */

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
  tree->Branch("CalibrationTime_1238", &Time_1238, "CalibrationTime_591/D");
  tree->Branch("Rate_511", &Rate_511, "Rate_511/D");
  tree->Branch("CalibrationTime_511", &Time_511, "CalibrationTime_591/D");
  tree->Branch("Rate_1771", &Rate_1771, "Rate_1771/D");
  tree->Branch("CalibrationTime_1771", &Time_1771, "CalibrationTime_1771/D");
  tree->Branch("Rate_1037", &Rate_1037, "Rate_1037/D");
  tree->Branch("CalibrationTime_1037", &Time_1037, "CalibrationTime_1037/D");
   tree->Branch("Rate_3254", &Rate_3254, "Rate_3254/D");
  tree->Branch("CalibrationTime_3254", &Time_3254, "CalibrationTime_3254/D");
  tree->Branch("Rate_2035", &Rate_2035, "Rate_2035/D");
  tree->Branch("CalibrationTime_2035", &Time_2035, "CalibrationTime_591/D");
  tree->Branch("Rate_1360", &Rate_1360, "Rate_1360/D");
  tree->Branch("CalibrationTime_1360", &Time_1360, "CalibrationTime_591/D");
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

      //Get channel
      Channel = TMath::CeilNint(Peak2599->GetBinCenter(n+1));

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

      /*
      cout << "*************" << endl;
      cout << "Events_2599: " << Events_2599 << endl;
      cout << "Events_847: " << Events_847 << endl;
      cout << "Events_1238: " << Events_1238 << endl;
      cout << "Events_511: " << Events_511 << endl;
      cout << "Events_1771: " << Events_1771 << endl;
      cout << "Events_1037: " << Events_1037 << endl;

      cout << "Events_Max: " << Events_Max << endl;
      cout << "Events_Two: " << Events_Two << endl;
      cout << "Events_Three: " << Events_Three << endl;
      cout << "Events_Four: " << Events_Four << endl;
      cout << "Events_Five: " << Events_Five << endl;
      cout << "Events_Min: " << Events_Min << endl;
      */
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
      Time_Four = eventsToCalibrate / (86 * Rate_Four);
      Time_Five = eventsToCalibrate / (86 * Rate_Five);
      Time_Six = eventsToCalibrate / (86.0 * Rate_Six);
      Time_Seven = eventsToCalibrate / (86.0 * Rate_Seven);
      Time_Eight = eventsToCalibrate / (86.0 * Rate_Eight);
      Time_Nine = eventsToCalibrate / (86 * Rate_Nine);
      Time_Ten = eventsToCalibrate / (86 * Rate_Ten);
      Time_Min = eventsToCalibrate / (86.0 * Rate_Max);

      /*if (Rate_511 <= 0.1 || Rate_2599 <= 0.1 || Rate_847 <= 0.1 || Rate_1238 <= 0.1 || Rate_1771 <= 0.1 || Rate_1037 <= 0.1 || Rate_3254 <=0.1 || Rate_2035 <=0.1 || Rate_1360 <=0.1 || Rate_3202 <=0.1 || Rate_3451 <=0.1) {
	cout << "Rate Error " << endl;
	cout << "Rate 511: " << Rate_511 << " Rate 2599: " << Rate_2599 << " Rate 847: " << Rate_847 << " Rate 1238: " << Rate_1238 << " Rate 1771: " << Rate_1771 << " Rate 1037: " << Rate_1037 << " Rate 3254: " << Rate_3254 << " Rate 2035: " << Rate_2035 << " Rate 1360: " << Rate_1360 << " Rate 3202: " << Rate_3202 << " Rate 3451: " << Rate_3451 << endl;
	}*/

      tree->Fill();
    }
    /*
      cout << "# of channels calibrated by peak" << endl;
      cout << "511: " << Calib_Peak511 << endl;
      cout << "2599: " << Calib_Peak2599 << endl;
      cout << "847: " << Calib_Peak847 << endl;
      cout << "1238: " << Calib_Peak1238 << endl;
      cout << "511 and 2599: " << Calib_Peak511and2599 << endl;
      cout << "511 and 847: " << Calib_Peak511and847 << endl;
      cout << "511 and 1238: " << Calib_Peak511and1238 << endl;
      cout << "2599 and 847: " << Calib_Peak2599and847 << endl;
      cout << "2599 and 1238: " << Calib_Peak2599and1238 << endl;
      cout << "847 and 1238: " << Calib_Peak847and1238 << endl;
      cout << "511 and 2599 and 847: " << Calib_Peak511and2599and847 << endl;
      cout << "511 and 2599 and 1238: " << Calib_Peak511and2599and1238 << endl;
      cout << "511 and 847 and 1238: " << Calib_Peak511and847and1238 << endl;
      cout << "2599 and 847 and 1238: " << Calib_Peak2599and847and1238 << endl;
    */


    TFile* output = new TFile("AllPeaks_rate.root","Recreate");
    tree->Write();
    output->Close();

}
