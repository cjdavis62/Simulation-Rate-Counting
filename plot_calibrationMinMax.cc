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


using namespace std;

void plot_calibrationMinMax() {

  int nbins = 30;
  double max = 3.0;

  TFile *f1 = new TFile("AllPeaks_rate.root");

  TTree* t1 = (TTree*)f1->Get("tree");
  TH1F* ThreeTime = new TH1F("ThreeTime", "Three Peaks", nbins, 0, max);
  TH1F* TwoTime = new TH1F("TwoTime", "Two Peaks", nbins, 0, max);
  TH1F* MaxTime = new TH1F("MaxTime", "Four Peaks", nbins, 0, max);
  TH1F* MinTime = new TH1F("MinTime", "One Peak", nbins, 0, max);

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  c2->cd();

  t1->Draw("CalibrationTime_Three >> ThreeTime");
  t1->Draw("CalibrationTime_Two >> TwoTime");
  t1->Draw("CalibrationTime_Max >> MaxTime");
  t1->Draw("CalibrationTime_Min >> MinTime");

  THStack *hs = new THStack("hs", "Calibration Times for peaks");



  ThreeTime->SetLineColor(kMagenta);
  ThreeTime->SetLineWidth(2);
  TwoTime->SetLineColor(kSpring);
  TwoTime->SetLineWidth(2);
  MaxTime->SetLineColor(kBlack);
  MaxTime->SetLineWidth(3);
  MinTime->SetLineColor(kCyan);
  MinTime->SetLineWidth(3);


  hs->Add(MinTime);
  hs->Add(TwoTime);
  hs->Add(ThreeTime);
  hs->Add(MaxTime);
  //hs->Add(Peak969);
  //hs->Add(Peak911);
  ThreeTime->GetXaxis()->SetTitle("Calibration Time [days]");
  ThreeTime->GetYaxis()->SetTitle("Channels per 0.1 days");

  hs->Draw("nostack");


}
