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

void plot_calibrationMinMax_channel() {

  int nbins = 988;
  int max = 988;

  TFile *f1 = new TFile("AllPeaks_rate.root");

  TTree* t1 = (TTree*)f1->Get("tree");
  TTree* t2 = (TTree*)f1->Get("tree");
  TTree* t3 = (TTree*)f1->Get("tree");
  TTree* t4 = (TTree*)f1->Get("tree");
  



  //TH1F* MaxTime = new TH1F("MaxTime", "Four Peaks", nbins, 0, max);
  //TH1F* MinTime = new TH1F("MinTime", "One Peak", nbins, 0, max);
  //TH1F* TwoTime = new TH1F("TwoTime", "Two Peaks", nbins, 0, max);
  //TH1F* ThreeTime = new TH1F("ThreeTime", "Three Peaks", nbins, 0, max);

 


  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  c2->cd();

  /* t1->Draw("CalibrationTime_Max:Channel >> MaxTime");
  t1->Draw("CalibrationTime_Three:Channel >> ThreeTime","","same");
  t1->Draw("CalibrationTime_Two:Channel >> TwoTime","","same");
  t1->Draw("CalibrationTime_Min:Channel >> MinTime","","same");
  */

  int n1 = t1->Draw("CalibrationTime_Max:Channel", "", "goff");
  TGraph *MaxTime = new TGraph(n1, tree->GetV2(), tree->GetV1());

  int n2 = t2->Draw("CalibrationTime_Three:Channel","","goff");
  TGraph *ThreeTime = new TGraph(n2, tree->GetV2(), tree->GetV1());

  int n3 = t3->Draw("CalibrationTime_Two:Channel","","goff");
  TGraph *TwoTime = new TGraph(n3, tree->GetV2(), tree->GetV1());

  int n4 = t4->Draw("CalibrationTime_Min:Channel","","goff");
  TGraph *MinTime = new TGraph(n4, tree->GetV2(), tree->GetV1());


  
  MaxTime->SetMarkerColor(kBlack);
  TwoTime->SetMarkerColor(kSpring);
  ThreeTime->SetMarkerColor(kMagenta);
  MinTime->SetMarkerColor(kCyan);

  MaxTime->SetFillColor(0);
  TwoTime->SetFillColor(0);
  ThreeTime->SetFillColor(0);
  MinTime->SetFillColor(0);

  MaxTime->SetLineColor(0);
  TwoTime->SetLineColor(0);
  ThreeTime->SetLineColor(0);
  MinTime->SetLineColor(0);

  MaxTime->SetMarkerStyle(8);
  MinTime->SetMarkerStyle(8);
  TwoTime->SetMarkerStyle(26);
  ThreeTime->SetMarkerStyle(32);


  MaxTime->SetTitle("Four Peaks");
  TwoTime->SetTitle("Two Peaks");
  ThreeTime->SetTitle("Three Peaks");
  MinTime->SetTitle("One Peak");



  TMultiGraph *mg = new TMultiGraph();

  mg->Add(MinTime,"p");
  mg->Add(TwoTime, "p");
  mg->Add(ThreeTime, "p");
  mg->Add(MaxTime, "p");


  mg->Draw("a");

  /*
  MinTime->Draw("ap");
  TwoTime->Draw("same");
  ThreeTime->Draw("same");
  MaxTime->Draw("same");
  */
 
}
