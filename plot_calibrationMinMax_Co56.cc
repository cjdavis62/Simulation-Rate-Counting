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

void plot_calibrationMinMax_Co56() {

  int nbins = 300;
  double max = 3.0;

  TFile *f1 = new TFile("AllPeaks_rate.root");

  TTree* t1 = (TTree*)f1->Get("tree");
  
  TH1F* Time2599 = new TH1F("Time2599", "Time2599", 90, 4, 13);
  TH1F* Time847 = new TH1F("Time847", "Time847", 140, 0, 3.5);
  TH1F* Time1238 = new TH1F("Time1238", "Time1238", 140, 0, 7);
  TH1F* Time511 = new TH1F("Time511", "Time511", 90, 1, 12);
  TH1F* Time1771 = new TH1F("Time1771", "Time1771", 250, 5, 30);
  TH1F* Time1037 = new TH1F("Time1037", "Time1037", 200, 5, 45);
  TH1F* Time3254 = new TH1F("Time3254", "Time3254", 250, 10, 35);
  TH1F* Time2035 = new TH1F("Time2035", "Time2035", 300, 5, 35);
  TH1F* Time1360 = new TH1F("Time1360", "Time1360", 150, 10, 110);
  TH1F* Time3202 = new TH1F("Time3202", "Time3202", 130, 25, 85);
  TH1F* Time3451 = new TH1F("Time3451", "Time3451", 100, 50, 500);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  TCanvas* c3 = new TCanvas("c3", "c3", 800, 800);
  TCanvas* c4 = new TCanvas("c4", "c4", 800, 800);
  TCanvas* c5 = new TCanvas("c5", "c5", 800, 800);
  TCanvas* c6 = new TCanvas("c6", "c6", 800, 800);
  TCanvas* c7 = new TCanvas("c7", "c7", 800, 800);
  TCanvas* c8 = new TCanvas("c8", "c8", 800, 800);
  TCanvas* c9 = new TCanvas("c9", "c9", 800, 800);
  TCanvas* c10 = new TCanvas("c10", "c10", 800, 800);
  TCanvas* c11 = new TCanvas("c11", "c11", 800, 800);

  c1->cd();
  t1->Draw("CalibrationTime_2599 >> Time2599");
  c2->cd();
  t1->Draw("CalibrationTime_847 >> Time847");
  c3->cd();
  t1->Draw("CalibrationTime_1238 >> Time1238");
  c4->cd();
  t1->Draw("CalibrationTime_511 >> Time511");
  c5->cd();
  t1->Draw("CalibrationTime_1771 >> Time1771");
  c6->cd();
  t1->Draw("CalibrationTime_1037 >> Time1037");
  c7->cd();
  t1->Draw("CalibrationTime_3254 >> Time3254");
  c8->cd();
  t1->Draw("CalibrationTime_2035 >> Time2035");
  c9->cd();
  t1->Draw("CalibrationTime_1360 >> Time1360");
  c10->cd();
  t1->Draw("CalibrationTime_3202 >> Time3202");
  c11->cd();
  t1->Draw("CalibrationTime_3451 >> Time3451");
 
}
