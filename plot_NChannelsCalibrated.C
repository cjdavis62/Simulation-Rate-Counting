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

void plot_NChannelsCalibrated() {

  int maxT = 25;
  const int nbins = 24 * maxT;


  double channels_2615[nbins+1];
  double channels_969[nbins+1];
  double channels_911[nbins+1];
  double channels_583[nbins+1];
  double channels_338[nbins+1];
  double channels_239[nbins+1];
  double x[nbins+1];

  TFile *f1 = new TFile("AllPeaks_rate.root");
  
  TTree* t1 = (TTree*)f1->Get("tree");
  
  TH1F * hist_2615 = new TH1F("hist_2615", "2615 Time", nbins, -1, maxT);
  TH1F * hist_969 = new TH1F("hist_969", "969 Time", nbins, -1, maxT);
  TH1F * hist_911 = new TH1F("hist_911", "911 Time", nbins, -1, maxT);
  TH1F * hist_583 = new TH1F("hist_583", "583 Time", nbins, -1, maxT);
  TH1F * hist_338 = new TH1F("hist_338", "338 Time", nbins, -1, maxT);
  TH1F * hist_239 = new TH1F("hist_239", "239 Time", nbins, -1, maxT);
  

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 800);
  c1->cd();

  t1->Draw("CalibrationTime_2615 >> hist_2615");
  t1->Draw("CalibrationTime_969 >> hist_969");
  t1->Draw("CalibrationTime_911 >> hist_911");
  t1->Draw("CalibrationTime_583 >> hist_583");
  t1->Draw("CalibrationTime_338 >> hist_338");
  t1->Draw("CalibrationTime_239  >> hist_239" );

  hist_2615->Draw();
  hist_969->Draw();
  hist_911->Draw();
  hist_583->Draw();
  hist_338->Draw();
  hist_239->Draw();

  TCanvas * c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->cd()

  for (int n = 0; n < nbins+1; n++)
    {
      x[n] = n * double(maxT)/double(nbins);
      channels_2615[n] = 987 - hist_2615->Integral(0, n);
      channels_969[n] = 987 - hist_969->Integral(0, n);
      channels_911[n] = 987 - hist_911->Integral(0, n);
      channels_583[n] = 987 - hist_583->Integral(0, n);
      channels_338[n] = 987 - hist_338->Integral(0, n);
      channels_239[n] = 987 - hist_239->Integral(0, n);
      std::cout << channels_583[n] << "\t" << x[n] << std::endl;
    }

  c2->SetLogy();

  TGraph* graph_2615 = new TGraph(nbins+1,x,channels_2615);
  graph_2615->SetTitle("2615 keV peak");
  TGraph* graph_969 = new TGraph(nbins+1,x,channels_969);
  graph_969->SetTitle("969 keV peak");
  TGraph* graph_911 = new TGraph(nbins+1,x,channels_911);
  graph_911->SetTitle("911 keV peak");
  TGraph* graph_583 = new TGraph(nbins+1,x,channels_583);
  graph_583->SetTitle("583 keV peak");
  TGraph* graph_338 = new TGraph(nbins+1,x,channels_338);
  graph_338->SetTitle("338 keV peak");
  TGraph* graph_239 = new TGraph(nbins+1,x,channels_239);
  graph_239->SetTitle("239 keV peak");



  //graph_2615->GetYaxis()->SetRangeUser(1e-1, 1e3);
  graph_2615->SetLineColor(kGreen);
  graph_969->SetLineColor(kMagenta);
  graph_911->SetLineColor(kOrange);
  graph_583->SetLineColor(kRed);
  graph_338->SetLineColor(kAzure);
  graph_239->SetLineColor(kCyan);

 
  TMultiGraph * mg = new TMultiGraph();
  mg->Add(graph_2615, "AC");
  mg->Add(graph_969, "AC");
  mg->Add(graph_911, "AC");
  mg->Add(graph_583, "AC");
  mg->Add(graph_338, "AC");
  mg->Add(graph_239, "AC");

  mg->Draw("a");

  mg->GetYaxis()->SetRangeUser(1,1000);
  mg->GetXaxis()->SetRangeUser(0, maxT);
  mg->GetYaxis()->SetTitle("Uncalibrated Channels");
  mg->GetXaxis()->SetTitle("Time [Days]");

  
}
