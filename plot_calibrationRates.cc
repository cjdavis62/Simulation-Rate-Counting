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

void plot_calibrationRates() {

  int nbins = 30;
  int max = 6;

  TFile *f1 = new TFile("AllPeaks_rate.root");

  TTree* t1 = (TTree*)f1->Get("tree");
  TH1F* Rate583 = new TH1F("Rate583", "583 Peak Rate", 200, 0, 5);
  TH1F* Rate2615 = new TH1F("Rate2615", "2615 Peak Rate", 200, 0, 5);
  TH1F* Rate969 = new TH1F("Rate969", "969 Peak Rate", 200, 0, 5);
  TH1F* Rate911 = new TH1F("Rate911", "911 Peak Rate", 200, 0, 5);

  TH1F* Rate239 = new TH1F("Rate239", "239 Peak Rate", 200,0,5);
  TH1F* Rate338 = new TH1F("Rate338", "338 Peak Rate", 200, 0, 5); // two peaks here

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  t1->Draw("Rate_239 >> Rate239","","goff");
  t1->Draw("Rate_338 >> Rate338","","goff");
  t1->Draw("Rate_583 >> Rate583","","goff");
  t1->Draw("Rate_2615 >> Rate2615","","goff");
  t1->Draw("Rate_969 >> Rate969","","goff");
  t1->Draw("Rate_911 >> Rate911","","goff");

  c2->Divide(3,2);
  c2->cd(1);
  Rate583->Draw();
  Rate583->GetXaxis()->SetRangeUser(0.,2.);
  Rate583->GetXaxis()->SetTitle("Rate [mHz]");
  Rate583->GetYaxis()->SetTitle("Channels / 0.025 mHz");

  c2->cd(2);
  Rate2615->Draw();
  Rate2615->GetXaxis()->SetRangeUser(0.,2.);
  Rate2615->GetXaxis()->SetTitle("Rate [mHz]");
  Rate2615->GetYaxis()->SetTitle("Channels / 0.025 mHz");

  c2->cd(3);
  Rate969->Draw();
  Rate969->GetXaxis()->SetRangeUser(0.,2.);
  Rate969->GetXaxis()->SetTitle("Rate [mHz]");
  Rate969->GetYaxis()->SetTitle("Channels / 0.025 mHz");

  c2->cd(4);
  Rate911->Draw();
  Rate911->GetXaxis()->SetRangeUser(0.,2.);
  Rate911->GetXaxis()->SetTitle("Rate [mHz]");
  Rate911->GetYaxis()->SetTitle("Channels / 0.025 mHz");

  c2->cd(5);
  Rate338->Draw();
  Rate338->GetXaxis()->SetRangeUser(0.,2.);
  Rate338->GetXaxis()->SetTitle("Rate [mHz]");
  Rate338->GetYaxis()->SetTitle("Channels / 0.025 mHz");
 
  c2->cd(6);
  Rate239->Draw();
  Rate239->GetXaxis()->SetRangeUser(0.,2.);
  Rate239->GetXaxis()->SetTitle("Rate [mHz]");
  Rate239->GetYaxis()->SetTitle("Channels / 0.025 mHz");

}
