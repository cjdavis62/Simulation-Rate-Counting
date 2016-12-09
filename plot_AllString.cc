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

void plot_AllString() {

  int nbins = 988; //988, 247, 19

  TFile* f1 = new TFile("AllString_g2root.root");
  
  TTree* t1 = (TTree*)f1->Get("outTree");
    
  TH1F* InteriorString = new TH1F("FullString", "FullString", nbins, 0, 988);
  
  //TCut cut1 = "Ener1 > 2605";
  //TCut cut2 = "Ener1 < 2625";
  //TCut cut3 = cut1 && cut2; 

  TCanvas* c2 = new TCanvas("c2", "c2", 1900, 1000);
  t1->Draw("Channel >> FullString");
  
  FullString->Scale(0.01115); // to per hour

  c2->cd();
  FullString->Scale(0.2778); // to mHz
  FullString->SetTitle("MC Events per Channel");
  FullString->GetXaxis()->SetTitle("Channel");
  FullString->GetYaxis()->SetTitle("Rate per Channel [mHz]");
  FullString->SetFillColor(kBlue);
  
  FullString->Draw();
}
