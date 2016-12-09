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

void plot_Energy() {

  int nbins = 10000; //988, 247, 19
  int max = 10000;

  TFile* f1 = new TFile("AllString_g2root.root");
  
  TTree* t1 = (TTree*)f1->Get("outTree");

  TH1F* string = new TH1F("string", "Energy", nbins, 0, max);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  t1->Draw("ESum2 >> string");


  string->Draw();

}
