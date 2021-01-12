#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include "AtlasLabels.h"
#include "AtlasUtils.h"

#include <vector>
#include <iostream>

#define useBTagOnly false

int main(int argc, char** argv) {

  double GeV = 1000.;
  gROOT->SetStyle("ATLAS");
  TCanvas* c1 = new TCanvas();
  
  TFile* input = TFile::Open("plots_data18.root", "update");
 
  TGraphErrors* _mc_scale = (TGraphErrors*) input->Get("mc_scale");
  TGraphErrors* _mc_resolution = (TGraphErrors*) input->Get("mc_resolution");

  TGraphErrors* mc_scale = new TGraphErrors();
  TGraphErrors* mc_resolution = new TGraphErrors();

  int npt = 0;
  for (int i=0; i<_mc_scale->GetN(); i++) {    
    //    if (i==1) continue;
    double yval = _mc_scale->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    mc_scale->SetPoint(npt, _mc_scale->GetX()[i], _mc_scale->GetY()[i]);
    mc_scale->SetPointError(npt, _mc_scale->GetEX()[i], std::max(_mc_scale->GetEY()[i],0.005));
    npt++;
  }
  npt = 0;
  for (int i=0; i<_mc_resolution->GetN(); i++) {    
    double yval = _mc_resolution->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    mc_resolution->SetPoint(npt, _mc_resolution->GetX()[i], _mc_resolution->GetY()[i]);
    mc_resolution->SetPointError(npt, _mc_resolution->GetEX()[i], _mc_resolution->GetEY()[i]);
    npt++;
  }
  mc_scale->Fit("pol1");
  mc_scale->SetLineColor(kRed);
  mc_scale->SetMarkerColor(kRed);
  mc_scale->SetMarkerStyle(20);
  mc_scale->Draw("AP");
  mc_scale->GetXaxis()->SetTitle("Jet reco p_{T} [GeV]");
  mc_scale->GetYaxis()->SetTitle("Z parallel energy scale");
  //  mc_scale->GetYaxis()->SetRangeUser(0.5, 1.0);
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  c1->Print("mc_zee_scale.png");
  c1->Print("mc_zee_scale.pdf");

  mc_resolution->SetLineColor(kRed);
  mc_resolution->SetMarkerColor(kRed);
  mc_resolution->SetMarkerStyle(20);
  mc_resolution->Draw("AP");
  mc_resolution->GetXaxis()->SetTitle("Jet reco p_{T} [GeV]");
  mc_resolution->GetYaxis()->SetTitle("Z parallel energy resolution");
  //  mc_resolution->GetYaxis()->SetRangeUser(0.8, 1.15);
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  c1->Print("mc_zee_resolution.png");
  c1->Print("mc_zee_resolution.pdf");

  mc_scale->Write("mc_scale_fit");
  input->Close();
}
