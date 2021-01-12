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
#include "TRatioPlot.h"
#include "AtlasLabels.h"
#include "AtlasUtils.h"

#include <vector>
#include <iostream>

#define useBTagOnly false

int main(int argc, char** argv) {

  double GeV = 1000.;
  gROOT->SetStyle("ATLAS");
  TCanvas* c1 = new TCanvas();
  std::cout << "1" << std::endl;
  TFile* data_file = TFile::Open("plots_data18.root");
  TFile* mc_file = TFile::Open("plots_zee.root");

  std::cout << "2" << std::endl; 
  TGraphErrors* zee_mc_scale = (TGraphErrors*) mc_file->Get("mc_scale_fit");
  TGraphErrors* zee_data_scale = (TGraphErrors*) data_file->Get("mc_scale_fit");

  std::cout << "3" << std::endl;
  TH1D* data_hist = new TH1D("data_hist", "data_hist", 10, 15., 25.);
  TH1D* mc_hist = new TH1D("mc_hist", "mc_hist", 10, 15., 25.);
  
  for (int i=0; i<10; i++) {
    data_hist->SetBinContent(i+1, zee_data_scale->GetY()[i]+1);
    mc_hist->SetBinContent(i+1, zee_mc_scale->GetY()[i]+1);
    data_hist->SetBinError(i+1, std::max(0.01,zee_data_scale->GetEY()[i]));
    mc_hist->SetBinError(i+1, std::max(0.01,zee_mc_scale->GetEY()[i]));
  }
  TLegend* leg = new TLegend(0.65,0.7,0.85,0.85);
  leg->SetBorderSize(0);

  mc_hist->SetLineColor(kBlue);
  data_hist->GetYaxis()->SetRangeUser(0.6, 1.1);
  mc_hist->GetYaxis()->SetRangeUser(0.6, 1.1);
  leg->AddEntry(mc_hist, "Z+jets MC", "l");
  leg->AddEntry(data_hist, "2018 Data", "pl");

  std::cout << "4" << std::endl;
  TRatioPlot* rp1 = new TRatioPlot(data_hist, mc_hist, "divsym");
  rp1->Draw();
  rp1->GetLowerRefYaxis()->SetTitle("Ratio");
  rp1->GetLowerRefYaxis()->SetRangeUser(0.7, 1.3);
  rp1->GetUpperRefYaxis()->SetTitle("p_{T,jet}^{#parallel}/p_{T,Z}");
  c1->Update();
  rp1->GetUpperPad()->SetFrameBorderSize(0);
  rp1->GetLowerPad()->SetFrameBorderSize(0);
  rp1->GetUpperPad()->SetLeftMargin(0.2);
  rp1->GetLowerPad()->SetLeftMargin(0.2);
  c1->Update();
  ATLASLabel(0.25, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("data_mc_zee_scale.png");
  c1->Print("data_mc_zee_scale.pdf");

  std::cout << "5" << std::endl;
  data_file->Close();
  mc_file->Close();
}
