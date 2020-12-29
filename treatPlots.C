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

int main(int argc, char** argv) {

  double GeV = 1000.;
  gROOT->SetStyle("ATLAS");
  TCanvas* c1 = new TCanvas();
  
  TFile* input = TFile::Open("plots.root");
 
  TGraphErrors* _old_jet_scale = (TGraphErrors*) input->Get("old_jet_scale");
  TGraphErrors* _old_jet_resolution = (TGraphErrors*) input->Get("old_jet_resolution");
  TGraphErrors* _old_jet_scale_pu = (TGraphErrors*) input->Get("old_jet_scale_pu");
  TGraphErrors* _old_jet_resolution_pu = (TGraphErrors*) input->Get("old_jet_resolution_pu");
  TGraphErrors* _new_jet_scale = (TGraphErrors*) input->Get("new_jet_scale");
  TGraphErrors* _new_jet_resolution = (TGraphErrors*) input->Get("new_jet_resolution");
  TGraphErrors* _new_jet_scale_pu = (TGraphErrors*) input->Get("new_jet_scale_pu");
  TGraphErrors* _new_jet_resolution_pu = (TGraphErrors*) input->Get("new_jet_resolution_pu");

  TGraphErrors* old_jet_scale = new TGraphErrors();
  TGraphErrors* old_jet_resolution = new TGraph();
  TGraphErrors* old_jet_scale_pu = new TGraphErrors();
  TGraphErrors* old_jet_resolution_pu = new TGraph();
  TGraphErrors* new_jet_scale = new TGraphErrors();
  TGraphErrors* new_jet_resolution = new TGraph();
  TGraphErrors* new_jet_scale_pu = new TGraphErrors();
  TGraphErrors* new_jet_resolution_pu = new TGraph();

  int npt = 0;
  for (int i=0; i<_old_jet_scale->GetN(); i++) {    
    double yval = _old_jet_scale->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    old_jet_scale->SetPoint(npt, _old_jet_scale->GetX()[i], _old_jet_scale->GetY()[i]);
    old_jet_scale->SetPointError(npt, _old_jet_scale->GetEX()[i], _old_jet_scale->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_new_jet_scale->GetN(); i++) {    
    double yval = _new_jet_scale->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    new_jet_scale->SetPoint(npt, _new_jet_scale->GetX()[i], _new_jet_scale->GetY()[i]);
    new_jet_scale->SetPointError(npt, _new_jet_scale->GetEX()[i], _new_jet_scale->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_old_jet_scale_pu->GetN(); i++) {    
    double yval = _old_jet_scale_pu->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    old_jet_scale_pu->SetPoint(npt, _old_jet_scale_pu->GetX()[i], _old_jet_scale_pu->GetY()[i]);
    old_jet_scale_pu->SetPointError(npt, _old_jet_scale_pu->GetEX()[i], _old_jet_scale_pu->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_new_jet_scale_pu->GetN(); i++) {    
    double yval = _new_jet_scale_pu->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    new_jet_scale_pu->SetPoint(npt, _new_jet_scale_pu->GetX()[i], _new_jet_scale_pu->GetY()[i]);
    new_jet_scale_pu->SetPointError(npt, _new_jet_scale_pu->GetEX()[i], _new_jet_scale_pu->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_old_jet_resolution->GetN(); i++) {    
    double yval = _old_jet_resolution->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    old_jet_resolution->SetPoint(npt, _old_jet_resolution->GetX()[i], _old_jet_resolution->GetY()[i]);
    old_jet_resolution->SetPointError(npt, _old_jet_resolution->GetEX()[i], _old_jet_resolution->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_new_jet_resolution->GetN(); i++) {    
    double yval = _new_jet_resolution->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    new_jet_resolution->SetPoint(npt, _new_jet_resolution->GetX()[i], _new_jet_resolution->GetY()[i]);
    new_jet_resolution->SetPointError(npt, _new_jet_resolution->GetEX()[i], _new_jet_resolution->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_old_jet_resolution_pu->GetN(); i++) {    
    double yval = _old_jet_resolution_pu->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    old_jet_resolution_pu->SetPoint(npt, _old_jet_resolution_pu->GetX()[i], _old_jet_resolution_pu->GetY()[i]);
    old_jet_resolution_pu->SetPointError(npt, _old_jet_resolution_pu->GetEX()[i], _old_jet_resolution_pu->GetEY()[i]);
    npt++;
  }
  npt = 0;
  for (int i=0; i<_new_jet_resolution_pu->GetN(); i++) {    
    double yval = _new_jet_resolution_pu->GetY()[i];
    if (yval != yval || yval > 2.) continue;
    new_jet_resolution_pu->SetPoint(npt, _new_jet_resolution_pu->GetX()[i], _new_jet_resolution_pu->GetY()[i]);
    new_jet_resolution_pu->SetPointError(npt, _new_jet_resolution_pu->GetEX()[i], _new_jet_resolution_pu->GetEY()[i]);
    npt++;
  }
  
  TF1* old_resolution_func = new TF1("old_resolution_function", "sqrt([0]+[1]/x+[2]/(x*x))", 15, 25);
  TF1* new_resolution_func = new TF1("new_resolution_function", "sqrt([0]+[1]/x+[2]/(x*x))", 15, 25);
  old_resolution_func->SetParameter(0, 0.1);
  old_resolution_func->SetParameter(1, 0.1);
  old_resolution_func->SetParameter(2, 0.1);
  new_resolution_func->SetParameter(0, 0.1);
  new_resolution_func->SetParameter(1, 0.1);
  new_resolution_func->SetParameter(2, 0.1);
  old_resolution_func->SetLineColor(kRed);
  new_resolution_func->SetLineColor(kBlue);
  old_jet_resolution->Fit(old_resolution_func, "V+");
  new_jet_resolution->Fit(new_resolution_func, "V+");

  TMultiGraph* jet_resolution = new TMultiGraph();
  old_jet_resolution->SetLineColor(kRed);
  new_jet_resolution->SetLineColor(kBlue);
  old_jet_resolution->SetMarkerColor(kRed);
  new_jet_resolution->SetMarkerColor(kBlue);
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(old_jet_resolution, "Option 0", "L");
  leg->AddEntry(new_jet_resolution, "Option 1", "L");
  jet_resolution->Add(old_jet_resolution, "P");
  jet_resolution->Add(new_jet_resolution, "P");
  jet_resolution->Draw("A");
  jet_resolution->GetXaxis()->SetTitle("Jet reco p_{T} [GeV]");
  jet_resolution->GetYaxis()->SetTitle("Jet energy resolution");
  jet_resolution->GetYaxis()->SetRangeUser(0.5, 1.1);
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_resolution.png");
  c1->Print("jet_resolution.pdf");

  TMultiGraph* jet_scale = new TMultiGraph();
  old_jet_scale->SetLineColor(kRed);
  old_jet_scale->SetMarkerColor(kRed);
  old_jet_scale->SetMarkerStyle(20);
  new_jet_scale->SetLineColor(kBlue);
  new_jet_scale->SetMarkerColor(kBlue);
  new_jet_scale->SetMarkerStyle(20);
  jet_scale->Add(old_jet_scale, "P");
  jet_scale->Add(new_jet_scale, "P");
  jet_scale->Draw("A");
  jet_scale->GetXaxis()->SetTitle("Jet reco p_{T} [GeV]");
  jet_scale->GetYaxis()->SetTitle("Jet energy scale");
  jet_scale->GetYaxis()->SetRangeUser(0.8, 1.15);
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_scale.png");
  c1->Print("jet_scale.pdf");

  TMultiGraph* jet_scale_pu = new TMultiGraph();
  old_jet_scale_pu->SetLineColor(kRed);
  old_jet_scale_pu->SetMarkerColor(kRed);
  old_jet_scale_pu->SetMarkerStyle(20);
  new_jet_scale_pu->SetLineColor(kBlue);
  new_jet_scale_pu->SetMarkerColor(kBlue);
  new_jet_scale_pu->SetMarkerStyle(20);
  jet_scale_pu->Add(old_jet_scale_pu, "P");
  jet_scale_pu->Add(new_jet_scale_pu, "P");
  jet_scale_pu->Draw("A");
  jet_scale_pu->GetXaxis()->SetTitle("Pile-up #mu (actual)");
  jet_scale_pu->GetYaxis()->SetTitle("Jet energy scale");
  jet_scale_pu->GetYaxis()->SetRangeUser(0.8, 1.15);
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_scale_pu.png");
  c1->Print("jet_scale_pu.pdf");

  input->Close();
}
