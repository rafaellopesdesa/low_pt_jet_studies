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
  TCanvas* c1 = new TCanvas("c1", "c1", 2000, 600);
  TCanvas* c2 = new TCanvas();
  
  TFile* input = TFile::Open("plots.root");

  TH1D* new_all_cutflow = (TH1D*) input->Get("new_all_cutflow");
  TH1D* new_selected_cutflow = (TH1D*) input->Get("new_selected_cutflow");
  TH1D* new_thiso_cutflow = (TH1D*) input->Get("new_thiso_cutflow");
  TH1D* new_expiso_cutflow = (TH1D*) input->Get("new_expiso_cutflow");
  TH1D* new_matched_cutflow = (TH1D*) input->Get("new_matched_cutflow");
  TH1D* old_all_cutflow = (TH1D*) input->Get("old_all_cutflow");
  TH1D* old_selected_cutflow = (TH1D*) input->Get("old_selected_cutflow");
  TH1D* old_thiso_cutflow = (TH1D*) input->Get("old_thiso_cutflow");
  TH1D* old_expiso_cutflow = (TH1D*) input->Get("old_expiso_cutflow");
  TH1D* old_matched_cutflow = (TH1D*) input->Get("old_matched_cutflow");


  TH1D* new_summary = new TH1D("new_summary", "", 5, 0., 5.);
  TH1D* old_summary = new TH1D("old_summary", "", 5, 0., 5.);

  new_all_cutflow->SetLineColor(kOrange);
  old_all_cutflow->SetLineColor(kOrange);
  old_all_cutflow->SetLineStyle(kDashed);

  new_selected_cutflow->SetLineColor(kRed);
  old_selected_cutflow->SetLineColor(kRed);
  old_selected_cutflow->SetLineStyle(kDashed);

  new_thiso_cutflow->SetLineColor(kGreen+3);
  old_thiso_cutflow->SetLineColor(kGreen+3);
  old_thiso_cutflow->SetLineStyle(kDashed);

  new_expiso_cutflow->SetLineColor(kBlue);
  old_expiso_cutflow->SetLineColor(kBlue);
  old_expiso_cutflow->SetLineStyle(kDashed);

  new_matched_cutflow->SetLineColor(kMagenta);
  old_matched_cutflow->SetLineColor(kMagenta);
  old_matched_cutflow->SetLineStyle(kDashed);

  new_summary->SetBinContent(1, new_all_cutflow->Integral());
  new_summary->SetBinContent(2, new_selected_cutflow->Integral());
  new_summary->SetBinContent(3, new_thiso_cutflow->Integral());
  new_summary->SetBinContent(4, new_expiso_cutflow->Integral());
  new_summary->SetBinContent(5, new_matched_cutflow->Integral());

  old_summary->SetBinContent(1, old_all_cutflow->Integral());
  old_summary->SetBinContent(2, old_selected_cutflow->Integral());
  old_summary->SetBinContent(3, old_thiso_cutflow->Integral());
  old_summary->SetBinContent(4, old_expiso_cutflow->Integral());
  old_summary->SetBinContent(5, old_matched_cutflow->Integral());

  old_summary->SetLineColor(kRed);
  new_summary->SetLineColor(kBlue);

  old_summary->GetXaxis()->SetBinLabel(1, "all");
  old_summary->GetXaxis()->SetBinLabel(2, "sel");
  old_summary->GetXaxis()->SetBinLabel(3, "thiso");
  old_summary->GetXaxis()->SetBinLabel(4, "exiso");
  old_summary->GetXaxis()->SetBinLabel(5, "match");
  old_summary->GetYaxis()->SetTitle("Yield");

  new_summary->GetXaxis()->SetBinLabel(1, "all");
  new_summary->GetXaxis()->SetBinLabel(2, "sel");
  new_summary->GetXaxis()->SetBinLabel(3, "thiso");
  new_summary->GetXaxis()->SetBinLabel(4, "exiso");
  new_summary->GetXaxis()->SetBinLabel(5, "match");
  new_summary->GetYaxis()->SetTitle("Yield");

  c1->cd();
  new_all_cutflow->GetYaxis()->SetRangeUser(0, 1.3*new_all_cutflow->GetMaximum());
  new_all_cutflow->SetTitle("cutflow per file;File number;Yield");
  new_all_cutflow->Draw("hist");
  old_all_cutflow->Draw("hist,same");
  new_selected_cutflow->Draw("hist,same");
  old_selected_cutflow->Draw("hist,same");  
  new_thiso_cutflow->Draw("hist,same");
  old_thiso_cutflow->Draw("hist,same");
  new_expiso_cutflow->Draw("hist,same");
  old_expiso_cutflow->Draw("hist,same");
  new_matched_cutflow->Draw("hist,same");
  old_matched_cutflow->Draw("hist,same");
  
  ATLASLabel(0.2, 0.8, "Work in Progress");
  c1->Print("cutflow.png");
  c1->Print("cutflow.pdf");

  c2->cd();
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(old_summary, "Option 0", "L");
  leg->AddEntry(new_summary, "Option 1", "L");  
  new_summary->Draw("hist");
  old_summary->Draw("hist,same");
  leg->Draw();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  c2->Print("summary.png");
  c2->Print("summary.pdf");
  c2->SetLogy();
  c2->Print("summary_log.png");
  c2->Print("summary_log.pdf");

  input->Close();

}
