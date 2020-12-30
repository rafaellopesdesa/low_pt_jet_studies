#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h" 
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooDoubleCBFast.h"
#include "RooGaussian.h"
#include "RooPlot.h" 
#include "RooDataHist.h"
#include "RooFitResult.h"

#include "AtlasLabels.h"
#include "AtlasUtils.h"

#include <vector>
#include <iostream>
#include <utility>

#define debug false
#define debug2 false
#define useBTagOnly false

void PutOverflowInLastBin(TH1D* h) {

  int nb = h->GetNbinsX();

  double bc_new = h->Integral(nb,nb+1);
  double be_new = TMath::Sqrt(h->GetBinError(nb)*h->GetBinError(nb) + h->GetBinError(nb+1)*h->GetBinError(nb+1));
  h->SetBinContent(nb,bc_new);
  h->SetBinError(nb,be_new);
  h->SetBinContent(nb+1,0);
  h->SetBinError(nb+1,0);
}

void PutUnderflowInFirstBin(TH1D* h) {

  int nb = 0;

  double bc_new = h->Integral(nb,nb+1);
  double be_new = TMath::Sqrt(h->GetBinError(nb)*h->GetBinError(nb) + h->GetBinError(nb+1)*h->GetBinError(nb+1));
  h->SetBinContent(nb+1,bc_new);
  h->SetBinError(nb+1,be_new);
  h->SetBinContent(nb,0);
  h->SetBinError(nb,0);
}

double calculateResolution(double sigma, double alpha1, double alpha2, double n1, double n2, double cl) {

  double leftval=0., rightval=0.;

  double N = (n1/(n1-1))*(1/alpha1)*TMath::Exp(-alpha1*alpha1*0.5)+(n2/(n2-1))*(1/alpha2)*TMath::Exp(-alpha2*alpha2*0.5) + TMath::Sqrt(TMath::Pi()/2)*(TMath::Erf(alpha2/TMath::Sqrt(2))+TMath::Erf(alpha1/TMath::Sqrt(2)));
  
  double right_gausarea = TMath::Sqrt(TMath::Pi()/2)*TMath::Erf(alpha2/TMath::Sqrt(2))/N;
  if (cl/2 < right_gausarea) {
    rightval = TMath::Sqrt(2)*TMath::ErfInverse(N*cl/TMath::Sqrt(2*TMath::Pi()));
  } else {
    double leftover = (cl/2 - right_gausarea) * N * TMath::Exp(alpha2*alpha2*0.5) * (alpha2/n2);
    double right_u = TMath::Power(1+(1-n2)*leftover,1/(1-n2));
    rightval = alpha2 - (n2/alpha2)*(1-right_u);
  }
  double left_gausarea = TMath::Sqrt(TMath::Pi()/2)*TMath::Erf(alpha1/TMath::Sqrt(2))/N;
  if (cl/2 < left_gausarea) {
    leftval = TMath::Sqrt(2)*TMath::ErfInverse(N*cl/TMath::Sqrt(2*TMath::Pi()));
  } else {
    double leftover = (cl/2 - left_gausarea) * N * TMath::Exp(alpha1*alpha1*0.5) * (alpha1/n1);
    double left_u = TMath::Power(1+(1-n1)*leftover,1/(1-n1));
    leftval = alpha1 - (n1/alpha1)*(1-left_u);
  }
  return (rightval + leftval)*sigma;

}
  
double doFit(RooDataSet& dataset, RooRealVar& response, RooRealVar& pt_reco, RooRealVar& option, TString cut, 
	     RooRealVar& mean, RooRealVar& sigma, RooRealVar& alpha1, RooRealVar& alpha2, RooRealVar& n1, RooRealVar& n2, 
	     RooGaussian& gaus, RooDoubleCBFast& cbfunc, double xmin, double xmax,
	     TCanvas* c1, TString name) {
  
  RooDataSet* reduced_dataset = (RooDataSet*) dataset.reduce(RooArgSet(response,pt_reco,option), cut.Data());
  std::cout << "Fitting " << reduced_dataset->numEntries() << " entries" << std::endl;
  mean.setVal(reduced_dataset->mean(response)-0.1);
  sigma.setVal(reduced_dataset->sigma(response));
  alpha1.setVal(1);
  alpha2.setVal(2);
  n1.setVal(20);
  n2.setVal(5);
  mean.setConstant(false);
  sigma.setConstant(false);
  alpha1.setConstant(false);
  alpha2.setConstant(false);
  n1.setConstant(true);
  n2.setConstant(true);
  gaus.fitTo(*reduced_dataset,RooFit::Range(0.2*xmin, 0.2*xmax), RooFit::NumCPU(4));
  cbfunc.fitTo(*reduced_dataset,RooFit::Range(xmin, xmax), RooFit::NumCPU(4));
  n1.setConstant(false);
  n2.setConstant(false);
  cbfunc.fitTo(*reduced_dataset,RooFit::Range(xmin, xmax), RooFit::NumCPU(4), RooFit::Save());
  RooFitResult *result = cbfunc.fitTo(*reduced_dataset,RooFit::Range(xmin, xmax), RooFit::NumCPU(4), RooFit::Save());
  if (result->status() > 2 && false) {
    n1.setVal(20);
    n2.setVal(5);
    n1.setConstant(true);
    n2.setConstant(true);
    cbfunc.fitTo(*reduced_dataset,RooFit::Range(xmin, xmax), RooFit::NumCPU(4));
  }

  RooPlot* frame = response.frame(RooFit::Range(xmin,xmax)) ; 
  reduced_dataset->plotOn(frame);
  cbfunc.plotOn(frame);

  c1->cd();
  frame->Draw();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, cut.Data());
  myText(0.2, 0.64, kBlack, TString::Format("scale = %.3f #pm %.3f", mean.getVal()+1., mean.getError()));
  double resolution = calculateResolution(sigma.getVal(), alpha1.getVal(), alpha2.getVal(), n1.getVal(), n2.getVal(), 0.68);
  myText(0.2, 0.56, kBlack, TString::Format("resolution = %.3f", resolution));
  c1->Print(name.Data());

  return resolution;
}

int main(int argc, char** argv) {

  double GeV = 1000.;
  int ntot = 10;
  
  TH1D* new_all_cutflow = new TH1D("new_all_cutflow", "new_all_cutflow", 450, 0, 450);
  TH1D* new_selected_cutflow = new TH1D("new_selected_cutflow", "new_selected_cutflow", 450, 0, 450);
  TH1D* new_thiso_cutflow = new TH1D("new_thiso_cutflow", "new_thiso_cutflow", 450, 0, 450);
  TH1D* new_expiso_cutflow = new TH1D("new_expiso_cutflow", "new_expiso_cutflow", 450, 0, 450);
  TH1D* new_matched_cutflow = new TH1D("new_matched_cutflow", "new_matched_cutflow", 450, 0, 450);

  TH1D* old_all_cutflow = new TH1D("old_all_cutflow", "old_all_cutflow", 450, 0, 450);
  TH1D* old_selected_cutflow = new TH1D("old_selected_cutflow", "old_selected_cutflow", 450, 0, 450);
  TH1D* old_thiso_cutflow = new TH1D("old_thiso_cutflow", "old_thiso_cutflow", 450, 0, 450);
  TH1D* old_expiso_cutflow = new TH1D("old_expiso_cutflow", "old_expiso_cutflow", 450, 0, 450);
  TH1D* old_matched_cutflow = new TH1D("old_matched_cutflow", "old_matched_cutflow", 450, 0, 450);

  TH1D* new_reco_pt = new TH1D("new_reco_pt", "new_reco_pt", 40, 15., 25.);
  TH1D* old_reco_pt = new TH1D("old_reco_pt", "old_reco_pt", 40, 15., 25.);

  TH1D* new_reco_all_pt = new TH1D("new_reco_all_pt", "new_reco_all_pt", 120, 5., 35.);
  TH1D* old_reco_all_pt = new TH1D("old_reco_all_pt", "old_reco_all_pt", 120, 5., 35.);

  TH1D* new_reco_iso_pt = new TH1D("new_reco_iso_pt", "new_reco_iso_pt", 25, 0., 50.);
  TH1D* old_reco_iso_pt = new TH1D("old_reco_iso_pt", "old_reco_iso_pt", 25, 0., 50.);

  TH1D* new_truth_pt = new TH1D("new_truth_pt", "new_truth_pt", 120, 5., 35.);
  TH1D* old_truth_pt = new TH1D("old_truth_pt", "old_truth_pt", 120, 5., 35.);

  TH1D* new_reco_eta = new TH1D("new_reco_eta", "new_reco_eta", 40, -2.5, 2.5);
  TH1D* old_reco_eta = new TH1D("old_reco_eta", "old_reco_eta", 40, -2.5, 2.5);

  TH1D* new_reco_all_eta = new TH1D("new_reco_all_eta", "new_reco_all_eta", 48, -3.0, 3.0);
  TH1D* old_reco_all_eta = new TH1D("old_reco_all_eta", "old_reco_all_eta", 48, -3.0, 3.0);

  TH1D* new_truth_eta = new TH1D("new_truth_eta", "new_truth_eta", 48, -3.0, 3.0);
  TH1D* old_truth_eta = new TH1D("old_truth_eta", "old_truth_eta", 48, -3.0, 3.0);

  TH1D* new_reco_truth_dist_low = new TH1D("new_reco_truth_dist_low", "", 30, 0., 3.);
  TH1D* old_reco_truth_dist_low = new TH1D("old_reco_truth_dist_low", "", 30, 0., 3.);

  TH1D* new_reco_truth_dist_high = new TH1D("new_reco_truth_dist_high", "", 30, 0., 3.);
  TH1D* old_reco_truth_dist_high = new TH1D("old_reco_truth_dist_high", "", 30, 0., 3.);
 
  TH1D* new_reco_response_low = new TH1D("new_reco_response_low", "", 40, -2., 2.);
  TH1D* old_reco_response_low = new TH1D("old_reco_response_low", "", 40, -2., 2.);

  TH1D* new_reco_response_high = new TH1D("new_reco_response_high", "", 40, -2., 2.);
  TH1D* old_reco_response_high = new TH1D("old_reco_response_high", "", 40, -2., 2.);

  TH1D* new_reco_eta_resolution_low = new TH1D("new_reco_eta_resolution_low", "", 40, -0.3, 0.3);
  TH1D* old_reco_eta_resolution_low = new TH1D("old_reco_eta_resolution_low", "", 40, -0.3, 0.3);

  TH1D* new_reco_eta_resolution_high = new TH1D("new_reco_eta_resolution_high", "", 40, -0.3, 0.3);
  TH1D* old_reco_eta_resolution_high = new TH1D("old_reco_eta_resolution_high", "", 40, -0.3, 0.3);

  TH1D* new_reco_phi_resolution_low = new TH1D("new_reco_phi_resolution_low", "", 40, -0.3, 0.3);
  TH1D* old_reco_phi_resolution_low = new TH1D("old_reco_phi_resolution_low", "", 40, -0.3, 0.3);

  TH1D* new_reco_phi_resolution_high = new TH1D("new_reco_phi_resolution_high", "", 40, -0.3, 0.3);
  TH1D* old_reco_phi_resolution_high = new TH1D("old_reco_phi_resolution_high", "", 40, -0.3, 0.3);

  TEfficiency* new_hs_jvt_efficiency = new TEfficiency("new_hs_jvt_efficiency", "", ntot, 15., 25.);
  TEfficiency* old_hs_jvt_efficiency = new TEfficiency("old_hs_jvt_efficiency", "", ntot, 15., 25.);

  TEfficiency* new_pu_jvt_efficiency = new TEfficiency("new_pu_jvt_efficiency", "", ntot, 15., 25.);
  TEfficiency* old_pu_jvt_efficiency = new TEfficiency("old_pu_jvt_efficiency", "", ntot, 15., 25.);

  TEfficiency* new_btag_b_efficiency = new  TEfficiency("new_btag_b_efficiency", "", ntot, 15., 25.);
  TEfficiency* old_btag_b_efficiency = new  TEfficiency("old_btag_b_efficiency", "", ntot, 15., 25.);

  TEfficiency* new_btag_c_efficiency = new  TEfficiency("new_btag_c_efficiency", "", ntot, 15., 25.);
  TEfficiency* old_btag_c_efficiency = new  TEfficiency("old_btag_c_efficiency", "", ntot, 15., 25.);

  TEfficiency* new_btag_light_efficiency = new  TEfficiency("new_btag_light_efficiency", "", ntot, 15., 25.);
  TEfficiency* old_btag_light_efficiency = new  TEfficiency("old_btag_light_efficiency", "", ntot, 15., 25.);
  
  RooRealVar _response("_response", "p_{T}^{reco}/p_{T}^{truth}-1", 0);
  RooRealVar _pt_reco("_pt_reco", "_pt_reco", 0);
  RooRealVar _eta_reco("_eta_reco", "_eta_reco", 0);
  RooRealVar _pt_truth("_pt_truth", "_pt_truth", 0);
  RooRealVar _eta_truth("_eta_truth", "_eta_truth", 0);
  RooRealVar _option("_option", "_option", 0);
  RooRealVar _mu_actual("_mu_actual", "_mu_actual", 0);
  RooRealVar _weight("_weight", "_weight", 0);
  RooArgSet _hs_jet("hs_jet");
  _hs_jet.add(_response);
  _hs_jet.add(_pt_reco);
  _hs_jet.add(_eta_reco);
  _hs_jet.add(_pt_truth);
  _hs_jet.add(_eta_truth);
  _hs_jet.add(_option);
  _hs_jet.add(_mu_actual);
  _hs_jet.add(_weight);
  RooDataSet _dataset("_dataset", "", _hs_jet, "_weight");

  float weight_mc, weight_pileup;
  ULong64_t eventNumber;
  float mu_actual;

  std::vector<float> reco_jet_pt, reco_jet_eta, reco_jet_phi, reco_jet_e, reco_jet_jvt ;
  std::vector<float> failJvt_reco_jet_pt, failJvt_reco_jet_eta, failJvt_reco_jet_phi, failJvt_reco_jet_e, failJvt_reco_jet_jvt ;
  std::vector<float> truth_jet_pt, truth_jet_eta, truth_jet_phi, truth_jet_e;
  std::vector<float> truth_mu_pt, truth_mu_eta, truth_mu_phi, truth_mu_e;
  std::vector<float> truth_el_pt, truth_el_eta, truth_el_phi, truth_el_e;
  std::vector<float> truth_tau_pt, truth_tau_eta, truth_tau_phi, truth_tau_e;
  std::vector<char> reco_jet_b;
  std::vector<int> truth_jet_b;
  std::vector<int> truth_jet_c;

  std::vector<float>* p_truth_jet_pt = &truth_jet_pt;
  std::vector<float>* p_truth_jet_eta = &truth_jet_eta;
  std::vector<float>* p_truth_jet_phi = &truth_jet_phi;
  std::vector<float>* p_truth_jet_e = &truth_jet_e;
  std::vector<int>* p_truth_jet_b = &truth_jet_b;
  std::vector<int>* p_truth_jet_c = &truth_jet_c;

  std::vector<float>* p_reco_jet_pt = &reco_jet_pt;
  std::vector<float>* p_reco_jet_eta = &reco_jet_eta;
  std::vector<float>* p_reco_jet_phi = &reco_jet_phi;
  std::vector<float>* p_reco_jet_e = &reco_jet_e;
  std::vector<float>* p_reco_jet_jvt = &reco_jet_jvt;
  std::vector<char>* p_reco_jet_b = &reco_jet_b;

  std::vector<float>* p_failJvt_reco_jet_pt = &failJvt_reco_jet_pt;
  std::vector<float>* p_failJvt_reco_jet_eta = &failJvt_reco_jet_eta;
  std::vector<float>* p_failJvt_reco_jet_phi = &failJvt_reco_jet_phi;
  std::vector<float>* p_failJvt_reco_jet_e = &failJvt_reco_jet_e;
  std::vector<float>* p_failJvt_reco_jet_jvt = &failJvt_reco_jet_jvt;

  std::vector<float>* p_truth_mu_pt = &truth_mu_pt;
  std::vector<float>* p_truth_mu_eta = &truth_mu_eta;
  std::vector<float>* p_truth_mu_phi = &truth_mu_phi;
  std::vector<float>* p_truth_mu_e = &truth_mu_e;

  std::vector<float>* p_truth_el_pt = &truth_el_pt;
  std::vector<float>* p_truth_el_eta = &truth_el_eta;
  std::vector<float>* p_truth_el_phi = &truth_el_phi;
  std::vector<float>* p_truth_el_e = &truth_el_e;

  std::vector<float>* p_truth_tau_pt = &truth_tau_pt;
  std::vector<float>* p_truth_tau_eta = &truth_tau_eta;
  std::vector<float>* p_truth_tau_phi = &truth_tau_phi;
  std::vector<float>* p_truth_tau_e = &truth_tau_e;

  for (int ifile = 1; ifile<=450; ifile++) {

    if (debug) {
      if (ifile == 2) break;
    }

    if (!debug2) std::cout << "Reading file " << ifile << std::endl;
    TFile* newFile = TFile::Open(TString::Format("/eos/user/r/rcoelhol/jet_studies/AnalysisTop/output_new_%d.root", ifile));  
    TTree* newReco = (TTree*) newFile->Get("nominal");
    TTree* newTruth = (TTree*) newFile->Get("particleLevel");
    newTruth->BuildIndex("eventNumber");
      
    newTruth->SetBranchAddress("jet_pt", &p_truth_jet_pt);
    newTruth->SetBranchAddress("jet_eta", &p_truth_jet_eta);
    newTruth->SetBranchAddress("jet_phi", &p_truth_jet_phi);
    newTruth->SetBranchAddress("jet_e", &p_truth_jet_e);
    newTruth->SetBranchAddress("jet_nGhosts_bHadron", &p_truth_jet_b);
    newTruth->SetBranchAddress("jet_nGhosts_cHadron", &p_truth_jet_c);
    newTruth->SetBranchAddress("mu_pt", &p_truth_mu_pt);
    newTruth->SetBranchAddress("mu_eta", &p_truth_mu_eta);
    newTruth->SetBranchAddress("mu_phi", &p_truth_mu_phi);
    newTruth->SetBranchAddress("mu_e", &p_truth_mu_e);
    newTruth->SetBranchAddress("el_pt", &p_truth_el_pt);
    newTruth->SetBranchAddress("el_eta", &p_truth_el_eta);
    newTruth->SetBranchAddress("el_phi", &p_truth_el_phi);
    newTruth->SetBranchAddress("el_e", &p_truth_el_e);
    newTruth->SetBranchAddress("tau_pt", &p_truth_tau_pt);
    newTruth->SetBranchAddress("tau_eta", &p_truth_tau_eta);
    newTruth->SetBranchAddress("tau_phi", &p_truth_tau_phi);
    newTruth->SetBranchAddress("tau_e", &p_truth_tau_e);

    newReco->SetBranchAddress("jet_pt", &p_reco_jet_pt);
    newReco->SetBranchAddress("jet_eta", &p_reco_jet_eta);
    newReco->SetBranchAddress("jet_phi", &p_reco_jet_phi);
    newReco->SetBranchAddress("jet_e", &p_reco_jet_e);
    newReco->SetBranchAddress("jet_jvt", &p_reco_jet_jvt);
    newReco->SetBranchAddress("jet_isbtagged_DL1r_77", &p_reco_jet_b);
    newReco->SetBranchAddress("failJvt_jet_pt", &p_failJvt_reco_jet_pt);
    newReco->SetBranchAddress("failJvt_jet_eta", &p_failJvt_reco_jet_eta);
    newReco->SetBranchAddress("failJvt_jet_phi", &p_failJvt_reco_jet_phi);
    newReco->SetBranchAddress("failJvt_jet_e", &p_failJvt_reco_jet_e);
    newReco->SetBranchAddress("failJvt_jet_jvt", &p_failJvt_reco_jet_jvt);
    newReco->SetBranchAddress("weight_mc", &weight_mc);
    newReco->SetBranchAddress("weight_pileup", &weight_pileup);
    newReco->SetBranchAddress("eventNumber", &eventNumber);
    newReco->SetBranchAddress("mu_actual", &mu_actual);
  
    // The study will be made with 25 > pT(reco) > 15 GeV
    // and |eta(reco)| < 2.5

    for (size_t ievent=0; ievent < newReco->GetEntries(); ievent++) {
      if (!debug2) {
	if (ievent % 1000 == 0) std::cout << "New correction, event " << ievent << "/" << newReco->GetEntries() << std::endl;
      }
      int recoEvent = newReco->GetEntry(ievent);
      int truthEvent = newTruth->GetEntryWithIndex((Int_t) eventNumber);
      if (debug) {
	std::cout << "Got event number " << eventNumber << std::endl;
	std::cout << "Trying to match " << truthEvent << " " << recoEvent << std::endl;
      }
      if (truthEvent < 0 || recoEvent < 0) continue;
    
      if (debug) {
	std::cout << "Got new correction, event number " << eventNumber << std::endl;
	if (ievent == 20) break;
      }

      double weight = weight_mc * weight_pileup;

      std::vector<float> all_jet_pt(reco_jet_pt);
      std::vector<float> all_jet_eta(reco_jet_eta);
      std::vector<float> all_jet_phi(reco_jet_phi);
      std::vector<float> all_jet_e(reco_jet_e);
      std::vector<float> all_jet_jvt(reco_jet_jvt);

      all_jet_pt.insert(all_jet_pt.end(), failJvt_reco_jet_pt.begin(), failJvt_reco_jet_pt.end());
      all_jet_eta.insert(all_jet_eta.end(), failJvt_reco_jet_eta.begin(), failJvt_reco_jet_eta.end());
      all_jet_phi.insert(all_jet_phi.end(), failJvt_reco_jet_phi.begin(), failJvt_reco_jet_phi.end());
      all_jet_e.insert(all_jet_e.end(), failJvt_reco_jet_e.begin(), failJvt_reco_jet_e.end());
      all_jet_jvt.insert(all_jet_jvt.end(), failJvt_reco_jet_jvt.begin(), failJvt_reco_jet_jvt.end());

      if (debug) {
	std::cout << "JVT size " << reco_jet_jvt.size() << " " << failJvt_reco_jet_jvt.size() << " " << all_jet_jvt.size() << std::endl;
      }

      // First, find the closest truth jet
      for (size_t irecoJet = 0; irecoJet < all_jet_pt.size(); irecoJet++) {

	new_all_cutflow->Fill(ifile-0.5,weight);
	new_reco_all_pt->Fill(all_jet_pt[irecoJet]/GeV, weight);
	new_reco_all_pt->Fill(all_jet_eta[irecoJet], weight);
            
	// Defines the relevant phase-space for this study
	if (all_jet_pt[irecoJet] < 15*GeV) continue;
	if (all_jet_pt[irecoJet] > 25*GeV) continue;
	if (std::abs(all_jet_eta[irecoJet]) > 2.5) continue;
	if (useBTagOnly && irecoJet > reco_jet_pt.size()) continue;
	if (useBTagOnly && !reco_jet_b[irecoJet]) continue;

	TLorentzVector reco_jet(0,0,0,0);
	reco_jet.SetPtEtaPhiE(all_jet_pt[irecoJet]/GeV, all_jet_eta[irecoJet], all_jet_phi[irecoJet], all_jet_e[irecoJet]/GeV);

	new_selected_cutflow->Fill(ifile-0.5,weight);

	// Check if jet is isolated https://arxiv.org/pdf/1703.09665.pdf
	int nIsolation = 0;
	for (size_t itruthMu2 = 0; itruthMu2 < truth_mu_pt.size(); itruthMu2++) {
	  if (truth_mu_pt[itruthMu2] < 7.0*GeV) continue;
	  TLorentzVector truth_mu2(0,0,0,0);
	  truth_mu2.SetPtEtaPhiE(truth_mu_pt[itruthMu2]/GeV, truth_mu_eta[itruthMu2], truth_mu_phi[itruthMu2], truth_mu_e[itruthMu2]/GeV);
	  if (truth_mu2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;
	for (size_t itruthEl2 = 0; itruthEl2 < truth_el_pt.size(); itruthEl2++) {
	  if (truth_el_pt[itruthEl2] < 7.0*GeV) continue;
	  TLorentzVector truth_el2(0,0,0,0);
	  truth_el2.SetPtEtaPhiE(truth_el_pt[itruthEl2]/GeV, truth_el_eta[itruthEl2], truth_el_phi[itruthEl2], truth_el_e[itruthEl2]/GeV);
	  if (truth_el2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;     
	for (size_t itruthTau2 = 0; itruthTau2 < truth_tau_pt.size(); itruthTau2++) {
	  if (truth_tau_pt[itruthTau2] < 7.0*GeV) continue;
	  TLorentzVector truth_tau2(0,0,0,0);
	  truth_tau2.SetPtEtaPhiE(truth_tau_pt[itruthTau2]/GeV, truth_tau_eta[itruthTau2], truth_tau_phi[itruthTau2], truth_tau_e[itruthTau2]/GeV);
	  if (truth_tau2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;     
	for (size_t itruthJet2 = 0; itruthJet2 < truth_jet_pt.size(); itruthJet2++) {
	  if (truth_jet_pt[itruthJet2] < 7.0*GeV) continue;
	  TLorentzVector truth_jet2(0,0,0,0);
	  truth_jet2.SetPtEtaPhiE(truth_jet_pt[itruthJet2]/GeV, truth_jet_eta[itruthJet2], truth_jet_phi[itruthJet2], truth_jet_e[itruthJet2]/GeV);
	  if (truth_jet2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	  }
	}
	if (nIsolation > 1) continue;
	
	new_thiso_cutflow->Fill(ifile-0.5,weight);

	int rIsolation = 0;
	TLorentzVector reco_isolation(0,0,0,0);
	for (size_t irecoJet2 = 0; irecoJet2 < all_jet_pt.size(); irecoJet2++) {
	  if (irecoJet == irecoJet2) continue;
	  if (all_jet_pt[irecoJet2] < 7.0*GeV) continue;
	  TLorentzVector reco_jet2(0,0,0,0);
	  reco_jet2.SetPtEtaPhiE(all_jet_pt[irecoJet2]/GeV, all_jet_eta[irecoJet2], all_jet_phi[irecoJet2], all_jet_e[irecoJet2]/GeV);
	  if (reco_jet2.DeltaR(reco_jet) < 0.6) {
	    reco_isolation += reco_jet2;
	    rIsolation++;
	  }
	}
	if (nIsolation == 1) new_reco_iso_pt->Fill(reco_isolation.Pt(), weight);
	if (rIsolation > 0) continue;
     
	new_expiso_cutflow->Fill(ifile-0.5,weight);

	// Find the closest truth jet
	double minDist = 99999;
	double minJet = -1;
	for (size_t itruthJet = 0; itruthJet < truth_jet_pt.size(); itruthJet++) {
	
	  TLorentzVector truth_jet(0,0,0,0);
	  truth_jet.SetPtEtaPhiE(truth_jet_pt[itruthJet]/GeV, truth_jet_eta[itruthJet], truth_jet_phi[itruthJet], truth_jet_e[itruthJet]/GeV);
	  if (truth_jet.DeltaR(reco_jet) < minDist) {
	    minJet = itruthJet;
	    minDist = truth_jet.DeltaR(reco_jet);
	  }
	}

	TLorentzVector truth_jet(0,0,0,0);
	truth_jet.SetPtEtaPhiE(truth_jet_pt[minJet]/GeV, truth_jet_eta[minJet], truth_jet_phi[minJet], truth_jet_e[minJet]/GeV);
	if (debug) {
	  std::cout << "Matched " << minDist << " " << nIsolation << " " << weight << std::endl;
	  std::cout << reco_jet.Pt() << " " << reco_jet.Eta() << " " << reco_jet.Phi() << " " << reco_jet.E() << std::endl;
	  std::cout << truth_jet.Pt() << " " << truth_jet.Eta() << " " << truth_jet.Phi() << " " << truth_jet.E() << std::endl;
	  std::cout << "------------------" << std::endl;
	}


	if (reco_jet.Pt() < 20) new_reco_truth_dist_low->Fill(minDist, weight);
	else new_reco_truth_dist_high->Fill(minDist, weight);

	if (minDist < 0.3) { // HS jet
	  new_matched_cutflow->Fill(ifile-0.5,weight);

	  new_reco_pt->Fill(reco_jet.Pt(), weight);
	  new_reco_eta->Fill(reco_jet.Eta(), weight);
	  new_truth_pt->Fill(truth_jet.Pt(), weight);
	  new_truth_eta->Fill(truth_jet.Eta(), weight);

	  if (reco_jet.Pt() < 20) {
	    new_reco_response_low->Fill((reco_jet.Pt()-truth_jet.Pt())/truth_jet.Pt(), weight);
	    new_reco_eta_resolution_low->Fill((reco_jet.Eta()-truth_jet.Eta()), weight);
	    new_reco_phi_resolution_low->Fill(TVector2::Phi_mpi_pi(reco_jet.Phi()-truth_jet.Phi()), weight);
	  } else {
	    new_reco_response_high->Fill((reco_jet.Pt()-truth_jet.Pt())/truth_jet.Pt(), weight);
	    new_reco_eta_resolution_high->Fill((reco_jet.Eta()-truth_jet.Eta()), weight);
	    new_reco_phi_resolution_high->Fill(TVector2::Phi_mpi_pi(reco_jet.Phi()-truth_jet.Phi()), weight);
	  }
	  new_hs_jvt_efficiency->FillWeighted(all_jet_jvt[irecoJet] > 0.50, weight, reco_jet.Pt());
	  if (irecoJet < reco_jet_pt.size()) {
	    if (truth_jet_b[minJet] > 0) new_btag_b_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	    else if (truth_jet_c[minJet] > 0) new_btag_c_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	    else new_btag_light_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	  }
	  _response.setVal((reco_jet.Pt()/truth_jet.Pt()) - 1.);
	  _pt_reco.setVal(reco_jet.Pt());
	  _pt_truth.setVal(truth_jet.Pt());
	  _eta_reco.setVal(reco_jet.Eta());
	  _eta_truth.setVal(truth_jet.Eta());
	  _mu_actual.setVal(mu_actual);
	  _option.setVal(1);
	  _weight.setVal(weight);
	  _dataset.add(_hs_jet);       
	} else { // PU jet
	  new_pu_jvt_efficiency->FillWeighted(all_jet_jvt[irecoJet] > 0.50, weight, reco_jet.Pt());
	}
      }
    }
    newFile->Close();

    TFile* oldFile = TFile::Open(TString::Format("/eos/user/r/rcoelhol/jet_studies/AnalysisTop/output_old_%d.root", ifile));
    TTree* oldReco = (TTree*) oldFile->Get("nominal");
    TTree* oldTruth = (TTree*) oldFile->Get("particleLevel");
    oldTruth->BuildIndex("eventNumber");

    oldTruth->SetBranchAddress("jet_pt", &p_truth_jet_pt);
    oldTruth->SetBranchAddress("jet_eta", &p_truth_jet_eta);
    oldTruth->SetBranchAddress("jet_phi", &p_truth_jet_phi);
    oldTruth->SetBranchAddress("jet_e", &p_truth_jet_e);
    oldTruth->SetBranchAddress("jet_nGhosts_bHadron", &p_truth_jet_b);
    oldTruth->SetBranchAddress("jet_nGhosts_cHadron", &p_truth_jet_c);
    oldTruth->SetBranchAddress("mu_pt", &p_truth_mu_pt);
    oldTruth->SetBranchAddress("mu_eta", &p_truth_mu_eta);
    oldTruth->SetBranchAddress("mu_phi", &p_truth_mu_phi);
    oldTruth->SetBranchAddress("mu_e", &p_truth_mu_e);
    oldTruth->SetBranchAddress("el_pt", &p_truth_el_pt);
    oldTruth->SetBranchAddress("el_eta", &p_truth_el_eta);
    oldTruth->SetBranchAddress("el_phi", &p_truth_el_phi);
    oldTruth->SetBranchAddress("el_e", &p_truth_el_e);
    oldTruth->SetBranchAddress("tau_pt", &p_truth_tau_pt);
    oldTruth->SetBranchAddress("tau_eta", &p_truth_tau_eta);
    oldTruth->SetBranchAddress("tau_phi", &p_truth_tau_phi);
    oldTruth->SetBranchAddress("tau_e", &p_truth_tau_e);

    oldReco->SetBranchAddress("jet_pt", &p_reco_jet_pt);
    oldReco->SetBranchAddress("jet_eta", &p_reco_jet_eta);
    oldReco->SetBranchAddress("jet_phi", &p_reco_jet_phi);
    oldReco->SetBranchAddress("jet_e", &p_reco_jet_e);
    oldReco->SetBranchAddress("jet_jvt", &p_reco_jet_jvt);
    oldReco->SetBranchAddress("jet_isbtagged_DL1r_77", &p_reco_jet_b);
    oldReco->SetBranchAddress("failJvt_jet_pt", &p_failJvt_reco_jet_pt);
    oldReco->SetBranchAddress("failJvt_jet_eta", &p_failJvt_reco_jet_eta);
    oldReco->SetBranchAddress("failJvt_jet_phi", &p_failJvt_reco_jet_phi);
    oldReco->SetBranchAddress("failJvt_jet_e", &p_failJvt_reco_jet_e);
    oldReco->SetBranchAddress("failJvt_jet_jvt", &p_failJvt_reco_jet_jvt);
    oldReco->SetBranchAddress("weight_mc", &weight_mc);
    oldReco->SetBranchAddress("weight_pileup", &weight_pileup);
    oldReco->SetBranchAddress("eventNumber", &eventNumber);
    oldReco->SetBranchAddress("mu_actual", &mu_actual);

    for (size_t ievent=0; ievent < oldReco->GetEntries(); ievent++) {
      if (!debug2) {
	if (ievent % 1000 == 0) std::cout << "Old correction, event " << ievent << "/" << oldReco->GetEntries() << std::endl;
      }
      int recoEvent = oldReco->GetEntry(ievent);
      int truthEvent = oldTruth->GetEntryWithIndex((Int_t) eventNumber);
      if (truthEvent < 0 || recoEvent < 0) continue;
    
      if (debug) {
	std::cout << "Got old correction, event number " << eventNumber << std::endl;
	if (ievent == 20) break;
      }

      double weight = weight_mc * weight_pileup;

      std::vector<float> all_jet_pt(reco_jet_pt);
      std::vector<float> all_jet_eta(reco_jet_eta);
      std::vector<float> all_jet_phi(reco_jet_phi);
      std::vector<float> all_jet_e(reco_jet_e);
      std::vector<float> all_jet_jvt(reco_jet_jvt);

      all_jet_pt.insert(all_jet_pt.end(), failJvt_reco_jet_pt.begin(), failJvt_reco_jet_pt.end());
      all_jet_eta.insert(all_jet_eta.end(), failJvt_reco_jet_eta.begin(), failJvt_reco_jet_eta.end());
      all_jet_phi.insert(all_jet_phi.end(), failJvt_reco_jet_phi.begin(), failJvt_reco_jet_phi.end());
      all_jet_e.insert(all_jet_e.end(), failJvt_reco_jet_e.begin(), failJvt_reco_jet_e.end());
      all_jet_jvt.insert(all_jet_jvt.end(), failJvt_reco_jet_jvt.begin(), failJvt_reco_jet_jvt.end());

      // First, find the closest truth jet
      for (size_t irecoJet = 0; irecoJet < all_jet_pt.size(); irecoJet++) {
      
      	old_all_cutflow->Fill(ifile-0.5,weight);
	old_reco_all_pt->Fill(all_jet_pt[irecoJet]/GeV, weight);
	old_reco_all_eta->Fill(all_jet_eta[irecoJet], weight);

	// Defines the relevant phase-space for this study
	if (all_jet_pt[irecoJet] < 15*GeV) continue;
	if (all_jet_pt[irecoJet] > 25*GeV) continue;
	if (std::abs(all_jet_eta[irecoJet]) > 2.5) continue;
	if (useBTagOnly && irecoJet > reco_jet_pt.size()) continue;
	if (useBTagOnly && !reco_jet_b[irecoJet]) continue;

	TLorentzVector reco_jet(0,0,0,0);
	reco_jet.SetPtEtaPhiE(all_jet_pt[irecoJet]/GeV, all_jet_eta[irecoJet], all_jet_phi[irecoJet], all_jet_e[irecoJet]/GeV);

	old_selected_cutflow->Fill(ifile-0.5,weight);

	// Check if jet is isolated https://arxiv.org/pdf/1703.09665.pdf
	int nIsolation = 0;
	for (size_t itruthMu2 = 0; itruthMu2 < truth_mu_pt.size(); itruthMu2++) {
	  if (truth_mu_pt[itruthMu2] < 7.0*GeV) continue;
	  TLorentzVector truth_mu2(0,0,0,0);
	  truth_mu2.SetPtEtaPhiE(truth_mu_pt[itruthMu2]/GeV, truth_mu_eta[itruthMu2], truth_mu_phi[itruthMu2], truth_mu_e[itruthMu2]/GeV);
	  if (truth_mu2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;
	for (size_t itruthEl2 = 0; itruthEl2 < truth_el_pt.size(); itruthEl2++) {
	  if (truth_el_pt[itruthEl2] < 7.0*GeV) continue;
	  TLorentzVector truth_el2(0,0,0,0);
	  truth_el2.SetPtEtaPhiE(truth_el_pt[itruthEl2]/GeV, truth_el_eta[itruthEl2], truth_el_phi[itruthEl2], truth_el_e[itruthEl2]/GeV);
	  if (truth_el2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;     
	for (size_t itruthTau2 = 0; itruthTau2 < truth_tau_pt.size(); itruthTau2++) {
	  if (truth_tau_pt[itruthTau2] < 7.0*GeV) continue;
	  TLorentzVector truth_tau2(0,0,0,0);
	  truth_tau2.SetPtEtaPhiE(truth_tau_pt[itruthTau2]/GeV, truth_tau_eta[itruthTau2], truth_tau_phi[itruthTau2], truth_tau_e[itruthTau2]/GeV);
	  if (truth_tau2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	    break;
	  }
	}
	if (nIsolation > 0) continue;     
	for (size_t itruthJet2 = 0; itruthJet2 < truth_jet_pt.size(); itruthJet2++) {
	  if (truth_jet_pt[itruthJet2] < 7.0*GeV) continue;
	  TLorentzVector truth_jet2(0,0,0,0);
	  truth_jet2.SetPtEtaPhiE(truth_jet_pt[itruthJet2]/GeV, truth_jet_eta[itruthJet2], truth_jet_phi[itruthJet2], truth_jet_e[itruthJet2]/GeV);
	  if (truth_jet2.DeltaR(reco_jet) < 1.0) {
	    nIsolation++;
	  }
	}
	if (nIsolation > 1) continue;

	old_thiso_cutflow->Fill(ifile-0.5,weight);

	int rIsolation = 0;
	TLorentzVector reco_isolation(0,0,0,0);
	for (size_t irecoJet2 = 0; irecoJet2 < all_jet_pt.size(); irecoJet2++) {
	  if (irecoJet == irecoJet2) continue;
	  if (all_jet_pt[irecoJet2] < 7.0*GeV) continue;
	  TLorentzVector reco_jet2(0,0,0,0);
	  reco_jet2.SetPtEtaPhiE(all_jet_pt[irecoJet2]/GeV, all_jet_eta[irecoJet2], all_jet_phi[irecoJet2], all_jet_e[irecoJet2]/GeV);
	  if (reco_jet2.DeltaR(reco_jet) < 0.6) {
	    reco_isolation += reco_jet2;
	    rIsolation++;
	  }
	}
	if (nIsolation == 1) old_reco_iso_pt->Fill(reco_isolation.Pt(), weight);
	if (rIsolation > 0) continue;

	old_expiso_cutflow->Fill(ifile-0.5,weight);
     
	// Find the closest truth jet
	double minDist = 99999;
	double minJet = -1;
	for (size_t itruthJet = 0; itruthJet < truth_jet_pt.size(); itruthJet++) {
	
	  TLorentzVector truth_jet(0,0,0,0);
	  truth_jet.SetPtEtaPhiE(truth_jet_pt[itruthJet]/GeV, truth_jet_eta[itruthJet], truth_jet_phi[itruthJet], truth_jet_e[itruthJet]/GeV);
	  if (truth_jet.DeltaR(reco_jet) < minDist) {
	    minJet = itruthJet;
	    minDist = truth_jet.DeltaR(reco_jet);
	  }
	}

	TLorentzVector truth_jet(0,0,0,0);
	truth_jet.SetPtEtaPhiE(truth_jet_pt[minJet]/GeV, truth_jet_eta[minJet], truth_jet_phi[minJet], truth_jet_e[minJet]/GeV);
	if (debug) {
	  std::cout << "Matched " << minDist << " " << nIsolation << " " << weight << std::endl;
	  std::cout << reco_jet.Pt() << " " << reco_jet.Eta() << " " << reco_jet.Phi() << " " << reco_jet.E() << std::endl;
	  std::cout << truth_jet.Pt() << " " << truth_jet.Eta() << " " << truth_jet.Phi() << " " << truth_jet.E() << std::endl;
	  std::cout << "------------------" << std::endl;
	}

	if (reco_jet.Pt() < 20) old_reco_truth_dist_low->Fill(minDist, weight);
	else old_reco_truth_dist_high->Fill(minDist, weight);

	if (minDist < 0.3) { // HS jet
	  old_matched_cutflow->Fill(ifile-0.5,weight);

	  old_reco_pt->Fill(reco_jet.Pt(), weight);
	  old_reco_eta->Fill(reco_jet.Eta(), weight);
	  old_truth_pt->Fill(truth_jet.Pt(), weight);
	  old_truth_eta->Fill(truth_jet.Eta(), weight);

	  if (reco_jet.Pt() < 20) {
	    old_reco_response_low->Fill((reco_jet.Pt()-truth_jet.Pt())/truth_jet.Pt(), weight);
	    old_reco_eta_resolution_low->Fill((reco_jet.Eta()-truth_jet.Eta()), weight);
	    old_reco_phi_resolution_low->Fill(TVector2::Phi_mpi_pi(reco_jet.Phi()-truth_jet.Phi()), weight);
	  } else {
	    old_reco_response_high->Fill((reco_jet.Pt()-truth_jet.Pt())/truth_jet.Pt(), weight);
	    old_reco_eta_resolution_high->Fill((reco_jet.Eta()-truth_jet.Eta()), weight);
	    old_reco_phi_resolution_high->Fill(TVector2::Phi_mpi_pi(reco_jet.Phi()-truth_jet.Phi()), weight);
	  }
	  old_hs_jvt_efficiency->FillWeighted(all_jet_jvt[irecoJet] > 0.50, weight, reco_jet.Pt());
	  if (irecoJet < reco_jet_pt.size()) {
	    if (truth_jet_b[minJet] > 0) old_btag_b_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	    else if (truth_jet_c[minJet] > 0) old_btag_c_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	    else old_btag_light_efficiency->FillWeighted(reco_jet_b[irecoJet], weight, reco_jet.Pt());
	  }
	  _response.setVal((reco_jet.Pt()/truth_jet.Pt()) - 1.);
	  _pt_reco.setVal(reco_jet.Pt());
	  _pt_truth.setVal(truth_jet.Pt());
	  _eta_reco.setVal(reco_jet.Eta());
	  _eta_truth.setVal(truth_jet.Eta());
	  _mu_actual.setVal(mu_actual);
	  _option.setVal(0);
	  _weight.setVal(weight);
	  _dataset.add(_hs_jet);       
	} else { // PU jet
	  old_pu_jvt_efficiency->FillWeighted(all_jet_jvt[irecoJet] > 0.50, weight, reco_jet.Pt());
	}
      }
    }
    oldFile->Close();
    if (debug2) {
      std::cout << "After file " << ifile << ", new correction: " << new_reco_pt->Integral() << "(" << new_reco_pt->GetEntries() << ")" << ", old correction: " << old_reco_pt->Integral() << "(" << old_reco_pt->GetEntries() << ")" << std::endl;
      std::cout << "After file " << ifile << ", new correction iso: " << new_reco_iso_pt->Integral() << " " << new_reco_iso_pt->Integral(0, 101) << " (" << new_reco_iso_pt->GetEntries() << ")" << ", old correction: iso " << old_reco_iso_pt->Integral() << " " << old_reco_iso_pt->Integral(0, 101) << " (" << old_reco_iso_pt->GetEntries() << ")" << std::endl;
    }
  }

  PutOverflowInLastBin(new_reco_pt); PutUnderflowInFirstBin(new_reco_pt);
  PutOverflowInLastBin(old_reco_pt); PutUnderflowInFirstBin(old_reco_pt);

  PutOverflowInLastBin(new_reco_iso_pt); PutUnderflowInFirstBin(new_reco_iso_pt);
  PutOverflowInLastBin(old_reco_iso_pt); PutUnderflowInFirstBin(old_reco_iso_pt);

  PutOverflowInLastBin(new_truth_pt); PutUnderflowInFirstBin(new_truth_pt);
  PutOverflowInLastBin(old_truth_pt); PutUnderflowInFirstBin(old_truth_pt);

  PutOverflowInLastBin(new_reco_eta); PutUnderflowInFirstBin(new_reco_eta);
  PutOverflowInLastBin(old_reco_eta); PutUnderflowInFirstBin(old_reco_eta);

  PutOverflowInLastBin(new_truth_eta); PutUnderflowInFirstBin(new_truth_eta);
  PutOverflowInLastBin(old_truth_eta); PutUnderflowInFirstBin(old_truth_eta);

  PutOverflowInLastBin(new_reco_truth_dist_low); PutUnderflowInFirstBin(new_reco_truth_dist_low);
  PutOverflowInLastBin(old_reco_truth_dist_low); PutUnderflowInFirstBin(old_reco_truth_dist_low);

  PutOverflowInLastBin(new_reco_truth_dist_high); PutUnderflowInFirstBin(new_reco_truth_dist_high);
  PutOverflowInLastBin(old_reco_truth_dist_high); PutUnderflowInFirstBin(old_reco_truth_dist_high);
 
  PutOverflowInLastBin(new_reco_response_low); PutUnderflowInFirstBin(new_reco_response_low);
  PutOverflowInLastBin(old_reco_response_low); PutUnderflowInFirstBin(old_reco_response_low);

  PutOverflowInLastBin(new_reco_response_high); PutUnderflowInFirstBin(new_reco_response_high);
  PutOverflowInLastBin(old_reco_response_high); PutUnderflowInFirstBin(old_reco_response_high);

  PutOverflowInLastBin(new_reco_eta_resolution_low); PutUnderflowInFirstBin(new_reco_eta_resolution_low);
  PutOverflowInLastBin(old_reco_eta_resolution_low); PutUnderflowInFirstBin(old_reco_eta_resolution_low);

  PutOverflowInLastBin(new_reco_phi_resolution_high); PutUnderflowInFirstBin(new_reco_phi_resolution_high);
  PutOverflowInLastBin(old_reco_phi_resolution_high); PutUnderflowInFirstBin(old_reco_phi_resolution_high);

  RooRealVar _mean("_mean", "_mean", -5., 5.);
  RooRealVar _sigma("_sigma", "_sigma", 0., 10.);
  RooRealVar _n1("_n1", "_n1", 1.01, 5000);
  RooRealVar _alpha1("_alpha1", "_alpha1", 0.5, 10);
  RooRealVar _n2("_n2", "_n2", 1.01, 5000);
  RooRealVar _alpha2("_alpha2", "_alpha2", 0.5, 10);
  RooDoubleCBFast _cbfunc("_cbfunc", "_cbfunc", _response, _mean, _sigma, _alpha1, _n1, _alpha2, _n2);
  RooGaussian _gaus("_gaus", "_gaus", _response, _mean, _sigma);

  gROOT->SetStyle("ATLAS");
  TCanvas* c1 = new TCanvas();

  double resolution;

  TGraphErrors* old_jet_scale = new TGraphErrors();
  TGraphErrors* old_jet_resolution = new TGraphErrors();  
  TGraphErrors* new_jet_scale = new TGraphErrors();
  TGraphErrors* new_jet_resolution = new TGraphErrors();  

  for (int i=0; i<ntot; i++) {
    double val = (25.-15.)/((double) ntot);

    resolution = doFit(_dataset, _response, _pt_reco, _option, TString::Format("_option==0 && _pt_reco > %.3f && _pt_reco < %.3f", 15.+val*((double) i), 15.+val*(((double) i)+1.)), 
		       _mean, _sigma, _alpha1, _alpha2, _n1, _n2, 
		       _gaus, _cbfunc, -1.0, 1.0,
		       c1, TString::Format("old_response_low_fit_%d.png", i));
    old_jet_scale->SetPoint(i, 15+val*(i+0.5), _mean.getVal()+1.);
    old_jet_scale->SetPointError(i, val*0.5, _mean.getError());
    old_jet_resolution->SetPoint(i, 15+val*(i+0.5), resolution);
    old_jet_resolution->SetPointError(i, val*0.5, resolution*_sigma.getError()/_sigma.getVal());
    resolution = doFit(_dataset, _response, _pt_reco, _option, TString::Format("_option==1 && _pt_reco > %.3f && _pt_reco < %.3f", 15.+val*((double) i), 15.+val*(((double) i)+1.)), 
		       _mean, _sigma, _alpha1, _alpha2, _n1, _n2, 
		       _gaus, _cbfunc, -1.0, 1.0,
		       c1, TString::Format("new_response_low_fit_%d.png", i));
    new_jet_scale->SetPoint(i, 15+val*(i+0.5), _mean.getVal()+1.);
    new_jet_scale->SetPointError(i, val*0.5, _mean.getError());
    new_jet_resolution->SetPoint(i, 15+val*(i+0.5), resolution);
    new_jet_resolution->SetPointError(i, val*0.5, resolution*_sigma.getError()/_sigma.getVal());
  }  

  TGraphErrors* old_jet_scale_pu = new TGraphErrors();
  TGraphErrors* old_jet_resolution_pu = new TGraphErrors();  
  TGraphErrors* new_jet_scale_pu = new TGraphErrors();
  TGraphErrors* new_jet_resolution_pu = new TGraphErrors();  

  for (int i=0; i<ntot; i++) {
    double val = (80.-0.)/((double) ntot);
    
    resolution = doFit(_dataset, _response, _pt_reco, _option, TString::Format("_option==0 && _mu_actual > %.3f && _mu_actual < %.3f && _pt_reco > 15 && _pt_reco < 20", val*((double) i), val*(((double) i)+1.)), 
		       _mean, _sigma, _alpha1, _alpha2, _n1, _n2, 
		       _gaus, _cbfunc, -1.0, 1.0,
		       c1, TString::Format("old_response_low_pu_fit_%d.png", i));
    old_jet_scale_pu->SetPoint(i, val*(i+0.5), _mean.getVal()+1.);
    old_jet_scale_pu->SetPointError(i, val*0.5, _mean.getError());
    old_jet_resolution_pu->SetPoint(i, val*(i+0.5), resolution);
    old_jet_resolution_pu->SetPointError(i, val*0.5, resolution*_sigma.getError()/_sigma.getVal());
    resolution = doFit(_dataset, _response, _pt_reco, _option, TString::Format("_option==1 && _mu_actual > %.3f && _mu_actual < %.3f && _pt_reco > 15 && _pt_reco < 20", val*((double) i), val*(((double) i)+1.)), 
		       _mean, _sigma, _alpha1, _alpha2, _n1, _n2, 
		       _gaus, _cbfunc, -1.0, 1.0,
		       c1, TString::Format("new_response_low_pu_fit_%d.png", i));
    new_jet_scale_pu->SetPoint(i, val*(i+0.5), _mean.getVal()+1.);
    new_jet_scale_pu->SetPointError(i, val*0.5, _mean.getError());
    new_jet_resolution_pu->SetPoint(i, val*(i+0.5), resolution);
    new_jet_resolution_pu->SetPointError(i, val*0.5, resolution*_sigma.getError()/_sigma.getVal());
  }  

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.85);
  leg->SetBorderSize(0);

  new_reco_truth_dist_low->SetLineColor(kBlue);
  old_reco_truth_dist_low->SetLineColor(kRed);
  new_reco_truth_dist_low->SetTitle("Truth matching 15 < p_{T} < 20 GeV; #Delta R(reco, truth); Entries");
  c1->SetLogy(true);
  new_reco_truth_dist_low->Draw("hist");
  old_reco_truth_dist_low->Draw("hist,same");  
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  // Do that only once
  leg->AddEntry(old_reco_truth_dist_low, "Option 0", "L");
  leg->AddEntry(new_reco_truth_dist_low, "Option 1", "L");
  //
  leg->Draw();
  c1->Print("reco_truth_dist_low.png");
  c1->Print("reco_truth_dist_low.pdf");
  c1->SetLogy(false);

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
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_scale.png");
  c1->Print("jet_scale.pdf");

  TMultiGraph* jet_resolution = new TMultiGraph();
  old_jet_resolution->SetLineColor(kRed);
  new_jet_resolution->SetLineColor(kBlue);
  jet_resolution->Add(old_jet_resolution, "L");
  jet_resolution->Add(new_jet_resolution, "L");
  jet_resolution->Draw("A");
  jet_resolution->GetXaxis()->SetTitle("Jet reco p_{T} [GeV]");
  jet_resolution->GetYaxis()->SetTitle("Jet energy resolution");
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_resolution.png");
  c1->Print("jet_resolution.pdf");

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
  gPad->Update();
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("jet_scale_pu.png");
  c1->Print("jet_scale_pu.pdf");

  new_reco_truth_dist_high->SetLineColor(kBlue);
  old_reco_truth_dist_high->SetLineColor(kRed);
  new_reco_truth_dist_high->SetTitle("Truth matching 20 < p_{T} < 25 GeV; #Delta R(reco, truth); Entries");
  c1->SetLogy(true);
  new_reco_truth_dist_high->Draw("hist");
  old_reco_truth_dist_high->Draw("hist,same");  
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  leg->Draw();
  c1->Print("reco_truth_dist_high.png");
  c1->Print("reco_truth_dist_high.pdf");
  c1->SetLogy(false);

  new_reco_response_low->SetLineColor(kBlue);
  old_reco_response_low->SetLineColor(kRed);
  new_reco_response_low->SetTitle("Jet response 15 < p_{T} < 20 GeV; p_{T}^{reco}/p_{T}^{truth} - 1; Entries");
  new_reco_response_low->SetNormFactor(1);
  old_reco_response_low->SetNormFactor(1);
  new_reco_response_low->Draw("hist");
  old_reco_response_low->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  leg->Draw();
  c1->Print("reco_response_low.png");
  c1->Print("reco_response_low.pdf");

  RooDataHist _old_response_low("_old_response_low", "_old_response_low", RooArgSet(_response), old_reco_response_low);
  RooDataHist _new_response_low("_new_response_low", "_new_response_low", RooArgSet(_response), new_reco_response_low);
  RooDoubleCBFast _cbfunc_old_response_low("_cbfunc_old_response_low", "_cbfunc_old_response_low", _response, _mean, _sigma, _alpha1, _n1, _alpha2, _n2);
  RooDoubleCBFast _cbfunc_new_response_low("_cbfunc_new_response_low", "_cbfunc_new_response_low", _response, _mean, _sigma, _alpha1, _n1, _alpha2, _n2);

  double _scale_old, _scale_old_err, _resolution_old;
  double _scale_new, _scale_new_err, _resolution_new;

  _mean.setVal(_old_response_low.mean(_response)-0.1);
  _sigma.setVal(_old_response_low.sigma(_response));
  _alpha1.setVal(1);
  _alpha2.setVal(2);
  _n1.setVal(20);
  _n2.setVal(5);

  _mean.setConstant(false);
  _sigma.setConstant(false);
  _alpha1.setConstant(false);
  _alpha2.setConstant(false);
  _n1.setConstant(true);
  _n2.setConstant(true);
  _gaus.fitTo(_old_response_low,RooFit::Range(-0.2,0.2), RooFit::NumCPU(4));
  _cbfunc_old_response_low.fitTo(_old_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _n1.setConstant(false);
  _n2.setConstant(false);
  _cbfunc_old_response_low.fitTo(_old_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _cbfunc_old_response_low.fitTo(_old_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _scale_old = _mean.getVal();
  _scale_old_err = _mean.getError();
  _resolution_old = calculateResolution(_sigma.getVal(), _alpha1.getVal(), _alpha2.getVal(), _n1.getVal(), _n2.getVal(), 0.68);

  _mean.setConstant(false);
  _sigma.setConstant(false);
  _alpha1.setConstant(false);
  _alpha2.setConstant(false);
  _n1.setConstant(true);
  _n2.setConstant(true);
  _gaus.fitTo(_new_response_low,RooFit::Range(-0.2,0.2), RooFit::NumCPU(4));
  _cbfunc_new_response_low.fitTo(_new_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _n1.setConstant(false);
  _n2.setConstant(false);
  _cbfunc_new_response_low.fitTo(_new_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _cbfunc_new_response_low.fitTo(_new_response_low,RooFit::Range(-1.,1.), RooFit::NumCPU(4));
  _scale_new = _mean.getVal();
  _scale_new_err = _mean.getError();
  _resolution_new = calculateResolution(_sigma.getVal(), _alpha1.getVal(), _alpha2.getVal(), _n1.getVal(), _n2.getVal(), 0.68);

  RooPlot* _new_response_low_frame = _response.frame(RooFit::Range(-1,1)) ; 
  _new_response_low.plotOn(_new_response_low_frame, RooFit::MarkerColor(kBlue), RooFit::LineColor(kBlue));
  _cbfunc_new_response_low.plotOn(_new_response_low_frame, RooFit::LineColor(kBlue));
  _new_response_low_frame->Draw();  
  RooPlot* _old_response_low_frame = _response.frame(RooFit::Range(-1,1)) ; 
  _old_response_low.plotOn(_old_response_low_frame, RooFit::MarkerColor(kRed), RooFit::LineColor(kRed));
  _cbfunc_old_response_low.plotOn(_old_response_low_frame, RooFit::LineColor(kRed));
  _old_response_low_frame->Draw("same");  
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  leg->Draw();
  myText(0.2, 0.64, kBlack, "Option 0");
  myText(0.2, 0.56, kBlack, TString::Format("scale = %.3f #pm %.3f", _scale_old+1, _scale_old_err));
  myText(0.2, 0.48, kBlack, TString::Format("resolution = %.3f", _resolution_old));
  myText(0.2, 0.40, kBlack, "Option 1");
  myText(0.2, 0.32, kBlack, TString::Format("scale = %.3f #pm %.3f", _scale_new+1, _scale_new_err));
  myText(0.2, 0.24, kBlack, TString::Format("resolution = %.3f", _resolution_new));

  c1->Print("reco_response_low_fit.png");
  c1->Print("reco_response_low_fit.pdf");

  new_reco_response_high->SetLineColor(kBlue);
  old_reco_response_high->SetLineColor(kRed);
  new_reco_response_high->SetTitle("Jet response 20 < p_{T} < 25 GeV; p_{T}^{reco}/p_{T}^{truth} - 1; Entries");
  new_reco_response_high->SetNormFactor(1);
  old_reco_response_high->SetNormFactor(1);
  new_reco_response_high->Draw("hist");
  old_reco_response_high->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  leg->Draw();
  c1->Print("reco_response_high.png");
  c1->Print("reco_response_high.pdf");

  new_reco_eta_resolution_low->SetLineColor(kBlue);
  old_reco_eta_resolution_low->SetLineColor(kRed);
  new_reco_eta_resolution_low->SetTitle("Jet eta resolution 15 < p_{T} < 20 GeV; #eta^{reco} - #eta^{truth}; Entries");
  new_reco_eta_resolution_low->SetNormFactor(1);
  old_reco_eta_resolution_low->SetNormFactor(1);
  new_reco_eta_resolution_low->Draw("hist");
  old_reco_eta_resolution_low->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  leg->Draw();
  c1->Print("reco_eta_resolution_low.png");
  c1->Print("reco_eta_resolution_low.pdf");

  new_reco_eta_resolution_high->SetLineColor(kBlue);
  old_reco_eta_resolution_high->SetLineColor(kRed);
  new_reco_eta_resolution_high->SetTitle("Jet eta resolution 20 < p_{T} < 25 GeV; #eta^{reco} - #eta^{truth}; Entries");
  new_reco_eta_resolution_high->SetNormFactor(1);
  old_reco_eta_resolution_high->SetNormFactor(1);
  new_reco_eta_resolution_high->Draw("hist");
  old_reco_eta_resolution_high->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  leg->Draw();
  c1->Print("reco_eta_resolution_high.png");
  c1->Print("reco_eta_resolution_high.pdf");

  new_reco_phi_resolution_low->SetLineColor(kBlue);
  old_reco_phi_resolution_low->SetLineColor(kRed);
  new_reco_phi_resolution_low->SetTitle("Jet phi resolution 15 < p_{T} < 20 GeV; #phi^{reco} - #phi^{truth}; Entries");
  new_reco_phi_resolution_low->SetNormFactor(1);
  old_reco_phi_resolution_low->SetNormFactor(1);
  new_reco_phi_resolution_low->Draw("hist");
  old_reco_phi_resolution_low->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  leg->Draw();
  c1->Print("reco_phi_resolution_low.png");
  c1->Print("reco_phi_resolution_low.pdf");

  new_reco_phi_resolution_high->SetLineColor(kBlue);
  old_reco_phi_resolution_high->SetLineColor(kRed);
  new_reco_phi_resolution_high->SetTitle("Jet phi resolution 20 < p_{T} < 25 GeV; #phi^{reco} - #phi^{truth}; Entries");
  new_reco_phi_resolution_high->SetNormFactor(1);
  old_reco_phi_resolution_high->SetNormFactor(1);
  new_reco_phi_resolution_high->Draw("hist");
  old_reco_phi_resolution_high->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  leg->Draw();
  c1->Print("reco_phi_resolution_high.png");
  c1->Print("reco_phi_resolution_high.pdf");

  new_reco_pt->SetMinimum(0);
  old_reco_pt->SetMinimum(0);
  new_reco_pt->SetMaximum(1.5*new_reco_pt->GetMaximum());
  old_reco_pt->SetMaximum(1.5*old_reco_pt->GetMaximum());
  new_reco_pt->SetLineColor(kBlue);
  old_reco_pt->SetLineColor(kRed);
  new_reco_pt->SetTitle("Jet p_{T}; p_{T}^{reco} [GeV]; Entries");
  new_reco_pt->Draw("hist");
  old_reco_pt->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("reco_pt.png");
  c1->Print("reco_pt.pdf");

  new_reco_all_pt->SetMinimum(0);
  old_reco_all_pt->SetMinimum(0);
  new_reco_all_pt->SetMaximum(1.5*new_reco_all_pt->GetMaximum());
  old_reco_all_pt->SetMaximum(1.5*old_reco_all_pt->GetMaximum());
  new_reco_all_pt->SetLineColor(kBlue);
  old_reco_all_pt->SetLineColor(kRed);
  new_reco_all_pt->SetTitle("Jet p_{T}; p_{T}^{reco} [GeV]; Entries");
  new_reco_all_pt->Draw("hist");
  old_reco_all_pt->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("reco_all_pt.png");
  c1->Print("reco_all_pt.pdf");

  new_reco_iso_pt->SetLineColor(kBlue);
  old_reco_iso_pt->SetLineColor(kRed);
  new_reco_iso_pt->SetNormFactor(1);
  old_reco_iso_pt->SetNormFactor(1);
  new_reco_iso_pt->SetTitle("Jet p_{T}; p_{T}^{iso} [GeV]; Entries");
  c1->SetLogy(true);
  new_reco_iso_pt->Draw("hist");
  old_reco_iso_pt->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("reco_iso_pt.png");
  c1->Print("reco_iso_pt.pdf");
  c1->SetLogy(false);

  new_truth_pt->SetMinimum(0);
  old_truth_pt->SetMinimum(0);
  new_truth_pt->SetLineColor(kBlue);
  old_truth_pt->SetLineColor(kRed);
  new_truth_pt->SetTitle("Jet p_{T}; p_{T}^{truth} [GeV]; Entries");
  new_truth_pt->Draw("hist");
  old_truth_pt->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("truth_pt.png");
  c1->Print("truth_pt.pdf");

  new_reco_eta->SetMinimum(0);
  old_reco_eta->SetMinimum(0);
  new_reco_eta->SetMaximum(1.5*new_reco_eta->GetMaximum());
  old_reco_eta->SetMaximum(1.5*old_reco_eta->GetMaximum());
  new_reco_eta->SetLineColor(kBlue);
  old_reco_eta->SetLineColor(kRed);
  new_reco_eta->SetTitle("Jet #eta; #eta^{reco} [GeV]; Entries");
  new_reco_eta->Draw("hist");
  old_reco_eta->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("reco_eta.png");
  c1->Print("reco_eta.pdf");

  new_reco_all_eta->SetMinimum(0);
  old_reco_all_eta->SetMinimum(0);
  new_reco_all_eta->SetMaximum(3*new_reco_all_eta->GetMaximum());
  old_reco_all_eta->SetMaximum(3*old_reco_all_eta->GetMaximum());
  new_reco_all_eta->SetLineColor(kBlue);
  old_reco_all_eta->SetLineColor(kRed);
  new_reco_all_eta->SetTitle("Jet p_{T}; p_{T}^{reco} [GeV]; Entries");
  new_reco_all_eta->Draw("hist");
  old_reco_all_eta->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("reco_all_eta.png");
  c1->Print("reco_all_eta.pdf");

  new_truth_eta->SetLineColor(kBlue);
  old_truth_eta->SetLineColor(kRed);
  new_truth_eta->SetTitle("Jet #eta; #eta^{truth} [GeV]; Entries");
  new_truth_eta->Draw("hist");
  old_truth_eta->Draw("hist,same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("truth_eta.png");
  c1->Print("truth_eta.pdf");

  new_hs_jvt_efficiency->SetLineColor(kBlue);
  old_hs_jvt_efficiency->SetLineColor(kRed);
  new_hs_jvt_efficiency->SetMarkerColor(kBlue);
  old_hs_jvt_efficiency->SetMarkerColor(kRed);
  new_hs_jvt_efficiency->SetMarkerStyle(20);
  old_hs_jvt_efficiency->SetMarkerStyle(20);
  new_hs_jvt_efficiency->SetTitle("JVT efficiency HS jets; Jet reco p_{T} [GeV]; JVT efficiency");
  TH1D* eff_hist = new TH1D("eff_hist", "", 1, 15, 25);
  eff_hist->SetTitle("JVT efficiency HS jets; Jet reco p_{T} [GeV]; JVT efficiency");
  eff_hist->GetYaxis()->SetRangeUser(0.6, 1.5);
  eff_hist->Draw();
  new_hs_jvt_efficiency->Draw("same");
  old_hs_jvt_efficiency->Draw("same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("hs_jvt_efficiency.png");
  c1->Print("hs_jvt_efficiency.pdf");

  new_pu_jvt_efficiency->SetLineColor(kBlue);
  old_pu_jvt_efficiency->SetLineColor(kRed);
  new_pu_jvt_efficiency->SetMarkerColor(kBlue);
  old_pu_jvt_efficiency->SetMarkerColor(kRed);
  new_pu_jvt_efficiency->SetMarkerStyle(20);
  old_pu_jvt_efficiency->SetMarkerStyle(20);
  new_pu_jvt_efficiency->SetTitle("JVT efficiency PU jets; Jet reco p_{T} [GeV]; JVT efficiency");
  eff_hist->SetTitle("JVT efficiency PU jets; Jet reco p_{T} [GeV]; JVT efficiency");
  eff_hist->GetYaxis()->SetRangeUser(0, 0.6);
  eff_hist->Draw();
  new_pu_jvt_efficiency->Draw("same");
  old_pu_jvt_efficiency->Draw("same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("pu_jvt_efficiency.png");
  c1->Print("pu_jvt_efficiency.pdf");

  new_btag_b_efficiency->SetLineColor(kBlue);
  old_btag_b_efficiency->SetLineColor(kRed);
  new_btag_b_efficiency->SetMarkerColor(kBlue);
  old_btag_b_efficiency->SetMarkerColor(kRed);
  new_btag_b_efficiency->SetMarkerStyle(20);
  old_btag_b_efficiency->SetMarkerStyle(20);
  new_btag_b_efficiency->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r b efficiency");
  eff_hist->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r efficiency");
  eff_hist->GetYaxis()->SetRangeUser(0.3, 1.6);
  eff_hist->Draw();
  new_btag_b_efficiency->Draw("same");
  old_btag_b_efficiency->Draw("same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("btag_b_efficiency.png");
  c1->Print("btag_b_efficiency.pdf");

  new_btag_c_efficiency->SetLineColor(kBlue);
  old_btag_c_efficiency->SetLineColor(kRed);
  new_btag_c_efficiency->SetMarkerColor(kBlue);
  old_btag_c_efficiency->SetMarkerColor(kRed);
  new_btag_c_efficiency->SetMarkerStyle(20);
  old_btag_c_efficiency->SetMarkerStyle(20);
  new_btag_c_efficiency->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r c efficiency");
  eff_hist->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r efficiency");
  eff_hist->GetYaxis()->SetRangeUser(0, 0.8);
  eff_hist->Draw();
  new_btag_c_efficiency->Draw("same");
  old_btag_c_efficiency->Draw("same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("btag_c_efficiency.png");
  c1->Print("btag_c_efficiency.pdf");

  new_btag_light_efficiency->SetLineColor(kBlue);
  old_btag_light_efficiency->SetLineColor(kRed);
  new_btag_light_efficiency->SetMarkerColor(kBlue);
  old_btag_light_efficiency->SetMarkerColor(kRed);
  new_btag_light_efficiency->SetMarkerStyle(20);
  old_btag_light_efficiency->SetMarkerStyle(20);
  new_btag_light_efficiency->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r light efficiency");
  eff_hist->SetTitle("DL1r efficiency (77% WP); Jet reco p_{T} [GeV]; DL1r efficiency");
  eff_hist->GetYaxis()->SetRangeUser(0, 0.04);
  eff_hist->Draw();
  new_btag_light_efficiency->Draw("same");
  old_btag_light_efficiency->Draw("same");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  leg->Draw();
  c1->Print("btag_light_efficiency.png");
  c1->Print("btag_light_efficiency.pdf");

  TFile* outputFile = TFile::Open("plots.root", "RECREATE");
  old_jet_scale->Write("old_jet_scale");
  old_jet_resolution->Write("old_jet_resolution");
  old_jet_scale_pu->Write("old_jet_scale_pu");
  old_jet_resolution_pu->Write("old_jet_resolution_pu");
  new_jet_scale->Write("new_jet_scale");
  new_jet_resolution->Write("new_jet_resolution");
  new_jet_scale_pu->Write("new_jet_scale_pu");
  new_jet_resolution_pu->Write("new_jet_resolution_pu");

  new_all_cutflow->Write();
  new_selected_cutflow->Write();
  new_thiso_cutflow->Write();
  new_expiso_cutflow->Write();
  new_matched_cutflow->Write();

  old_all_cutflow->Write();
  old_selected_cutflow->Write();
  old_thiso_cutflow->Write();
  old_expiso_cutflow->Write();
  old_matched_cutflow->Write();
  outputFile->Close();

}

