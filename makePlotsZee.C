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
#include "TVector2.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
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
#define doOnlyMuons false
#define doOnlyElectrons false

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
  
double doFit(RooDataSet& dataset, RooRealVar& response, RooRealVar& pt_jet, RooRealVar& weight, TString cut, 
	     RooRealVar& mean, RooRealVar& sigma, RooRealVar& alpha1, RooRealVar& alpha2, RooRealVar& n1, RooRealVar& n2, 
	     RooGaussian& gaus, RooDoubleCBFast& cbfunc, double xmin, double xmax,
	     TCanvas* c1, TString name) {
  
  RooDataSet* reduced_dataset = (RooDataSet*) dataset.reduce(RooArgSet(response,pt_jet,weight), cut.Data());
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
  
  TH1D* zee_mass_low = new TH1D("zee_mass_low", "zee_mass_low", 40, 81, 101);
  TH1D* zee_mass_high = new TH1D("zee_mass_high", "zee_mass_high", 40, 81, 101);
  TH1D* zee_pt_low = new TH1D("zee_pt_low", "zee_pt_low", 50, 0, 50);
  TH1D* zee_pt_high = new TH1D("zee_pt_high", "zee_pt_high", 50, 0, 50);
  TH1D* response_low = new TH1D("response_low", "response_low", 50, -1, 1);
  TH1D* response_high = new TH1D("response_high", "response_high", 50, -1, 1);

  RooRealVar _response("_response", "p_{T}^{#parallel}/p_{T}^{Z}-1", 0);
  RooRealVar _pt_jet("_pt_jet", "_pt_jet", 0);
  RooRealVar _pt_zee("_pt_zee", "_pt_zee", 0);
  RooRealVar _eta_jet("_eta_jet", "_eta_jet", 0);
  RooRealVar _eta_zee("_eta_zee", "_eta_zee", 0);
  RooRealVar _phi_jet("_phi_jet", "_phi_jet", 0);
  RooRealVar _phi_zee("_phi_zee", "_phi_zee", 0);
  RooRealVar _weight("_weight", "_weight", 0);
  RooArgSet _hs_jet("hs_jet");
  _hs_jet.add(_response);
  _hs_jet.add(_pt_jet);
  _hs_jet.add(_pt_zee);
  _hs_jet.add(_eta_jet);
  _hs_jet.add(_eta_zee);
  _hs_jet.add(_phi_jet);
  _hs_jet.add(_phi_zee);
  _hs_jet.add(_weight);
  RooDataSet _dataset("_dataset", "", _hs_jet, "_weight");

  int mm, ee;
  float weight_mc, weight_pileup, weight_leptonSF, weight_globalLeptonTriggerSF, weight_jvt;
  float weight, response;

  std::vector<float> jet_pt, jet_eta, jet_phi, jet_e, jet_jvt ;
  std::vector<char> jet_b;

  std::vector<float> mu_pt, mu_eta, mu_phi, mu_e;
  std::vector<float> el_pt, el_eta, el_phi, el_e;

  std::vector<float>* p_jet_pt = &jet_pt;
  std::vector<float>* p_jet_eta = &jet_eta;
  std::vector<float>* p_jet_phi = &jet_phi;
  std::vector<float>* p_jet_e = &jet_e;
  std::vector<float>* p_jet_jvt = &jet_jvt;
  std::vector<char>* p_jet_b = &jet_b;

  std::vector<float>* p_mu_pt = &mu_pt;
  std::vector<float>* p_mu_eta = &mu_eta;
  std::vector<float>* p_mu_phi = &mu_phi;
  std::vector<float>* p_mu_e = &mu_e;

  std::vector<float>* p_el_pt = &el_pt;
  std::vector<float>* p_el_eta = &el_eta;
  std::vector<float>* p_el_phi = &el_phi;
  std::vector<float>* p_el_e = &el_e;

  for (int ifile = 1; ifile<=110; ifile++) {

    if (!debug2) std::cout << "Reading file " << ifile << std::endl;
    TFile* zeeFile = TFile::Open(TString::Format("/eos/user/r/rcoelhol/jet_studies/AnalysisTop/Zee/output_old_%d.root", ifile));  
    if (!zeeFile || zeeFile->IsZombie()) continue;
    TTree* zeeReco = (TTree*) zeeFile->Get("nominal");

    zeeReco->SetBranchAddress("mm_data", &mm);
    zeeReco->SetBranchAddress("ee_data", &ee);

    zeeReco->SetBranchAddress("jet_pt", &p_jet_pt);
    zeeReco->SetBranchAddress("jet_eta", &p_jet_eta);
    zeeReco->SetBranchAddress("jet_phi", &p_jet_phi);
    zeeReco->SetBranchAddress("jet_e", &p_jet_e);
    zeeReco->SetBranchAddress("jet_jvt", &p_jet_jvt);
    zeeReco->SetBranchAddress("jet_isbtagged_DL1r_77", &p_jet_b);

    zeeReco->SetBranchAddress("mu_pt", &p_mu_pt);
    zeeReco->SetBranchAddress("mu_eta", &p_mu_eta);
    zeeReco->SetBranchAddress("mu_phi", &p_mu_phi);
    zeeReco->SetBranchAddress("mu_e", &p_mu_e);

    zeeReco->SetBranchAddress("el_pt", &p_el_pt);
    zeeReco->SetBranchAddress("el_eta", &p_el_eta);
    zeeReco->SetBranchAddress("el_phi", &p_el_phi);
    zeeReco->SetBranchAddress("el_e", &p_el_e);

    zeeReco->SetBranchAddress("weight_mc", &weight_mc);
    zeeReco->SetBranchAddress("weight_pileup", &weight_pileup);
    zeeReco->SetBranchAddress("weight_leptonSF", &weight_leptonSF);
    zeeReco->SetBranchAddress("weight_globalLeptonTriggerSF", &weight_globalLeptonTriggerSF);
  
    for (size_t ievent=0; ievent < zeeReco->GetEntries(); ievent++) {
      if (!debug2) {
	if (ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << zeeReco->GetEntries() << std::endl;
      }
      int recoEvent = zeeReco->GetEntry(ievent);

      TLorentzVector lep0(0,0,0,0);
      TLorentzVector lep1(0,0,0,0);
      if (doOnlyMuons && mm == 0) continue;
      if (doOnlyElectrons && ee == 0) continue;

      if (mm == 1) {
	lep0.SetPtEtaPhiE(mu_pt[0]/GeV, mu_eta[0], mu_phi[0], mu_e[0]/GeV);
	lep1.SetPtEtaPhiE(mu_pt[1]/GeV, mu_eta[1], mu_phi[1], mu_e[1]/GeV);
      } else if (ee == 1) {
	lep0.SetPtEtaPhiE(el_pt[0]/GeV, el_eta[0], el_phi[0], el_e[0]/GeV);
	lep1.SetPtEtaPhiE(el_pt[1]/GeV, el_eta[1], el_phi[1], el_e[1]/GeV);
      } else continue;

      TLorentzVector z = lep0+lep1;

      if (jet_pt[0]/GeV < 15. || jet_pt[0]/GeV > 25.) continue;
      if (std::abs(jet_eta[0]) > 2.5) continue;      
      if (jet_pt.size() > 1) continue;

      TLorentzVector jet(0,0,0,0);
      jet.SetPtEtaPhiE(jet_pt[0]/GeV, jet_eta[0], jet_phi[0], jet_e[0]/GeV);

      if (z.DeltaPhi(jet) < 2.9) continue;
    
      TVector2 ptz = z.Vect().XYvector();
      TVector2 ptj = jet.Vect().XYvector();    
      TVector2 ptj_ref = ptj.Proj(ptz);

      weight = weight_mc * weight_pileup * weight_leptonSF * weight_globalLeptonTriggerSF;
      response = ptj_ref.Mod()/ptz.Mod();

      if (jet.Pt() < 20) {
	zee_mass_low->Fill(z.M(), weight);
	zee_pt_low->Fill(z.Pt(), weight);
	response_low->Fill(response-1, weight);
      } else {
	zee_mass_high->Fill(z.M(), weight);
	zee_pt_high->Fill(z.Pt(), weight);
	response_high->Fill(response-1, weight);
      }

      _response.setVal(response - 1.);
      _pt_jet.setVal(jet.Pt());
      _eta_jet.setVal(jet.Eta());
      _phi_jet.setVal(jet.Phi());
      _pt_zee.setVal(z.Pt());
      _eta_zee.setVal(z.Eta());
      _phi_zee.setVal(z.Phi());
      _weight.setVal(weight);
      _dataset.add(_hs_jet);       
    }
    zeeFile->Close();
  }
  
  PutOverflowInLastBin(zee_mass_low); PutUnderflowInFirstBin(zee_mass_low);
  PutOverflowInLastBin(zee_pt_low); PutUnderflowInFirstBin(zee_pt_low);
  //  PutOverflowInLastBin(response_low); PutUnderflowInFirstBin(response_low);

  PutOverflowInLastBin(zee_mass_high); PutUnderflowInFirstBin(zee_mass_high);
  PutOverflowInLastBin(zee_pt_high); PutUnderflowInFirstBin(zee_pt_high);
  //  PutOverflowInLastBin(response_high); PutUnderflowInFirstBin(response_high);

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

  TGraphErrors* mc_scale = new TGraphErrors();
  TGraphErrors* mc_resolution = new TGraphErrors();  

  for (int i=0; i<ntot; i++) {
    double val = (25.-15.)/((double) ntot);
    double meanEstimate = -0.627476+2.00548e-02*(15+val*(i+0.5)) + 1.;
    _mean.setRange(0.8*meanEstimate-1., 1.2*meanEstimate-1.);
    resolution = doFit(_dataset, _response, _pt_jet, _weight, TString::Format("_pt_jet > %.3f && _pt_jet < %.3f", 15.+val*((double) i), 15.+val*(((double) i)+1.)), 
		       _mean, _sigma, _alpha1, _alpha2, _n1, _n2, 
		       _gaus, _cbfunc, -1.0, 1.0,
		       c1, TString::Format("zee_response_fit_%d.png", i));
    mc_scale->SetPoint(i, 15+val*(i+0.5), _mean.getVal());
    mc_scale->SetPointError(i, val*0.5, _mean.getError());
    mc_resolution->SetPoint(i, 15+val*(i+0.5), resolution);
    mc_resolution->SetPointError(i, val*0.5, resolution*_sigma.getError()/_sigma.getVal());
  }  
  
  RooDataSet* _reduced_dataset_1 = (RooDataSet*) _dataset.reduce(RooArgSet(_response,_pt_jet,_weight), "_pt_jet > 15.0 && _pt_jet < 17.5");
  RooDataSet* _reduced_dataset_2 = (RooDataSet*) _dataset.reduce(RooArgSet(_response,_pt_jet,_weight), "_pt_jet > 17.5 && _pt_jet < 20.0");
  RooDataSet* _reduced_dataset_3 = (RooDataSet*) _dataset.reduce(RooArgSet(_response,_pt_jet,_weight), "_pt_jet > 20.0 && _pt_jet < 22.5");
  RooDataSet* _reduced_dataset_4 = (RooDataSet*) _dataset.reduce(RooArgSet(_response,_pt_jet,_weight), "_pt_jet > 22.5 && _pt_jet < 25.0");

  TH1D* _binned_response_1 = (TH1D*) _reduced_dataset_1->createHistogram("binned_response_1", _response, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_response_2 = (TH1D*) _reduced_dataset_2->createHistogram("binned_response_2", _response, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_response_3 = (TH1D*) _reduced_dataset_3->createHistogram("binned_response_3", _response, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_response_4 = (TH1D*) _reduced_dataset_4->createHistogram("binned_response_4", _response, RooFit::Binning(12,-0.6,0.6));

  zee_mass_low->SetTitle("Z mass 15 < p_{T} < 20 GeV; Dilepton mass [GeV]; Entries");
  zee_mass_low->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  c1->Print("zee_mass_low.png");
  c1->Print("zee_mass_low.pdf");

  zee_pt_low->SetTitle("Z pt 15 < p_{T} < 20 GeV; Dilepton p_{T} [GeV]; Entries");
  zee_pt_low->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  c1->Print("zee_pt_low.png");
  c1->Print("zee_pt_low.pdf");

  response_low->SetTitle("Response 15 < p_{T} < 20 GeV; p_{T}^{#parallel}/p_T^{Z}; Entries");
  response_low->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "15 < p_{T}^{jet} < 20 GeV");
  c1->Print("response_low.png");
  c1->Print("response_low.pdf");

  zee_mass_high->SetTitle("Z mass 20 < p_{T} < 25 GeV; Dilepton mass [GeV]; Entries");
  zee_mass_high->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  c1->Print("zee_mass_high.png");
  c1->Print("zee_mass_high.pdf");

  zee_pt_high->SetTitle("Z pt 20 < p_{T} < 25 GeV; Dilepton p_{T} [GeV]; Entries");
  zee_pt_high->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  c1->Print("zee_pt_high.png");
  c1->Print("zee_pt_high.pdf");

  response_high->SetTitle("Response 20 < p_{T} < 25 GeV; p_{T}^{#parallel}/p_T^{Z}; Entries");
  response_high->Draw("hist");
  ATLASLabel(0.2, 0.8, "Work in Progress");
  myText(0.2, 0.72, kBlack, "20 < p_{T}^{jet} < 25 GeV");
  c1->Print("response_high.png");
  c1->Print("response_high.pdf");

  TFile* outputFile = TFile::Open("plots_zee.root", "RECREATE");
  mc_scale->Write("mc_scale");
  mc_resolution->Write("mc_resolution");
  _binned_response_1->Write();
  _binned_response_2->Write();
  _binned_response_3->Write();
  _binned_response_4->Write();
  outputFile->Close();

}

