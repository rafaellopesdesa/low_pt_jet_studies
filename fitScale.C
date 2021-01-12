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
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TVector2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDoubleCBFast.h"
#include "RooGaussian.h"
#include "RooPlot.h" 
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"

#include "AtlasLabels.h"
#include "AtlasUtils.h"

#include <vector>
#include <iostream>
#include <utility>

#define debug false
#define doOnlyMuons false
#define doOnlyElectrons false

TF1* _mc_scale_func;
double _p0_func;
double _p1_func;
int _fit_bin;

TH1D* _binned_response_1;
TH1D* _binned_response_2;
TH1D* _binned_response_3;
TH1D* _binned_response_4;

RooRealVar pt_jet("pt_jet", "pt_jet", 0);
RooRealVar pt_zee("pt_zee", "pt_zee", 0);
RooRealVar eta_jet("eta_jet", "eta_jet", 0);
RooRealVar eta_zee("eta_zee", "eta_zee", 0);
RooRealVar phi_jet("phi_jet", "phi_jet", 0);
RooRealVar phi_zee("phi_zee", "phi_zee", 0);
RooArgSet hs_jet("hs_jet");

RooDataSet* _dataset;
TFile* _outputFile_1;
TFile* _outputFile_2;
TFile* _outputFile_3;
TFile* _outputFile_4;

bool _save_file;
int _save_i;
int _save_j;

double fcn(const double* par) {

  RooFormulaVar response_form("response_form", "p_{T}^{#parallel}/p_{T}^{Z}-1", "((-1.*@0*cos(@1-@2)/@3)-1.)", RooArgList(pt_jet, phi_jet, phi_zee, pt_zee));    
  RooDataSet* _reduced_dataset_1 = (RooDataSet*) _dataset->reduce(RooArgSet(pt_jet, pt_zee, phi_jet, phi_zee), TString::Format("%f*pt_jet > 15.0 && %f*pt_jet <= 17.5", par[0], par[0]));
  RooDataSet* _reduced_dataset_2 = (RooDataSet*) _dataset->reduce(RooArgSet(pt_jet, pt_zee, phi_jet, phi_zee), TString::Format("%f*pt_jet > 17.5 && %f*pt_jet <= 20.0", par[2], par[2]));
  RooDataSet* _reduced_dataset_3 = (RooDataSet*) _dataset->reduce(RooArgSet(pt_jet, pt_zee, phi_jet, phi_zee), TString::Format("%f*pt_jet > 20.0 && %f*pt_jet <= 22.5", par[4], par[4]));
  RooDataSet* _reduced_dataset_4 = (RooDataSet*) _dataset->reduce(RooArgSet(pt_jet, pt_zee, phi_jet, phi_zee), TString::Format("%f*pt_jet > 22.5 && %f*pt_jet <= 25.0", par[6], par[6]));

  RooRealVar* _response_1 = (RooRealVar*) _reduced_dataset_1->addColumn(response_form);
  RooFormulaVar calib_form_1("calib_form_1", "calib", TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[0], _p1_func, _p0_func, par[1], _p1_func, _p0_func), RooArgList(*_response_1, pt_jet));
  RooRealVar* _calib_response_1 = (RooRealVar*) _reduced_dataset_1->addColumn(calib_form_1);

  RooRealVar* _response_2 = (RooRealVar*) _reduced_dataset_2->addColumn(response_form);
  RooFormulaVar calib_form_2("calib_form_2", "calib", TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[2], _p1_func, _p0_func, par[3], _p1_func, _p0_func), RooArgList(*_response_2, pt_jet));
  RooRealVar* _calib_response_2 = (RooRealVar*) _reduced_dataset_2->addColumn(calib_form_2);

  RooRealVar* _response_3 = (RooRealVar*) _reduced_dataset_3->addColumn(response_form);
  RooFormulaVar calib_form_3("calib_form_3", "calib", TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[4], _p1_func, _p0_func, par[5], _p1_func, _p0_func), RooArgList(*_response_3, pt_jet));
  RooRealVar* _calib_response_3 = (RooRealVar*) _reduced_dataset_3->addColumn(calib_form_3);

  RooRealVar* _response_4 = (RooRealVar*) _reduced_dataset_4->addColumn(response_form);
  RooFormulaVar calib_form_4("calib_form_4", "calib", TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[6], _p1_func, _p0_func, par[7], _p1_func, _p0_func), RooArgList(*_response_4, pt_jet));
  RooRealVar* _calib_response_4 = (RooRealVar*) _reduced_dataset_4->addColumn(calib_form_4);

  TH1D* _binned_data_1 = (TH1D*) _reduced_dataset_1->createHistogram("binned_data_1", *_calib_response_1, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_data_2 = (TH1D*) _reduced_dataset_2->createHistogram("binned_data_2", *_calib_response_2, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_data_3 = (TH1D*) _reduced_dataset_3->createHistogram("binned_data_3", *_calib_response_3, RooFit::Binning(12,-0.6,0.6));
  TH1D* _binned_data_4 = (TH1D*) _reduced_dataset_4->createHistogram("binned_data_4", *_calib_response_4, RooFit::Binning(12,-0.6,0.6));

  double chi2 = 0;
  
  double scale_1 = _binned_data_1->Integral();
  double scale_2 = _binned_data_2->Integral();
  double scale_3 = _binned_data_3->Integral();
  double scale_4 = _binned_data_4->Integral();
  
  if (_fit_bin == 1) {
    for (int i=1; i<_binned_data_1->GetNbinsX(); i++) {
      if (_binned_data_1->GetBinContent(i) > 0) 
	chi2 += (_binned_data_1->GetBinContent(i) - _binned_response_1->GetBinContent(i)*scale_1)*(_binned_data_1->GetBinContent(i) - _binned_response_1->GetBinContent(i)*scale_1)/_binned_data_1->GetBinContent(i);
    }
  } else if (_fit_bin == 2) {
    for (int i=1; i<_binned_data_2->GetNbinsX(); i++) {
      if (_binned_data_2->GetBinContent(i) > 0) 
	chi2 += (_binned_data_2->GetBinContent(i) - _binned_response_2->GetBinContent(i)*scale_2)*(_binned_data_2->GetBinContent(i) - _binned_response_2->GetBinContent(i)*scale_2)/_binned_data_2->GetBinContent(i);
    }
  } else if (_fit_bin == 3) {
    for (int i=1; i<_binned_data_3->GetNbinsX(); i++) {
      if (_binned_data_3->GetBinContent(i) > 0) 
	chi2 += (_binned_data_3->GetBinContent(i) - _binned_response_3->GetBinContent(i)*scale_3)*(_binned_data_3->GetBinContent(i) - _binned_response_3->GetBinContent(i)*scale_3)/_binned_data_3->GetBinContent(i);
    }
  } else if (_fit_bin == 4) {
    for (int i=1; i<_binned_data_4->GetNbinsX(); i++) {
      if (_binned_data_4->GetBinContent(i) > 0) 
	chi2 += (_binned_data_4->GetBinContent(i) - _binned_response_4->GetBinContent(i)*scale_4)*(_binned_data_4->GetBinContent(i) - _binned_response_4->GetBinContent(i)*scale_4)/_binned_data_4->GetBinContent(i);
    }
  } 

  if (debug) {
    std::cout << "RCLSA " <<  TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[0], _p1_func, _p0_func, par[1], _p1_func, _p0_func).Data() << std::endl;
    std::cout << "RCLSA " <<  TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[2], _p1_func, _p0_func, par[3], _p1_func, _p0_func).Data() << std::endl;
    std::cout << "RCLSA " <<  TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[4], _p1_func, _p0_func, par[5], _p1_func, _p0_func).Data() << std::endl;
    std::cout << "RCLSA " <<  TString::Format("(%f*((%f*@1)+(%f)+1)+%f*(@0-(%f*@1)-(%f))-1)", par[6], _p1_func, _p0_func, par[7], _p1_func, _p0_func).Data() << std::endl;
    std::cout << "RCLSA " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << " " << par[7] << " " << chi2 << std::endl;
    std::cout << "RCLSA " << _binned_data_1->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_data_1->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_data_1->GetBinLowEdge(i) << " " << _binned_data_1->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_data_2->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_data_2->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_data_2->GetBinLowEdge(i) << " " << _binned_data_2->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_data_3->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_data_3->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_data_3->GetBinLowEdge(i) << " " << _binned_data_3->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_data_4->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_data_4->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_data_4->GetBinLowEdge(i) << " " << _binned_data_4->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_response_1->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_response_1->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_response_1->GetBinLowEdge(i) << " " << _binned_response_1->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_response_2->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_response_2->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_response_2->GetBinLowEdge(i) << " " << _binned_response_2->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_response_3->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_response_3->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_response_3->GetBinLowEdge(i) << " " << _binned_response_3->GetBinContent(i) << std::endl;
    std::cout << "RCLSA " << _binned_response_4->GetNbinsX() << std::endl;
    for (int i=1; i<= _binned_response_4->GetNbinsX(); i++)
      std::cout << "RCLSA    " << _binned_response_4->GetBinLowEdge(i) << " " << _binned_response_4->GetBinContent(i) << std::endl;
  }

  if (_fit_bin==1 && _save_file) {
    _outputFile_1->cd();
    _binned_data_1->Write(TString::Format("%d_%d", _save_i, _save_j));
  } else if (_fit_bin==2 && _save_file) {
    _outputFile_2->cd();
    _binned_data_2->Write(TString::Format("%d_%d", _save_i, _save_j));
  } else if (_fit_bin==3 && _save_file) {
    _outputFile_3->cd();
    _binned_data_3->Write(TString::Format("%d_%d", _save_i, _save_j));
  } else if (_fit_bin==4 && _save_file) {
    _outputFile_4->cd();
    _binned_data_4->Write(TString::Format("%d_%d", _save_i, _save_j));
  }
  
  delete _binned_data_1;
  delete _binned_data_2;
  delete _binned_data_3;
  delete _binned_data_4;
  delete _reduced_dataset_1;
  delete _reduced_dataset_2;
  delete _reduced_dataset_3;
  delete _reduced_dataset_4;
  /*
  delete _response_1;
  delete _response_2;
  delete _response_3;
  delete _response_4;
  delete _calib_response_1;
  delete _calib_response_2;
  delete _calib_response_3;
  delete _calib_response_4;
  */
  return chi2;

}

int main(int argc, char** argv) {

  double GeV = 1000.;

  hs_jet.add(pt_jet);
  hs_jet.add(pt_zee);
  hs_jet.add(eta_jet);
  hs_jet.add(eta_zee);
  hs_jet.add(phi_jet);
  hs_jet.add(phi_zee);
  _dataset = new RooDataSet("dataset", "", hs_jet);

  int mm, ee;

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

  for (int ifile = 1; ifile<=437; ifile++) {
    //  for (int ifile = 1; ifile<=40; ifile++) {

    std::cout << "Reading file " << ifile << std::endl;
    TFile* zeeFile = TFile::Open(TString::Format("/eos/user/r/rcoelhol/jet_studies/AnalysisTop/data18/output_old_%d.root", ifile));  
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
  
    for (size_t ievent=0; ievent < zeeReco->GetEntries(); ievent++) {

      if (ievent % 1000 == 0) std::cout << "Event " << ievent << "/" << zeeReco->GetEntries() << std::endl;
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

      if (jet_pt[0]/GeV < 10. || jet_pt[0]/GeV > 30.) continue;
      if (std::abs(jet_eta[0]) > 2.5) continue;      
      if (jet_pt.size() > 1) continue;

      TLorentzVector jet(0,0,0,0);
      jet.SetPtEtaPhiE(jet_pt[0]/GeV, jet_eta[0], jet_phi[0], jet_e[0]/GeV);

      if (z.DeltaPhi(jet) < 2.9) continue;
    
      TVector2 ptz = z.Vect().XYvector();
      TVector2 ptj = jet.Vect().XYvector();    
      TVector2 ptj_ref = ptj.Proj(ptz);

      pt_jet.setVal(jet.Pt());
      eta_jet.setVal(jet.Eta());
      phi_jet.setVal(jet.Phi());
      pt_zee.setVal(z.Pt());
      eta_zee.setVal(z.Eta());
      phi_zee.setVal(z.Phi());
      _dataset->add(hs_jet);       
      
      if (debug) std::cout << ptj_ref.Mod()/ptz.Mod() - 1. << " " << jet.Pt() << " " << jet.Phi() << " " << z.Pt() << " "  << z.Phi() << " " << ptz.Mod() << " " << -jet.Pt()*cos(jet.Phi()-z.Phi()) << " " << ptj_ref.Mod() << " " << (-jet.Pt()*cos(jet.Phi()-z.Phi())/z.Pt()) - 1 << std::endl;
    }
    zeeFile->Close();
  }

  TFile* modelFile = TFile::Open("plots_zee.root");
  TGraphErrors* mc_scale_fit = (TGraphErrors*) modelFile->Get("mc_scale_fit");
  _mc_scale_func = (TF1*) (mc_scale_fit->GetListOfFunctions()->At(0));
  _p0_func = _mc_scale_func->GetParameter(0);
  _p1_func = _mc_scale_func->GetParameter(1);
  _binned_response_1 = (TH1D*) modelFile->Get("binned_response_1___response");
  _binned_response_2 = (TH1D*) modelFile->Get("binned_response_2___response");
  _binned_response_3 = (TH1D*) modelFile->Get("binned_response_3___response");
  _binned_response_4 = (TH1D*) modelFile->Get("binned_response_4___response");

  double scale_1 = _binned_response_1->Integral();
  _binned_response_1->Scale(1./scale_1);
  double scale_2 = _binned_response_2->Integral();
  _binned_response_2->Scale(1./scale_2);
  double scale_3 = _binned_response_3->Integral();
  _binned_response_3->Scale(1./scale_3);
  double scale_4 = _binned_response_4->Integral();
  _binned_response_4->Scale(1./scale_4);

  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(1000000);  // for GSL
  minimum->SetErrorDef(1);
  minimum->SetTolerance(50);
  minimum->SetPrintLevel(1);
 
  ROOT::Math::Functor f(&fcn,8);
  double step[8] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  double variable[8] = {1.,1.,1.,1.,1.,1.,1.,1.}; 
  double var_down[8] = {0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8}; 
  double var_up[8] = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2}; 
  minimum->SetFunction(f);
 
  minimum->SetVariable(0,"scale 15",variable[0], step[0]);
  minimum->SetVariableLimits(0,var_down[0], var_up[0]);
  minimum->SetVariable(1,"smear 15",variable[1], step[1]);
  minimum->SetVariableLimits(1,var_down[1], var_up[1]);
  minimum->SetVariable(2,"scale 17.5",variable[2], step[2]);
  minimum->SetVariableLimits(2,var_down[2], var_up[2]);
  minimum->SetVariable(3,"smear 17.5",variable[3], step[3]);
  minimum->SetVariableLimits(3,var_down[3], var_up[3]);
  minimum->SetVariable(4,"scale 20",variable[4], step[4]);
  minimum->SetVariableLimits(4,var_down[4], var_up[4]);
  minimum->SetVariable(5,"smear 20",variable[5], step[5]);
  minimum->SetVariableLimits(5,var_down[5], var_up[5]);
  minimum->SetVariable(6,"scale 22.5",variable[6], step[6]);
  minimum->SetVariableLimits(6,var_down[6], var_up[6]);
  minimum->SetVariable(7,"smear 22.5",variable[7], step[7]);
  minimum->SetVariableLimits(7,var_down[7], var_up[7]);

  _save_file = false;
  _fit_bin = 1;
  minimum->ReleaseVariable(0);
  minimum->ReleaseVariable(1);
  minimum->FixVariable(2);
  minimum->FixVariable(3);
  minimum->FixVariable(4);
  minimum->FixVariable(5);
  minimum->FixVariable(6);
  minimum->FixVariable(7);
  minimum->Minimize();

  _fit_bin = 2;
  minimum->ReleaseVariable(2);
  minimum->ReleaseVariable(3);
  minimum->FixVariable(0);
  minimum->FixVariable(1);
  minimum->FixVariable(4);
  minimum->FixVariable(5);
  minimum->FixVariable(6);
  minimum->FixVariable(7);
  minimum->Minimize();

  _fit_bin = 3;
  minimum->ReleaseVariable(4);
  minimum->ReleaseVariable(5);
  minimum->FixVariable(0);
  minimum->FixVariable(1);
  minimum->FixVariable(2);
  minimum->FixVariable(3);
  minimum->FixVariable(6);
  minimum->FixVariable(7);
  minimum->Minimize();

  _fit_bin = 4;
  minimum->ReleaseVariable(6);
  minimum->ReleaseVariable(7);
  minimum->FixVariable(0);
  minimum->FixVariable(1);
  minimum->FixVariable(2);
  minimum->FixVariable(3);
  minimum->FixVariable(4);
  minimum->FixVariable(5);
  minimum->Minimize();

  const double *xs = minimum->X();
  const double *es = minimum->Errors();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "," << xs[3] << "," << xs[4] << "," << xs[5] << "," << xs[6] <<"," << xs[7] << "): "
	    << minimum->MinValue()  << std::endl;
  std::cout << "Error  : f(" << es[0] << "," << es[1] << "," << es[2] << "," << es[3] << "," << es[4] << "," << es[5] << "," << es[6] <<"," << es[7] << ")" << std::endl;
  
  _outputFile_1 = TFile::Open("fit_hists_1.root", "RECREATE");
  _outputFile_2 = TFile::Open("fit_hists_2.root", "RECREATE");
  _outputFile_3 = TFile::Open("fit_hists_3.root", "RECREATE");
  _outputFile_4 = TFile::Open("fit_hists_4.root", "RECREATE");
  _save_file = true;
  for (int i=0; i<50; i++) {
    for (int j=0; j<50; j++) {
      _save_i = i;
      _save_j = j;
      std::cout << "Saving " << i << " " << j << std::endl;
      double p0 = 0.99 + i*0.0004;
      double p1 = 0.99 + j*0.0004;
      double par_1[8] = {p0,p1,1.,1.,1.,1.,1.,1.};
      double par_2[8] = {1.,1.,p0,p1,1.,1.,1.,1.};
      double par_3[8] = {1.,1.,1.,1.,p0,p1,1.,1.};
      double par_4[8] = {1.,1.,1.,1.,1.,1.,p0,p1};
      _fit_bin = 1;
      fcn(par_1);
      _fit_bin = 2;
      fcn(par_2);
      _fit_bin = 3;
      fcn(par_3);
      _fit_bin = 4;
      fcn(par_4);
    }
  }
  _outputFile_1->cd();
  _binned_response_1->Write("mc");
  _outputFile_2->cd();
  _binned_response_2->Write("mc");
  _outputFile_3->cd();
  _binned_response_3->Write("mc");
  _outputFile_4->cd();
  _binned_response_4->Write("mc");
  
  _outputFile_1->Close();
  _outputFile_2->Close();
  _outputFile_3->Close();
  _outputFile_4->Close();

  modelFile->Close();

}  
