void drawChi2() {

  TFile* output = TFile::Open("chi2_file.root", "RECREATE");
  for (int m=1; m<=4; m++) {
    TFile* input = TFile::Open(TString::Format("fit_hists_%d.root", m));
    TH1D* mc = (TH1D*) input->Get("mc");

    TH2D* chi2 = new TH2D(TString::Format("chi2_%d", m), "chi2", 50, 0.99, 1.01, 50, 0.99, 1.01);
    double minchi2 = 999999999999999;
    for (int i=0; i<50; i++) {
      for (int j=0; j<50; j++) {
	TString histname = TString::Format("%d_%d", i, j);
	TH1D* data = (TH1D*) input->Get(histname);
	std::cout << histname.Data() << " " << data << " " << i << " " << j << " " << m << std::endl;
	double scale = data->Integral();
	double c2 = 0;
	for (int k=1; k<=data->GetNbinsX(); k++) {
	  if (data->GetBinContent(i) >= 0 && mc->GetBinContent(i) > 0) 
	    c2 += 2*mc->GetBinContent(i)*scale - 2*data->GetBinContent(i)*TMath::Log(mc->GetBinContent(i)*scale);
	}
	if (c2 < minchi2) minchi2 = c2;
	chi2->SetBinContent(i,j,c2);
      }
    }
    for (int i=-50; i<50; i++) {
      for (int j=-50; j<50; j++) {
	chi2->SetBinContent(i,j,chi2->GetBinContent(i,j)-minchi2);
      }
    }
    output->cd();
    chi2->Write();
    input->Close();
  }
  output->Close();
}
