#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TFileCollection.h>
#include <TChain.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLatex.h>
#include <sstream>
#include <iostream>
#include <iomanip>

void printProgBar( int percent );

void pileupCorrelation()
{
  
  // *** 0. Set some overall style options
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  // *** 1. Declaring some histograms and tree variables

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TCanvas *c3 = new TCanvas("c3","c3",800,800);

  TH2D* h_lep_lowmu = new TH2D("h_lep_lowmu","h_lep_lowmu",15,-1,1,41,0,40);
  TH2D* h_lep_highmu = new TH2D("h_lep_highmu","h_lep_highmu",15,-1,1,41,0,40);
  TH2D* h_had_lowmu = new TH2D("h_had_lowmu","h_had_lowmu",15,-1,1,41,0,40);
  TH2D* h_had_highmu = new TH2D("h_had_highmu","h_had_highmu",15,-1,1,41,0,40);
  TH1D* h_mu = new TH1D("h_mu","h_mu",41,0,40);
  TH1D* h_costheta_lep_lowmu = new TH1D("h_costheta_lep_lowmu","h_costheta_lep_lowmu",15,-1,1);
  TH1D* h_costheta_lep_highmu = new TH1D("h_costheta_lep_highmu","h_costheta_lep_highmu",15,-1,1);
  TH1D* h_costheta_had_lowmu = new TH1D("h_costheta_had_lowmu","h_costheta_had_lowmu",15,-1,1);
  TH1D* h_costheta_had_highmu = new TH1D("h_costheta_had_highmu","h_costheta_had_highmu",15,-1,1);
 
  Float_t mu;
  Float_t cos_theta_had_reco;
  Float_t cos_theta_lep_reco;
  Float_t nBtags;

  // *** 2. Read in files from .txt and add to chain
  TChain* fChain = new TChain("nominal");
  //TString inputFile = "eos_ttbar_el.txt";
  TString inputFile = "eos_ttbar_mu.txt";

  if(gSystem->AccessPathName(inputFile))
    {
      std::cout<<"Input list file /"<<inputFile<<" DNE."<<std::endl;
      gSystem->Exit();
    }
  TFileCollection *fc= new TFileCollection("fc","files",inputFile);
  fChain->AddFileInfoList(fc->GetList());

  fChain->SetBranchAddress("mu",&mu);
  fChain->SetBranchAddress("cos_theta_had_reco",&cos_theta_had_reco);
  fChain->SetBranchAddress("cos_theta_lep_reco",&cos_theta_lep_reco);
  fChain->SetBranchAddress("nBtags",&nBtags);

  // *** 3. Loop over events in file and fill histograms
  Long64_t nentries = fChain->GetEntries();
  //Long64_t nentries = 50e3;
  std::cout<<"Nentries (ttbar): "<<nentries<<std::endl;
  
  for (Long64_t i=0;i<nentries;i++) {
    
    // keep track with progress bar!
    if (nentries > 100) {
      if ((i+1)%(5*nentries/100)==0)	printProgBar(100*i/nentries +1); }
    if (i == nentries-1) {printProgBar(100); std::cout << std::endl;}
    
    //read entry
    fChain->GetEntry(i);
    if(nBtags>=2){
    //if(nBtags==1){
      h_mu->Fill(mu);
      if(mu < 20){
	h_lep_lowmu->Fill(cos_theta_lep_reco,mu);
	h_had_lowmu->Fill(cos_theta_had_reco,mu);
	h_costheta_lep_lowmu->Fill(cos_theta_lep_reco);
	h_costheta_had_lowmu->Fill(cos_theta_had_reco);
      }
      else{
	h_lep_highmu->Fill(cos_theta_lep_reco,mu);
	h_had_highmu->Fill(cos_theta_had_reco,mu);
	h_costheta_lep_highmu->Fill(cos_theta_lep_reco);
	h_costheta_had_highmu->Fill(cos_theta_had_reco);
      }
    }
  }

  //make latex labels
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextFont(72);
  l1.SetTextSize(0.035);
  l1.SetNDC();

  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextSize(0.025);
  l2.SetNDC();
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextSize(0.025);
  l3.SetNDC();

  TLegend* leg1=new TLegend(.2,.51,.45,.81);
  //leg1-> SetBorderSize(0);
  leg1->AddEntry(h_costheta_lep_lowmu,"#mu < 20","lf");
  leg1->AddEntry(h_costheta_lep_highmu,"#mu #geq 20","lf");
  
  c1->cd();
  h_mu->GetYaxis()->SetRangeUser(0,1.2*h_mu->GetMaximum());
  h_mu->SetXTitle("#mu");
  h_mu->SetYTitle("Entries / 1");

  h_mu->Draw();
  l1.DrawLatex(0.20, 0.890, "ATLAS simulation");
  l2.DrawLatex(0.20, 0.850, "work in progress");
  l3.DrawLatex(0.20, 0.810, "#mu+jets, #geq 2 b-tags");
  //l3.DrawLatex(0.20, 0.810, "#mu+jets, =1 b-tags");

  c2->cd();
  h_costheta_lep_lowmu->SetXTitle("Cos #theta * (Leptonic)");
  h_costheta_lep_lowmu->SetYTitle("Normalized Entries / 0.133");
  h_costheta_lep_lowmu->GetYaxis()->SetRangeUser(0,1.3*h_costheta_lep_lowmu->GetMaximum());
  h_costheta_lep_lowmu->SetLineColor(kBlue);
  h_costheta_lep_lowmu->SetMarkerStyle(0);
  h_costheta_lep_highmu->SetLineColor(kBlack);
  h_costheta_lep_highmu->SetMarkerStyle(0);
  h_costheta_lep_lowmu->Sumw2();
  h_costheta_lep_highmu->Sumw2();
  h_costheta_lep_lowmu->DrawNormalized("el");
  h_costheta_lep_highmu->DrawNormalized("el same");
  l1.DrawLatex(0.20, 0.890, "ATLAS simulation");
  l2.DrawLatex(0.20, 0.850, "work in progress");
  l3.DrawLatex(0.20, 0.810, "#mu+jets, #geq 2 b-tags");
  //l3.DrawLatex(0.20, 0.810, "#mu+jets, =1 b-tags");
  TLegend* leg1=new TLegend(.75,.71,.9,.91);
  //leg1-> SetBorderSize(0);
  leg1->AddEntry(h_costheta_lep_lowmu,"#mu < 20","lf");
  leg1->AddEntry(h_costheta_lep_highmu,"#mu #geq 20","lf");
  leg1->Draw("same");

  c3->cd();
  h_costheta_had_lowmu->SetXTitle("Cos #theta * (Hadronic)");
  h_costheta_had_lowmu->SetYTitle("Normalized Entries / 0.133");
  h_costheta_had_lowmu->GetYaxis()->SetRangeUser(0,1.3*h_costheta_had_lowmu->GetMaximum());
  h_costheta_had_lowmu->SetLineColor(kBlue);
  h_costheta_had_lowmu->SetMarkerStyle(0);
  h_costheta_had_highmu->SetLineColor(kBlack);
  h_costheta_had_highmu->SetMarkerStyle(0);
  h_costheta_had_lowmu->Sumw2();
  h_costheta_had_highmu->Sumw2();
  h_costheta_had_lowmu->DrawNormalized("el");
  h_costheta_had_highmu->DrawNormalized("el same");
  l1.DrawLatex(0.20, 0.890, "ATLAS simulation");
  l2.DrawLatex(0.20, 0.850, "work in progress");
  l3.DrawLatex(0.20, 0.810, "#mu+jets, #geq 2 b-tags");
  //l3.DrawLatex(0.20, 0.810, "#mu+jets, =1 b-tags");
  leg1->Draw("same");

  std::cout<<"Leptonic, low mu correlation factor: "<< std::setprecision(2)<<h_lep_lowmu->GetCorrelationFactor()<<std::endl;
  std::cout<<"Leptonic, high mu correlation factor: "<< std::setprecision(2)<<h_lep_highmu->GetCorrelationFactor()<<std::endl;
  std::cout<<"Hadronic, low mu correlation factor: "<< std::setprecision(2)<<h_had_lowmu->GetCorrelationFactor()<<std::endl;
  std::cout<<"Hadronic, high mu correlation factor: "<< std::setprecision(2)<<h_had_highmu->GetCorrelationFactor()<<std::endl;

  return;
}


void printProgBar( int percent )
{
  string bar;
  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  
  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}
