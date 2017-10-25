#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TLatex.h>
#include <sstream>
#include <iomanip>

void addHistoLabel(std::string channel, TH2D* h2, std::string x_name, std::string y_name, Long64_t nPEs);
void printProgBar( int percent );

void makeCorrelationPlots()
{
  // *** 0. Set some overall style options
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  // *** 1. Declare histograms and canvasses
  // lep+had
  //TH2D* fl_vs_f0=new TH2D("fl_vs_f0","fl_vs_f0",100,.26,.34,100,.64,.75);
  //TH2D* fr_vs_f0=new TH2D("fr_vs_f0","fr_vs_f0",100,-0.03,0.03,100,.64,.75);
  //TH2D* fr_vs_fl=new TH2D("fr_vs_fl","fr_vs_fl",100,-0.03,0.03,100,.26,.34);
  
  // leptonic
  TH2D* fl_vs_f0=new TH2D("fl_vs_f0","fl_vs_f0",100,.26,.34,100,.64,.75);
  TH2D* fr_vs_f0=new TH2D("fr_vs_f0","fr_vs_f0",100,-0.03,0.03,100,.64,.75);
  TH2D* fr_vs_fl=new TH2D("fr_vs_fl","fr_vs_fl",100,-0.03,0.03,100,.27,.35);

  // hadronic
  //TH2D* fl_vs_f0=new TH2D("fl_vs_f0","fl_vs_f0",100,.2,.4,100,.64,.75);
  //TH2D* fr_vs_f0=new TH2D("fr_vs_f0","fr_vs_f0",100,-0.1,0.1,100,.64,.75);
  //TH2D* fr_vs_fl=new TH2D("fr_vs_fl","fr_vs_fl",100,-0.1,0.1,100,.22,.42);
  
  TCanvas *c1= new TCanvas("c1","c_1",200,10,700,500);
  TCanvas *c2= new TCanvas("c2","c2",200,10,700,500);
  TCanvas *c3= new TCanvas("c3","c3",200,10,700,500);

  // *** 2. Declare TFile and TTree
    TFile* fFile = new TFile("/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr26_allCF/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep/Syst_Nominal_el_mu_bTag.root","READ"); //leptonic
  //TFile* fFile = new TFile("/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr26_allCF/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/Syst_Nominal_el_mu_bTag.root","READ"); //hadronic
  //TFile* fFile = new TFile("/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr26_allCF/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  //TFile* fFile = new TFile("/afs/cern.ch/user/m/mkareem/public/for_Ben/output_RadTest_20Feb/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  TTree* fTree= (TTree*)fFile->Get("EnsembleTree");

  // *** 3. Declare and set variables to read from TTree
  Double_t f0;
  Double_t fL;
  Double_t fR;
  fTree->SetBranchAddress("F0", &f0);
  fTree->SetBranchAddress("FL", &fL);
  fTree->SetBranchAddress("FR", &fR);

  // *** 4. Loop over events in file and fill histograms
  Long64_t nentries = fTree->GetEntries();
  cout<<"Nentries (ttbar): "<<nentries<<endl;
  
  for (Long64_t i=0;i<nentries;i++) {
    
    // keep track with progress bar!
    if (nentries > 100) {
      if ((i+1)%(5*nentries/100)==0)	printProgBar(100*i/nentries +1); }
    if (i == nentries-1) {printProgBar(100); cout << endl;}
    
    //read entry
    fTree->GetEntry(i);

    fl_vs_f0->Fill(fL,f0);
    fr_vs_f0->Fill(fR,f0);
    fr_vs_fl->Fill(fR,fL);
  }

  c1->cd();
  //addHistoLabel("Lep+Had e+#mu", fl_vs_f0, "FL","F0", nentries);
  addHistoLabel("Lep e+#mu", fl_vs_f0, "FL","F0", nentries);
  c1->Print("may31_fractionCorrelation/lep_fl_vs_f0.eps");
  c2->cd();
  //addHistoLabel("Lep+Had e+#mu", fr_vs_f0, "FR","F0", nentries);
  addHistoLabel("Lep e+#mu", fr_vs_f0, "FR","F0", nentries);
  c2->Print("may31_fractionCorrelation/lep_fr_vs_f0.eps");
  c3->cd();
  //addHistoLabel("Lep+Had e+#mu", fr_vs_fl, "FR","FL", nentries);
  addHistoLabel("Lep e+#mu", fr_vs_fl, "FR","FL", nentries);
  c3->Print("may31_fractionCorrelation/lep_fr_vs_fl.eps");

  return;
}

void addHistoLabel(std::string channel, TH2D* h2, std::string x_name, std::string y_name, Long64_t nPEs)
//void addHistoLabel(std::string channel, TH1D* h1, TH1D* h2, TH1D* h3, std::string h1_name, std::string h2_name, std::string h3_name)
{
  TPad *pad0 =new TPad("pad0","pad0",0,0,1,1);
  pad0->SetRightMargin(0.105);
  pad0->SetLeftMargin(0.11);
  pad0->Draw();
  pad0->cd();

  h2->SetXTitle(x_name.c_str());
  h2->SetYTitle(y_name.c_str());
  h2->SetTitle(" ");
  //h2->GetYaxis()->SetTitleOffset(1.55);
  h2->GetYaxis()->SetTitleOffset(1.1);
  h2->Draw("colz");

  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextFont(72);
  l1.SetTextSize(0.035);
  l1.SetNDC();
  l1.DrawLatex(0.20, 0.850, "ATLAS simulation");

  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextSize(0.025);
  l2.SetNDC();
  l2.DrawLatex(0.20, 0.810, "Internal");

  std::string Ssqrt ="   (#sqrt{s}=8 TeV)"; 
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextFont(72);
  l3.SetTextSize(0.035);
  l3.SetNDC();
  l3.DrawLatex(0.6, 0.77,"#sqrt{s}=8 TeV");

  TLatex l4;
  l4.SetTextAlign(9);
  l4.SetTextSize(0.03);
  l4.SetNDC();
  l4.DrawLatex(0.55, 0.83,(channel+" Channel").c_str());

  stringstream ss;
  ss<<nPEs;
  std::string sPE=ss.str();
  TLatex l5;
  l5.SetTextAlign(9);
  l5.SetTextSize(0.025);
  l5.SetNDC();
  l5.DrawLatex(0.20, 0.765, (sPE+" PEs").c_str());

  stringstream ss2;
  ss2<< std::setprecision(2)<<h2->GetCorrelationFactor();
  std::string sCF=ss2.str();

  TLatex l6;
  l6.SetTextAlign(9);
  l6.SetTextSize(0.045);
  l6.SetNDC();
  l6.DrawLatex(0.27, 0.23, ("Correlation Factor: " + sCF).c_str());

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
  
  cout<< "\r" "[" << bar << "] ";
  cout.width( 3 );
  cout<< percent << "%     " << std::flush;
}
