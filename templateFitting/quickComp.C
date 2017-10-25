#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TLatex.h>
#include <sstream>
#include <iomanip>

void dumpNumbers(TFile* nom, TFile* var, string name);

void quickComp()
{
  TFile* fNom = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  TFile* fAFII = new TFile("../13_Jun_AFIInom_results_with_STunc17/AFIInom_results_with_STunc17/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_AFII_el_mu_lephad_bTag.root","READ");
  //TFile* fAFII = new TFile("../06_Jun_SingleTop_17percent_unc/nomAFIIresults_with_STunc17/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  TFile* fbTagUp = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_BTAG_bTagVarUp_5_el_mu_lephad_bTag.root","READ");
  TFile* fbTagDown = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_BTAG_bTagVarDown_5_el_mu_lephad_bTag.root","READ");
  TFile* fRadHi = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_mu=0.5 hdamp=2mtop_el_mu_lephad_bTag.root","READ");
  TFile* fRadLo = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_mu=2 hdamp=mtop_el_mu_lephad_bTag.root","READ");
  TFile* fAnom = new TFile("../06_Jun_SingleTopT_Protos_Anomalous/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  TFile* fSM = new TFile("../06_Jun_SingleTopT_Protos_SM/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");

  TFile* fNom_current = new TFile("ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_el_mu_lephad_bTag.root","READ");
  TFile* fAFII_current = new TFile("ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/Syst_Nominal_AFII_el_mu_lephad_bTag.root","READ");
  
  TFile* fNom_had = new TFile("../06_Jun_SingleTop_17percent_unc/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/Syst_Nominal_el_mu_bTag.root","READ");
  TFile* fAFII_had = new TFile("../06_Jun_SingleTop_17percent_unc/nomAFIIresults_with_STunc17/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/Syst_Nominal_AFII_el_mu_bTag.root","READ");



  dumpNumbers(fNom, fbTagUp, "bTagUp");
  dumpNumbers(fNom, fbTagDown, "bTagDown");
  dumpNumbers(fAFII, fRadHi, "rad_mu_2-0");
  dumpNumbers(fAFII, fRadLo, "rad_mu_0-5fg");
  dumpNumbers(fNom, fAnom, "Protos_Anomalous");
  dumpNumbers(fNom, fSM, "Protos_SM");
  dumpNumbers(fNom, fAFII, "17percent_Sanity");
  dumpNumbers(fNom_current, fAFII_current, "Current_Configuration_Sanity");
  dumpNumbers(fNom_had, fAFII_had, "17percent_Hadronic_Sanity");
  
}

void dumpNumbers(TFile* nom, TFile* var, string name)
{
  gStyle->SetOptFit();


  TTree* tNom = (TTree*)nom->Get("EnsembleTree");
  TTree* tVar = (TTree*)var->Get("EnsembleTree");

  TF1* fgaus0 = new TF1("fgaus0","gaus", 0.65, 0.75);
  TF1* fgausL = new TF1("fgausL","gaus", 0.25, 0.35);
  TF1* fgausR = new TF1("fgausR","gaus", -0.05, 0.05);

  TH1D* f0 = new TH1D("f0","f0",100, 0.65, 0.75);
  TH1D* fL = new TH1D("fL","fL",100, 0.25, 0.35);
  TH1D* fR = new TH1D("fR","fR",100, -0.05, 0.05);
  tNom->Draw("F0>>+f0"); 
  tNom->Draw("FL>>+fL");
  tNom->Draw("FR>>+fR");
   //tNom->Fit("fgaus0","F0","","Q");
  //tNom->Fit("fgausL","FL","","Q");
  //tNom->Fit("fgausR","FR","","Q");
 
  f0->Fit("fgaus0","F0","Q");
  fL->Fit("fgausL","FL","Q");
  fR->Fit("fgausR","FR","Q");
  double nomF0 = fgaus0->GetParameter(1);
  double nomFL = fgausL->GetParameter(1);
  double nomFR = fgausR->GetParameter(1);

  TH1D* f0_var = new TH1D("f0_var","f0_var",100, 0.65, 0.75);
  TH1D* fL_var = new TH1D("fL_var","fL_var",100, 0.25, 0.35);
  TH1D* fR_var = new TH1D("fR_var","fR_var",100, -0.05, 0.05);
  tVar->Draw("F0>>+f0_var"); 
  tVar->Draw("FL>>+fL_var");
  tVar->Draw("FR>>+fR_var");

  f0_var->Fit("fgaus0","F0","Q");
  fL_var->Fit("fgausL","FL","Q");
  fR_var->Fit("fgausR","FR","Q");
  //tVar->Fit("fgaus0","F0","","Q");
  //tVar->Fit("fgausL","FL","","Q");
  //tVar->Fit("fgausR","FR","","Q");
  double varF0 = fgaus0->GetParameter(1);
  double varFL = fgausL->GetParameter(1);
  double varFR = fgausR->GetParameter(1);

  std::cout<<"================================================================="<<std::endl;
  std::cout<<"\tVariation: "<<name<<std::endl;
  std::cout<<"Nom F0: "<<nomF0<<"\tNom FL: "<<nomFL<<"\tNom FR: "<<nomFR<<std::endl;
  std::cout<<"Var F0: "<<varF0<<"\tVar FL: "<<varFL<<"\tVar FR: "<<varFR<<std::endl;
  std::cout<<"Diff F0: "<<nomF0-varF0<<"\tDiff FL: "<<nomFL-varFL<<"\tDiff FR: "<<nomFR-varFR<<std::endl;
  std::cout<<"================================================================="<<std::endl;

   
  TCanvas *c0 = new TCanvas ("f0","f0",800,800);
  TCanvas *cL = new TCanvas ("fL","fL",800,800);
  TCanvas *cR = new TCanvas ("fR","fR",800,800);

  gStyle->SetOptFit();
  //gStyle->SetOptStat(1111);
  
  c0->cd();
  //tVar->Draw("F0");
  //fgaus0->Draw("same");
  //f0->Draw();
  f0_var->SetXTitle("F_{0}");
  f0_var->SetYTitle("Pseudo-Experiments / Bin");
  f0_var->Draw();
  fgaus0->Draw("same");
  c0->Print(("jun12_quickComp/"+name+"_F0.png").c_str());

  cL->cd();
  //tVar->Draw("FL");
  //fgausL->Draw("same");
  fL_var->SetXTitle("F_{L}");
  fL_var->SetYTitle("Pseudo-Experiments / Bin");
  fL_var->Draw();
  fgausL->Draw("same");
  cL->Print(("jun12_quickComp/"+name+"_FL.png").c_str());

  cR->cd();
  //tVar->Draw("FR");
  //fgausR->Draw("same");
  fR_var->SetXTitle("F_{R}");
  fR_var->SetYTitle("Pseudo-Experiments / Bin");
  fR_var->Draw();
  fgausR->Draw("same");
  cR->Print(("jun12_quickComp/"+name+"_FR.png").c_str());

  return;
}
