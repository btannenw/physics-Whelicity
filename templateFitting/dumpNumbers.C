#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>

void dumpNumbers()
{
  //std::string regionDir = "../syst_outputs_8ch_03Aug/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/"; //8-channel
  //std::string regionDir = "../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/"; // had 1excl2incl
  //std::string regionDir = "../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_2incl_3D_3W_had/"; // had 2incl
  //std::string regionDir = "../syst_outputs_lep_2incl_aug02/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/"; //lep 2incl
  //std::string regionDir = "ExternalSystematicsOutput_el_1excl_3D_3W_had/"; // had 1excl, el
  //std::string regionDir = "ExternalSystematicsOutput_mu_1excl_3D_3W_had/"; // had 1excl, mu
  std::string regionDir = "ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/"; // had 1excl2incl
  
  TFile* fNom = new TFile( (regionDir + "Syst_Nominal_el_mu_bTag.root").c_str(), "READ");
  //TFile* fAFII = new TFile( (regionDir + "Syst_Nominal_AFII_el.root").c_str(), "READ");
  //TFile* fPS_pythia = new TFile( (regionDir + "Syst_Powheg+Pythia6(PS)_el.root").c_str(), "READ");
  //TFile* fPS_herwig = new TFile( (regionDir + "Syst_Powheg+fHerwig(PS)_el.root").c_str(), "READ");
  //TFile* fRad_hi = new TFile( (regionDir + "Syst_mu\=2\ hdamp\=mtop_el.root").c_str(), "READ");
  //TFile* fRad_lo = new TFile( (regionDir + "Syst_mu\=0.5\ hdamp\=2mtop_el.root").c_str(), "READ");
  //TFile* fME_powheg = new TFile( (regionDir + "Syst_Powheg+fHerwig_el.root").c_str(), "READ");
  //TFile* fME_nlo = new TFile( (regionDir + "Syst_MC@NLO+fHerwig_el.root").c_str(), "READ");
  TFile* fbtag5_up = new TFile( (regionDir + "Syst_BTAG_bTagVarUp_5_el_mu_bTag.root ").c_str(), "READ");
  TFile* fbtag5_down = new TFile( (regionDir + "Syst_BTAG_bTagVarDown_5_el_mu_bTag.root ").c_str(), "READ");

  cout<<"############ Nominal ############"<<endl;
  fitNumbers(fNom);
  /*cout<<"############ AFII Nominal ############"<<endl;
  fitNumbers(fAFII);
  cout<<"########## PS(Pythia) ##########"<<endl;
  fitNumbers(fPS_pythia);
  cout<<"########## PS(Herwig) [HERE] ##########"<<endl;
  fitNumbers(fPS_herwig);
  cout<<"########## Radiation(mu=2, hdamp=mtop) ##########"<<endl;
  fitNumbers(fRad_hi);
  cout<<"########## Radiation(mu=0.5, hdamp=2*mtop) ##########"<<endl;
  fitNumbers(fRad_lo);
  cout<<"##########  ME Gen (Powheg+Herwig)  [HERE]  ##########"<<endl;
  fitNumbers(fME_powheg);
  cout<<"##########  ME Gen (MC@NLO+Herwig)  ##########"<<endl;
  fitNumbers(fME_nlo);*/
  cout<<"##########  bTag5 Up  ##########"<<endl;
  fitNumbers(fbtag5_up);
  cout<<"##########  bTag5 Down  ##########"<<endl;
  fitNumbers(fbtag5_down);

}

void fitNumbers(TFile* fIn)
{
  Double_t m_QCD_4incl_2incl_el_nom;
  Double_t m_QCD_4incl_2incl_mu_nom;
  Double_t m_RemBkg_4incl_2incl_el_nom;
  Double_t m_RemBkg_4incl_2incl_mu_nom;
  
  Double_t m_QCD_4incl_1excl_el_nom;
  Double_t m_QCD_4incl_1excl_mu_nom;
  Double_t m_RemBkg_4incl_1excl_el_nom;
  Double_t m_RemBkg_4incl_1excl_mu_nom;

  Double_t m_WLight_4incl_nom;
  Double_t m_Wc_4incl_nom;
  Double_t m_Wbbcc_4incl_nom;
  Double_t m_RemBkg_4incl_nom;
  
  TTree* tree = (TTree*)fIn->Get("EnsembleTree");
  
  //tree->SetBranchAddress("WLight_nom", &m_WLight_4incl_nom);
  //tree->SetBranchAddress("Wc_nom", &m_Wc_4incl_nom);
  //tree->SetBranchAddress("Wbbcc_nom", &m_Wbbcc_4incl_nom); 
  //tree->SetBranchAddress("WLight_4incl_2incl_nom", &m_WLight_4incl_nom);
  //tree->SetBranchAddress("Wc_4incl_2incl_nom", &m_Wc_4incl_nom);
  //tree->SetBranchAddress("Wbbcc_4incl_2incl_nom", &m_Wbbcc_4incl_nom); 
  tree->SetBranchAddress("WLight_4incl_nom", &m_WLight_4incl_nom);
  tree->SetBranchAddress("Wc_4incl_nom", &m_Wc_4incl_nom);
  tree->SetBranchAddress("Wbbcc_4incl_nom", &m_Wbbcc_4incl_nom); 
  tree->SetBranchAddress("QCD_4incl_2incl_el_nom", &m_QCD_4incl_2incl_el_nom);
  tree->SetBranchAddress("QCD_4incl_2incl_mu_nom", &m_QCD_4incl_2incl_mu_nom);
  tree->SetBranchAddress("QCD_4incl_1excl_el_nom", &m_QCD_4incl_1excl_el_nom);
  tree->SetBranchAddress("QCD_4incl_1excl_mu_nom", &m_QCD_4incl_1excl_mu_nom);
  //tree->SetBranchAddress("QCD_nom", &m_QCD_4incl_1excl_mu_nom);
  tree->SetBranchAddress("RemBkg_4incl_nom", &m_RemBkg_4incl_nom);
  //tree->SetBranchAddress("RemBkg_4incl_2incl_nom", &m_RemBkg_4incl_nom);
  //tree->SetBranchAddress("RemBkg_nom", &m_RemBkg_4incl_nom);
  tree->GetEntry(0);
   
  TF1* n0_gaus = new TF1("n0_gaus","gaus");
  TF1* nL_gaus = new TF1("nL_gaus","gaus");
  TF1* nR_gaus = new TF1("nR_gaus","gaus");
  TF1* f0_gaus = new TF1("f0_gaus","gaus");
  TF1* fL_gaus = new TF1("fL_gaus","gaus");
  TF1* fR_gaus = new TF1("fR_gaus","gaus");
  TF1* wlight_gaus_2incl = new TF1("wlight_gaus_2incl","gaus");
  TF1* wc_gaus_2incl = new TF1("wc_gaus_2incl","gaus");
  TF1* wbbcc_gaus_2incl = new TF1("wbbcc_gaus_2incl","gaus");
  TF1* wlight_gaus_1excl = new TF1("wlight_gaus_1excl","gaus");
  TF1* wc_gaus_1excl = new TF1("wc_gaus_1excl","gaus");
  TF1* wbbcc_gaus_1excl = new TF1("wbbcc_gaus_1excl","gaus");
  TF1* qcd_gaus_el_2incl = new TF1("qcd_gaus_el_2incl","gaus");
  TF1* qcd_gaus_mu_2incl = new TF1("qcd_gaus_mu_2incl","gaus");
  TF1* qcd_gaus_el_1excl = new TF1("qcd_gaus_el_1excl","gaus");
  TF1* qcd_gaus_mu_1excl = new TF1("qcd_gaus_mu_1excl","gaus");
  TF1* rem_gaus = new TF1("rem_gaus","gaus");

  tree->Fit("n0_gaus","N0","","Q");
  tree->Fit("nL_gaus","NL","","Q");
  tree->Fit("nR_gaus","NR","","Q");
  tree->Fit("f0_gaus","F0","","Q");
  tree->Fit("fL_gaus","FL","","Q");
  tree->Fit("fR_gaus","FR","","Q");
  //tree->Fit("wlight_gaus_2incl","WLight","","Q");
  //tree->Fit("wc_gaus_2incl","Wc","","Q");
  //tree->Fit("wbbcc_gaus_2incl","Wbbcc","","Q");
  //tree->Fit("wlight_gaus_2incl","WLight_4incl_2incl","","Q");
  //tree->Fit("wc_gaus_2incl","Wc_4incl_2incl","","Q");
  //tree->Fit("wbbcc_gaus_2incl","Wbbcc_4incl_2incl","","Q");
  tree->Fit("wlight_gaus_2incl","WLight_4incl","","Q");
  tree->Fit("wc_gaus_2incl","Wc_4incl","","Q");
  tree->Fit("wbbcc_gaus_2incl","Wbbcc_4incl","","Q");
  tree->Fit("qcd_gaus_el_2incl","QCD_4incl_2incl_el","","Q");
  tree->Fit("qcd_gaus_mu_2incl","QCD_4incl_2incl_mu","","Q");
  tree->Fit("qcd_gaus_el_1excl","QCD_4incl_1excl_el","","Q");
  tree->Fit("qcd_gaus_mu_1excl","QCD_4incl_1excl_mu","","Q");
  //tree->Fit("qcd_gaus_el_1excl","QCD","","Q");
  //tree->Fit("rem_gaus","RemBkg_4incl_2incl","","Q");
  tree->Fit("rem_gaus","RemBkg_4incl","","Q");
  //tree->Fit("rem_gaus","RemBkg","","Q");

  cout<<"Param \t\t\t Fitted \t Fit Err \t Nominal \t % Diff"<<endl;
  cout<<"N0:\t\t\t"<<n0_gaus->GetParameter(1)<<"\t+/- "<<n0_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
  cout<<"NL:\t\t\t"<<nL_gaus->GetParameter(1)<<"\t\t+/- "<<nL_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
  cout<<"NR:\t\t\t"<<nR_gaus->GetParameter(1)<<"\t\t+/- "<<nR_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
  cout<<"F0:\t\t\t"<<f0_gaus->GetParameter(1)<<"\t+/- "<<f0_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
  cout<<"FL:\t\t\t"<<fL_gaus->GetParameter(1)<<"\t\t+/- "<<fL_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
  cout<<"FR:\t\t\t"<<fR_gaus->GetParameter(1)<<"\t\t+/- "<<fR_gaus->GetParameter(2)<<"\t\t\t    -"<<endl;
 
  //cout<<"Wlight:\t\t\t"<<wlight_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wlight_gaus_2incl->GetParameter(2)<<"\t"<<m_WLight_4incl_nom<<"\t\t"<<100*(1 - wlight_gaus_2incl->GetParameter(1)/m_WLight_4incl_nom)<<endl;
  //cout<<"Wc:\t\t\t"<<wc_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wc_gaus_2incl->GetParameter(2)<<"\t"<<m_Wc_4incl_nom<<"\t\t"<<100*(1 - wc_gaus_2incl->GetParameter(1)/m_Wc_4incl_nom)<<endl;
  //cout<<"Wbbcc:\t\t\t"<<wbbcc_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wbbcc_gaus_2incl->GetParameter(2)<<"\t"<<m_Wbbcc_4incl_nom<<"\t\t"<<100*(1 - wbbcc_gaus_2incl->GetParameter(1)/m_Wbbcc_4incl_nom)<<endl;
 
  cout<<"Wlight:\t\t\t"<<wlight_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wlight_gaus_2incl->GetParameter(2)<<"\t"<<m_WLight_4incl_nom<<"\t\t"<<100*(1 - wlight_gaus_2incl->GetParameter(1)/m_WLight_4incl_nom)<<endl;
  cout<<"Wc:\t\t\t"<<wc_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wc_gaus_2incl->GetParameter(2)<<"\t"<<m_Wc_4incl_nom<<"\t\t"<<100*(1 - wc_gaus_2incl->GetParameter(1)/m_Wc_4incl_nom)<<endl;
  cout<<"Wbbcc:\t\t\t"<<wbbcc_gaus_2incl->GetParameter(1)<<"\t\t+/- "<<wbbcc_gaus_2incl->GetParameter(2)<<"\t"<<m_Wbbcc_4incl_nom<<"\t\t"<<100*(1 - wbbcc_gaus_2incl->GetParameter(1)/m_Wbbcc_4incl_nom)<<endl;
  cout<<"QCD (2incl, el):\t"<<qcd_gaus_el_2incl->GetParameter(1)<<"\t\t+/- "<<qcd_gaus_el_2incl->GetParameter(2)<<"\t"<<m_QCD_4incl_2incl_el_nom<<"\t\t"<<100*(1 - qcd_gaus_el_2incl->GetParameter(1)/m_QCD_4incl_2incl_el_nom)<<endl;
  cout<<"QCD (2incl, mu):\t"<<qcd_gaus_mu_2incl->GetParameter(1)<<"\t\t+/- "<<qcd_gaus_mu_2incl->GetParameter(2)<<"\t"<<m_QCD_4incl_2incl_mu_nom<<"\t\t"<<100*(1 - qcd_gaus_mu_2incl->GetParameter(1)/m_QCD_4incl_2incl_mu_nom)<<endl;
  cout<<"QCD (1excl, el):\t"<<qcd_gaus_el_1excl->GetParameter(1)<<"\t\t+/- "<<qcd_gaus_el_1excl->GetParameter(2)<<"\t"<<m_QCD_4incl_1excl_el_nom<<"\t\t"<<100*(1 - qcd_gaus_el_1excl->GetParameter(1)/m_QCD_4incl_1excl_el_nom)<<endl;
  cout<<"QCD (1excl, mu):\t"<<qcd_gaus_mu_1excl->GetParameter(1)<<"\t\t+/- "<<qcd_gaus_mu_1excl->GetParameter(2)<<"\t"<<m_QCD_4incl_1excl_mu_nom<<"\t\t"<<100*(1 - qcd_gaus_mu_1excl->GetParameter(1)/m_QCD_4incl_1excl_mu_nom)<<endl;
  //cout<<"QCD (1excl, el):\t"<<qcd_gaus_el_1excl->GetParameter(1)<<"\t\t+/- "<<qcd_gaus_el_1excl->GetParameter(2)<<"\t"<<m_QCD_4incl_1excl_mu_nom<<"\t\t"<<100*(1 - qcd_gaus_el_1excl->GetParameter(1)/m_QCD_4incl_1excl_mu_nom)<<endl;
  cout<<"RemBkg:\t\t\t"<<rem_gaus->GetParameter(1)<<"\t\t+/- "<<rem_gaus->GetParameter(2)<<"\t"<<m_RemBkg_4incl_nom<<"\t\t"<<100*(1 - rem_gaus->GetParameter(1)/m_RemBkg_4incl_nom)<<endl;
}
