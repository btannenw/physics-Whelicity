#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <math.h>
#include "TImage.h"
#include "TLatex.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TF1.h>

void addHistoLabel(std::string);
void addYieldLine(ofstream& file, std::string param, Double_t nom, Double_t fit);

void FracFit()
{
  gROOT->SetBatch();
  Double_t *m_F0, *m_FL, *m_FR;
  Double_t *m_N0, *m_NL, *m_NR;
  Double_t m_Wjets_4incl_2incl_el_nom;
  Double_t m_Wjets_4incl_2incl_mu_nom;
  Double_t m_Wjets_4incl_1excl_el_nom;
  Double_t m_Wjets_4incl_1excl_mu_nom;
  Double_t m_Wjets_4incl_el_nom;
  Double_t m_Wjets_4incl_mu_nom;
  
  Double_t m_WLight_4incl_2incl_el_nom;
  Double_t m_WLight_4incl_2incl_mu_nom;
  Double_t m_Wc_4incl_2incl_el_nom;
  Double_t m_Wc_4incl_2incl_mu_nom;
  Double_t m_Wbbcc_4incl_2incl_el_nom;
  Double_t m_Wbbcc_4incl_2incl_mu_nom;

  Double_t m_QCD_4incl_2incl_el_nom;
  Double_t m_QCD_4incl_2incl_mu_nom;
  Double_t m_RemBkg_4incl_2incl_el_nom;
  Double_t m_RemBkg_4incl_2incl_mu_nom;

  Double_t m_QCD_4incl_1excl_el_nom;
  Double_t m_QCD_4incl_1excl_mu_nom;
  Double_t m_RemBkg_4incl_1excl_el_nom;
  Double_t m_RemBkg_4incl_1excl_mu_nom;

  Double_t m_WLight_4incl_el_nom;
  Double_t m_WLight_4incl_mu_nom;
  Double_t m_Wc_4incl_el_nom;
  Double_t m_Wc_4incl_mu_nom;
  Double_t m_Wbbcc_4incl_el_nom;
  Double_t m_Wbbcc_4incl_mu_nom;

  Double_t m_RemBkg_4incl_el_nom;
  Double_t m_RemBkg_4incl_mu_nom;


  std::string workDir = "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne/ExternalSystematicNoBkgConstraint/";
  //std::string OutDir = "/afs/cern.ch/work/m/mkareem/FitPackage_allinOne/output_RadTest_14Feb_Ben";
  std::string bTag="2incl"; // 2incl , 1excl2incl
  std::string channel="el_mu"; //el_mu, el_mu_bTag
  std::string wmode="3W"; // 1W , 3W
  std::string angleMode= "had"; // lep, had, lephad
  std::string bkgMode= ""; // fixedBkg_ or "" ==> "" means float

  //TFile *f_file = TFile::Open((workDir+"/ExternalCalibrationOutput_"+channel+"_"+bTag+"_3D_"+wmode+"_lephad/Syst_F0=0.7_"+channel+".root").c_str(),"READ");
  TFile *f_file = TFile::Open((workDir+"/ExternalSystematicsOutput_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"/Syst_Nominal_"+channel+".root").c_str(),"READ");

  ofstream myfile;
  //myfile.open ("example.txt",ios::out);
  std::string fname=(workDir+"/"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"_Nominal_YIELDS.txt").c_str();
  cout<<fname<<endl;
  myfile.open ( fname.c_str() ,ios::out);
  myfile<<"Parameter\t\t\t|\tNominal\t\t|\tFit\t\t|\t% Diff\t\t|"<<endl;
  myfile<<"========================================================================================================="<<endl;

  TTree *t_tree= (TTree*) f_file->Get("EnsembleTree");
  TF1 *fgaus= new TF1("fgaus","gaus");

  t_tree->SetBranchAddress("F0", &m_F0);
	t_tree->SetBranchAddress("FL", &m_FL);
	t_tree->SetBranchAddress("FR", &m_FR);
  if(bTag=="2incl"){

    if(wmode=="1W"){
      t_tree->SetBranchAddress("Wjets_4incl_2incl_el_nom", &m_Wjets_4incl_2incl_el_nom);
      t_tree->SetBranchAddress("Wjets_4incl_2incl_mu_nom", &m_Wjets_4incl_2incl_mu_nom);
    }
    else
    {
      t_tree->SetBranchAddress("WLight_4incl_2incl_el_nom", &m_WLight_4incl_2incl_el_nom);
      t_tree->SetBranchAddress("WLight_4incl_2incl_mu_nom", &m_WLight_4incl_2incl_mu_nom);
      t_tree->SetBranchAddress("Wc_4incl_2incl_el_nom", &m_Wc_4incl_2incl_el_nom);
      t_tree->SetBranchAddress("Wc_4incl_2incl_mu_nom", &m_Wc_4incl_2incl_mu_nom);
      t_tree->SetBranchAddress("Wbbcc_4incl_2incl_el_nom", &m_Wbbcc_4incl_2incl_el_nom);
      t_tree->SetBranchAddress("Wbbcc_4incl_2incl_mu_nom", &m_Wbbcc_4incl_2incl_mu_nom); 
    }
    
    t_tree->SetBranchAddress("QCD_4incl_2incl_el_nom", &m_QCD_4incl_2incl_el_nom);
    t_tree->SetBranchAddress("QCD_4incl_2incl_mu_nom", &m_QCD_4incl_2incl_mu_nom);
    t_tree->SetBranchAddress("RemBkg_4incl_2incl_el_nom", &m_RemBkg_4incl_2incl_el_nom);
    t_tree->SetBranchAddress("RemBkg_4incl_2incl_mu_nom", &m_RemBkg_4incl_2incl_mu_nom);
  
  	
    TH1D* h_FL =new TH1D("FL", "FL", 100, 0.25,    0.35); h_FL->SetLineWidth(2); h_FL->SetXTitle("F_{L}"); h_FL->SetYTitle("# Pseu. Exp.");
    TH1D* h_F0 =new TH1D("F0", "F0", 100, 0.60,    0.8); h_F0->SetLineWidth(2); h_F0->SetXTitle("F_{0}"); h_F0->SetYTitle("# Pseu. Exp.");
    TH1D* h_FR =new TH1D("FR", "FR", 100, -0.04,    0.04); h_FR->SetLineWidth(2); h_FR->SetXTitle("F_{R}"); h_FR->SetYTitle("# Pseu. Exp.");
    
    h_FL->Sumw2();
    h_F0->Sumw2();
    h_FR->Sumw2();
    
    std::cout<<"t_tree->GetEntries()= " << t_tree->GetEntries() << std::endl;
    gStyle->SetOptFit();
    
    /*
      t_tree->Draw("F0>>+F0","");
      h_F0->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/F0_"+bTag+"_"+channel+".png").c_str());
      
      t_tree->Draw("FL>>+FL","");
      h_FL->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/FL_"+bTag+"_"+channel+".png").c_str());
      
      t_tree->Draw("FR>>+FR","");
      h_FR->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/FR_"+bTag+"_"+channel+".png").c_str());
      
    */
    std::cout<< channel+"_"+bTag <<"\n"<<std::endl;
    t_tree->GetEntry(0);
    
    std::cout<< "\nN0"<< "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","N0");
    std::cout<< "\nNL" << "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","NL");
    std::cout<< "\nNR" << "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","NR");
  
    if(wmode=="1W"){
      std::cout<< "\nWjets_4incl_2incl_el_nom: "<< m_Wjets_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_2incl_el");
      t_tree->Fit("fgaus","Wjets_4incl_2incl_el","","Q");
      addYieldLine(myfile, "Wjets_4incl_2incl_el", m_Wjets_4incl_2incl_el_nom, fgaus->GetParameter(1) );
    
      std::cout<< "\nWjets_4incl_2incl_mu_nom: "<< m_Wjets_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_2incl_mu");  
      t_tree->Fit("fgaus","Wjets_4incl_2incl_mu","","Q");
      addYieldLine(myfile, "Wjets_4incl_2incl_mu", m_Wjets_4incl_2incl_mu_nom, fgaus->GetParameter(1) );
    }
    else
    {
      std::cout<< "\nWLight_4incl_2incl_el_nom: "<< m_WLight_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","WLight_4incl_2incl_el");
      t_tree->Fit("fgaus","WLight_4incl_2incl_el","","Q");
      addYieldLine(myfile, "WLight_4incl_2incl_el", m_WLight_4incl_2incl_el_nom, fgaus->GetParameter(1) );

      std::cout<< "\nWLight_4incl_2incl_mu_nom: "<< m_WLight_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","WLight_4incl_2incl_mu");
      t_tree->Fit("fgaus","WLight_4incl_2incl_mu","","Q");
      addYieldLine(myfile, "WLight_4incl_2incl_mu", m_WLight_4incl_2incl_mu_nom, fgaus->GetParameter(1) );

      std::cout<< "\nWc_4incl_2incl_el_nom: "<< m_Wc_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wc_4incl_2incl_el");
      t_tree->Fit("fgaus","Wc_4incl_2incl_el","","Q");
      addYieldLine(myfile, "Wc_4incl_2incl_el", m_Wc_4incl_2incl_el_nom, fgaus->GetParameter(1) );

      std::cout<< "\nWc_4incl_2incl_mu_nom: "<< m_Wc_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wc_4incl_2incl_mu");   
      t_tree->Fit("fgaus","Wc_4incl_2incl_mu","","Q");
      addYieldLine(myfile, "Wc_4incl_2incl_mu", m_Wc_4incl_2incl_mu_nom, fgaus->GetParameter(1) );
  
      std::cout<< "\nWbbcc_4incl_2incl_el_nom: "<< m_Wbbcc_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wbbcc_4incl_2incl_el");
      t_tree->Fit("fgaus","Wbbcc_4incl_2incl_el","","Q");
      addYieldLine(myfile, "Wbbcc_4incl_2incl_el", m_Wbbcc_4incl_2incl_el_nom, fgaus->GetParameter(1) );

      std::cout<< "\nWbbcc_4incl_2incl_mu_nom: "<< m_Wbbcc_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wbbcc_4incl_2incl_mu");   
      t_tree->Fit("fgaus","Wbbcc_4incl_2incl_mu","","Q");
      addYieldLine(myfile, "Wbbcc_4incl_2incl_mu", m_Wbbcc_4incl_2incl_mu_nom, fgaus->GetParameter(1) );
      
    }
  
    std::cout<< "\nQCD_4incl_2incl_el_nom: "<< m_QCD_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_2incl_el");
    t_tree->Fit("fgaus","QCD_4incl_2incl_el","","Q");
    addYieldLine(myfile, "QCD_4incl_2incl_el", m_QCD_4incl_2incl_el_nom, fgaus->GetParameter(1) );
  
    std::cout<< "\nQCD_4incl_2incl_mu_nom: "<< m_QCD_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_2incl_mu");
    t_tree->Fit("fgaus","QCD_4incl_2incl_mu","","Q");
    addYieldLine(myfile, "QCD_4incl_2incl_mu", m_QCD_4incl_2incl_mu_nom, fgaus->GetParameter(1) );

    std::cout<< "\nRemBkg_4incl_2incl_el_nom: "<< m_RemBkg_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","RemBkg_4incl_2incl_el");
    t_tree->Fit("fgaus","RemBkg_4incl_2incl_el","","Q");
    addYieldLine(myfile, "RemBkg_4incl_2incl_el", m_RemBkg_4incl_2incl_el_nom, fgaus->GetParameter(1) );
  
    std::cout<< "\nRemBkg_4incl_2incl_mu_nom: "<< m_RemBkg_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","RemBkg_4incl_2incl_mu");
    t_tree->Fit("fgaus","RemBkg_4incl_2incl_mu","","Q");
    addYieldLine(myfile, "RemBkg_4incl_2incl_mu", m_RemBkg_4incl_2incl_mu_nom, fgaus->GetParameter(1) );

  }
  if(bTag=="1excl2incl"){
    
    if(wmode=="1W"){
      t_tree->SetBranchAddress("Wjets_4incl_2incl_el_nom", &m_Wjets_4incl_2incl_el_nom);
      t_tree->SetBranchAddress("Wjets_4incl_2incl_mu_nom", &m_Wjets_4incl_2incl_mu_nom);
      t_tree->SetBranchAddress("Wjets_4incl_1excl_el_nom", &m_Wjets_4incl_1excl_el_nom);
      t_tree->SetBranchAddress("Wjets_4incl_1excl_mu_nom", &m_Wjets_4incl_1excl_mu_nom);
      
    }
    else
      {
	t_tree->SetBranchAddress("WLight_4incl_el_nom", &m_WLight_4incl_el_nom);
	t_tree->SetBranchAddress("WLight_4incl_mu_nom", &m_WLight_4incl_mu_nom);
	t_tree->SetBranchAddress("Wc_4incl_el_nom", &m_Wc_4incl_el_nom);
	t_tree->SetBranchAddress("Wc_4incl_mu_nom", &m_Wc_4incl_mu_nom);
	t_tree->SetBranchAddress("Wbbcc_4incl_el_nom", &m_Wbbcc_4incl_el_nom);
	t_tree->SetBranchAddress("Wbbcc_4incl_mu_nom", &m_Wbbcc_4incl_mu_nom); 
      }
    
    t_tree->SetBranchAddress("QCD_4incl_2incl_el_nom", &m_QCD_4incl_2incl_el_nom);
    t_tree->SetBranchAddress("QCD_4incl_2incl_mu_nom", &m_QCD_4incl_2incl_mu_nom);
    t_tree->SetBranchAddress("QCD_4incl_1excl_el_nom", &m_QCD_4incl_1excl_el_nom);
    t_tree->SetBranchAddress("QCD_4incl_1excl_mu_nom", &m_QCD_4incl_1excl_mu_nom);
    t_tree->SetBranchAddress("RemBkg_4incl_el_nom", &m_RemBkg_4incl_el_nom);
    t_tree->SetBranchAddress("RemBkg_4incl_mu_nom", &m_RemBkg_4incl_mu_nom);
  
    TH1D* h_FL =new TH1D("FL", "FL", 100, 0.25,    0.35); h_FL->SetLineWidth(2); h_FL->SetXTitle("F_{L}"); h_FL->SetYTitle("# Pseu. Exp.");
    TH1D* h_F0 =new TH1D("F0", "F0", 100, 0.60,    0.8); h_F0->SetLineWidth(2); h_F0->SetXTitle("F_{0}"); h_F0->SetYTitle("# Pseu. Exp.");
    TH1D* h_FR =new TH1D("FR", "FR", 100, -0.04,    0.04); h_FR->SetLineWidth(2); h_FR->SetXTitle("F_{R}"); h_FR->SetYTitle("# Pseu. Exp.");
    h_FL->Sumw2();
    h_F0->Sumw2();
    h_FR->Sumw2();
    
    std::cout<<"t_tree->GetEntries()= " << t_tree->GetEntries() << std::endl;
    gStyle->SetOptFit();
    /*
      t_tree->Draw("F0>>+F0","");
      h_F0->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/F0_"+bTag+"_"+channel+".png").c_str());
      t_tree->Draw("FL>>+FL","");
      h_FL->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/FL_"+bTag+"_"+channel+".png").c_str());
      
      t_tree->Draw("FR>>+FR","");
      h_FR->Fit("gaus");
      addHistoLabel(channel);
      c1->SaveAs((OutDir+"/FR_"+bTag+"_"+channel+".png").c_str());
    */
    std::cout<< channel+"_"+bTag <<"\n"<<std::endl;
    t_tree->GetEntry(0);
    std::cout<< "\nN0"<< "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","N0");
    std::cout<< "\nNL" << "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","NL");
    std::cout<< "\nNR" << "\n----------------------------------------------"<<std::endl;
    t_tree->Fit("gaus","NR");
    if(wmode=="1W"){
      std::cout<< "\nWjets_4incl_2incl_el_nom: "<< m_Wjets_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_2incl_el");
      t_tree->Fit("fgaus","Wjets_4incl_2incl_el","","Q");
      addYieldLine(myfile, "Wjets_4incl_2incl_el", m_Wjets_4incl_2incl_el_nom, fgaus->GetParameter(1) );
      
      std::cout<< "\nWjets_4incl_2incl_mu_nom: "<< m_Wjets_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_2incl_mu","","Q");  
      t_tree->Fit("fgaus","Wjets_4incl_2incl_mu","","Q");  
      addYieldLine(myfile, "Wjets_4incl_2incl_mu", m_Wjets_4incl_2incl_mu_nom, fgaus->GetParameter(1) );
      
      std::cout<< "\nWjets_4incl_1excl_el_nom: "<< m_Wjets_4incl_1excl_el_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_1excl_el","","Q");
      t_tree->Fit("fgaus","Wjets_4incl_1excl_el","","Q");
      addYieldLine(myfile, "Wjets_4incl_1excl_el", m_Wjets_4incl_1excl_el_nom, fgaus->GetParameter(1) );
      
      std::cout<< "\nWjets_4incl_1excl_mu_nom: "<< m_Wjets_4incl_1excl_mu_nom << "\n----------------------------------------------"<<std::endl;
      //t_tree->Fit("gaus","Wjets_4incl_1excl_mu","","Q");  
      t_tree->Fit("fgaus","Wjets_4incl_1excl_mu","","Q");  
      addYieldLine(myfile, "Wjets_4incl_1excl_mu", m_Wjets_4incl_1excl_mu_nom, fgaus->GetParameter(1) );
    }
    else
      {
	std::cout<< "\nWLight_4incl_el_nom: "<< m_WLight_4incl_el_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","WLight_4incl_el","","Q");
	t_tree->Fit("fgaus","WLight_4incl_el","","Q");
	addYieldLine(myfile, "WLight_4incl_el", m_WLight_4incl_el_nom, fgaus->GetParameter(1) );

	std::cout<< "\nWLight_4incl_mu_nom: "<< m_WLight_4incl_mu_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","WLight_4incl_mu","","Q");
	t_tree->Fit("fgaus","WLight_4incl_mu","","Q");
	addYieldLine(myfile, "WLight_4incl_mu", m_WLight_4incl_mu_nom, fgaus->GetParameter(1) );
		
	std::cout<< "\nWc_4incl_el_nom: "<< m_Wc_4incl_el_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","Wc_4incl_el","","Q");
	t_tree->Fit("fgaus","Wc_4incl_el","","Q");
	addYieldLine(myfile, "Wc_4incl_el", m_Wc_4incl_el_nom, fgaus->GetParameter(1) );
  
	std::cout<< "\nWc_4incl_mu_nom: "<< m_Wc_4incl_mu_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","Wc_4incl_mu","","Q");   
	t_tree->Fit("fgaus","Wc_4incl_mu","","Q");   
	addYieldLine(myfile, "Wc_4incl_mu", m_Wc_4incl_mu_nom, fgaus->GetParameter(1) );

	std::cout<< "\nWbbcc_4incl_el_nom: "<< m_Wbbcc_4incl_el_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","Wbbcc_4incl_el","","Q");
	t_tree->Fit("fgaus","Wbbcc_4incl_el","","Q");
	addYieldLine(myfile, "Wbbcc_4incl_el", m_Wbbcc_4incl_el_nom, fgaus->GetParameter(1) );

	std::cout<< "\nWbbcc_4incl_mu_nom: "<< m_Wbbcc_4incl_mu_nom << "\n----------------------------------------------"<<std::endl;
	//t_tree->Fit("gaus","Wbbcc_4incl_mu","","Q");   
	t_tree->Fit("fgaus","Wbbcc_4incl_mu","","Q");   
	addYieldLine(myfile, "Wbbcc_4incl_mu", m_Wbbcc_4incl_mu_nom, fgaus->GetParameter(1) );

      }
    
    std::cout<< "\nQCD_4incl_2incl_el_nom: "<< m_QCD_4incl_2incl_el_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_2incl_el","","Q");
    t_tree->Fit("fgaus","QCD_4incl_2incl_el","","Q");
    addYieldLine(myfile, "QCD_4incl_2incl_el", m_QCD_4incl_2incl_el_nom, fgaus->GetParameter(1) );

    std::cout<< "\nQCD_4incl_2incl_mu_nom: "<< m_QCD_4incl_2incl_mu_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_2incl_mu","","Q");
    t_tree->Fit("fgaus","QCD_4incl_2incl_mu","","Q");
    addYieldLine(myfile, "QCD_4incl_2incl_mu", m_QCD_4incl_2incl_mu_nom, fgaus->GetParameter(1) );

    std::cout<< "\nQCD_4incl_1excl_el_nom: "<< m_QCD_4incl_1excl_el_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_1excl_el","","Q");
    t_tree->Fit("fgaus","QCD_4incl_1excl_el","","Q");
    addYieldLine(myfile, "QCD_4incl_1excl_el", m_QCD_4incl_1excl_el_nom, fgaus->GetParameter(1) );

    std::cout<< "\nQCD_4incl_1excl_mu_nom: "<< m_QCD_4incl_1excl_mu_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","QCD_4incl_1excl_mu","","Q");
    t_tree->Fit("fgaus","QCD_4incl_1excl_mu","","Q");
    addYieldLine(myfile, "QCD_4incl_1excl_mu", m_QCD_4incl_1excl_mu_nom, fgaus->GetParameter(1) );

    std::cout<< "\nRemBkg_4incl_el_nom: "<< m_RemBkg_4incl_el_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","RemBkg_4incl_el","","Q");
    t_tree->Fit("fgaus","RemBkg_4incl_el","","Q");
    addYieldLine(myfile, "RemBkg_4incl_el", m_RemBkg_4incl_el_nom, fgaus->GetParameter(1) );

    std::cout<< "\nRemBkg_4incl_mu_nom: "<< m_RemBkg_4incl_mu_nom << "\n----------------------------------------------"<<std::endl;
    //t_tree->Fit("gaus","RemBkg_4incl_mu","","Q");
    t_tree->Fit("fgaus","RemBkg_4incl_mu","","Q");
    addYieldLine(myfile, "RemBkg_4incl_mu", m_RemBkg_4incl_mu_nom, fgaus->GetParameter(1) );

  }
  
  
  
  std::cout<< "\n----------------END------------------"<<std::endl; 
  

  //myfile.close();
  
}

void addYieldLine(ofstream& file, std::string param, Double_t nom, Double_t fit)
{
  double diff= (1. - nom/fit)*100;
  file<<param<<setprecision(5)<<"\t\t|\t"<<nom<<"\t\t|\t"<<fit<<"\t\t|\t"<<setprecision(3)<<diff<<"\t\t|"<<endl;

  return;
}

void addHistoLabel(std::string channel)
{
    // leg = new TLegend(0.70,0.7,0.89,0.89,"");
    // leg-> SetBorderSize(0);
    // leg->AddEntry(h1,h1_name.c_str(),"l");
    // leg->AddEntry(h2,h2_name.c_str(),"l");
    // leg->AddEntry(h3,h3_name.c_str(),"l");
    // leg->AddEntry(h4,h4_name.c_str(),"l");
    // leg->Draw();

    std::string Ssqrt ="   (#sqrt{s}=8 TeV)";

    std::string m_channel = channel;
    if(channel=="mu") channel = "#mu";

    TLatex l1;
    l1.SetTextAlign(9);
    l1.SetTextFont(72);
    l1.SetTextSize(0.035);
    l1.SetNDC();
    l1.DrawLatex(0.16, 0.840, "ATLAS simulation");
    TLatex l2;
    l2.SetTextAlign(9);
    l2.SetTextSize(0.025);
    l2.SetNDC();
    l2.DrawLatex(0.16, 0.800, "work in progress");

    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextSize(0.025);
    l3.SetNDC();
    l3.DrawLatex(0.16, 0.760, (channel+"+jets"+Ssqrt).c_str());

}



