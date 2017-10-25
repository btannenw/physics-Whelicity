#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
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
void MassFit_plot(std::string Frac);

void FracFit_all()
{
  //std::string workDir = "/afs/cern.ch/work/m/mkareem/FitPackage_allinOneCF/output_CF_2CF";
  //std::string workDir = "/afs/cern.ch/work/m/mkareem/FitPackage_allinOne";
  std::string workDir = "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr26_fixWlight/";

  //std::string fitMod="Systematics"; // Datafit , Systematics, PE
  //std::string fitMod="PE"; // Datafit , Systematics, PE
  std::string fitMod="Datafit"; // Datafit , Systematics, PE

  std::string angleMode= "lephad"; // lep, had, lephad
  std::string bTag="1excl2incl"; // 2incl , 1excl2incl
  std::string channel="el_mu_lephad_bTag"; //el_mu, el_mu_lephad, el_mu_bTag
  std::string wmode="3W"; // 1W , 3W
  std::string bkgMode= ""; // fixedBkg_ or "" ==> "" means float
  
  
  //gROOT->SetBatch();
  //gStyle->SetOptFit();
  Double_t m_F0, m_FL, m_FR;
  Double_t m_N0, m_NL, m_NR;
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

  int nbins;
  std::vector<std::string> labels;
  std::vector<double> normVec;
  std::vector<double> fittedVec;
  

  if(fitMod=="Datafit")
    TFile *f_file = TFile::Open((workDir+"/External"+fitMod+"Output_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"/Syst_Data_"+channel+".root").c_str(),"READ");
  else if(fitMod=="PE")
    TFile *f_file = TFile::Open((workDir+"/ExternalSystematicsOutput_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"/Syst_Nominal_"+channel+".root").c_str(),"READ");
  else if(fitMod=="Systematics")
    TFile *f_file = TFile::Open((workDir+"/External"+fitMod+"Output_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"/Syst_mu=2 hdamp=mtop_"+channel+".root").c_str(),"READ");
    //TFile *f_file = TFile::Open((workDir+"/External"+fitMod+"Output_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"/Syst_mu=0.5 hdamp=2mtop_"+channel+".root").c_str(),"READ");
  
  TTree *t_tree= (TTree*) f_file->Get("EnsembleTree");
  std::cout<<"t_tree->GetEntries()= " << t_tree->GetEntries() << std::endl;
  
  TObjArray *branchList = t_tree->GetListOfBranches()->Clone();
  std::cout<<"branchList->GetEntries()= " << branchList->GetEntries() << std::endl;
  std::vector<double> var(branchList->GetEntries(),0.0);
  
  ofstream myfile;
  std::string fname=(workDir+"/"+fitMod+"_"+channel+"_"+bTag+"_3D_"+bkgMode+wmode+"_"+angleMode+"_YIELDS.txt").c_str();
  cout<<fname<<endl;
  myfile.open ( fname.c_str() ,ios::out);
  
  for(int i=0; i<branchList->GetEntries(); i++){
      t_tree->SetBranchAddress(branchList->At(i)->GetName(), &var.at(i));
      std::cout<<branchList->At(i)->GetName()<< std::endl;
  }
    

  t_tree->GetEntry(0);

  myfile<<"N_total= " << var.at(0)+var.at(3)+var.at(6) << " ± "<< sqrt(var.at(1)*var.at(1)+var.at(4)*var.at(4)+var.at(7)*var.at(7)) <<std::endl; 
  myfile<<std::setw(20)<<"Parameter"<< "\t|\t"<<std::setw(10)<<"pre-fit \t|\t"<<std::setw(10)<<"Post-fit \t\t|\t" <<std::setw(15)<<"% Diff \t|\t"<< std::endl;
  myfile<<"========================================================================================================="<<endl;
  for(int i=0; i<branchList->GetEntries(); i+=3)
    addYieldLine(myfile, branchList->At(i)->GetName(), var.at(i), var.at(i+1), var.at(i+2));

  
   //MassFit_plot("F0");
  // MassFit_plot("FL");
  // MassFit_plot("FR");

  std::cout<< "\n----------------DONE------------------"<<std::endl; 
  
}


void addYieldLine(ofstream& file, std::string param, double fit, double err, double nom)
{
  double diff= (1. - nom/fit)*100;
  
  file<< std::setw(20)<<param<<"\t|\t"<<std::setw(10)<<setprecision(5)<<nom<< "\t|\t"<<std::setw(10)<<fit<<" ± "<< setprecision(3)<< err << "\t|\t"<<setprecision(3)<<diff<<"\t\t|"<<endl;
  //std::cout<< std::setw(20)<<param<<"\t|\t"<<std::setw(10)<<setprecision(5)<<nom<< "\t|\t"<<std::setw(10)<<fit<<" ± "<< setprecision(3)<< err << "\t|\t"<<setprecision(3)<<diff<<"\t\t|"<<endl;

  return;
}

void PlotDiffHisto(int nbins)
{
  TH1F *h1 = new TH1F("","",nbins,0,nbins);
  TH1F *h2 = h1->Clone();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->cd();

}

void addHistoLabel()
{
    TLatex l1;
    l1.SetTextAlign(9);
    l1.SetTextFont(72);
    l1.SetTextSize(0.035);
    l1.SetNDC();
    l1.DrawLatex(0.12, 0.85, "ATLAS simulation");
    TLatex l2;
    l2.SetTextAlign(9);
    l2.SetTextSize(0.033);
    l2.SetNDC();
    l2.DrawLatex(0.12, 0.80, "work in progress");
    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextSize(0.030);
    l3.SetNDC();
    l3.DrawLatex(0.12, 0.75, "e+#mu Lep+had combined");

    TLatex l4;
    l4.SetTextAlign(9);
    l4.SetTextSize(0.028);
    l4.SetNDC();
    l4.DrawLatex(0.12, 0.70, "#sqrt{s}=8 TeV");

}

void MassFit_plot(std::string Frac){

  std::string m_frac= Frac;
   
  TCanvas *c = new TCanvas("c", "c", 1000, 1000);
  c->SetGrid();
  
  TH1 *frame;
  if(m_frac == "F0"){
    frame= c->DrawFrame(165,0.65,180,0.75);
    frame->GetYaxis()->SetTitle("F_{0}");
  }
    
  else if(m_frac == "FL"){
    frame= c->DrawFrame(165,0.25,180,0.35);
    frame->GetYaxis()->SetTitle("F_{L}");
  }

  else if(m_frac == "FR"){
    frame= c->DrawFrame(165,-0.02,180,0.02);
    frame->GetYaxis()->SetTitle("F_{R}");
  }
  frame->GetXaxis()->SetTitle("m_{top-quark} [GeV]");
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->SetTitleOffset(1.1);
  
  //the top-quark mass world avarage 173.34±0.27(stat)±0.71(syst) 
  Double_t m0 = 173.34;
  Double_t m0_unc = sqrt(0.27*0.27+0.71*0.71); 
  
  
   const int n =2;
   TBox *bv = new TBox(m0-m0_unc, frame->GetMinimum() , m0+m0_unc, frame->GetMaximum());
   bv->SetFillColor(3); bv->SetFillStyle(3001);
   bv->Draw();

   // Double_t x[n][1] = {{170}, {175} }; // 170,175
   // Double_t ex[n][1] = {{0}, {0} }; // 170,175
   // Double_t y_f0[n][1] = {{0.69163}, {0.68869} }; // 170,175
   // Double_t ey_f0[n][1] = {{0.00901},{0.00872}};
   Double_t x[n] = {170, 175 }; // 170,175
   Double_t ex[n] = {0, 0 }; // 170,175
   Double_t y_f0[n] = {0.70163, 0.68869 }; // 170,175
   Double_t ey_f0[n] = {0.00901,0.00872};

  if(m_frac == "F0"){  
    //gr_m170 = new TGraphErrors(1,x[0],y_f0[0],ex[0],ey_f0[0]);
    //gr_m175 = new TGraphErrors(1,x[1],y_f0[1],ex[1],ey_f0[1]);
    gr = new TGraphErrors(n,x,y_f0,ex,ey_f0);
  }

  //TLine *line = new TLine(x[0][0]-4,y_f0[0][0],x[0][1]+4,y_f0[0][1]);
  //line->SetLineColor(kGray); line->SetLineWidth(3); 
  //line->Draw(); 
    
   //gr_m170->SetMarkerColor(1); gr_m170->SetMarkerStyle(20);  gr_m170->SetMarkerSize(1.5); 
   //gr_m175->SetMarkerColor(2); gr_m175->SetMarkerStyle(20);  gr_m175->SetMarkerSize(1.5);
   gr->SetMarkerStyle(20);  gr->SetMarkerSize(1.5); gr->SetLineWidth(3); 
   

    //gr_m170->SetLineColor(1); gr_m170->SetLineWidth(3); 
    //gr_m175->SetLineColor(2); gr_m175->SetLineWidth(3); 
     
    
   
   //gr_m170->Draw("Psame");
   //gr_m175->Draw("Psame");
   gr->Draw("Psame");
   gr->Fit("pol1");
   

    leg = new TLegend(0.65,0.75,0.84,0.80,"");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(bv,"m_{top} world avarage","f");
    //leg->AddEntry(gr_m170,"m_{top}=170 GeV","l");
    //leg->AddEntry(gr_m175,"m_{top}=175 GeV","l");
    //leg->AddEntry(gr,"Fraction","l");
    
    leg->Draw();

    addHistoLabel();

  // if(m_frac=="F0")
  //   {
  //     c->SaveAs("DataFit_F0_lephad.png");
  //     c->SaveAs("DataFit_F0_lephad.png");
  // }
  // else if(m_frac=="FL")
  // {
  //   c->SaveAs("DataFit_FL_lephad.png");
  //   c->SaveAs("DataFit_FL_lephad.png");
  // }
    
  // else if(m_frac=="FR")
  // {
  //  c->SaveAs("DataFit_FR_lephad.png");
  //  c->SaveAs("DataFit_FR_lephad.png"); 
  // }
  
}

