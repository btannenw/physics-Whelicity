#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <TPad.h>
#include <TLegend.h>

void calcStatError(string name, string upPS, string downPS)
{

  bool dumpAll=false;

  TCanvas *c1= new TCanvas("c1","c_1",200,10,700,500);

  //TFile* f0=new TFile("../plotMaker/new_22Jan2016_LHcut/Leptonic/Templates_110404_syst_2incl_KLF5jOPT_el_mu.root","READ");
  //TFile* f0=new TFile("/afs/cern.ch/work/m/mkareem/Wpol_mc/templateMaker/TemplateFiles/systTemplates/systTemplates_110404/new_24Jan2016_LHcut/Hadronic/Templates_110404_syst_2incl_KLF5jOPT_el_mu.root","READ");
  TFile* f0=new TFile("/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/new_25Feb2016_LHcut/LepHad/Templates_110404_syst_1excl2incl_KLF5jOPT_el_mu_lephad_bTag.root","READ");
  TH1D* nom=(TH1D*)f0->Get("PseudoData");
  //TH1D* varup=(TH1D*)f0->Get("PseudoData_jes_up_FlavourComp");
  //TH1D* vardown=(TH1D*)f0->Get("PseudoData_jes_down_FlavourComp");
  TH1D* varup=(TH1D*)f0->Get(upPS.c_str());
  TH1D* vardown=(TH1D*)f0->Get(downPS.c_str());

  TH1D* nomstat=(TH1D*)nom->Clone();
  TH1D* nomstat2=(TH1D*)nom->Clone();
  TH1D* updiff=(TH1D*)varup->Clone();
  TH1D* downdiff=(TH1D*)vardown->Clone();
  
  double sum_nom = 0;
  double sum_varup = 0;
  double sum_vardown = 0;
  double err_nom = 0;
  double diffup = 0;
  double diffdown = 0;
  double offbins = 0;

  for (int i = 1; i<=nom->GetNbinsX(); i++){
    sum_nom+= nom->GetBinContent(i);
    sum_varup+= varup->GetBinContent(i);
    sum_vardown+= vardown->GetBinContent(i);
    err_nom+= nom->GetBinError(i)*nom->GetBinError(i);
    
    //nomstat->SetBinContent(i,nom->GetBinError(i));
    //nomstat2->SetBinContent(i,-nom->GetBinError(i));
    nomstat->SetBinContent(i,nom->GetBinError(i)/nom->GetBinContent(i));
    nomstat2->SetBinContent(i,-nom->GetBinError(i)/nom->GetBinContent(i));
    
    updiff->SetBinContent(i,(varup->GetBinContent(i)-nom->GetBinContent(i))/nom->GetBinContent(i) ) ;
    downdiff->SetBinContent(i,(vardown->GetBinContent(i)-nom->GetBinContent(i))/nom->GetBinContent(i) );
    
    //updiff->SetBinContent(i,varup->GetBinContent(i)-nom->GetBinContent(i) );
    //downdiff->SetBinContent(i,vardown->GetBinContent(i)-nom->GetBinContent(i) );
  }
  diffup = sum_varup - sum_nom;
  diffdown = sum_vardown - sum_nom;
  err_nom= sqrt(err_nom);
  
  for(int i=1; i<= nomstat->GetNbinsX(); i++){
    if( fabs(updiff->GetBinContent(i)) > fabs(nomstat->GetBinContent(i)) )   {offbins++;}
    if( fabs(downdiff->GetBinContent(i)) > fabs(nomstat->GetBinContent(i)) ) {offbins++;}
  }

  if(dumpAll){
    std::cout<< "sum nominal: "<< sum_nom << std::endl
	     << "sum sum_varup: "<< sum_varup << std::endl
	     << "sum sum_vardown: "<< sum_vardown << std::endl
	     << "sum err_nom: "<< err_nom << std::endl
	     << "sum diffup: "<< diffup << std::endl
	     << "sum diffdown: "<< diffdown << std::endl
	     << std::endl;
    if (fabs(diffup)<fabs(err_nom)) std::cout<<"Up variation OK!"<<std::endl;
    else std::cout<<"Up variation SHIT!"<<std::endl;
    if (fabs(diffdown)<fabs(err_nom)) std::cout<<"Down variation OK!"<<std::endl;
    else std::cout<<"Down variation SHIT!"<<std::endl;
    if ((fabs(diffdown)<fabs(err_nom))&&(fabs(diffdown)<fabs(err_nom))) std::cout<<"Both variations OK!"<<std::endl;
    else std::cout<<"At least one variation SHIT!"<<std::endl;
    
    if(offbins>1)
      std::cout<<"Two or more bins outside statistical error! ... SHIT!"<<endl;
    else
      std::cout<<"Fewer than two bins outside statistical error! ... OK!"<<endl;
  }
  else{
    if(offbins>1 || !(fabs(diffdown)<fabs(err_nom))&&(fabs(diffdown)<fabs(err_nom)) )
      std::cout<<"======================================================================\n"
	       <<"================= SHIT! Keep syst variation: "<<name<<" =================\n"
	       <<"======================================================================"<<std::endl;
}
    if()
  c1->cd();
  nom->SetLineColor(kBlack);
  varup->SetLineColor(kRed);
  vardown->SetLineColor(kBlue);
  updiff->SetLineColor(kRed);
  downdiff->SetLineColor(kBlue);
  nomstat->SetLineColor(kGreen);
  nomstat->SetFillColor(kGreen);
  nomstat2->SetLineColor(kGreen);
  nomstat2->SetFillColor(kGreen);
  
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextFont(72);
  l1.SetTextSize(0.05);
  l1.SetNDC();


  TLegend* leg1=new TLegend(.77,.65,.95,.95);
  //leg1->SetLineColor(0);
  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);
  leg1->AddEntry(nom,"Nominal","l");
  leg1->AddEntry(varup,"Up Variation","l");
  leg1->AddEntry(vardown,"Down Variation","l");
  leg1->AddEntry(nomstat,"Nominal Stat. Unc.","f");

  TPad *pad1 =new TPad("pad1","pad1",0,0.2,1,1);
  pad1->Draw();
  pad1->cd();
  nom->SetXTitle("");
  nom->GetYaxis()->SetTitleOffset(1);
  //nom->SetYTitle("Entries / 0.067");
  nom->SetYTitle("Entries / 0.067");
  //nom->GetXaxis()->SetLabelOffset(999);
  nom->GetXaxis()->SetLabelSize(0);
  nom->GetYaxis()->SetRangeUser(0,1.3*nom->GetMaximum());
  nom->Draw("HIST");
  varup->Draw("HIST same");
  vardown->Draw("HIST same");
  leg1->Draw();
  l1.DrawLatex(0.4, 0.15, name.c_str());

  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextFont(72);
  l2.SetTextSize(0.035);
  l2.SetNDC();
  l2.DrawLatex(0.15, 0.85, "ATLAS simulation");
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextFont(72);
  l3.SetTextSize(0.030);
  l3.SetNDC();
  l3.DrawLatex(0.15, 0.80, "work in progress");
  l3.DrawLatex(0.15, 0.75, "#sqrt{s}=8 TeV, e+#mu combination");

  c1->cd();
  TPad *pad2=new TPad("pad2","pad2",0,0,1,0.2*1.55);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.5);
  pad2->Draw();
  pad2->cd();

  if( nomstat->GetMaximum() < updiff->GetMaximum() )
    nomstat->SetMaximum(1.3*updiff->GetMaximum());
  else
    nomstat->SetMaximum(1.3*nomstat->GetMaximum());

  if( nomstat2->GetMinimum() > downdiff->GetMinimum() )
    nomstat->SetMinimum(1.3*downdiff->GetMinimum());
  else
    nomstat->SetMinimum(1.3*nomstat2->GetMinimum());

  nomstat->GetYaxis()->SetNdivisions(5);
  nomstat->GetYaxis()->SetTitleSize(0.12);
  nomstat->GetYaxis()->SetTitleOffset(0.40);
  nomstat->GetYaxis()->SetLabelSize(0.1);
  nomstat->GetYaxis()->SetTitle("Rel. Diff.");
  nomstat->GetXaxis()->SetNdivisions(507);
  nomstat->GetXaxis()->SetLabelSize(0.14);
  nomstat->GetXaxis()->SetTitleSize(0.16);
  nomstat->GetXaxis()->SetTitleOffset(1.2);

  nomstat->SetXTitle("cos #theta^{*}");
  nomstat->Draw("HIST");  
  nomstat2->Draw("HIST same");  
  updiff->Draw("HIST same");
  downdiff->Draw("HIST same");

  c1->Print(("./feb26_testSignificantPlots/"+name+"_systVsStat_Eval.eps").c_str());

  //c1->cd();
  //h_ll->SetXTitle("Log Likelihood");
  //h_ll->SetYTitle("Entries / 5");
  //h_ll->Draw();

  return;
}

void calcModellingError(string name, string upPS, string downPS)
{
  std::cout<<"%%%%%%%%%%%%%%%%%%%% HI &&&&&&&&&&&&&&&"<<std::endl;

  return;
}
