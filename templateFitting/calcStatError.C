#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <TPad.h>
#include <TLegend.h>
#include "TLine.h"

std:string AngleSide;
std:string channel;
std:string tagmode;
std:string ratioMode;
std:string tempSubDir;

void calcStatError(string name, string upPS, string downPS, string sModelling)
{

  bool isModelling = sModelling == "true" ? true : false;

  
  TCanvas *c1= new TCanvas("c1","c_1",200,10,1300,500);
  //std:string AngleSide = "Leptonic"; // Hadronic , LepHad, Leptonic
  //std:string AngleSide = "Hadronic"; // Hadronic , LepHad, Leptonic
  AngleSide = "LepHad"; // Hadronic , LepHad, Leptonic
  
  channel = "el_mu_lephad_bTag"; //el_mu, el_mu_bTag
  //std:string channel = "el_mu_bTag"; //el_mu, el_mu_bTag
 //std:string channel = "el_mu"; //el_mu, el_mu_bTag
  
  //std:string tagmode = "2incl"; // 1excl, 2incl, 1excl2incl
  // std:string tagmode = "1excl"; // 1excl, 2incl, 1excl2incl
  tagmode = "1excl2incl"; // 1excl, 2incl, 1excl2incl
  
  ratioMode = "rel"; // abs, rel
  //std:string ratioMode = "abs"; // abs, rel

  tempSubDir = "new_26Apr2016_LHcut";
  //std:string tempSubDir = "new_11Mar2016_LHcut";

  //Leptonic
  //TFile*f0 =new TFile( ("/afs/cern.ch/work/m/mkareem/Wpol_mc/templateMaker/TemplateFiles/systTemplates/systTemplates_110404/"+tempSubDir+"/"+AngleSide+"/Templates_110404_syst_"+tagmode+"_KLF5jOPT_"+channel+".root").c_str(),"READ" );
  TFile*f0 =new TFile( ("/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/"+tempSubDir+"/"+AngleSide+"/Templates_110404_syst_"+tagmode+"_KLF5jOPT_"+channel+".root").c_str(),"READ" );
  
  if (isModelling){
    calcModellingError(name, upPS, downPS, f0, c1);
    return;
  }
  
  bool dumpAll=false;
  
  TH1D* nom=(TH1D*)f0->Get("PseudoData3W");
  TH1D* varup=(TH1D*)f0->Get(upPS.c_str());
  TH1D* vardown=(TH1D*)f0->Get(downPS.c_str());

  TH1D* nomstat=(TH1D*)nom->Clone();
  TH1D* nomstat2=(TH1D*)nom->Clone();
  TH1D* updiff=(TH1D*)varup->Clone();
  TH1D* downdiff=(TH1D*)vardown->Clone();
  
  //nom->Rebin(4);
  //varup->Rebin(4);
  //vardown->Rebin(4);
  //nomstat->Rebin(4);
  //nomstat2->Rebin(4);
  //updiff->Rebin(4);
  //downdiff->Rebin(4);
  
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
    
    
    if(ratioMode == "abs")
    {
      updiff->SetBinContent(i,varup->GetBinContent(i)-nom->GetBinContent(i) );
      downdiff->SetBinContent(i,vardown->GetBinContent(i)-nom->GetBinContent(i) );
      nomstat->SetBinContent(i,nom->GetBinError(i));
      nomstat2->SetBinContent(i,-nom->GetBinError(i));
    }
    else if(ratioMode == "rel"){
      updiff->SetBinContent(i,(varup->GetBinContent(i)-nom->GetBinContent(i))/nom->GetBinContent(i) ) ;
      downdiff->SetBinContent(i,(vardown->GetBinContent(i)-nom->GetBinContent(i))/nom->GetBinContent(i) );
      nomstat->SetBinContent(i,nom->GetBinError(i)/nom->GetBinContent(i));
      nomstat2->SetBinContent(i,-nom->GetBinError(i)/nom->GetBinContent(i));
      }
      
  }
  //updiff = (TH1D*)varup->Clone();
  //updiff->Add(nom, -1);
  //updiff->Divide(nom);
  //downdiff = (TH1D*)vardown->Clone();
  //downdiff->Add(nom, -1);
  //downdiff->Divide(nom);

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
  pad1->SetLeftMargin(1);

  nom->SetXTitle("");
  nom->GetYaxis()->SetTitleOffset(1);
  //nom->SetYTitle("Entries / 0.067");
  nom->SetYTitle("Entries / 0.067");
  //nom->GetXaxis()->SetLabelOffset(999);
  nom->GetYaxis()->SetTitleOffset(0.75);
  nom->GetXaxis()->SetLabelSize(0);
  nom->GetYaxis()->SetRangeUser(0,1.3*nom->GetMaximum());
  nom->Draw("HIST");
  varup->Sumw2();
  vardown->Sumw2();
  varup->SetMarkerStyle(1);
  vardown->SetMarkerStyle(1);
  varup->Draw("HIST E same");
  vardown->Draw("HIST E same");
  leg1->Draw();
  l1.DrawLatex(0.4, 0.85, name.c_str());

  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextFont(72);
  l2.SetTextSize(0.042);
  l2.SetNDC();
  l2.DrawLatex(0.125, 0.85, "ATLAS simulation");
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextFont(72);
  l3.SetTextSize(0.040);
  l3.SetNDC();
  l3.DrawLatex(0.125, 0.80, "Internal");
  l3.DrawLatex(0.125, 0.72, "#sqrt{s}=8 TeV, e+#mu combination");

  TLatex l4;
  l4.SetTextAlign(9);
  l4.SetTextFont(72);
  l4.SetTextSize(0.038);
  l4.SetNDC();
  l4.DrawLatex(0.14, 0.24, "lep_e");
  l4.DrawLatex(0.24, 0.24, "had_e");
  l4.DrawLatex(0.35, 0.24, "lep_e");
  l4.DrawLatex(0.45, 0.24, "had_e");
  l4.DrawLatex(0.57, 0.24, "lep_#mu");
  l4.DrawLatex(0.67, 0.24, "had_#mu");
  l4.DrawLatex(0.78, 0.24, "lep_#mu");
  l4.DrawLatex(0.88, 0.24, "had_#mu");
  TLatex l5;
  l5.SetTextAlign(9);
  l5.SetTextFont(72);
  l5.SetTextSize(0.038);
  l5.SetNDC();
  l5.DrawLatex(0.19, 0.20, "1excl.");
  l5.DrawLatex(0.41, 0.20, "2incl.");
  l5.DrawLatex(0.61, 0.20, "1excl.");
  l5.DrawLatex(0.83, 0.20, "2incl.");


  double max = varup->GetBinContent(varup->GetMaximumBin())*1.65;
  TLine line1, line2, line3;
  line1 = TLine(0.0, 0.0, 0.0, max/2.0);
  line1.SetLineStyle(2);
  line1.SetLineWidth(1);
  line1.Draw("SAME");
  line2 = TLine(-4, 0.0, -4.0, max/2.0);
  line2.SetLineStyle(2);
  line2.SetLineWidth(1);
  line2.Draw("SAME");
  line3 = TLine(4, 0.0, 4.0, max/2.0);
  line3.SetLineStyle(2);
  line3.SetLineWidth(1);
  line3.Draw("SAME");
  c1->cd();
  TPad *pad2=new TPad("pad2","pad2",0,0,1,0.2*1.55);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.5);
  pad2->SetLeftMargin(1);
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
  nomstat->GetYaxis()->SetTitleOffset(0.25);
  //nomstat->GetYaxis()->SetTitleOffset(0.40);
  nomstat->GetYaxis()->SetLabelSize(0.1);
  if(ratioMode=="rel")
    nomstat->GetYaxis()->SetTitle("Rel. Diff.");
  else
    nomstat->GetYaxis()->SetTitle("Abs. Diff.");
  nomstat->GetXaxis()->SetNdivisions(507);
  nomstat->GetXaxis()->SetLabelSize(0.14);
  nomstat->GetXaxis()->SetTitleSize(0.16);
  nomstat->GetXaxis()->SetTitleOffset(1.2);

  nomstat->SetXTitle("cos #theta*");
  nomstat->Draw("HIST");  
  nomstat2->Draw("HIST same");
  updiff->SetMarkerStyle(1);
  downdiff->SetMarkerStyle(1);
  updiff->Draw("HIST E same");
  downdiff->Draw("HIST E same");

 std::string outputDir = "calcStatError_"+tempSubDir+"_"+tagmode+"_"+channel+"_"+AngleSide+"_3W_"+ratioMode;
 
 gSystem->Exec( ("mkdir "+outputDir).c_str() );
  
  c1->Print( ("./"+outputDir+"/"+name+"_systVsStat_Eval.eps").c_str() );

  //c1->cd();
  //h_ll->SetXTitle("Log Likelihood");
  //h_ll->SetYTitle("Entries / 5");
  //h_ll->Draw();

  return;
}

void calcModellingError(string name, string upPS, string downPS, TFile* f0, TCanvas *c0)
{
  TH1D* varup=(TH1D*)f0->Get(upPS.c_str());
  TH1D* vardown=(TH1D*)f0->Get(downPS.c_str());

  TH1D* diff=(TH1D*)varup->Clone();
  
  double sum_varup = 0;
  double sum_vardown = 0;
  double sum_diff = 0;
  double offbins = 0;
  bool dumpAll = false;

  varup->Sumw2();
  vardown->Sumw2();
  diff->Sumw2();

  for (int i = 1; i<=varup->GetNbinsX(); i++){
    sum_varup+= varup->GetBinContent(i);
    sum_vardown+= vardown->GetBinContent(i);
    //diff->SetBinContent(i,varup->GetBinContent(i)-vardown->GetBinContent(i));
  }
  diff->Add(vardown,-1);
  sum_diff = sum_varup - sum_vardown;
  
  for(int i=1; i<= varup->GetNbinsX(); i++){
    if( fabs(diff->GetBinContent(i)) > fabs(diff->GetBinError(i)) )   {offbins++;}
    //if( fabs(diff->GetBinContent(i)) > fabs(diff->GetBinError(i)) ) {offbins++;}
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
    if(offbins>1 )//|| !(fabs(diffdown)<fabs(err_nom))&&(fabs(diffdown)<fabs(err_nom)) )
      std::cout<<"======================================================================\n"
         <<"================= SHIT! Keep syst variation: "<<name<<" =================\n"
         <<"======================================================================"<<std::endl;
  }
  
  c1->cd();
  //nom->SetLineColor(kBlack);
  varup->SetLineColor(kRed);
  vardown->SetLineColor(kBlue);
  diff->SetLineColor(kRed);
  //downdiff->SetLineColor(kBlue);
  //nomstat->SetLineColor(kGreen);
  //nomstat->SetFillColor(kGreen);
  //nomstat2->SetLineColor(kGreen);
  //nomstat2->SetFillColor(kGreen);
  
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextFont(72);
  l1.SetTextSize(0.05);
  l1.SetNDC();


  TLegend* leg1=new TLegend(.73,.62,.93,.92);
  leg1->SetLineColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  //leg1->AddEntry(nom,"Nominal","l");
  leg1->AddEntry(varup,upPS.c_str(),"l");
  leg1->AddEntry(vardown,downPS.c_str(),"l");
  //leg1->AddEntry(nomstat,"Nominal Stat. Unc.","f");

  TPad *pad1 =new TPad("pad1","pad1",0,0.2,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(1);

  varup->SetXTitle("");
  varup->SetYTitle("Entries / 0.067");
  varup->GetXaxis()->SetLabelOffset(999);
  varup->GetXaxis()->SetLabelSize(0);
  varup->GetYaxis()->SetTitleOffset(0.75);
  varup->GetYaxis()->SetRangeUser(0,1.3*varup->GetMaximum());
  varup->SetMarkerStyle(1);
  vardown->SetMarkerStyle(1);

  varup->Draw("HIST E");
  vardown->Draw("HIST E same");
  leg1->Draw();
  l1.DrawLatex(0.4, 0.85, name.c_str());

  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextFont(72);
  l2.SetTextSize(0.042);
  l2.SetNDC();
  l2.DrawLatex(0.125, 0.85, "ATLAS simulation");
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextFont(72);
  l3.SetTextSize(0.025);
  l3.SetNDC();
  l3.DrawLatex(0.125, 0.80, "Internal");
  l3.DrawLatex(0.125, 0.72, "#sqrt{s}=8 TeV, e+#mu combination");

  TLatex l4;
  l4.SetTextAlign(9);
  l4.SetTextFont(72);
  l4.SetTextSize(0.038);
  l4.SetNDC();
  l4.DrawLatex(0.14, 0.24, "lep_e");
  l4.DrawLatex(0.24, 0.24, "had_e");
  l4.DrawLatex(0.35, 0.24, "lep_e");
  l4.DrawLatex(0.45, 0.24, "had_e");
  l4.DrawLatex(0.57, 0.24, "lep_#mu");
  l4.DrawLatex(0.67, 0.24, "had_#mu");
  l4.DrawLatex(0.78, 0.24, "lep_#mu");
  l4.DrawLatex(0.88, 0.24, "had_#mu");
  TLatex l5;
  l5.SetTextAlign(9);
  l5.SetTextFont(72);
  l5.SetTextSize(0.038);
  l5.SetNDC();
  l5.DrawLatex(0.19, 0.20, "1excl.");
  l5.DrawLatex(0.41, 0.20, "2incl.");
  l5.DrawLatex(0.61, 0.20, "1excl.");
  l5.DrawLatex(0.83, 0.20, "2incl.");

  double max = varup->GetBinContent(varup->GetMaximumBin())*1.65;
  TLine line1, line2, line3;
  line1 = TLine(0.0, 0.0, 0.0, max/2.0);
  line1.SetLineStyle(2);
  line1.SetLineWidth(1);
  line1.Draw("SAME");
  line2 = TLine(-4, 0.0, -4.0, max/2.0);
  line2.SetLineStyle(2);
  line2.SetLineWidth(1);
  line2.Draw("SAME");
  line3 = TLine(4, 0.0, 4.0, max/2.0);
  line3.SetLineStyle(2);
  line3.SetLineWidth(1);
  line3.Draw("SAME");

  c1->cd();
  TPad *pad2=new TPad("pad2","pad2",0,0,1,0.2*1.55);
  pad2->SetLeftMargin(1);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.5);
  pad2->Draw();
  pad2->cd();

  /*if( nomstat->GetMaximum() < updiff->GetMaximum() )
    nomstat->SetMaximum(1.3*updiff->GetMaximum());
  else
    nomstat->SetMaximum(1.3*nomstat->GetMaximum());

  if( nomstat2->GetMinimum() > downdiff->GetMaximum() )
    nomstat->SetMinimum(1.3*downdiff->GetMinimum());
  else
    nomstat->SetMinimum(1.3*nomstat2->GetMinimum());
  */
  diff->GetYaxis()->SetNdivisions(5);
  diff->GetYaxis()->SetTitleSize(0.12);
  //diff->GetYaxis()->SetTitleOffset(0.55);
  diff->GetYaxis()->SetTitleOffset(0.25);
  diff->GetYaxis()->SetLabelSize(0.1);
  //if(ratioMode=="rel")
  //  diff->GetYaxis()->SetTitle("Rel. Diff.");
  //else
    diff->GetYaxis()->SetTitle("Abs. Diff.");
  diff->GetXaxis()->SetNdivisions(507);
  diff->GetXaxis()->SetTitle("cos #theta*");
  diff->GetXaxis()->SetLabelSize(0.14);
  diff->GetXaxis()->SetTitleSize(0.16);
  diff->GetXaxis()->SetTitleOffset(1.2);
  diff->SetLineColor(kBlack);
  diff->SetMarkerStyle(0);
  diff->Draw("e");  

  TLine *line4=new TLine(diff->GetXaxis()->GetXmin(),1,diff->GetXaxis()->GetXmax(),1);
  line4->SetLineStyle(2);
  line4->Draw("same");

  //nomstat2->Draw("HIST same");  
  //updiff->Draw("HIST same");
  //downdiff->Draw("HIST same");

  std::string outputDir = "calcStatError_"+tempSubDir+"_"+tagmode+"_"+channel+"_"+AngleSide+"_3W_"+ratioMode;
  gSystem->Exec( ("mkdir "+outputDir).c_str() );  
  c1->Print( ("./"+outputDir+"/"+name+"_systVsStat_Eval.eps").c_str() );

  return;
}
