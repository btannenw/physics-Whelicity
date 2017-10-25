/*
 plots.cxx
*/

#include <iostream>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <sstream>
#include <iomanip>

#include <TROOT.h>
#include <TSystem.h>
#include "TPad.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSpline.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TF1.h"

#include "plots.h"
#include "fit.h"

#include "AtlasStyle.h"

 using namespace std;  

///___________________________PLOTS_______________________________



void PlotDataFit(TH1D hist_sum0, TH1D hist_bkg0, TH1D hist_sum_SM0, TGraphErrors datahist0, TH1D hist_unc, std::string input_channel, std::string fOutputFolder)
{
	SetAtlasStyle();

	TH1D         *hist_sum    = &hist_sum0;
	TH1D         *hist_bkg    = &hist_bkg0;
	TH1D         *hist_sum_SM = &hist_sum_SM0;
	TGraphErrors *datahist    = &datahist0;
	TH1D         *hist_ratio  = (TH1D*) hist_sum -> Clone("");
	TH1D         *datahist2   = (TH1D*) hist_sum -> Clone("");

	//--------new ratio plot (added by MJK: 19.08.2016)
	//TH1D* hLineInOne(0);

	int nBins = hist_unc.GetNbinsX();

	TH1D HistSumUnc = *(TH1D*) hist_sum -> Clone("");

	for(int iBin = 1; iBin <= nBins; ++iBin){

	  double x,y,xe,ye;

	  double val       = hist_sum -> GetBinContent(iBin);
	  double err       = hist_unc.GetBinContent(iBin);
	  
	  datahist -> GetPoint(iBin-1, x, y);
	  ye = datahist -> GetErrorY(iBin-1);
	  
	  double ratio     = y/val;
	  double ratio_unc = ratio*sqrt(ye*ye/y/y + err*err/val/val); 

	  //	  std::cout << iBin << "\t" << val << "\t" << y << "\t" << err << "\t" << ye << "\t" << ratio << "\t" << ratio_unc << std::endl;

	  hist_ratio -> SetBinContent(iBin, ratio);
	  
	  hist_ratio -> SetBinError(iBin, ratio_unc);
	  //hist_ratio -> SetBinError(iBin, 0); //---TEST: for new ratio plot (take into account only fit unc. from hLineInOne)

	  HistSumUnc.SetBinContent(iBin, val);
	  HistSumUnc.SetBinError(iBin,   err);

	  hist_sum    -> SetBinError(iBin,  0.0);
	  hist_sum_SM -> SetBinError(iBin,  0.0);
	  hist_bkg    -> SetBinError(iBin,  0.0);

	  double data_x,data_y,data_ye;

	  datahist -> GetPoint(iBin-1, data_x, data_y);

	  data_ye = datahist -> GetErrorY(iBin-1);

	  datahist2->SetBinContent(iBin, data_y);
	  datahist2->SetBinError(iBin, data_ye);

	}
	//--------new ratio plot (added by MJK: 19.08.2016)
	// hLineInOne = (TH1D*)datahist2->Clone();
 //  	hLineInOne->Divide(datahist2);
 //  	for (int i=0; i<datahist2->GetNbinsX(); i++) {
 //  	    //line in 1 with hPrediction unc.
 //      hLineInOne->SetBinContent(i+1,0.);      hLineInOne->SetBinError(i+1,0.); 
 //      hLineInOne->SetBinContent(i+1, 1.);
 //      if (hist_sum->GetBinContent(i+1)!=0.)
 //        hLineInOne->SetBinError(i+1, float(hist_unc.GetBinContent(i+1)*datahist2->GetBinContent(i+1)/TMath::Power(hist_sum->GetBinContent(i+1),2)));  //take into account only hPrediction unc.
 //      else hLineInOne->SetBinContent(i+1, -1000.);
 //    }

	cout << "Plot fit results... " << endl;
	cout << endl;
	
	//TCanvas
	TCanvas * c0;
	if(input_channel == "el_mu")
	  c0 = new TCanvas("c0","c0", 1200, 800);
	else if(input_channel == "el" || input_channel == "mu")
	  c0 = new TCanvas("c0","c0",  600, 600);
	else if(input_channel == "el_mu_lephad_bTag" )
		c0 = new TCanvas("c0","c0",  2000, 800);
	else
	  c0 = new TCanvas("c0","c0", 1800, 800);

	c0->SetFillColor(0);

	TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
	pad1->SetBottomMargin(0.0);
	//pad1->SetRightMargin(-0.1);
	if(input_channel == "el_mu" || input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag")
	  pad1->SetLeftMargin(0.105);
	else
	  pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	//pad1->SetGrid();

	double max = HistSumUnc.GetBinContent(HistSumUnc.GetMaximumBin())*1.75;


	//HistSumUnc.GetXaxis()->SetTitle("cos #theta*");
	HistSumUnc.GetYaxis()->SetTitle("Events");
	//HistSumUnc.GetXaxis()->SetTitleSize(0.055);
	HistSumUnc.GetYaxis()->SetTitleSize(0.062);
	//HistSumUnc.GetXaxis()->SetTitleOffset(0.9);
	HistSumUnc.GetYaxis()->SetTitleOffset(1.);

	if(input_channel == "el_BTag" || input_channel == "mu_BTag"){
	  HistSumUnc.GetYaxis()->SetTitleOffset(0.80);
	  HistSumUnc.GetYaxis()->SetLabelSize(0.060);
	}
	if(input_channel == "el" || input_channel == "mu"){
          HistSumUnc.GetYaxis()->SetTitleOffset(1.025);
          HistSumUnc.GetYaxis()->SetLabelSize(0.050);
        }
	else{
	  HistSumUnc.GetYaxis()->SetTitleOffset(0.85);
	  HistSumUnc.GetYaxis()->SetLabelSize(0.055);
	}

	//	HistSumUnc.GetXaxis()->SetLabelSize(0.042);

	/*	if(input_channel == "el_mu")
	  HistSumUnc.GetYaxis()->SetLabelSize(0.055);
	else
	HistSumUnc.GetYaxis()->SetLabelSize(0.048); */

 	//hist_bkg->SetFillColor(kGreen);
	HistSumUnc.SetStats(kFALSE);
	HistSumUnc.SetFillStyle(0);
	HistSumUnc.SetLineColor(kGray+2);
	if(input_channel == "el_mu" || input_channel == "el_mu_BTag" || input_channel == "el_BTag" || input_channel == "mu_BTag" || input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag")
	  HistSumUnc.SetLineWidth(1.25);
	else
	  HistSumUnc.SetLineWidth(2);

	HistSumUnc.SetLineStyle(2);
	HistSumUnc.SetMarkerStyle(1);
	//HistSumUnc.SetMarkerColor(kGray+2);
	HistSumUnc.SetMinimum(1.0);
	HistSumUnc.SetMaximum(max);
	//	hist_sum_SM->Draw();	


	//HistSumUnc.SetMarkerStyle(1);
	HistSumUnc.SetMarkerColor(kWhite); //-------------- FIXME (for paper only: kGreen to kWhite)
	HistSumUnc.SetLineColor(kWhite);
	HistSumUnc.SetFillColor(kWhite);
	//HistSumUnc.SetFillStyle(3144);
	HistSumUnc.Draw("E2"); //-------------- FIXME (for paper only)

	hist_sum->SetStats(kFALSE);
        hist_sum->SetFillStyle(0);
        //hist_sum->SetFillColor(kOrange);
        hist_sum->SetLineColor(kRed);
	if(input_channel == "el_mu" || input_channel == "el_mu_BTag" || input_channel == "el_BTag" ||input_channel == "mu_BTag" || input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag")
	  hist_sum->SetLineWidth(1.25);
	else
	  hist_sum->SetLineWidth(2.0);
        hist_sum->SetMarkerStyle(1);
        hist_sum->SetMarkerColor(kRed);
        hist_sum->SetMaximum(max); //-------------- FIXME (for paper only)
        hist_sum->Draw("SAME");

	hist_bkg->SetFillStyle(1001);
	hist_bkg->SetFillColor(kAzure+1);
	hist_bkg->SetLineColor(kAzure+1);
	hist_bkg->SetLineWidth(2);
	hist_bkg->SetMarkerStyle(1);
	//hist_bkg->SetMarkerColor(kGreen);
	hist_bkg->Draw("SAME");

	//	hist_sum->Draw("sameaxis");

	// datahist->SetFillStyle(0);
	datahist->SetLineColor(kBlack); 
	datahist->SetLineWidth(2);
	datahist->SetMarkerStyle(20);
	datahist->SetMarkerSize(1.1);
	datahist->SetMarkerColor(kBlack);
	datahist->Draw("PSAME");
	
	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.046); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	if(input_channel == "el_mu" || input_channel == "el_mu_bTag" || input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag"){
	  
	  l.DrawLatex(0.140, 0.855, "ATLAS");
	  l.SetTextFont(42);
	  l.DrawLatex(0.220, 0.855, "Preliminary");
	  
	}
	
	else if(input_channel == "el_BTag" || input_channel == "mu_BTag"){
	  
	  l.SetTextSize(0.05);
	  l.DrawLatex(0.160, 0.855, "ATLAS");
          l.SetTextFont(42);
          l.DrawLatex(0.237, 0.855, "Preliminary");
	  
	}
	else if(input_channel == "el" || input_channel == "mu"){
	  
	  l.DrawLatex(0.160, 0.860, "ATLAS");
          l.SetTextFont(42);  
	  
	  l.DrawLatex(0.275, 0.860, "Preliminary");
	  
        }
        else{
	  
          l.DrawLatex(0.160, 0.855, "ATLAS");
          l.SetTextFont(42);
          l.DrawLatex(0.265, 0.855, "Preliminary");
	  
        }

	TLine line1,line2,line3,line4,line5;
	
	if(nBins == 30){
	  
	  line1 = TLine(0.0, 0.0, 0.0, max/2.0);
	  line1.SetLineStyle(2);
	  
	  line1.Draw("SAME");
	  
	}
	if(nBins == 45){
	  
	  line1 = TLine(-1.0, 0.0, -1.0, max/2.0);
          line1.SetLineStyle(2);
	  
          line1.Draw("SAME");
	  
	  line2 = TLine(1.0,  0.0,  1.0, max/2.0);
          line2.SetLineStyle(2);
	  
          line2.Draw("SAME"); 
	}

	if(nBins == 60){
	  line1 = TLine(0.0, 0.0, 0.0, max/2.0);
      line1.SetLineStyle(2);
      line1.Draw("SAME");
      line2 = TLine(-2, 0.0, -2.0, max/2.0);
      line2.SetLineStyle(2);
      line2.Draw("SAME");
      line3 = TLine(2, 0.0, 2.0, max/2.0);
      line3.SetLineStyle(2);
      line3.Draw("SAME");
	}

	if(nBins == 120){
	  line1 = TLine(0.0, 0.0, 0.0, max/2.0);
      line1.SetLineStyle(2);
      line1.Draw("SAME");
      line2 = TLine(-4, 0.0, -4.0, max/2.0);
      line2.SetLineStyle(2);
      line2.Draw("SAME");
      line3 = TLine(4, 0.0, 4.0, max/2.0);
      line3.SetLineStyle(2);
      line3.Draw("SAME");
	}
	
	if(nBins == 90){
	  
          line1 = TLine(-4.0, 0.0, -4.0, max/2.0);
          line1.SetLineStyle(2);
	  
          line1.Draw("SAME");
	  
          line2 = TLine(-2.0,  0.0, -2.0, max/2.0);
          line2.SetLineStyle(2);
	  
          line2.Draw("SAME");
	  
	  line3 = TLine(0.0,   0.0,  0.0, max/2.0);
          line2.SetLineStyle(2);
	  
          line3.Draw("SAME");
	  
	  line4 = TLine(2.0,  0.0,  2.0, max/2.0);
          line4.SetLineStyle(2);
	  
          line4.Draw("SAME");

	  line5 = TLine(4.0,  0.0,  4.0, max/2.0);
          line5.SetLineStyle(2);
	  
          line5.Draw("SAME");
	  
	  
        }
	
	//TLatex latex0;
	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.040);
	if (input_channel == "el")    l.DrawLatex(0.160,0.81, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.160,0.81, "#mu+jets channel");
	
	if (input_channel == "el_mu"){
	  l.DrawLatex(0.140, 0.800,  "Leptonic analyser");
	  l.DrawLatex(0.140, 0.730,  "#int L dt = 20.2 fb^{-1}, #sqrt{s} = 8 TeV ");
	  
	  l.SetTextFont(62);
	  l.DrawLatex(0.260, 0.580, "e+jets(#geq 2 b-tags)");
	  l.DrawLatex(0.695, 0.580, "#mu+jets(#geq 2 b-tags)");
	}
	if (input_channel == "el_BTag"){
          l.DrawLatex(0.160, 0.795,  "Combined e+jets channel");
          l.DrawLatex(0.215, 0.125, "- 0 b-tag -");
          l.DrawLatex(0.500, 0.125, "- 1 b-tag -");
	  l.DrawLatex(0.760, 0.125, "- #geq 2 b-tags -");
        }
	if (input_channel == "mu_BTag"){
          l.DrawLatex(0.160, 0.795,  "Combined #mu+jets channel");
          l.DrawLatex(0.215, 0.125, "- 0 B-tag -");
          l.DrawLatex(0.500, 0.125, "- 1 B-tag -");
          l.DrawLatex(0.760, 0.125, "- #geq 2 b-tags -");
        }
	if (input_channel == "el_mu_lephad"){
          l.DrawLatex(0.140, 0.795,  "Combined channel");
          l.DrawLatex(0.140, 0.120, "- e+jets(lep.) -");
	  l.DrawLatex(0.350, 0.120, "- #mu+jets(lep.) -");
	  l.DrawLatex(0.560, 0.120, "- e+jets(had.) -");
	  l.DrawLatex(0.760, 0.120, "- #mu+jets(had.) -");
        }
	if (input_channel == "el_mu_bTag"){
          l.DrawLatex(0.140, 0.795,  "Hadronic analyser");
          l.DrawLatex(0.140, 0.730,  "#int L dt = 20.2 fb^{-1}, #sqrt{s} = 8 TeV");
          l.SetTextFont(62);
          l.DrawLatex(0.200, 0.580, "e+jets(1 b-tag)");
	  l.DrawLatex(0.380, 0.580, "e+jets(#geq 2 b-tags)"); 
	  l.DrawLatex(0.600, 0.580, "#mu+jets(1 b-tag)");
	  l.DrawLatex(0.800, 0.580, "#mu+jets(#geq 2 b-tags)");
	}
	
	if (input_channel == "el_mu_lephad_bTag"){
	  l.DrawLatex(0.140, 0.795,  "Combined channel");
	  l.DrawLatex(0.140, 0.730,  "#int L dt = 20.2 fb^{-1}, #sqrt{s} = 8 TeV");
	  l.SetTextFont(62);
	  l.DrawLatex(0.150, 0.550, "e+jets(1 b-tag)");
	  l.DrawLatex(0.360, 0.550, "e+jets(#geq 2 b-tags)");
	  l.DrawLatex(0.570, 0.550, "#mu+jets(1 b-tag)");
	  l.DrawLatex(0.770, 0.550, "#mu+jets(#geq 2 b-tags)");
	}
	
	std::stringstream oss, oss2, oss3;
	std::string chi2_ndf;
	
	double chi2 = datahist2 -> Chi2Test(hist_sum, "CHI2/NDF");
	double prob = datahist2 -> Chi2Test(hist_sum, "UW");
	double prob_k = datahist2 -> KolmogorovTest(hist_sum);


	oss  << std::setprecision(3) << chi2;
	oss2 << std::setprecision(3) << prob;
	oss3 << std::setprecision(3) << prob_k;
	/*
	if(input_channel == "el_BTag" || input_channel == "mu_BTag"){
	 
	  l.DrawLatex(0.160, 0.725, ("#chi^{2}/ndf = "+oss.str()).c_str());
	  l.DrawLatex(0.160, 0.670, ("Prob.  = "+oss2.str()).c_str());

	}
	else if(input_channel == "el" || input_channel == "mu"){

	  l.SetTextSize(0.040);
          l.DrawLatex(0.160, 0.730, ("#chi^{2}/ndf = "+oss.str()).c_str());
          l.DrawLatex(0.160, 0.685, ("Prob.  = "+oss2.str()).c_str());

        }
	else if(input_channel == "el_mu"){

          l.SetTextSize(0.0425);
          l.DrawLatex(0.140, 0.670, ("#chi^{2}/ndf = "+oss.str()).c_str());
          //l.DrawLatex(0.140, 0.675, ("Prob.  = "+oss2.str()).c_str());
          //l.DrawLatex(0.140, 0.675, ("Kolmogorov Prob.  = "+oss3.str()).c_str());

        }
	else if(input_channel == "el_mu_bTag" ||input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag"){
	  
          l.SetTextSize(0.0425);
          l.DrawLatex(0.140, 0.670, ("#chi^{2}/ndf = "+oss.str()).c_str());
          //l.DrawLatex(0.140, 0.675, ("Prob.  = "+oss2.str()).c_str());
          //l.DrawLatex(0.140, 0.675, ("Kolmogorov Prob.  = "+oss3.str()).c_str());
	  
        }
	*/

	TLegend * legend0;
	if(input_channel == "el_mu")
	  legend0 = new TLegend(0.68,0.64,0.92,0.90);
	else if(input_channel == "el_BTag" || input_channel == "mu_BTag")
          legend0 = new TLegend(0.71,0.62,0.94,0.90);
	else
	  legend0 = new TLegend(0.65,0.61,0.93,0.91);
	
	legend0->SetTextFont(42);
	legend0->AddEntry(hist_sum,     "Best Fit",                "l");
	//legend0->AddEntry(&HistSumUnc,  "Best Fit Unc. (Stat.)",   "f"); //--------- FIXME (for paper only)
	legend0->AddEntry(hist_bkg,     "Background",              "f");
	//	legend0->AddEntry(hist_sum_SM,  "SM Expectation",          "l");
	//legend0->AddEntry(datahist,     "Data (#sqrt{s} = 8 TeV)", "p");	
	legend0->AddEntry(datahist,     "Data", "p");	
	legend0->SetFillColor(0);
	legend0->SetBorderSize(0);
	legend0->Draw("same");

	c0->cd();
//---------------------------- Ratio Plot-------------------------------//
	TPad *pad2 = new TPad("pad2","pad2",0,0.00,1,0.22);
	pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.38);
	if(input_channel == "el_mu" || input_channel == "el_mu_lephad")
	  pad2->SetLeftMargin(0.105);
	else if(input_channel == "el_mu_lephad_bTag")
		pad2->SetLeftMargin(0.105);
	else
	  pad2->SetLeftMargin(0.12);
	pad2->Draw();
	pad2->cd();
	
	// //-------new 
	// //hLineInOne->SetMarkerSize(0.);
 //  	//hLineInOne->SetMarkerStyle(0);
 //  	//hLineInOne->SetFillColor(0);
 //  	hLineInOne->SetLineColor(2);
 //  	hLineInOne->SetFillStyle(3256);
 //  	hLineInOne->Draw("L SAME"); //"L E2 SAME");
	//------------
	
	hist_ratio->GetXaxis()->SetTitle("cos #theta*");
	hist_ratio->GetYaxis()->SetTitle("Data/Fit");
	hist_ratio->GetYaxis()->CenterTitle();
	hist_ratio->GetXaxis()->SetTitleSize(0.21);
	if(input_channel == "el_mu" || input_channel == "el_mu_bTag" ||input_channel == "el_mu_lephad" || input_channel == "el_mu_lephad_bTag" ){
	  hist_ratio->GetYaxis()->SetTitleSize(0.15);
	  hist_ratio->GetXaxis()->SetTitleOffset(0.90);
	  hist_ratio->GetXaxis()->SetLabelOffset(1000);
	  hist_ratio->GetYaxis()->SetTitleOffset(0.325);
	  hist_ratio->GetYaxis()->SetLabelSize(0.140);
	}
	else if(input_channel == "el" || input_channel == "mu" ){
	  hist_ratio->GetYaxis()->SetTitleSize(0.15);
	  hist_ratio->GetYaxis()->SetTitleOffset(0.415);
	  hist_ratio->GetXaxis()->SetTitleOffset(0.80);
	  hist_ratio->GetYaxis()->SetLabelSize(0.125);
	}
	else if(input_channel == "el_BTag" || input_channel == "mu_BTag"){
          hist_ratio->GetYaxis()->SetTitleSize(0.15);
          hist_ratio->GetXaxis()->SetTitleOffset(0.90);
	  hist_ratio->GetYaxis()->SetTitleOffset(0.35);
          hist_ratio->GetYaxis()->SetLabelSize(0.165);
	  //	  hist_ratio->GetXaxis()->SetTitleSize(0.2);
        }
	
    
	/*	if(input_channel == "el_mu")
	  hist_ratio->GetYaxis()->SetTitleOffset(0.25);
	else
	hist_ratio->GetYaxis()->SetTitleOffset(0.45); */

	//if(input_channel != "el" && input_channel != "mu")
	 // hist_ratio->GetXaxis()->SetLabelSize(0.00);
	//else
	hist_ratio->GetXaxis()->SetLabelSize(0.165);

	hist_ratio->SetMinimum(0.86);
	hist_ratio->SetMaximum(1.14);
	hist_ratio->SetMarkerStyle(20);
	hist_ratio->SetMarkerColor(kBlack);
	hist_ratio->SetStats(kFALSE);
	hist_ratio->SetLineColor(kWhite);
	hist_ratio->SetLineWidth(0);
	hist_ratio->GetXaxis()->SetNdivisions(505);
	hist_ratio->GetYaxis()->SetNdivisions(505);

	//hist_ratio->DrawCopy("hist");
	hist_ratio->SetFillColor(kOrange);
	hist_ratio->SetFillStyle(3002);
	//hist_ratio->Draw("E2same");
	
	//--------------------------------------------Normal line
	TF1 *norm1 = new TF1("fa1","1", -100.0, 500.0);
        norm1 -> SetLineColor(kRed);
        norm1 -> SetLineStyle(1);
        norm1 -> SetLineWidth(1.4);
        //norm1 -> Draw("Same");

    double min = 0.94;

    ShowRatio(datahist2, hist_sum, hist_unc, min, 1.07);

    TLine line7,line8,line9,line10,line11;
    
    if(nBins == 30){
      
      line7 = TLine(0.0, min, 0.0, 2.0-min);
      line7.SetLineStyle(2);
      
      line7.Draw();
      
      TGaxis *axis2 = new TGaxis(-2.0, min, 0.0, min, -1.0, 1.0, 505, "");
      axis2->SetTextAlign(9);
      axis2->SetLabelSize(0.17);
      axis2->SetLabelFont(42);
      axis2->SetLabelOffset(0.025);
      
      TGaxis *axis3 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
      axis3->SetTextAlign(9);
      axis3->SetLabelSize(0.17);
      axis3->SetLabelFont(42);
      axis3->SetLabelOffset(0.025);
      
      axis2->Draw();  
      axis3->Draw();
      
    }
    
    else if(nBins == 45){
      
      line7 = TLine(-1.0, min, -1.0, 2.0-min);
      line7.SetLineStyle(2);
      
      line7.Draw();
      
      line8 = TLine( 1.0, min,  1.0, 2.0-min);
      line8.SetLineStyle(2);
      
      line8.Draw();
      
      TGaxis *axis2 = new TGaxis(-3.0, min, -1.0, min, -1.0, 1.0, 505, "");
      axis2->SetTextAlign(9);
      axis2->SetLabelSize(0.17);
      axis2->SetLabelFont(42);
      axis2->SetLabelOffset(0.025);
      
      TGaxis *axis3 = new TGaxis(-0.9, min, 1.0, min, -0.9, 1.0, 505, "");
      axis3->SetTextAlign(9);
      axis3->SetLabelSize(0.17);
      axis3->SetLabelFont(42);
      axis3->SetLabelOffset(0.025);
      
      TGaxis *axis4 = new TGaxis(1.1, min, 3.0, min, -0.9, 1.0, 505, "");
      axis4->SetTextAlign(9);
      axis4->SetLabelSize(0.17);
      axis4->SetLabelFont(42);
      axis4->SetLabelOffset(0.025);
      
      axis2->Draw();
      axis3->Draw();
      axis4->Draw();
      
    }
    
    if(nBins == 60){
      
      line7 = TLine(-4.0, min, -4.0, 2.0-min);
      // line7.SetLineStyle(2);
      
      // line7.Draw();
      
      // line8 = TLine(-2.0, min, -2.0, 2.0-min);
      // line8.SetLineStyle(2);
      
      // line8.Draw();
      
      // line9 = TLine(0.0,  min,  0.0, 2.0-min);
      // line9.SetLineStyle(2);
      
      // line9.Draw();
      
      // line10 = TLine(2.0,  min,  2.0, 2.0-min);
      // line10.SetLineStyle(2);
      
      // line10.Draw();
      
      // line11 = TLine(4.0,  min,  4.0, 2.0-min);
      // line11.SetLineStyle(2);
      
      // line11.Draw();
      
      TGaxis *axis3 = new TGaxis(-4.0, min, -2.0, min, -0.9, 1.0, 505, "");
      axis3->SetTextAlign(9);
      axis3->SetLabelSize(0.17);
      axis3->SetLabelFont(42);
      axis3->SetLabelOffset(0.025);
      
      TGaxis *axis4 = new TGaxis(-1.9, min, 0.0, min, -0.9, 1.0, 505, "");
      axis4->SetTextAlign(9);
      axis4->SetLabelSize(0.17);
      axis4->SetLabelFont(42);
      axis4->SetLabelOffset(0.025);
      
      TGaxis *axis5 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
      axis5->SetTextAlign(9);
      axis5->SetLabelSize(0.17);
      axis5->SetLabelFont(42);
      axis5->SetLabelOffset(0.025);
      
      TGaxis *axis6 = new TGaxis(2.1, min, 4.0, min, -0.9, 1.0, 505, "");
      axis6->SetTextAlign(9);
      axis6->SetLabelSize(0.17);
      axis6->SetLabelFont(42);
      axis6->SetLabelOffset(0.025);
      
      
      axis3->Draw();
      axis4->Draw();
      axis5->Draw();
      axis6->Draw();
    }
    else if(nBins == 90){
      
      line7 = TLine(-4.0, min, -4.0, 2.0-min);
      line7.SetLineStyle(2);
      
      line7.Draw();
      
      line8 = TLine(-2.0, min, -2.0, 2.0-min);
      line8.SetLineStyle(2);
      
      line8.Draw();
      
      line9 = TLine(0.0,  min,  0.0, 2.0-min);
      line9.SetLineStyle(2);
      
      line9.Draw();
      
      line10 = TLine(2.0,  min,  2.0, 2.0-min);
      line10.SetLineStyle(2);
      
      line10.Draw();
      
      line11 = TLine(4.0,  min,  4.0, 2.0-min);
      line11.SetLineStyle(2);
      
      line11.Draw();
      
      TGaxis *axis2 = new TGaxis(-6.0, min, -4.0, min, -1.0, 1.0, 505, "");
      axis2->SetTextAlign(9);
      axis2->SetLabelSize(0.17);
      axis2->SetLabelFont(42);
      axis2->SetLabelOffset(0.025);
      
      TGaxis *axis3 = new TGaxis(-3.9, min, -2.0, min, -0.9, 1.0, 505, "");
      axis3->SetTextAlign(9);
      axis3->SetLabelSize(0.17);
      axis3->SetLabelFont(42);
      axis3->SetLabelOffset(0.025);
      
      TGaxis *axis4 = new TGaxis(-1.9, min, 0.0, min, -0.9, 1.0, 505, "");
      axis4->SetTextAlign(9);
      axis4->SetLabelSize(0.17);
      axis4->SetLabelFont(42);
      axis4->SetLabelOffset(0.025);
      
      TGaxis *axis5 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
      axis5->SetTextAlign(9);
      axis5->SetLabelSize(0.17);
      axis5->SetLabelFont(42);
      axis5->SetLabelOffset(0.025);
      
      TGaxis *axis6 = new TGaxis(2.1, min, 4.0, min, -0.9, 1.0, 505, "");
      axis6->SetTextAlign(9);
      axis6->SetLabelSize(0.17);
      axis6->SetLabelFont(42);
      axis6->SetLabelOffset(0.025);
      
      TGaxis *axis7 = new TGaxis(4.1, min, 6.0, min, -0.9, 1.0, 505, "");
      axis7->SetTextAlign(9);
      axis7->SetLabelSize(0.17);
      axis7->SetLabelFont(42);
      axis7->SetLabelOffset(0.025);
      
      axis2->Draw();
      axis3->Draw();
      axis4->Draw();
      axis5->Draw();
      axis6->Draw();
      axis7->Draw();
      
    }
    
    
    if(nBins == 120){
      
      TGaxis *axis7 = new TGaxis(-8.0, min, -6.0, min, -0.9, 1.0, 505, "");
      axis7->SetTextAlign(9);
      axis7->SetLabelSize(0.17);
      axis7->SetLabelFont(42);
      axis7->SetLabelOffset(0.025);
      
      TGaxis *axis8 = new TGaxis(-5.9, min, -4.0, min, -0.9, 1.0, 505, "");
      axis8->SetTextAlign(9);
      axis8->SetLabelSize(0.17);
      axis8->SetLabelFont(42);
      axis8->SetLabelOffset(0.025);
      
      TGaxis *axis3 = new TGaxis(-3.9, min, -2.0, min, -0.9, 1.0, 505, "");
      axis3->SetTextAlign(9);
      axis3->SetLabelSize(0.17);
      axis3->SetLabelFont(42);
      axis3->SetLabelOffset(0.025);
      
      TGaxis *axis4 = new TGaxis(-1.9, min, 0.0, min, -0.9, 1.0, 505, "");
      axis4->SetTextAlign(9);
      axis4->SetLabelSize(0.17);
      axis4->SetLabelFont(42);
      axis4->SetLabelOffset(0.025);
      
      TGaxis *axis5 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
      axis5->SetTextAlign(9);
      axis5->SetLabelSize(0.17);
      axis5->SetLabelFont(42);
      axis5->SetLabelOffset(0.025);
      
      TGaxis *axis6 = new TGaxis(2.1, min, 4.0, min, -0.9, 1.0, 505, "");
      axis6->SetTextAlign(9);
      axis6->SetLabelSize(0.17);
      axis6->SetLabelFont(42);
      axis6->SetLabelOffset(0.025);
      
      //----------
      TGaxis *axis9 = new TGaxis(4.1, min, 6.0, min, -0.9, 1.0, 505, "");
      axis9->SetTextAlign(9);
      axis9->SetLabelSize(0.12);
      axis9->SetLabelFont(42);
      axis9->SetLabelOffset(0.025);
      
      TGaxis *axis10 = new TGaxis(6.1, min, 8.0, min, -0.9, 1.0, 505, "");
      axis10->SetTextAlign(9);
      axis10->SetLabelSize(0.12);
      axis10->SetLabelFont(42);
      axis10->SetLabelOffset(0.025);
      
      axis3->Draw();
      axis4->Draw();
      axis5->Draw();
      axis6->Draw();
      axis7->Draw();
      axis8->Draw();
      axis9->Draw();
      axis10->Draw();
    }
    
    
    std::string printtitle1 = fOutputFolder+"/Datafit_StatUnc_" + input_channel + ".png";
    std::string printtitle2 = fOutputFolder+"/Datafit_StatUnc_" + input_channel + ".gif";
    std::string printtitle3 = fOutputFolder+"/Datafit_StatUnc_" + input_channel + ".eps";
    std::string printtitle4 = fOutputFolder+"/Datafit_StatUnc_" + input_channel + ".pdf";
    
    const char *print1, *print2, *print3, *print4;
    print1=printtitle1.c_str();
    print2=printtitle2.c_str();
    print3=printtitle3.c_str();
    print4=printtitle4.c_str();
    
    c0->Print(print1);
    c0->Print(print3); //for eps
    c0->Print(print4); //for pdf
    
    legend0->Delete();
    pad1->Delete();
    pad2->Delete();
    
}


void PlotDataFitSimple(TH1D hist_sum, TH1D hist_bkg, TH1D hist_sum_SM, TGraphErrors datahist, std::string input_channel, std::string fOutputFolder)
{
  //	SetAtlasStyle();

	cout << "Plot fit results simple... " << endl;
	cout << endl;
	
	//TCanvas
	TCanvas * c0 = new TCanvas("c0","c0", 600, 600);
	c0->SetFillColor(0);
	//c0->SetLogy(1);

	SetAtlasStyle();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1,1);
        pad1->SetBottomMargin(0.12);
        //pad1->SetRightMargin(-0.1);
        pad1->SetLeftMargin(0.135);
	pad1->SetBorderMode(0);
        pad1->Draw();
        pad1->cd();


	double max = hist_sum.GetBinContent(hist_sum.GetMaximumBin())*1.65;

	hist_sum.GetXaxis()->SetTitle("cos #theta*");
	hist_sum.GetYaxis()->SetTitle("Events");
	hist_sum.GetYaxis()->SetTitleSize(0.049);
	hist_sum.GetXaxis()->SetTitleSize(0.052);
	hist_sum.GetXaxis()->SetTitleOffset(1.03);
	hist_sum.GetYaxis()->SetTitleOffset(1.45);
	hist_sum.GetXaxis()->SetLabelSize(0.043);
	hist_sum.GetYaxis()->SetLabelSize(0.043);
	hist_sum.SetMinimum(1.0);
	hist_sum.SetMaximum(max);

	hist_sum.SetStats(kFALSE);
	hist_sum.SetFillStyle(0);
	//hist_sum.SetFillColor(kOrange);
	hist_sum.SetLineColor(kGreen);
	hist_sum.SetLineWidth(2);
	hist_sum.SetMarkerStyle(1);
	hist_sum.SetMarkerColor(kGreen);
	hist_sum.Draw("hist");	


 	//hist_bkg->SetFillColor(kGreen);
	hist_sum_SM.SetStats(kFALSE);
	hist_sum_SM.SetFillStyle(0);
	hist_sum_SM.SetLineColor(kRed);
	hist_sum_SM.SetLineWidth(2);
	hist_sum_SM.SetLineStyle(2);
	hist_sum_SM.SetMarkerStyle(1);
	hist_sum_SM.SetMarkerColor(kRed);
	hist_sum_SM.Draw("histsame");	


	hist_bkg.SetFillStyle(1001);
	hist_bkg.SetFillColor(kAzure+1);
	hist_bkg.SetLineColor(kAzure+1);
	hist_bkg.SetLineWidth(2);
	hist_bkg.SetMarkerStyle(1);
	//hist_bkg.SetMarkerColor(kGreen);
	hist_bkg.Draw("histsame");

	hist_sum.Draw("sameaxis");

	datahist.SetLineColor(kBlack); 
	datahist.SetLineWidth(2);
	datahist.SetMarkerStyle(20);
	datahist.SetMarkerSize(1.2);
	datahist.SetMarkerColor(kBlack);
	datahist.Draw("PSAME");
	
	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.038); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.160,0.87,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.290,0.87,"Preliminary");
	
	//TLatex latex0;
	l.SetNDC(1);
	l.SetTextFont(42);
	// l.DrawLatex(0.245,0.775, "TMinuit");
	// l.DrawLatex(0.155,0.79, "TMinuit");
	l.SetTextSize(0.038);
	if (input_channel == "el")    l.DrawLatex(0.160,0.82, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.160,0.82, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.160,0.82, "Combined channel");

	/*if (input_channel == "el") filename_output = "e+jets channel";
	else if (input_channel == "mu") filename_output = "#mu+jets channel";
	else filename_output = " ";*/
	//const char *filename_output;
	//filename_output=output.c_str();
	
	l.SetTextSize(0.037);
	//l.DrawLatex(0.245, 0.735, "Rebinned Sample");
	//l.DrawLatex(0.155, 0.66, "Profile likelihood fit");
	//l.DrawLatex(0.155, 0.66, "Template fit");

	
	TLegend * legend0 = new TLegend(0.70,0.65,0.92,0.92);
	//legend0->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	legend0->SetTextFont(42);
	legend0->AddEntry(&hist_sum,"Best Fit","l");
	legend0->AddEntry(&hist_bkg,"Background","f");
	legend0->AddEntry(&hist_sum_SM,"SM Expectation","l");
	legend0->AddEntry(&datahist,"Data","p");	
	legend0->SetFillColor(0);
	legend0->SetBorderSize(0);
	legend0->Draw("same");
	
	//	c0->cd();
	
	std::string printtitle1 = fOutputFolder+"/Datafit_" + input_channel + ".png";
	std::string printtitle2 = fOutputFolder+"/Datafit_" + input_channel + ".gif";
	std::string printtitle3 = fOutputFolder+"/Datafit_" + input_channel + ".eps";

	//char filename_output[sizeof output + 100];
	//sprintf(filename_output, output);
	
	const char *print1, *print2, *print3;
	print1=printtitle1.c_str();
	print2=printtitle2.c_str();
	print3=printtitle3.c_str();

	/*char printtitle[100] = "plots/first_test.png";// use %f for double, %d for int
	char print[sizeof printtitle + 100];
	sprintf(print, printtitle);*/

	c0->Print(print1);
	//c0->Print(print2); //for gif
	c0->Print(print3); //for eps

	legend0->Delete();

}

void PlotDataFit(TH1D hist_sum0, TH1D hist_bkg0, TH1D hist_sum_SM0, TH1D datahist0, TGraphAsymmErrors* hist_unc, TH1D hist_ratio0, std::string input_channel, std::string fOutputFolder)
{
	SetAtlasStyle();

	TH1D * hist_sum = &hist_sum0;
	TH1D * hist_bkg = &hist_bkg0;
	TH1D * hist_sum_SM = &hist_sum_SM0;
	TH1D * datahist = &datahist0;
	TH1D * hist_ratio = &hist_ratio0;
	//TH1D * hist_ratio_err = &hist_ratio_err0;

	cout << "Plot fit results... " << endl;
	cout << endl;
	
	//TCanvas
	TCanvas * c0 = new TCanvas("c0","c0", 1000, 700);
	c0->SetFillColor(0);
	//c0->SetLogy(1);

	TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
	pad1->SetBottomMargin(0);
	//pad1->SetRightMargin(-0.1);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	//pad1->SetGrid();

	double max = hist_sum->GetBinContent(hist_sum->GetMaximumBin())*1.55;

	hist_sum->GetXaxis()->SetTitle("#scale[1.1]{cos #theta}");
	hist_sum->GetYaxis()->SetTitle("#scale[1.23]{Events}");
	hist_sum->GetYaxis()->SetTitleSize(0.05);
	hist_sum->GetXaxis()->SetTitleOffset(0.9);
	hist_sum->GetYaxis()->SetTitleOffset(1.03);
	hist_sum->GetXaxis()->SetLabelSize(0.042);
	hist_sum->GetYaxis()->SetLabelSize(0.048);
	hist_sum->SetMinimum(1.0);
	hist_sum->SetMaximum(max);

	hist_sum->SetStats(kFALSE);
	hist_sum->SetFillStyle(0);
	//hist_sum->SetFillColor(kOrange);
	hist_sum->SetLineColor(kGreen);
	hist_sum->SetLineWidth(2);
	hist_sum->SetMarkerStyle(1);
	hist_sum->SetMarkerColor(kGreen);
	hist_sum->Draw("hist");	


 	//hist_bkg->SetFillColor(kGreen);
	hist_sum_SM->SetStats(kFALSE);
	hist_sum_SM->SetFillStyle(0);
	hist_sum_SM->SetLineColor(kRed);
	hist_sum_SM->SetLineWidth(2);
	hist_sum_SM->SetLineStyle(2);
	hist_sum_SM->SetMarkerStyle(1);
	hist_sum_SM->SetMarkerColor(kRed);
	hist_sum_SM->Draw("histsame");	


	hist_unc->SetMarkerStyle(1);
	hist_unc->SetMarkerColor(kGreen);
	//hist_unc->SetLineWidth(2);
	hist_unc->SetLineColor(kGreen);
	hist_unc->SetFillColor(kGreen);
	hist_unc->SetFillStyle(3001);
	hist_unc->Draw("pz2same");


	hist_bkg->SetFillStyle(1001);
	hist_bkg->SetFillColor(kAzure+1);
	hist_bkg->SetLineColor(kAzure+1);
	hist_bkg->SetLineWidth(2);
	hist_bkg->SetMarkerStyle(1);
	//hist_bkg->SetMarkerColor(kGreen);
	hist_bkg->Draw("histsame");

	hist_sum->Draw("sameaxis");

	datahist->SetFillStyle(0);
	datahist->SetLineColor(kBlack); 
	datahist->SetLineWidth(2);
	datahist->SetMarkerStyle(20);
	datahist->SetMarkerColor(kBlack);
	datahist->Draw("Esame");
	
	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.06); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.155,0.85,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.27,0.85,"Preliminary");
	
	//TLatex latex0;
	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.06);
	//l.DrawLatex(0.245,0.775, "TMinuit");
	l.DrawLatex(0.155,0.79, "TMinuit");
	l.SetTextSize(0.046);
	if (input_channel == "el")    l.DrawLatex(0.155,0.78, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.155,0.78, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.155,0.78, "Combined channel");

	/*if (input_channel == "el") filename_output = "e+jets channel";
	else if (input_channel == "mu") filename_output = "#mu+jets channel";
	else filename_output = " ";*/
	//const char *filename_output;
	//filename_output=output.c_str();
	
	l.SetTextSize(0.037);
	//l.DrawLatex(0.245, 0.735, "Rebinned Sample");
	//l.DrawLatex(0.155, 0.66, "Profile likelihood fit");
	l.DrawLatex(0.155, 0.66, "Template fit");

	
	TLegend * legend0 = new TLegend(0.65,0.63,0.89,0.89);
	//legend0->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	legend0->SetTextFont(42);
	legend0->AddEntry(hist_sum,"Best Fit","l");
	legend0->AddEntry(hist_unc,"Best Fit Unc.","f");
	legend0->AddEntry(hist_bkg,"Background","f");
	legend0->AddEntry(hist_sum_SM,"SM Expectation","l");
	legend0->AddEntry(datahist,"Data","lep");	
	//legend0->AddEntry(hb_c,"Reconstructed as a b-jet","f");
	legend0->SetFillColor(0);
	legend0->SetBorderSize(0);
	legend0->Draw("same");
	

	c0->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.25);
	pad2->SetTopMargin(0);
	//pad2->SetBottomMargin(0.55);
	//pad2->SetRightMargin(-0.1);
        pad2->SetLeftMargin(0.12);
	pad2->Draw();
	pad2->cd();

	//TGraph *graph_line = new TGraph(2);

	hist_ratio->GetXaxis()->SetTitle("#scale[5.1]{cos #theta*}");
	hist_ratio->GetYaxis()->SetTitle("#scale[3.9]{Ratio Data/Fit}");
	hist_ratio ->GetXaxis()->SetTitleSize(0.047);
	hist_ratio ->GetYaxis()->SetTitleSize(0.05);
	hist_ratio->GetXaxis()->SetTitleOffset(4.0);
	hist_ratio->GetYaxis()->SetTitleOffset(1.0);
	hist_ratio->GetXaxis()->SetLabelSize(0.18);
	hist_ratio->GetYaxis()->SetLabelSize(0.14);
	hist_ratio->SetMinimum(0.825);
	hist_ratio->SetMaximum(1.175);
	hist_ratio->SetMarkerStyle(20);
	hist_ratio->SetMarkerColor(kBlack);
	hist_ratio->SetStats(kFALSE);
	hist_ratio->SetLineColor(kWhite);
	hist_ratio->SetLineWidth(0);
	
	hist_ratio->SetFillStyle(0);
	
	hist_ratio->DrawCopy("hist");
	hist_ratio->SetFillColor(kBlack);
	hist_ratio->SetFillStyle(3002);
	//hist_ratio->SetFillStyle(3018);
	//hist_ratio->Draw("hist");
	hist_ratio->Draw("E2same");
	
	TF1 *norm1 = new TF1("fa1","1", -100.0, 500.0);
        norm1 -> SetLineColor(kRed);
        norm1 -> SetLineStyle(1);
        norm1 -> SetLineWidth(2);
        norm1 -> Draw("Same");

	std::string printtitle1 = fOutputFolder+"/Datafit_" + input_channel + ".png";
	std::string printtitle2 = fOutputFolder+"/Datafit_" + input_channel + ".gif";
	std::string printtitle3 = fOutputFolder+"/Datafit_" + input_channel + ".eps";

	//char filename_output[sizeof output + 100];
	//sprintf(filename_output, output);
	
	const char *print1, *print2, *print3;
	print1=printtitle1.c_str();
	print2=printtitle2.c_str();
	print3=printtitle3.c_str();

	/*char printtitle[100] = "plots/first_test.png";// use %f for double, %d for int
	char print[sizeof printtitle + 100];
	sprintf(print, printtitle);*/

	c0->Print(print1);
	//c0->Print(print2); //for gif
	c0->Print(print3); //for eps

	legend0->Delete();
	pad1->Delete();
	pad2->Delete();
}


void PlotDistribution (TH1D* hist, double F_true, std::string xlabel, std::string xtitle, std::string input_channel, std::string sFi_number, TF1* fit, int mode, int mode_gauss, std::string OutputFolder)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 700, 700);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	//Define all necessary labels/titles
	const char *label;
	std::string slabel = xlabel;
	label = slabel.c_str();
	const char *title1, *title2;
	std::string stitle1 = OutputFolder+"/" + xtitle + ".eps";
	title1 = stitle1.c_str();
	std::string stitle2 = OutputFolder+"/" + xtitle + ".gif";
	title2 = stitle2.c_str();
	//	char filename_title1[sizeof title1 + 100];
	//	sprintf(filename_title1, title1, F_true);
	//	char filename_title2[sizeof title2 + 100];
	//	sprintf(filename_title2, title2, F_true);

	//	hist->SetRightMargin(0.10);
	hist->GetXaxis()->SetTitle(label);
	hist->GetYaxis()->SetTitle("Events");
	hist->GetXaxis()->SetTitleOffset(0.3);
	hist->GetYaxis()->SetTitleOffset(0.3);
	hist->GetXaxis()->SetTitleSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->SetMinimum(0.001);
       
	//	hist->SetMaximum(max);
	hist->SetStats(kFALSE);
	//hist->GetXaxis()->SetRange(low,up);
	//cout << low << endl;
	//	if (low < 0) hist->SetAxisRange(4*low,4*up);
	//	else hist->SetAxisRange(0.25*low,2.1*up);

	hist->SetFillStyle(1001);
	hist->SetFillColor(kGray+1);
	hist->SetLineColor(kCyan+2);
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(1);	
	hist->SetMarkerColor(kCyan+2);
	hist->Draw("hist");

	fit->SetLineColor(kViolet+2);
	fit->SetLineWidth(2);
	if (mode_gauss == 10) fit->Draw("same");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.04); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.145,0.87,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.28,0.87,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.035);
	//	l.DrawLatex(0.125,0.79, "TMinuit");

	if (input_channel == "el")    l.DrawLatex(0.145,0.75, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.145,0.75, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.145,0.75, "Combined channel");

// 	cout << "F_true old " << F_true << endl;
// 	if (abs(F_true) < 0.001) {
// 	
// 		 F_true = 0.0;
// 		cout << "F_true " << F_true << endl;
// 	}

	if (mode_gauss == 0) {
	double number = atof(sFi_number.c_str());
	if (abs(number) < 0.001 )  sFi_number = "0";
	const char *Fi_number;
	Fi_number = sFi_number.c_str();
	char file_Fi_number[sizeof Fi_number + 100];
	sprintf(file_Fi_number, Fi_number, F_true);
	//	l.DrawLatex(0.28,0.79, file_Fi_number);
	}

	std::string smean;
	if (mode == 0) smean = "Mean:  %.3f #pm %.3f";// use %f for double, %d for int
	else if (mode == 7) smean = "Mean:  %.4f #pm %.4f";
	else smean = "Mean:  %.0f #pm %.0f";
	const char *mean;
	mean = smean.c_str();
	char latex_mean[sizeof mean + 100];
	sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	std::string ssigma;
	if (mode == 0) ssigma = "Sigma: %.3f #pm %.3f";// use %f for double, %d for int
	else if (mode == 7) ssigma = "Sigma: %.4f #pm %.4f";
	else ssigma = "Sigma: %.0f #pm %.0f";
	const char *sigma;
	sigma = ssigma.c_str();
	char latex_sigma[sizeof sigma + 100];
	sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );

	l.SetTextSize(0.034);
	l.DrawLatex(0.58, 0.85, latex_mean);
	l.DrawLatex(0.58, 0.80, latex_sigma);

	std::string sgauss_mean;
	if (mode == 0) sgauss_mean = "Gaussian Mean:   %.3f #pm %.3f";// use %f for double, %d for int
	else sgauss_mean = "Gaussian Mean:   %.0f #pm %.0f";
	const char *gauss_mean;
	gauss_mean = sgauss_mean.c_str();
	char latex_gauss_mean[sizeof gauss_mean + 100];
	sprintf(latex_gauss_mean, gauss_mean, fit->GetParameter(1), fit->GetParError(1) );

	std::string sgauss_sigma;
	if (mode == 0) sgauss_sigma = "Gaussian Sigma:   %.3f #pm %.3f";// use %f for double, %d for int
	else sgauss_sigma = "Gaussian Sigma:   %.0f #pm %.0f";
	const char *gauss_sigma;
	gauss_sigma = sgauss_sigma.c_str();
	char latex_gauss_sigma[sizeof gauss_sigma + 100];
	sprintf(latex_gauss_sigma, gauss_sigma, fit->GetParameter(2), fit->GetParError(2) );

	l.SetTextSize(0.03);
	//if (mode_gauss == 10) l.DrawLatex(0.55, 0.55, latex_gauss_mean);
	//if (mode_gauss == 10) l.DrawLatex(0.55, 0.50, latex_gauss_sigma);
	
	//legend
	TLegend *legend;
	if (mode_gauss == 10) legend = new TLegend(0.66,0.75,0.89,0.89);
 	else legend = new TLegend(0.70,0.845,0.89,0.89);
	legend->SetTextFont(42);
	legend->AddEntry(hist,"Histogram","f");
	if (mode_gauss == 10) legend->AddEntry(fit,"Gaussian Exp.","l");
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	//legend->SetHeader("#font[132]{legend}");
	legend->Draw("same");

	c_F->Print(stitle1.c_str());
	//c_F->Print(filename_title2);

	//c_F->Delete();
}


void PlotDistribution_N (TH1D* hist, double F_true, double xsec, double lumi, std::string xlabel, std::string xtitle, std::string input_channel, std::string sFi_number, int mode, std::string OutputFolder, std::string mode_mean)
{

    gErrorIgnoreLevel = kError;
 
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);
	
	SetAtlasStyle();

	gStyle -> SetPadRightMargin(0.08);
	gStyle -> SetPadLeftMargin(0.10);
	gStyle -> SetPadBottomMargin(0.10);

	//Define all necessary labels/titles
	const char *label;
	std::string slabel = xlabel;
	label = slabel.c_str();
	const char *title1, *title2;
	std::string stitle1 = OutputFolder+"/" + xtitle + ".eps";
	title1 = stitle1.c_str();
	std::string stitle2 = OutputFolder+"/" + xtitle + ".png";
	title2 = stitle2.c_str();
	char filename_title1[sizeof title1 + 100];
	sprintf(filename_title1, title1, F_true);
	char filename_title2[sizeof title2 + 100];
	sprintf(filename_title2, title2, F_true);

	hist->GetXaxis()->SetTitle(label);
	//	hist->GetXaxis()->SetTitleSize(0.042);
	hist->GetYaxis()->SetTitle("Events");
	hist->GetXaxis()->SetTitleOffset(1.0);
	hist->GetYaxis()->SetTitleOffset(1.15);
	hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetXaxis()->SetNdivisions(505);
	hist->SetMinimum(0.001);

	double max = hist->GetBinContent(hist->GetMaximumBin());

	hist->SetMaximum(max*1.4);
	hist->SetStats(kFALSE);

	hist->SetFillStyle(1001);
	hist->SetFillColor(kGray+1);
	hist->SetLineColor(kCyan+2);
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(1);	
	hist->SetMarkerColor(kCyan+2);
	hist->Draw("hist");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.032); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.135,0.88,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.24,0.88,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.032);
	//	l.DrawLatex(0.125,0.79, "TMinuit    Normal Constraint");

	if (input_channel == "el") l.DrawLatex(0.135,0.84, "e+jets channel");
	if (input_channel == "mu") l.DrawLatex(0.135,0.84, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.135,0.84, "Combined channel");

	const char *Fi_number;
	Fi_number = sFi_number.c_str();
	char file_Fi_number[sizeof Fi_number + 100];
	if (mode == 1) sprintf(file_Fi_number, Fi_number, xsec*lumi);
	else sprintf(file_Fi_number, Fi_number, F_true*xsec*lumi);
	l.DrawLatex(0.135,0.80, file_Fi_number);
	
	l.SetTextSize(0.030);
	
	if(mode_mean == "Pull"){

	  char mean[100] = "Mean:  %.2f #pm %.2f";// use %f for double, %d for int
	  char latex_mean[sizeof mean + 100];
	  sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	  char sigma[100] = "Sigma: %.2f #pm %.2f";// use %f for double, %d for int
	  char latex_sigma[sizeof sigma + 100];
	  sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );

	  l.SetTextSize(0.030);

	  l.DrawLatex(0.63, 0.88, latex_mean);
	  l.DrawLatex(0.63, 0.83, latex_sigma);

	}
	else{

	  char mean[100] = "Mean:  %.1f #pm %.1f";// use %f for double, %d for int
	  char latex_mean[sizeof mean + 100];
	  sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	  char sigma[100] = "Sigma: %.1f #pm %.1f";// use %f for double, %d for int
	  char latex_sigma[sizeof sigma + 100];
	  sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );

	  l.SetTextSize(0.030);

	  l.DrawLatex(0.58, 0.88, latex_mean);
          l.DrawLatex(0.58, 0.83, latex_sigma);

	}

	c_F->Print(filename_title1);
	c_F->Print(filename_title2);

	//c_F->Delete();
}


void PlotErrorDistribution (TH1D* hist, double F_true, std::string xlabel, std::string xtitle, double max, std::string input_channel, std::string sFi_number, int mode)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 800, 600);
	c_F->SetFillColor(0);

	//Define all necessary labels/titles
	const char *label;
	std::string slabel = "#scale[1.1]{" + xlabel + " [-]}";
	label = slabel.c_str();
	const char *title1, *title2;
	std::string stitle1 = "plots/" + xtitle + ".png";
	title1 = stitle1.c_str();
	std::string stitle2 = "plots/" + xtitle + ".gif";
	title2 = stitle2.c_str();
	char filename_title1[sizeof title1 + 100];
	sprintf(filename_title1, title1, F_true);
	char filename_title2[sizeof title2 + 100];
	sprintf(filename_title2, title2, F_true);

	hist->GetXaxis()->SetTitle(label);
	hist->GetYaxis()->SetTitle("#scale[1.1]{Events [-]}");
	hist->GetXaxis()->SetTitleOffset(1.25);
	hist->GetYaxis()->SetTitleOffset(1.4);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->SetMinimum(0.001);
	hist->SetMaximum(max);
	hist->SetStats(kFALSE);

	hist->SetFillStyle(1001);
	hist->SetFillColor(kGray+1);
	hist->SetLineColor(kCyan+2);
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(1);	
	hist->SetMarkerColor(kCyan+2);
	hist->Draw("hist");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.04); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.125,0.84,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.24,0.84,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.035);
	l.DrawLatex(0.125,0.79, "TMinuit");

	if (input_channel == "el") l.DrawLatex(0.125,0.75, "e channel");
	if (input_channel == "mu") l.DrawLatex(0.125,0.75, "#mu channel");
	
	const char *Fi_number;
	Fi_number = sFi_number.c_str();
	char file_Fi_number[sizeof Fi_number + 100];
	sprintf(file_Fi_number, Fi_number, F_true);
	l.DrawLatex(0.25,0.75, file_Fi_number);
	
	char mean[100] = "Mean:  %.4f #pm %.4f";// use %f for double, %d for int
	char latex_mean[sizeof mean + 100];
	sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	char sigma[100] = "Sigma: %.4f #pm %.4f";// use %f for double, %d for int
	char latex_sigma[sizeof sigma + 100];
	sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );

	l.SetTextSize(0.03);
	if (mode == 1) l.DrawLatex(0.65, 0.85, latex_mean);
	if (mode == 1) l.DrawLatex(0.65, 0.80, latex_sigma);

	c_F->Print(filename_title1);
	//c_F->Print(filename_title2);

	//c_F->Delete();
}



void PlotErrorDataComp (TH1D* hist, double F_true, double F_data, std::string xlabel, std::string input_channel, std::string sFi_number, std::string filename, int mode)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	c_F->SetBottomMargin(0.12);
	//pad1->SetRightMargin(-0.1);
	c_F->SetLeftMargin(0.13);

	double max = hist -> GetBinContent(hist->GetMaximumBin())*1.65;

	hist->GetXaxis()->SetTitle(xlabel.c_str());
	hist->GetYaxis()->SetTitle("Events");
	hist->GetXaxis()->SetTitleOffset(1.2);
	hist->GetYaxis()->SetTitleOffset(1.45);
	hist->GetXaxis()->SetTitleSize(0.045);
        hist->GetYaxis()->SetTitleSize(0.047);
	hist->GetXaxis()->SetLabelSize(0.045);
	hist->GetYaxis()->SetLabelSize(0.047);
	hist->SetMinimum(0.001);
	hist->SetMaximum(max);
	hist->SetStats(kFALSE);
	hist->GetXaxis()->SetNdivisions(504);
	hist->GetXaxis()->SetNoExponent();

	hist->SetFillStyle(1001);
	hist->SetFillColor(kGray+1);
	hist->SetLineColor(kCyan+2);
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(1);	
	hist->SetMarkerColor(kCyan+2);
	hist->Draw("hist");

	//	fit->SetLineColor(kViolet+2);
	//	fit->SetLineWidth(2);
	//fit->Draw("same");

	TGraph *graph_line = new TGraph(2);
	graph_line->SetPoint(1, F_data, 0);
	graph_line->SetPoint(2, F_data, max*0.7);


	std::stringstream oss;
	oss << std::setprecision(3) << F_data;

	graph_line->RemovePoint(0);
	graph_line->SetLineStyle(2);
	graph_line->SetLineColor(kRed);
	graph_line->SetLineWidth(3);
	graph_line->SetMarkerStyle(20);
	graph_line->SetMarkerColor(kRed);
	graph_line->Draw("Lsame");
	//graph_line->Draw("same");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.038); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.155,0.88,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.285,0.88,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	//	l.DrawLatex(0.15,0.79, "TMinuit");

	if (input_channel == "el")    l.DrawLatex(0.155, 0.84, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.155, 0.84, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.155, 0.84, "Combined channel");

// 	const char *Fi_number;
// 	Fi_number = sFi_number.c_str();
// 	char file_Fi_number[sizeof Fi_number + 100];
// 	sprintf(file_Fi_number, Fi_number, F_true);
// 	l.DrawLatex(0.25,0.75, file_Fi_number);
	
	char mean[100] = "Mean:  %.5f #pm %.5f";// use %f for double, %d for int
	char latex_mean[sizeof mean + 100];
	sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	char sigma[100] = "Sigma: %.5f #pm %.5f";// use %f for double, %d for int
	char latex_sigma[sizeof sigma + 100];
	sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );

	l.SetTextSize(0.028);
	if (mode == 1) l.DrawLatex(0.60, 0.74, latex_mean);
	if (mode == 1) l.DrawLatex(0.60, 0.70, latex_sigma);

	TLegend * legend0 = new TLegend(0.58,0.78,0.92,0.92);
	//legend0->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	legend0->SetTextFont(42);
	legend0->AddEntry(hist,"Pseudo-experiments","f");
	//legend0->AddEntry(fit,"Gaussian Fit","l");
	legend0->AddEntry(graph_line, ("Observed: "+oss.str()).c_str(),"l");
	legend0->SetFillColor(0);
	legend0->SetBorderSize(0);
	legend0->Draw("same");

	c_F->Print(filename.c_str());
	//c_F->Print(filename_title2);

	//c_F->Delete();

	graph_line->Delete();
	legend0->Delete();
}




void Plot2Dcorrelation_F (TH2D* hist2, double Fa_true, double Fb_true, std::string xlabel_x, std::string xlabel_y, std::string filename, std::string input_channel)
{

	TCanvas *c_F2 = new TCanvas("c_F2","c_F2", 800, 800);
	c_F2->SetFillColor(0);
	
	const char *label_x;
	std::string slabel_x = "#scale[1.05]{" + xlabel_x + " [-]}";
	label_x = slabel_x.c_str();
	const char *label_y;
	std::string slabel_y = "#scale[1.05]{" + xlabel_y + " [-]}";
	label_y = slabel_y.c_str();
	const char *title1, *title2;
	/*std::string stitle1 = "plots/" + xtitle + ".png";
	title1 = stitle1.c_str();
	char filename_title1[sizeof title1 + 100];
	sprintf(filename_title1, title1, Fa_true, Fb_true);
	std::string stitle2 = "plots/" + xtitle + ".gif";
	title2 = stitle2.c_str();
	char filename_title2[sizeof title2 + 100];
	sprintf(filename_title2, title2, Fa_true, Fb_true);*/

	c_F2->SetGrid();

	hist2->GetXaxis()->SetTitle(label_x);
	hist2->GetYaxis()->SetTitle(label_y);
	hist2->GetXaxis()->SetTitleOffset(1.25);
	hist2->GetYaxis()->SetTitleOffset(1.4);
	hist2->GetXaxis()->SetLabelSize(0.032);
	hist2->GetXaxis()->SetLabelOffset(0.0);
	//hist2->GetXaxis()->SetNoExponent();
	//hist2->GetXaxis()->SetNdivisions(6);
	hist2->GetYaxis()->SetLabelSize(0.032);
	hist2->SetStats(kFALSE);

	hist2->SetFillStyle(1001);
	hist2->SetFillColor(kGray+3);
	hist2->SetLineColor(kCyan+2);
	hist2->SetLineWidth(2);
	hist2->SetMarkerStyle(1);
	hist2->SetMarkerColor(kCyan+2);
	hist2->Draw("COLZ");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.04); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.125,0.84,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.26,0.84,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.035);
	l.DrawLatex(0.125,0.79, "Correlation");
	
		if (input_channel == "el") l.DrawLatex(0.125,0.75, "e channel");
		if (input_channel == "mu") l.DrawLatex(0.125,0.75, "#mu channel");
	
	char corr[100] = "Corr:  %.4f";// use %f for double, %d for int
	char latex_corr[sizeof corr + 100];
	sprintf(latex_corr, corr, hist2->GetCorrelationFactor());

	l.SetTextSize(0.035);
	l.DrawLatex(0.67, 0.17, latex_corr);

	//c_F2->Print(filename_title1);
	c_F2->Print(filename.c_str());
	//c_F2->Print(filename_title2);
}



void Plot2Dcorrelation_N (TH2D* hist2, double F0_true, std::string xlabel_x, std::string xlabel_y, std::string filename, std::string input_channel)
{

	TCanvas *c_F2 = new TCanvas("c_F2","c_F2", 600, 600);
	c_F2->SetFillColor(0);

	TGaxis::SetMaxDigits(3);
	
        TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1,1);
	pad1->SetBottomMargin(0.12);
        pad1->SetRightMargin(0.10);
        pad1->SetLeftMargin(0.14);
        pad1->Draw();
        pad1->cd();

	const char *label_x;
	std::string slabel_x = "#scale[1.05]{" + xlabel_x + "}";
	label_x = slabel_x.c_str();
	const char *label_y;
	std::string slabel_y = "#scale[1.05]{" + xlabel_y + "}";
	label_y = slabel_y.c_str();
	const char *title1, *title2;
	/*std::string stitle1 = "plots/" + xtitle + "_F0=%.1f.png";
	title1 = stitle1.c_str();
	char filename_title1[sizeof title1 + 100];
	sprintf(filename_title1, title1, F0_true);
	std::string stitle2 = "plots/" + xtitle + "_F0=%.1f.gif";
	title2 = stitle2.c_str();
	char filename_title2[sizeof title2 + 100];
	sprintf(filename_title2, title2, F0_true);*/

	c_F2->SetGrid();

	hist2->GetXaxis()->SetTitle(label_x);
	hist2->GetYaxis()->SetTitle(label_y);
	hist2->GetXaxis()->SetTitleOffset(1.05);
	hist2->GetYaxis()->SetTitleOffset(1.375);
	hist2->GetXaxis()->SetLabelSize(0.045);
	hist2->GetZaxis()->SetLabelSize(0.038);
	hist2->GetXaxis()->SetLabelOffset(0.0);
	hist2->GetXaxis()->SetNoExponent();
	hist2->GetXaxis()->SetNdivisions(505);
	hist2->GetYaxis()->SetLabelSize(0.045);
	hist2->SetStats(kFALSE);

	hist2->SetFillStyle(1001);
	hist2->SetFillColor(kGray+3);
	hist2->SetLineColor(kCyan+2);
	hist2->SetLineWidth(2);
	hist2->SetMarkerStyle(1);
	hist2->SetMarkerColor(kCyan+2);
	hist2->Draw("COLZ");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.036); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.175,0.88,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.301,0.88,"Preliminary");

	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.035);
	//	l.DrawLatex(0.175,0.84, "Correlation");
	
		if(input_channel == "el")    l.DrawLatex(0.175,0.84, "e+jets channel");
		if(input_channel == "mu")    l.DrawLatex(0.175,0.84, "#mu+jets channel");
		if(input_channel == "el_mu") l.DrawLatex(0.175,0.84, "Combined channel");
		
	char corr[100] = "Corr:  %.4f";// use %f for double, %d for int
	char latex_corr[sizeof corr + 100];
	sprintf(latex_corr, corr, hist2->GetCorrelationFactor());
	
	l.SetTextSize(0.035);
	l.DrawLatex(0.67, 0.19, latex_corr);
	
	//c_F2->Print(filename_title1);
	c_F2->Print(filename.c_str());
	//c_F2->Print(filename_title2);
}




void PlotPullDistribution (TH1D* hist, double F_true, std::string xlabel, std::string xtitle, double max, std::string input_channel, std::string sFi_number, TF1* fit)
{

	TCanvas *c_F = new TCanvas("c_F","c_F", 900, 600);
	c_F->SetFillColor(0);

	const char *label;
	std::string slabel = "#scale[1.1]{" + xlabel + " [-]}";
	label = slabel.c_str();
	const char *title1, *title2;
	std::string stitle1 = "plots/" + xtitle + ".png";
	title1 = stitle1.c_str();
	char filename_title1[sizeof title1 + 100];
	sprintf(filename_title1, title1, F_true);
	std::string stitle2 = "plots/" + xtitle + ".gif";
	title2 = stitle2.c_str();
	char filename_title2[sizeof title2 + 100];
	sprintf(filename_title2, title2, F_true);

	hist->GetXaxis()->SetTitle(label);
	hist->GetYaxis()->SetTitle("#scale[1.1]{Events [-]}");
	hist->GetXaxis()->SetTitleOffset(1.25);
	hist->GetYaxis()->SetTitleOffset(1.4);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->SetMinimum(0.001);
	hist->SetMaximum(max);
	hist->SetStats(kFALSE);	

	hist->SetFillStyle(1001);
	hist->SetFillColor(kGray+1);
	hist->SetLineColor(kCyan+2);
	hist->SetLineWidth(2);
	hist->SetMarkerStyle(1);
	hist->SetMarkerColor(kCyan+2);
	hist->Draw("hist");

	fit->SetLineColor(kViolet+2);
	fit->SetLineWidth(2);
	fit->Draw("same");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.04); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.125,0.84,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.24,0.84,"Preliminary");
	
	l.SetNDC(1);
	l.SetTextFont(42);
	l.SetTextSize(0.035);
	l.DrawLatex(0.125,0.79, "TMinuit   Pull distribution");
	
		if (input_channel == "el") l.DrawLatex(0.125,0.75, "e channel");
		if (input_channel == "mu") l.DrawLatex(0.125,0.75, "#mu channel");

	const char *Fi_number;
	Fi_number = sFi_number.c_str();
	char file_Fi_number[sizeof Fi_number + 100];
	sprintf(file_Fi_number, Fi_number, F_true);
	l.DrawLatex(0.25,0.75, file_Fi_number);

	char mean[100] = "Mean:  %.4f #pm %.4f";// use %f for double, %d for int
	char latex_mean[sizeof mean + 100];
	sprintf(latex_mean, mean, hist->GetMean(), hist->GetMeanError() );
	char sigma[100] = "Sigma: %.4f #pm %.4f";// use %f for double, %d for int
	char latex_sigma[sizeof sigma + 100];
	sprintf(latex_sigma, sigma, hist->GetRMS(), hist->GetRMSError() );
	
	l.SetTextSize(0.03);
	l.DrawLatex(0.65, 0.65, latex_mean);
	l.DrawLatex(0.65, 0.60, latex_sigma);
	
	std::string sgauss_mean = "Gaussian mean:   %.4f #pm %.4f";// use %f for double, %d for int
	const char *gauss_mean;
	gauss_mean = sgauss_mean.c_str();
	char latex_gauss_mean[sizeof gauss_mean + 100];
	sprintf(latex_gauss_mean, gauss_mean, fit->GetParameter(1), fit->GetParError(1) );

	std::string sgauss_sigma = "Gaussian Sigma:   %.4f #pm %.4f";// use %f for double, %d for int
	const char *gauss_sigma;
	gauss_sigma = sgauss_sigma.c_str();
	char latex_gauss_sigma[sizeof gauss_sigma + 100];
	sprintf(latex_gauss_sigma, gauss_sigma, fit->GetParameter(2), fit->GetParError(2) );

	l.SetTextSize(0.03);
	l.DrawLatex(0.55, 0.55, latex_gauss_mean);
	l.DrawLatex(0.55, 0.50, latex_gauss_sigma);
	
	//legend
	TLegend *legend = new TLegend(0.68,0.75,0.89,0.89);
	legend->SetTextFont(42);
	legend->AddEntry(hist,"Histogram","f");
	legend->AddEntry(fit,"Gaussian Fit","l");
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	//legend->SetHeader("#font[132]{legend}");
	legend->Draw("same");

	c_F->Print(filename_title1);
	//c_F->Print(filename_title2);
		
	//c_F->Delete();
}


void PlotCalDistribution (TGraphErrors* graph, TGraphErrors* diff, std::string xlabel, std::string ylabel, std::string xtitle, double min_x, double max_x, double min, double max, std::string input_channel, int mode, std::string output_folder)
{

    gErrorIgnoreLevel = kError;
 
    SetAtlasStyle();

	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 800, 700);
	c_F->SetFillColor(0);

	//Define all necessary labels/titles
	const char *label_x;
	const char *label_y;
	std::string slabel_x;
	if (mode == 5) slabel_x = "{Counter}";
	else slabel_x = "{Input " + xlabel + "}";
	
	std::string slabel_y = "#{Output " + ylabel + "}";
	label_x = slabel_x.c_str();
	label_y = slabel_y.c_str();

	const char *title1, *title2;
	//std::string stitle1 = "Output_mu_Cal_Plots/" + xtitle + ".png";
	std::string stitle1 = output_folder + "/" + xtitle + ".png";
	title1 = stitle1.c_str();
	//std::string stitle2 = "Output_mu_Cal_Plots/" + xtitle + ".eps";
	std::string stitle2 = output_folder + "/" + xtitle + ".eps";
	title2 = stitle2.c_str();
	//char filename_title[sizeof title + 100];
	//sprintf(filename_title, title, F0_true);
	
	const char *label_diff;
	std::string slabel_diff;
	if (mode==0) slabel_diff = "#{#Delta" + xlabel + "}";
	else if (mode==2) slabel_diff = "#{#Delta " + ylabel + " [-]}";
	else slabel_diff = "#{#Delta N/N}";
	//std::string slabel_diff = "#scale[5.3]{Ratio " + xlabel + " [-]}";
	//std::string slabel_diff = "#Delta/Nom";
	label_diff = slabel_diff.c_str();

	TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
	pad1->SetBottomMargin(0);
	pad1->SetRightMargin(-0.10);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	pad1->SetGrid();
	graph->RemovePoint(0);
	graph->SetTitle("");
	//graph->GetXaxis()->SetTitle(label_x);
	graph->GetYaxis()->SetTitle(label_y);
	graph->GetXaxis()->SetTitleOffset(1.2);
	graph->GetYaxis()->SetTitleOffset(1.7);
	graph->GetXaxis()->SetLabelSize(0.04);
	graph->GetYaxis()->SetLabelSize(0.04);
	//graph->GetYaxis()->SetTitleSize(0.04);
	if (mode == 1) graph->GetXaxis()->SetNoExponent();
	graph->SetMinimum(min);
	graph->SetMaximum(max);
	if (mode == 5) graph->GetXaxis()->SetLimits(0.5,8.5);
	else graph->GetXaxis()->SetLimits(min_x,max_x);

	//graph->SetStats(kFALSE);
	graph->SetFillStyle(0);
	graph->SetLineColor(kCyan+2);
	graph->SetLineWidth(2);
	graph->SetMarkerStyle(20);
	graph->SetMarkerColor(kCyan+2);
	graph->Draw("APZ");

	TF1 *fit = new TF1("fit", "[0]*x + [1]",min_x, max_x);

	fit->SetLineColor(kCyan+2);
	fit->SetLineWidth(2);
	graph->Fit("fit", "Q");
	fit->Draw("same");

	//cout << "Fit result " << setprecision(4) << fit->GetParameter(0) << " +/- " << setprecision(4) << fit->GetParError(0) << endl;

	TLegend *legend = new TLegend(0.69,0.83,0.89,0.89);
	//legend->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	legend->SetTextFont(42);
	if (mode == 0) legend->AddEntry(graph,"#scale[1.5]{Helicity fraction}","lep");
	else legend->AddEntry(graph,"#scale[1.5]{Mean values}","lep");
	legend->AddEntry(fit,"#scale[1.5]{Linear fit}","l");
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->Draw("same");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.042); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.14,0.84,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.25,0.84,"Preliminary");
	l.SetTextSize(0.04);
	l.DrawLatex(0.14,0.79,"TMinuit");
	l.SetTextSize(0.04); 
	//l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	l.DrawLatex(0.25,0.79,"Calibration curves");

	if (input_channel == "el") l.DrawLatex(0.14,0.74,    "e+jets channel");
	if (input_channel == "mu") l.DrawLatex(0.14,0.74,    "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.14,0.74, "Combined channel");

	std::string sslope;
	if (mode == 1) sslope = "Slope:      %.3f #pm %.3f";// use %f for double, %d for int
	else sslope = "Slope:      %.4f #pm %.4f";// use %f for double, %d for int
	const char *slope;
	slope = sslope.c_str();
	char latex_slope[sizeof slope + 100];
	sprintf(latex_slope, slope, fit->GetParameter(0), fit->GetParError(0) );
	
	std::string sintercept;
	if (mode == 1) sintercept = "Intercept: %.0f #pm %.0f";// use %f for double, %d for int
	else sintercept = "Intercept: %.4f #pm %.4f";// use %f for double, %d for int
	const char *intercept;
	intercept = sintercept.c_str();
	char latex_intercept[sizeof intercept + 100];
	sprintf(latex_intercept, intercept,  fit->GetParameter(1), fit->GetParError(1));

	l.SetTextSize(0.042);
	l.SetTextFont(42);
	l.DrawLatex(0.545, 0.095, latex_slope);
	l.DrawLatex(0.545, 0.045, latex_intercept);

	c_F->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0.06,1,0.2);
	pad2->SetTopMargin(0);
	//pad2->SetBottomMargin(0.55);
	//pad2->SetRightMargin(-0.1);
        pad2->SetLeftMargin(0.12);
	pad2->Draw();
	pad2->cd();

	TGraph *graph_line = new TGraph(2);

	diff->RemovePoint(0);
	diff->SetTitle("");
	diff->GetXaxis()->SetTitle(label_x);
	diff->GetYaxis()->SetTitle(label_diff);
	diff->GetXaxis()->SetTitleOffset(5.6);
	diff->GetYaxis()->SetTitleOffset(1.7);
	diff->GetXaxis()->SetLabelSize(0.225);
	diff->GetYaxis()->SetLabelSize(0.225);
	//diff->GetXaxis()->SetTitleSize(0.05);
        //diff->GetYaxis()->SetTitleSize(0.04);
	diff->GetYaxis()->SetNdivisions(6);
	if (mode == 0){
	diff->SetMinimum(-0.029);
	diff->SetMaximum(0.029);
	//diff->SetMinimum(0.901);
	//diff->SetMaximum(1.099);
	diff->GetXaxis()->SetLimits(min_x, max_x);
	}
	else {
	diff->SetMinimum(-0.059);
	diff->SetMaximum(0.059);
	//diff->SetMinimum(0.901);
	//diff->SetMaximum(1.099);
	diff->GetXaxis()->SetLimits(min_x,max_x);
	}

	//diff->SetStats(kFALSE);
	diff->SetFillStyle(0);
	diff->SetLineColor(kBlue-7);
	diff->SetLineWidth(2);
	diff->SetMarkerStyle(20);
	diff->SetMarkerColor(kBlue-7);
	diff->Draw("APZ");
		
	graph_line->SetPoint(1,-10000000,0);
	graph_line->SetPoint(2,1000000,0);

	graph_line->RemovePoint(0);
	graph_line->SetFillStyle(0);
	graph_line->SetLineColor(kGray);
	graph_line->SetLineWidth(2);
	graph_line->SetMarkerStyle(20);
	graph_line->SetMarkerColor(kGray);
	graph_line->Draw("Lsame");
	
	c_F->Print(title1);
	c_F->Print(title2);

	fit->Delete();
	graph_line->Delete();
	pad1->Delete();
	pad2->Delete();
}



void PlotCalDistribution_N (TGraphErrors* graph, TGraphErrors* diff, std::string xlabel, std::string ylabel, std::string xtitle, double min, double max, std::string input_channel, std::string output_folder, std::string ValueLabel)
{
        gErrorIgnoreLevel = kError;

	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	//Define all necessary labels/titles
	const char *label_x;
	const char *label_y;
	//	std::string slabel_x = "#scale[5.7]{Input " + xlabel + " [-]}";
	//	std::string slabel_y = "#scale[1.3]{Output " + xlabel + " [-]}";

	std::string slabel_x = "Input " + xlabel;
	std::string slabel_y = "Output " + ylabel;

	label_x = slabel_x.c_str();
	label_y = slabel_y.c_str();
	
	const char *title1, *title2;
	std::string stitle1 = output_folder+"/" + xtitle + ".png";
	title1 = stitle1.c_str();
	std::string stitle2 = output_folder+"/" + xtitle + ".eps";
	title2 = stitle2.c_str();
	//char filename_title[sizeof title + 100];
	//sprintf(filename_title, title, F0_true);
	
	const char *label_diff;
	std::string slabel_diff;

	double maxy = TMath::MaxElement(graph -> GetN(), (graph->GetY()));
    double miny = TMath::MinElement(graph -> GetN(), (graph->GetY()));

	if(maxy > 2){
	  //slabel_diff = "#Delta N/N";
	  slabel_diff = "#Delta N"; // modified by MJK
	  double scale = maxy/10.0;
	  maxy = maxy+3*scale;
	  miny = miny-3*scale;
	}
	else
	  //slabel_diff = "#Delta F/F";
	slabel_diff = "#Delta F"; // modified by MJK


	//std::string slabel_diff = "#scale[5.3]{Ratio " + xlabel + " [-]}";
	label_diff = slabel_diff.c_str();

	TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
	pad1->SetBottomMargin(0);
	pad1->SetRightMargin(0.08);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	//	pad1->SetGrid();
	//	graph->RemovePoint(0);
	graph->SetTitle("");
	//graph->GetXaxis()->SetTitle(label_x);
	graph->GetYaxis()->SetTitle(label_y);
	graph->GetYaxis()->SetTitleSize(0.062);
	graph->GetXaxis()->SetTitleOffset(1.2);
	graph->GetYaxis()->SetTitleOffset(0.875);
	//	graph->GetXaxis()->SetLabelSize(0.04);
       	graph->GetYaxis()->SetLabelSize(0.051);
	graph->GetXaxis()->SetNdivisions(505);
	graph->GetYaxis()->SetNdivisions(505);

	graph->GetXaxis()->SetNoExponent();

	if(maxy > 2 && max < 2){
	  graph->SetMinimum(miny);
	  graph->SetMaximum(maxy);
	}
	graph->GetXaxis()->SetLimits(min,max);

	//graph->SetStats(kFALSE);
	graph->SetFillStyle(0);
	graph->SetLineColor(kOrange+9);
	graph->SetLineWidth(2);
	graph->SetMarkerStyle(20);
	graph->SetMarkerColor(kOrange+9);
	graph->Draw("APZ");

	TF1 *fit = new TF1("fit", "[0]*x + [1]",min,max);

	fit->SetLineColor(kCyan+2);
	fit->SetLineWidth(2);
	graph->Fit("fit", "Q");
	fit->Draw("same");

	std::cout << "Slope:  " << fit->GetParameter(0) << " $\\pm$ " << fit->GetParError(0) << std::endl;
	std::cout << "Offset: " << fit->GetParameter(1) << " $\\pm$ " << fit->GetParError(1) << std::endl;

	//cout << "Fit result " << setprecision(4) << fit->GetParameter(0) << " +/- " << setprecision(4) << fit->GetParError(0) << endl;

	TLegend *legend = new TLegend(0.58,0.81,0.78,0.89);
	//legend->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	legend->SetTextFont(42);
	legend->AddEntry(graph,"#scale[1.25]{Mean Signal}","lep");
	legend->AddEntry(fit,"#scale[1.25]{Linear fit}","l");
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->Draw("same");

	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.042); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.170,0.87,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.280,0.87,"Preliminary");
	//	l.SetTextSize(0.036);
	//	l.DrawLatex(0.14,0.79,"TMinuit");
	l.SetTextSize(0.041); 
	//l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	//l.DrawLatex(0.25,0.79,"Calibration");

	if (input_channel == "el")    l.DrawLatex(0.170,0.82, "e+jets channel");
	if (input_channel == "mu")    l.DrawLatex(0.170,0.82, "#mu+jets channel");
	if (input_channel == "el_mu") l.DrawLatex(0.170,0.82, "Combined channel");

	std::string sslope;
	sslope = "Slope:   %.3f #pm %.3f";// use %f for double, %d for int
	const char *slope;
	slope = sslope.c_str();
	char latex_slope[sizeof slope + 100];
	sprintf(latex_slope, slope, fit->GetParameter(0), fit->GetParError(0) );
	
	std::string sintercept;
	sintercept = "Offset:   %.3f #pm %.3f";// use %f for double, %d for int
	const char *intercept;
	intercept = sintercept.c_str();
	char latex_intercept[sizeof intercept + 100];
	sprintf(latex_intercept, intercept,  fit->GetParameter(1), fit->GetParError(1));

	std::string schi2;
        schi2 = "#chi^{2}/ndof: %.2f / %d";// use %f for double, %d for int
        const char *chi2;
        chi2 = schi2.c_str();
        char latex_chi2[sizeof chi2 + 100];
        sprintf(latex_chi2, chi2,  fit->GetChisquare(), fit->GetNDF());

	l.SetTextSize(0.041);
	l.DrawLatex(0.52, 0.18, latex_chi2);
	l.DrawLatex(0.52, 0.11, latex_slope);
	l.DrawLatex(0.52, 0.06, latex_intercept);

	const char *Fi_number;
        Fi_number = ValueLabel.c_str();
       	l.DrawLatex(0.170, 0.77, ValueLabel.c_str());
	
	c_F->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0.003,1,0.25); // ========================= The lower part of plots
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.45);
        pad2->SetRightMargin(0.08);
	pad2->SetLeftMargin(0.12);
	pad2->Draw();
	pad2->cd();

	TGraph *graph_line = new TGraph(2);

	//	diff->RemovePoint(0);
	diff->SetTitle("");
	diff->GetXaxis()->SetTitle(label_x);
	diff->GetYaxis()->SetTitle(label_diff);
	diff->GetXaxis()->SetTitleSize(0.21);
    diff->GetYaxis()->SetTitleSize(0.18);
	diff->GetXaxis()->SetTitleOffset(0.90);
	diff->GetYaxis()->SetTitleOffset(0.34);
	diff->GetXaxis()->SetLabelSize(0.15);
	diff->GetYaxis()->SetLabelSize(0.14);
	diff->GetXaxis()->SetNdivisions(505);
	diff->GetYaxis()->SetNdivisions(505);
	diff->SetMinimum(-0.010); //FIXME
	diff->SetMaximum(0.010);  //FIXME
	//diff->SetMinimum(0.901);
	//diff->SetMaximum(1.099);
	diff->GetXaxis()->SetLimits(min,max);
	

	//diff->SetStats(kFALSE);
	diff->SetFillStyle(0);
	diff->SetLineColor(kOrange+9);
	diff->SetLineWidth(2);
	diff->SetMarkerStyle(20);
	diff->SetMarkerSize(0.7);
	diff->SetMarkerColor(kOrange+9);
	diff->GetYaxis()->CenterTitle();
	diff->Draw("APZ");
		
	graph_line->SetPoint(1,-200000,0);
	graph_line->SetPoint(2,1000000,0);

	graph_line->RemovePoint(0);
	graph_line->SetFillStyle(0);
	graph_line->SetLineColor(kGray);
	graph_line->SetLineWidth(2);
	graph_line->SetMarkerStyle(20);
	graph_line->SetMarkerColor(kGray);
	graph_line->Draw("Lsame");
      
	//	c_F->Print(title1);
	c_F->Print(title2);

	fit->Delete();
	graph_line->Delete();
	pad1->Delete();
	pad2->Delete();
}


void PlotCalPullDistribution (TGraphErrors* graph, std::string xlabel, std::string ylabel, std::string xtitle, double xmin, double xmax, double ymin, double ymax, std::string input_channel, std::string smode, std::string output_folder)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
	pad1->SetBottomMargin(0.12);
        pad1->SetRightMargin(0.08);
        pad1->SetLeftMargin(0.125);
	pad1->SetTopMargin(0.045);
        pad1->Draw();
        pad1->cd();

	//Define all necessary labels/titles
	const char *label_x;
	const char *label_y;
	std::string slabel_x = "Input " + xlabel;
	std::string slabel_y = "Pull " + ylabel;
	label_x = slabel_x.c_str();
	label_y = slabel_y.c_str();
	
	const char *title1, *title2;
	std::string stitle2 = output_folder+  "/" + xtitle + ".png";
	std::string stitle1 = output_folder + "/" + xtitle + ".eps";
	title1 = stitle1.c_str();
	title2 = stitle2.c_str();
	
	/*TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
	pad1->SetBottomMargin(0);
	pad1->SetRightMargin(-0.10);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	pad1->SetGrid();*/

        c_F->SetLeftMargin(0.12);
	//	c_F->SetGrid();
	//	graph->RemovePoint(0);
	graph->SetTitle("");
	graph->GetXaxis()->SetTitle(label_x);
	graph->GetYaxis()->SetTitle(label_y);
	graph->GetXaxis()->SetTitleSize(0.0505);
        graph->GetYaxis()->SetTitleSize(0.0505);
	graph->GetXaxis()->SetTitleOffset(1.05);
	graph->GetYaxis()->SetTitleOffset(1.175);
	graph->GetXaxis()->SetLabelSize(0.040);
	graph->GetYaxis()->SetLabelSize(0.040);
	graph->SetMinimum(ymin);
	graph->SetMaximum(ymax);
	graph->GetXaxis()->SetNdivisions(505);

	//graph->SetStats(kFALSE);
	graph->SetFillStyle(0);
	graph->SetLineColor(kOrange+9);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(1);
	graph->SetMarkerColor(kOrange+9);
	//	graph->GetXaxis()->SetLimits(xmin,xmax);
	graph->Draw("APZ");

	TF1 *fit = new TF1("fit", "[0]*x + [1]",xmin,xmax);
        fit->SetLineColor(kCyan+2);
        fit->SetLineWidth(2);
        graph->Fit("fit", "Q");
        fit->Draw("same");

        TLegend *legend = new TLegend(0.64,0.86,0.84,0.93);
        legend->SetTextFont(42);
        legend->AddEntry(graph,"#scale[1.24]{Pull values}","lep");
        legend->AddEntry(fit,"#scale[1.24]{Linear fit}","l");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->Draw("same");

	std::string sslope;
        sslope = "Slope:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *slope;
        slope = sslope.c_str();
        char latex_slope[sizeof slope + 100];
        sprintf(latex_slope, slope, fit->GetParameter(0), fit->GetParError(0) );

	std::cout << "Slope:  " << fit->GetParameter(0) << " $\\pm$ " << fit->GetParError(0) << std::endl;
	std::cout << "Offset: " << fit->GetParameter(1) << " $\\pm$ " << fit->GetParError(1) << std::endl;

	std::string sintercept;
        sintercept = "Offset:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *intercept;
        intercept = sintercept.c_str();
        char latex_intercept[sizeof intercept + 100];
        sprintf(latex_intercept, intercept,  fit->GetParameter(1), fit->GetParError(1));

	std::string schi2;
        schi2 = "#chi^{2}/ndof: %.2f / %d";// use %f for double, %d for int
        const char *chi2;
        chi2 = schi2.c_str();
        char latex_chi2[sizeof chi2 + 100];
        sprintf(latex_chi2, chi2,  fit->GetChisquare(), fit->GetNDF());

	TLatex l; //l.SetTextAlign(12);
        l.SetTextSize(0.032);
        l.SetNDC();
	l.SetTextFont(42);
	l.DrawLatex(0.54, 0.27, latex_chi2);
        l.DrawLatex(0.54, 0.215, latex_slope);
        l.DrawLatex(0.54, 0.175, latex_intercept);

	//	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.033); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.165,0.90,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.285,0.90,"Preliminary");
	l.SetTextSize(0.031);
	//	l.DrawLatex(0.14,0.79,"TMinuit");
	//	l.SetTextSize(0.038); 
	
	const char *mode;
	mode = smode.c_str();
	l.DrawLatex(0.26,0.79,mode);	
	//l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	//l.DrawLatex(0.14,0.74,"Calibration");

	if(input_channel == "el")    l.DrawLatex(0.165,0.855, "e+jets channel");
	if(input_channel == "mu")    l.DrawLatex(0.165,0.855, "#mu+jets channel");
	if(input_channel == "el_mu") l.DrawLatex(0.165,0.855, "Combined channel");

	c_F->Print(title1);
	c_F->Print(title2);

	//fit->Delete();
}


void PlotCalRMSDistribution (TGraphErrors* graph, std::string xlabel, std::string ylabel, std::string xtitle, double xmin, double xmax, double ymin, double ymax, std::string input_channel, std::string smode, std::string output_folder)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
	pad1->SetBottomMargin(0.12);
        pad1->SetRightMargin(0.08);
        pad1->SetLeftMargin(0.125);
	pad1->SetTopMargin(0.045);
        pad1->Draw();
        pad1->cd();

	//Define all necessary labels/titles
	const char *label_x;
	const char *label_y;
	std::string slabel_x = "Input " + xlabel;
	std::string slabel_y = "RMS " + ylabel;
	label_x = slabel_x.c_str();
	label_y = slabel_y.c_str();
	
	const char *title1, *title2;
	std::string stitle2 = output_folder+  "/" + xtitle + ".png";
	std::string stitle1 = output_folder + "/" + xtitle + ".eps";
	title1 = stitle1.c_str();
	title2 = stitle2.c_str();
	
	/*TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
	pad1->SetBottomMargin(0);
	pad1->SetRightMargin(-0.10);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	pad1->SetGrid();*/

        c_F->SetLeftMargin(0.12);
	//	c_F->SetGrid();
	//	graph->RemovePoint(0);
	graph->SetTitle("");
	graph->GetXaxis()->SetTitle(label_x);
	graph->GetYaxis()->SetTitle(label_y);
	graph->GetXaxis()->SetTitleSize(0.0505);
        graph->GetYaxis()->SetTitleSize(0.0505);
	graph->GetXaxis()->SetTitleOffset(1.05);
	graph->GetYaxis()->SetTitleOffset(1.175);
	graph->GetXaxis()->SetLabelSize(0.040);
	graph->GetYaxis()->SetLabelSize(0.040);
	graph->SetMinimum(ymin);
	graph->SetMaximum(ymax);
	graph->GetXaxis()->SetNdivisions(505);

	//graph->SetStats(kFALSE);
	graph->SetFillStyle(0);
	graph->SetLineColor(kOrange+9);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(1);
	graph->SetMarkerColor(kOrange+9);
	//	graph->GetXaxis()->SetLimits(xmin,xmax);
	graph->Draw("APZ");

	TF1 *fit = new TF1("fit", "[0]*x + [1]",xmin,xmax);
        fit->SetLineColor(kCyan+2);
        fit->SetLineWidth(2);
        graph->Fit("fit", "Q");
        fit->Draw("same");

        TLegend *legend = new TLegend(0.64,0.86,0.84,0.93);
        legend->SetTextFont(42);
        legend->AddEntry(graph,"#scale[1.24]{RMS values}","lep");
        legend->AddEntry(fit,"#scale[1.24]{Linear fit}","l");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->Draw("same");

	std::string sslope;
        sslope = "Slope:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *slope;
        slope = sslope.c_str();
        char latex_slope[sizeof slope + 100];
        sprintf(latex_slope, slope, fit->GetParameter(0), fit->GetParError(0) );

	std::cout << "Slope:  " << fit->GetParameter(0) << " $\\pm$ " << fit->GetParError(0) << std::endl;
	std::cout << "Offset: " << fit->GetParameter(1) << " $\\pm$ " << fit->GetParError(1) << std::endl;

	std::string sintercept;
        sintercept = "Offset:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *intercept;
        intercept = sintercept.c_str();
        char latex_intercept[sizeof intercept + 100];
        sprintf(latex_intercept, intercept,  fit->GetParameter(1), fit->GetParError(1));

	std::string schi2;
        schi2 = "#chi^{2}/ndof: %.2f / %d";// use %f for double, %d for int
        const char *chi2;
        chi2 = schi2.c_str();
        char latex_chi2[sizeof chi2 + 100];
        sprintf(latex_chi2, chi2,  fit->GetChisquare(), fit->GetNDF());

	TLatex l; //l.SetTextAlign(12);
        l.SetTextSize(0.032);
        l.SetNDC();
	l.SetTextFont(42);
	l.DrawLatex(0.54, 0.27, latex_chi2);
        l.DrawLatex(0.54, 0.215, latex_slope);
        l.DrawLatex(0.54, 0.175, latex_intercept);

	//	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.033); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.165,0.90,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.285,0.90,"Preliminary");
	l.SetTextSize(0.031);
	//	l.DrawLatex(0.14,0.79,"TMinuit");
	//	l.SetTextSize(0.038); 
	
	const char *mode;
	mode = smode.c_str();
	l.DrawLatex(0.26,0.79,mode);	
	//l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	//l.DrawLatex(0.14,0.74,"Calibration");

	if(input_channel == "el")    l.DrawLatex(0.165,0.855, "e+jets channel");
	if(input_channel == "mu")    l.DrawLatex(0.165,0.855, "#mu+jets channel");
	if(input_channel == "el_mu") l.DrawLatex(0.165,0.855, "Combined channel");

	c_F->Print(title1);
	c_F->Print(title2);

	//fit->Delete();
}


void PlotPileupDistribution (TGraphErrors* graph, std::string xlabel, std::string ylabel, std::string xtitle, double xmin, double xmax, double ymin, double ymax, std::string input_channel, std::string smode, std::string output_file)
{
	//Define TCanvas
	TCanvas *c_F = new TCanvas("c_F","c_F", 600, 600);
	c_F->SetFillColor(0);

	SetAtlasStyle();

	TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
	pad1->SetBottomMargin(0.12);
        pad1->SetRightMargin(0.08);
        pad1->SetLeftMargin(0.125);
	pad1->SetTopMargin(0.045);
        pad1->Draw();
        pad1->cd();

	//Define all necessary labels/titles
	const char *label_x;
	const char *label_y;
	std::string slabel_x = xlabel;
	std::string slabel_y = "Fraction " + ylabel;
	label_x = slabel_x.c_str();
	label_y = slabel_y.c_str();
	
	const char *title1, *title2;
	std::string stitle2 = output_file;
	std::string stitle1 = output_file;
	title1 = stitle1.c_str();
	title2 = stitle2.c_str();
	
	/*TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
	pad1->SetBottomMargin(0);
	pad1->SetRightMargin(-0.10);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy(1);

	pad1->SetGrid();*/

        c_F->SetLeftMargin(0.12);
	//	c_F->SetGrid();
	//	graph->RemovePoint(0);
	graph->SetTitle("");
	graph->GetXaxis()->SetTitle(label_x);
	graph->GetYaxis()->SetTitle(label_y);
	graph->GetXaxis()->SetTitleSize(0.0505);
        graph->GetYaxis()->SetTitleSize(0.0505);
	graph->GetXaxis()->SetTitleOffset(1.05);
	graph->GetYaxis()->SetTitleOffset(1.175);
	graph->GetXaxis()->SetLabelSize(0.040);
	graph->GetYaxis()->SetLabelSize(0.040);
	graph->SetMinimum(ymin);
	graph->SetMaximum(ymax);
	graph->GetXaxis()->SetNdivisions(505);

	//graph->SetStats(kFALSE);
	graph->SetFillStyle(0);
	graph->SetLineColor(kOrange+9);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(1);
	graph->SetMarkerColor(kOrange+9);
	//	graph->GetXaxis()->SetLimits(xmin,xmax);
	graph->Draw("APZ");

	TF1 *fit = new TF1("fit", "[0]*x + [1]",xmin,xmax);
        fit->SetLineColor(kCyan+2);
        fit->SetLineWidth(2);
        graph->Fit("fit", "Q");
        fit->Draw("same");

        TLegend *legend = new TLegend(0.64,0.86,0.84,0.93);
        legend->SetTextFont(42);
        legend->AddEntry(graph,"#scale[1.24]{Fractions}","lep");
        legend->AddEntry(fit,"#scale[1.24]{Linear fit}","l");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->Draw("same");

	std::string sslope;
        sslope = "Slope:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *slope;
        slope = sslope.c_str();
        char latex_slope[sizeof slope + 100];
        sprintf(latex_slope, slope, fit->GetParameter(0), fit->GetParError(0) );

	std::cout << xlabel.c_str() << "\t  Slope: " << fit->GetParameter(0) << "\t" << fit->GetParError(0) << std::endl;
	std::cout << xlabel.c_str() << "\t  Offset:" << fit->GetParameter(1) << "\t" << fit->GetParError(1) << std::endl;

	std::string sintercept;
        sintercept = "Offset:   %.3f #pm %.3f";// use %f for double, %d for int
        const char *intercept;
        intercept = sintercept.c_str();
        char latex_intercept[sizeof intercept + 100];
        sprintf(latex_intercept, intercept,  fit->GetParameter(1), fit->GetParError(1));

	std::string schi2;
        schi2 = "#chi^{2}/ndof: %.2f / %d";// use %f for double, %d for int
        const char *chi2;
        chi2 = schi2.c_str();
        char latex_chi2[sizeof chi2 + 100];
        sprintf(latex_chi2, chi2,  fit->GetChisquare(), fit->GetNDF());

	TLatex l; //l.SetTextAlign(12);
        l.SetTextSize(0.032);
        l.SetNDC();
	l.SetTextFont(42);
	l.DrawLatex(0.54, 0.27, latex_chi2);
        l.DrawLatex(0.54, 0.215, latex_slope);
        l.DrawLatex(0.54, 0.175, latex_intercept);

	//	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(0.033); 
	l.SetNDC();
	l.SetTextFont(72);
	//l.SetTextColor();
	l.DrawLatex(0.165,0.90,"ATLAS");
	l.SetTextFont(42);
	l.DrawLatex(0.285,0.90,"Preliminary");
	l.SetTextSize(0.031);
	//	l.DrawLatex(0.14,0.79,"TMinuit");
	//	l.SetTextSize(0.038); 
	
	const char *mode;
	mode = smode.c_str();
	l.DrawLatex(0.26,0.79,mode);	
	//l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	//l.DrawLatex(0.14,0.74,"Calibration");

	if(input_channel == "el")    l.DrawLatex(0.165,0.855, "e+jets channel");
	if(input_channel == "mu")    l.DrawLatex(0.165,0.855, "#mu+jets channel");
	if(input_channel == "el_mu") l.DrawLatex(0.165,0.855, "Combined channel");

	c_F->Print(title1);
	c_F->Print(title2);

	//fit->Delete();
}

void ShowRatio(TH1D* hData, TH1D* hExp, TH1D hExpUnc, Float_t lowerLimit, Float_t upperLimit)
{
  TH1D* hRatio(0);
  TH1D* hLineInOne(0);
  hRatio = (TH1D*)hData->Clone();
  hRatio->Divide(hExp);

  //for (int ibin=0; ibin<hRatio->GetNbinsX(); ibin++)
  //  cout<<"hRatio: "<<hRatio->GetBinContent(ibin+1)<<" +- "<<hRatio->GetBinError(ibin+1)<<endl;
  
  hLineInOne = (TH1D*)hData->Clone();
  hLineInOne->Divide(hData);
   for (int i=0; i<hData->GetNbinsX(); i++) {
  //     //ratio: data/exp. with hData unc.
      hRatio->SetBinContent(i+1,0.);   hRatio->SetBinError(i+1, 0.);
      if (hExp->GetBinContent(i+1)!=0.) {
        hRatio->SetBinContent(i+1, hData->GetBinContent(i+1)/hExp->GetBinContent(i+1));
        hRatio->SetBinError(i+1, float(hData->GetBinError(i+1)/hExp->GetBinContent(i+1))); //take into account only hData unc.
      }
      else { hRatio->SetBinContent(i+1, -1000.); }
  
      //line in 1 with hPrediction unc.
      hLineInOne->SetBinContent(i+1,0.);      hLineInOne->SetBinError(i+1,0.); 
      hLineInOne->SetBinContent(i+1, 1.);
      if (hExp->GetBinContent(i+1)!=0.)
        //hLineInOne->SetBinError(i+1, float(hExp->GetBinError(i+1)*hData->GetBinContent(i+1)/TMath::Power(hExp->GetBinContent(i+1),2)));  //take into account only hPrediction unc.
        hLineInOne->SetBinError(i+1, float(hExpUnc.GetBinContent(i+1)*hData->GetBinContent(i+1)/TMath::Power(hExp->GetBinContent(i+1),2)));  //take into account only hPrediction unc.
      else hLineInOne->SetBinContent(i+1, -1000.);
     }

  Float_t myFontSize = 0.18;
  
  hRatio->SetMarkerSize(1.2);  hRatio->SetMarkerStyle(8);    hRatio->SetMarkerColor(1);
  hRatio->SetLineWidth(1);
  
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetYaxis()->SetTitleOffset(0.25);
  hRatio->GetYaxis()->SetTitleSize(0.19);
  hRatio->GetYaxis()->SetLabelSize(myFontSize);

  hRatio->GetXaxis()->SetTitle("cos #theta*");
  hRatio->GetXaxis()->SetTitleOffset(0.97);
  hRatio->GetXaxis()->SetTitleSize(0.22);
  hRatio->GetXaxis()->SetLabelSize(myFontSize);

  
  
  
  
  //hRatio->GetYaxis()->SetRangeUser(lowerLimit+0.001, upperLimit-0.001);
  hRatio->GetYaxis()->SetRangeUser(lowerLimit, upperLimit);
  
  hRatio->GetYaxis()->SetNdivisions(804);
  //hRatio->GetXaxis()->SetNdivisions(0);
  //hRatio->GetYaxis()->SetTickLength(0.05);
  hRatio->GetXaxis()->SetTickLength(0.08);
  hRatio->GetXaxis()->SetLabelSize(0);
  hRatio->SetLineColor(1);
  hRatio->SetLineStyle(1);
  //hRatio->SetAxisRange(lowerLimit, upperLimit, "Y");
  //hRatio->SetAxisRange(lowerLimit+0.01, upperLimit-0.01, "Y");
  hRatio->Draw("P E");
  hLineInOne->SetMarkerSize(0.);
  hLineInOne->SetMarkerStyle(0);
  hLineInOne->SetFillColor(1);
  hLineInOne->SetLineColor(2);
  hLineInOne->SetFillStyle(3256);//3256);
  hLineInOne->Draw("L E2 SAME"); //"L E2 SAME");

  TLine* myLine = new TLine();
  myLine->SetLineStyle(2);
  //myLine->SetLineColor(kOrange);
  myLine->DrawLine(hRatio->GetXaxis()->GetXmin(), 1., hRatio->GetXaxis()->GetXmax(), 1.);
}



























