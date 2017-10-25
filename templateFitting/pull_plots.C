#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <math.h>
#include "TImage.h"
#include "TLatex.h"

#include <fstream>
#include <iostream>
#include <sstream>
void pull_plots()
{
	
	TGraphErrors *graph_pull = new TGraphErrors(3);
	TGraphErrors *graph_pull_RMS = new TGraphErrors(3);
  std::string lepton = "el";
  std::string frac = "FR"; // F0, FL , FR 


  if(lepton=="el")
  {
    float F0_mean[7]=    {0.00, -0.01, 0.01, 0.02, -0.02, 0.02, 0.03};
    float F0_mean_err[7]= {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.00};
    float F0_sigma[7]=    {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
    float F0_sigma_err[7]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0};


  float FL_mean[7]= {0.00, 0.01, -0.01, -0.02, 0.01, -0.02,-0.03 };
  float FL_mean_err[7]= {0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01};
  float FL_sigma[7]=  {1.00, 0.99,  1.00, 1.00, 1.00, 0.99, 0.99};
  float FL_sigma_err[7]={0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01};

  float FR_mean[7]=     {0.00,  0.00,  -0.01,   -0.02, 0.03, -0.01, -0.02};
  float FR_mean_err[7]= {0.01,   0.01,  0.01,   0.01, 0.01, 0.01,   0.01};
  float FR_sigma[7]=    {1.0,    1.00,  1.00,   0.99, 1.00, 1.00,   1.00};
  float FR_sigma_err[7]={0.01,   0.01,  0.01,   0.01, 0.01, 0.01,   0.01};  
  }

	

  if(lepton=="mu")
  {
  float F0_mean[7]=      {-0.01, -0.01, 0.02, 0.01, 0.00, 0.03, 0.04};
  float F0_mean_err[7]=  {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  float F0_sigma[7]=     {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
  float F0_sigma_err[7]= {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};


  float FL_mean[7]=     {0.00, 0.01, -0.01, -0.01, 0.01, -0.02,-0.04 };
  float FL_mean_err[7]= {0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01};
  float FL_sigma[7]=    {1.00, 1.00,  1.00, 1.00, 1.00, 1.00, 1.00};
  float FL_sigma_err[7]={0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01};

  float FR_mean[7]=     {0.02,  0.00,  -0.03,   -0.01, 0.00, -0.02, -0.03};
  float FR_mean_err[7]= {0.01,   0.01,  0.01,   0.01, 0.01, 0.01,   0.01};
  float FR_sigma[7]=    {1.0,    1.00,  1.00,   1.00, 1.00, 1.00,   1.00};
  float FR_sigma_err[7]={0.01,   0.01,  0.01,   0.01, 0.01, 0.01,   0.01};
  }


//-------------------
	float F0frac[7]={0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85};
	float FLfrac[7]={0.447,0.398 ,0.349 ,0.30 ,0.251, 0.202, 0.153};
	float FRfrac[7]={0.003, 0.002, 0.001, 0.0, -0.001, -0.002, -0.003};

	double min_y =-0.19 ;
  double max_y =1.35 ;

  	double min_x ;
  	double max_x ;


	
	char lfrac;
	

	if(frac == "F0"){
  		min_x = 0.5;
  		max_x = 0.9;

		for(int i=0;i<7;i++)
		{
		graph_pull->SetPoint(i+1,F0frac[i],F0_mean[i]);
		graph_pull->SetPointError(i+1,0,F0_mean_err[i]);
		
		graph_pull_RMS->SetPoint(i+1,F0frac[i],F0_sigma[i] );
		graph_pull_RMS->SetPointError(i+1,0,F0_sigma_err[i]);

		}
	}
	
	else if(frac == "FL"){
		min_x = 0.15;
		max_x = 0.5;

		for(int i=0;i<7;i++)
		{	
		graph_pull->SetPoint(i+1,FLfrac[i],FL_mean[i]);
		graph_pull->SetPointError(i+1,0,FL_mean_err[i]);
		
		graph_pull_RMS->SetPoint(i+1,FLfrac[i],FL_sigma[i] );
		graph_pull_RMS->SetPointError(i+1,0,FL_sigma_err[i]);

		}	
	}
	
	else if(frac=="FR")
	{
		min_x = -0.004;
		max_x = 0.004;

		for(int i=0;i<7;i++)
		{	
		graph_pull->SetPoint(i+1,FRfrac[i],FR_mean[i]);
		graph_pull->SetPointError(i+1,0,FR_mean_err[i]);
		
		graph_pull_RMS->SetPoint(i+1,FRfrac[i],FR_sigma[i] );
		graph_pull_RMS->SetPointError(i+1,0,FR_sigma_err[i]);

		}	
	}

	
	//Pull plots:
  
  	TCanvas *c0 = new TCanvas("c0","c0",800,600);
  	c0->cd();

  graph_pull->RemovePoint(0);
  graph_pull->SetTitle("");
  
  if(frac=="F0")
  	graph_pull->GetXaxis()->SetTitle("#scale[1.0]{F_{0}}");
  else if(frac=="FL")
  	graph_pull->GetXaxis()->SetTitle("#scale[1.0]{F_{L}}");
  else if(frac=="FR")
  	graph_pull->GetXaxis()->SetTitle("#scale[1.0]{F_{R}}");

  graph_pull->GetYaxis()->SetTitle("#scale[0.8]{ [a.u.] }");
  graph_pull->GetXaxis()->SetTitleOffset(1.0);
  graph_pull->GetYaxis()->SetTitleOffset(1.0);
  graph_pull->GetXaxis()->SetLabelSize(0.045);
  graph_pull->GetYaxis()->SetLabelSize(0.045);
  graph_pull->SetMinimum(min_y);
  graph_pull->SetMaximum(max_y);
  graph_pull->GetXaxis()->SetLimits(min_x,max_x);
  graph_pull->SetFillStyle(0);
  graph_pull->SetLineColor(kOrange+1);
  graph_pull->SetLineWidth(2);
  graph_pull->SetMarkerStyle(20);
  graph_pull->SetMarkerColor(kOrange+1);
  graph_pull->Draw("APZ");

  graph_pull_RMS->RemovePoint(0);
  graph_pull_RMS->SetFillStyle(0);
  graph_pull_RMS->SetLineColor(kOrange+9);
  graph_pull_RMS->SetLineWidth(2);
  graph_pull_RMS->SetMarkerStyle(20);
  graph_pull_RMS->SetMarkerColor(kOrange+9);
  graph_pull_RMS->Draw("PZsame");
  

  TGraph *graph_line_pull = new TGraph(2);
  graph_line_pull->SetPoint(1,-1.0,0);
  graph_line_pull->SetPoint(2,80,0);
  graph_line_pull->RemovePoint(0);
  graph_line_pull->SetFillStyle(0);
  graph_line_pull->SetLineColor(kOrange+1);
  graph_line_pull->SetLineWidth(2);
  graph_line_pull->SetMarkerStyle(20);
  graph_line_pull->SetMarkerColor(kOrange+1);
  graph_line_pull->SetLineStyle(2);
  graph_line_pull->Draw("Lsame");
  
 
  TGraph *graph_line_RMS = new TGraph(2);
  graph_line_RMS->SetPoint(1,-1.0,1);
  graph_line_RMS->SetPoint(2,80,1);
  graph_line_RMS->RemovePoint(0);
  graph_line_RMS->SetFillStyle(0);
  graph_line_RMS->SetLineColor(kOrange+9);
  graph_line_RMS->SetLineWidth(2);
  graph_line_RMS->SetMarkerStyle(20);
  graph_line_RMS->SetMarkerColor(kOrange+9);
  graph_line_RMS->SetLineStyle(2);
  graph_line_RMS->Draw("Lsame");
 
  TLegend *legend_pull = new TLegend(0.69,0.5,0.90,0.6);
  legend_pull->SetTextFont(42);
  legend_pull->AddEntry(graph_pull,"#scale[1.0]{Pull mean}","lep");
  legend_pull->AddEntry(graph_pull_RMS,"#scale[1.0]{Pull sigma}","lep");
  legend_pull->SetFillColor(0);
  legend_pull->SetBorderSize(0);
  legend_pull->Draw("same");
 
  TLatex l;
  l.SetNDC(1);
  l.SetTextFont(42);
  
  TLatex atlas;
  atlas.SetNDC();
  atlas.SetTextFont(72);

  TLatex atlastext; 
  atlastext.SetNDC();
  atlastext.SetTextFont(42);

  atlas.SetTextSize(0.040);
  atlas.DrawLatex(0.15,0.85, "ATLAS");
  atlastext.SetTextSize(0.040);
  atlastext.DrawLatex(0.25,0.85, "Simulation work in progress");
  l.SetTextSize(0.040);
  l.DrawLatex(0.15,0.79,"#sqrt{s} = 8 TeV");
  if(lepton=="el")
    l.DrawLatex(0.62,0.79,"e+jets");
  else
    l.DrawLatex(0.62,0.79,"#mu+jets");
  if(lepton == "el")
  {
    if(frac=="F0")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_el_2incl_F0.png");
  else if(frac=="FL")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_el_2incl_FL.png");
  else if(frac=="FR")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_el_2incl_FR.png");  
  }
  else
  {
    if(frac=="F0")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_mu_2incl_F0.png");
  else if(frac=="FL")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_mu_2incl_FL.png");
  else if(frac=="FR")
    c0->Print("/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/pull_curve_mu_2incl_FR.png");   
  }

}