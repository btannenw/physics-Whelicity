/*
 functions.cxx
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
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TRandom3.h"

#include <Math/PdfFuncMathCore.h>


#include "functions.h"
#include "TemplateInfo.h"
#include "SystematicInfo.h"

using namespace std;

///______________________Global variables/functions_______________________

std::vector<std::vector<TemplateInfo::SigInterpolObject> > fSignalInterpolation;
std::vector<std::vector<TemplateInfo::BkgInterpolObject> > fBkgInterpolation;

std::string fInterpolMethod;

int fNSignal;
int fNBkg;
int fNNui;

void AddSignalSystematic(std::vector<TemplateInfo::MySignalTemplate> SigTempl, std::vector<SystematicInfo::MySystematic> Syst)
{

  //  fSignalInterpolation.SystematicsType = Syst[0].SystematicsType;

  //  fNParameters = ].ParamName.size();

  //  fNSignal2 = SigTempl.size();

  // create for each parameter a new InterpolationObject, store k-Values and corresponding Histos
  std::vector<TemplateInfo::SigInterpolObject> HelpVecObject;

  TemplateInfo::SigInterpolObject HelpObject;
  
  TH1D* HelpHist;
  double HelpK;

  for(int k = 0; k < SigTempl.size(); ++k){

    HelpObject.hist.push_back(SigTempl[k].hist);
    HelpObject.k.push_back(0.0);
    HelpObject.SelEff = SigTempl[k].SelEff;

    for(int l = 0; l < Syst.size(); ++l){

      if(l == 0) HelpObject.SystematicType = Syst[0].SystematicType;

      HelpObject.hist.push_back(Syst[l].HistUp[k]);
      HelpObject.k.push_back(Syst[l].kup);
      HelpObject.hist.push_back(Syst[l].HistDown[k]);
      HelpObject.k.push_back(Syst[l].kdown);

    }

    HelpVecObject.push_back(HelpObject);

    HelpObject.hist.clear();
    HelpObject.k.clear();
    
  }

  fSignalInterpolation.push_back(HelpVecObject);

}


void AddBkgSystematic(std::vector<TemplateInfo::MyBkgTemplate> BkgTempl, std::vector<SystematicInfo::MySystematic> Syst)
{

  //  fNBkg = BkgTempl.size();

  // create for each parameter a new InterpolationObject, store k-Values and corresponding Histos
  std::vector<TemplateInfo::BkgInterpolObject> HelpVecObject;

  TemplateInfo::BkgInterpolObject HelpObject;

  TH1D* HelpHist;
  double HelpK;

  for(int k = 0; k < BkgTempl.size(); ++k){

    HelpObject.hist.push_back(BkgTempl[k].hist);
    HelpObject.k.push_back(0.0);
    HelpObject.BkgExp = BkgTempl[k].Exp;
    HelpObject.BkgUnc = BkgTempl[k].Unc;

    for(int l = 0; l < Syst.size(); ++l){
      
      if(l == 0) HelpObject.SystematicType = Syst[0].SystematicType;

      HelpObject.hist.push_back(Syst[l].HistUp[k+fNSignal]);
      HelpObject.k.push_back(Syst[l].kup);
      HelpObject.hist.push_back(Syst[l].HistDown[k+fNSignal]);
      HelpObject.k.push_back(Syst[l].kdown);
      
    }

    HelpVecObject.push_back(HelpObject);

    HelpObject.hist.clear();
    HelpObject.k.clear();

  }

  fBkgInterpolation.push_back(HelpVecObject);

}



void SetInterpolationMethod(std::string interpol)
{

  fInterpolMethod = interpol;

};

std::string GetInterpolationMethod()
{
  return fInterpolMethod;
}


//Normalise histograms during fit                                                                                                                                    
TH1D * normalise (TH1D* hist)
{
  TH1D * norm_hist = (TH1D*) hist->Clone("hist");
  double norm_factor = hist->Integral();//or GetSumOfWeights                                                                                                   

  //cout << "norm old " << norm_hist << endl;                                                                                                                  

  norm_hist->Scale(1.0/norm_factor);

  //cout << "norm new " << hist->GetSumOfWeights() << endl;                                                                                                    

  return norm_hist;

}




//Normalise histograms during fit
TH1D normalise (TH1D hist)
{
	TH1D  norm_hist = hist;
	double norm_factor = hist.Integral();//or GetSumOfWeights
	
	//cout << "norm old " << norm_hist << endl;
	
	norm_hist.Scale(1.0/norm_factor);
	
	//cout << "norm new " << hist->GetSumOfWeights() << endl; 
	
	return norm_hist; 

}



//Normalise histograms during fit - without clone (much faster fit!!!)
TH1D normalise_fast (TH1D hist)
{
	double norm_factor = hist.Integral();//or GetSumOfWeights
	hist.Scale(1.0/norm_factor);

	return hist; 
}


void CreateDummyFiles (std::vector<TH1D*> &hist_nom_vec)
{

	std::vector <TH1D*> hist_nom_vec_clone_down;	
	std::vector <TH1D*> hist_nom_vec_clone_up;
	
	for (int i = 0; i != hist_nom_vec.size(); i++)
	{
		TH1D * hist_dummy_down = (TH1D*) hist_nom_vec[i]->Clone("hist");
		TH1D * hist_dummy_up = (TH1D*) hist_nom_vec[i]->Clone("hist");
		hist_nom_vec_clone_down.push_back(hist_dummy_down);
		hist_nom_vec_clone_up.push_back(hist_dummy_up);
	}
	
	for (int i = 0; i != hist_nom_vec.size(); i++)
	{
		//hist_nom_vec_clone_down[i]->Scale(0.8);
		//hist_nom_vec_clone_up[i]->Scale(1.0);
	
		for (int k =1; k != hist_nom_vec[0]->GetNbinsX()+1; k++)
		{
		double value_down= hist_nom_vec_clone_down[i]->GetBinContent(k);
		double value_up= hist_nom_vec_clone_up[i]->GetBinContent(k);
		
		if (k < 6) {
		hist_nom_vec_clone_down[i]->SetBinContent(k,0.99*value_down);
		hist_nom_vec_clone_up[i]->SetBinContent(k,1.01*value_up);
		}
		else if (k > 5 && k < 11) {
		hist_nom_vec_clone_down[i]->SetBinContent(k,0.995*value_down);
		hist_nom_vec_clone_up[i]->SetBinContent(k,1.005*value_up);
		}
		else {
		hist_nom_vec_clone_down[i]->SetBinContent(k,0.99*value_down);
		hist_nom_vec_clone_up[i]->SetBinContent(k,1.01*value_up);
		}
		}

	}

	//Create output file down
	TFile output_file_down("root_files/Templates_dummyd_mu.root", "recreate");

	hist_nom_vec_clone_down[0]->Write("CosTheta_F0");
	hist_nom_vec_clone_down[1]->Write("CosTheta_FL");
	hist_nom_vec_clone_down[2]->Write("CosTheta_FR");

	//	hist_nom_vec_clone_down[3]->Write("Wjets");/
	//	hist_nom_vec_clone_down[4]->Write("QCD");
	//	hist_nom_vec_clone_down[5]->Write("RemBkg");	

	//	hist_nom_vec_clone_down[3]->Write("Wjets");
        hist_nom_vec_clone_down[3]->Write("QCD");
        hist_nom_vec_clone_down[4]->Write("RemBkg");

	
	output_file_down.Close();

	//Create output file up
	TFile output_file_up("root_files/Templates_dummyu_mu.root", "recreate");

	hist_nom_vec_clone_up[0]->Write("CosTheta_F0");
	hist_nom_vec_clone_up[1]->Write("CosTheta_FL");
	hist_nom_vec_clone_up[2]->Write("CosTheta_FR");
	//	hist_nom_vec_clone_up[3]->Write("Wjets");
	//	hist_nom_vec_clone_up[4]->Write("QCD");
	//	hist_nom_vec_clone_up[5]->Write("RemBkg");	

	//	hist_nom_vec_clone_up[3]->Write("Wjets");
        hist_nom_vec_clone_up[3]->Write("QCD");
        hist_nom_vec_clone_up[4]->Write("RemBkg");

	
	output_file_up.Close();
}
