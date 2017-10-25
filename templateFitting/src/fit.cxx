/*
 fit.cxx
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

#include "interpolation.h"
#include "functions.h"
#include "ProfilingClass.h"
#include "ValidationClass.h"
#include "fit_validation.h"
#include "fit.h"
#include "BAT/BCMath.h"
//#include "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/RootCoreBin/obj/x86_64-slc6-gcc48-opt/Asg_BAT/lib/local/include/BAT/BCMath.h"

using namespace std;  

///_____________________________FIT FUNCTION_______________________________



//minuit will call eval_chisqrt and it will return the chi-sqrt between the TInterpolation object (with parameters k=par[...]) and the datahist
void eval_chisqrt(int &npar, double *gin, double &f, double *par, int iflag)
{
		
	std::vector<TH1D> fParameterHistNorm;
	std::vector<TH1D> fInterpolHist;
	std::vector<TH1D> fInterpolHistNorm;

	fParameterHistNorm.clear();
	fInterpolHist.clear();
	fInterpolHistNorm.clear();

	//Fit parameters	
	std::vector<double> k;
	k.clear();
	if(fBkgMode == "BkgFit"){
	    
	  for (int i = fNBkg+fNSignal; i != fNNui+fNBkg+fNSignal; i++)
	    {
	      k.push_back(par[i]);
	    }
	}
	else{

	  for (int i = fNSignal; i != fNNui+fNSignal; i++)
            {
              k.push_back(par[i]);
            }


	}

	const int ksize = k.size();

	std::string interpol_method = GetInterpolationMethod();

	if(ksize > 0){

	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam){
	    fInterpolHist.push_back(inter_obj[iParam].interpolate(interpol_method, k));	

	    //   std::cout << fInterpolHist[iParam].Integral() << std::endl;
	    

	    
	    
	  }


	}
	else{
	  
	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam){
            fInterpolHist.push_back(fParameterHist[iParam]);
	    
	    //	    std::cout << iParam << "\t" << fParameterHist[iParam].GetBinLowEdge(1) << std::endl; 
	    
          }
	  
	  
	}

	std::vector<double> norm_factors;
	norm_factors.clear();

	for(int iParam = 0; iParam < fNSignal; ++iParam){
	  
	  norm_factors.push_back(fParameterHist[iParam].Integral());

	}

	if(ksize > 0){ 

	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam)
	    fInterpolHistNorm.push_back(normalise_fast(fInterpolHist[iParam]));
	}
	else{


	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam)
            fParameterHistNorm.push_back(normalise_fast(fParameterHist[iParam]));


	}

	//	std::cout << "\t" << std::endl;


	//for(int i = 0; i < fNSignal+fNBkg+fNNui; ++i)
	//  std::cout << i << "\t" << par[i] << std::endl;
	  
	TH1D hist_sum = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

	int gNBins = hist_sum.GetNbinsX();

	if(ksize > 0){

	  for(int iParam = 0; iParam < fNSignal; ++iParam){
	    	
	    double new_sel_eff = SelEff[iParam];

	    //	    double new_sel_eff = fInterpolHist[iParam].Integral()/415307.97;
		//cout << fInterpolHist[iParam].Integral() << endl;
	    if(iParam == 0){
	      
	      for(int l = 1; l <= fNBins; ++l)
	       	hist_sum.SetBinContent(l, 0);
	      
	    }
	    
	    hist_help = fInterpolHistNorm[iParam];

	    hist_help.Scale(par[iParam]*new_sel_eff);
	    
	    //cout << "Hier 1" << endl;

	    hist_sum = hist_sum + hist_help;

	    //cout << "Hier 2" << endl;


	//cout << hist_sum.GetBinContent(1) << endl;
	    
	  }
	  for(int iParam = 0; iParam < fNBkg; ++iParam){

	    if(fBkgMode == "BkgFit"){
	          
	      hist_help = fInterpolHistNorm[iParam+fNSignal];
	      hist_help.Scale(par[iParam+fNSignal]);
	          
	    }
	    else{

	      hist_help = fInterpolHist[iParam+fNSignal];

	    }
	    	    
	    hist_sum = hist_sum + hist_help;
    	//	cout << hist_sum.GetBinContent(1) << endl;
	  }
	}
	else{

	  for(int iParam = 0; iParam < fNSignal; ++iParam){
	    
            if(iParam == 0){
	      
              for(int l = 1; l <= fNBins; ++l)
		hist_sum.SetBinContent(l, 0.0);
	      
	    }
	    
	    //std::cout << "Sel Eff!!!!" << "\t" << SelEff[iParam] << std::endl;


	    hist_help = fParameterHistNorm[iParam];
            hist_help.Scale(par[iParam]*SelEff[iParam]);
	    
	    //	    hist_help.Scale(par[iParam]);

	    //cout << "Hier 1a" << endl;


	    hist_sum = hist_sum + hist_help;

	    //cout << "Hier 1b" << endl;


          }
          for(int iParam = 0; iParam < fNBkg; ++iParam){

            if(fBkgMode == "BkgFit"){
	          
	      //	      std::cout << fParameterHistNorm[iParam+fNSignal].Integral() << std::endl;
	      
	      hist_help = fParameterHistNorm[iParam+fNSignal];
	      hist_help.Scale(par[iParam+fNSignal]);
	      
	      //	      std::cout << hist_sum.Integral() << "\t" <<  hist_help.Integral() << std::endl;

	      //  cout << hist_help.GetNbinsX() << "\t" << hist_help.GetXaxis()->GetXmin() << "\t" << hist_help.GetXaxis()->GetXmax()<< endl;
	     	      
		
	    }
	    else{

	      hist_help = fInterpolHist[iParam+fNSignal];

	    }
	 
	    /////	    cout << "Hier 1c   " << iParam << endl;
   
            hist_sum = hist_sum + hist_help;
    
	    //	    cout << "Hier 1d   " <<  iParam << endl;


	  }
	}

	//FIT
	long double like = 0.0; //like corresponds to -2*ln(L), error poisson sqrt

	//	std::cout << "Data integral = " << hist_data.Integral() << "\t" << hist_sum.Integral() << std::endl;

	//Calculation of likelihood
	for (int i = 1; i != fNBins+1; i++){
	  
	  double val   = hist_sum.GetBinContent(i);
	  double data  = hist_data.GetBinContent(i);
	  double sigma = sqrt(data);
		//cout << i << " " << val << endl;


	  //	  std::cout << hist_sum.GetBinLowEdge(1) << "\t" << hist_data.GetBinLowEdge(1) << std::endl;

	  if(val > 0)
	    like -= BCMath::LogPoisson(data, val);
	  else
	    like += 10e5;
	
	}	
	


	std::vector<double> num_bkg;
	num_bkg.clear();
	for(int i = fNSignal; i < fNSignal+fNBkg; ++i)
	  num_bkg.push_back(par[i]);

	if(fBkgMode == "BkgFit"){
	
	  //Gaussian constraints on background
	  for (int i = 0; i != fNBkg; i++)
	    {
	      like -= BCMath::LogGaus(num_bkg[i], fBkgNorm[i], fBkgUnc[i], true);
	      
	    }
	
	}
	

	for (int i = 0; i != ksize; i++)//start with 0 since not sure whether cos_theta without penalty
	  {
	    
	    like -= BCMath::LogGaus(k[i], 0.0, 1.0, true);

	  }

	//cout << "par : " << par[0] << " : " << par[3] << " : " << par[6] << endl; 

	f = 2.0*like;
	
}


