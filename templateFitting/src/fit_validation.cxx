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
#include "ExternalSystematics.h"
#include "fit_validation.h"
//#include "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/RootCoreBin/obj/x86_64-slc6-gcc48-opt/Asg_BAT/lib/local/include/BAT/BCMath.h"
#include "BAT/BCMath.h"

using namespace std;  

///_____________________________FIT FUNCTION_______________________________


//Pseudo-histogram used for fitting
TH1D hist_ensemble; // extern TH1D hist_ensemble

//TH1D hist      = TH1D("hist", "hist", gNBins, -1.0, 1.0);
//TH1D hist_sum  = TH1D("hist", "hist", gNBins, -1.0, 1.0);
//TH1D hist_help = TH1D("hist", "hist", gNBins, -1.0, 1.0);

//cout << "Number of bins " << gNBins << endl;

//minuit will call eval_chisqrt and it will return the chi-sqrt between the TInterpolation object (with parameters k=par[...]) and the datahist
void eval_chisqrt_vali(int &npar, double *gin, double &f, double *par, int iflag)
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
	
	if(fBkgMode == "BkgFit" || fBkgMode == ""){
	  
	  for (int i = fNBkg+fNSignal; i != fNNui+fNBkg+fNSignal; i++) // fNNui  = fSystematicFiles.size()+fNBkg; (Number of Nuisance Parameter)
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


	const int ksize = k.size(); // mjk: ksize=0 by using EvaluateExternalUpdate

	std::string interpol_method = GetInterpolationMethod();
	
	

	if(ksize > 0){

	  for(int iParam = 0; iParam < inter_obj.size(); ++iParam){
	    fInterpolHist.push_back(inter_obj[iParam].interpolate(interpol_method, k));	

	  }


	}
	else{
		//fParameterHist.size() = 6 (3 signal + 3 Bkg)
	  for(int iParam = 0; iParam < fParameterHist.size(); ++iParam){ // validationClass.h: extern std::vector<TH1D> fParameterHist (outside of the class)
            fInterpolHist.push_back(fParameterHist[iParam]);

	    //	    std::cout << iParam << "\t" << fParameterHist[iParam].Integral() << std::endl;

          }


	}

	std::vector<double> norm_factors;
	norm_factors.clear();

	for(int iParam = 0; iParam < fNSignal; ++iParam){
	  
	  norm_factors.push_back(fParameterHist[iParam].Integral());

	}

	if(ksize > 0){ 

	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam)
	    fInterpolHistNorm.push_back(normalise_fast(fInterpolHist[iParam])); // normalise to 1
	}
	else{

	  for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam){
            
	    fParameterHistNorm.push_back(normalise_fast(fParameterHist[iParam]));
	    
	    //	    std::cout << iParam << "\t" << fParameterHistNorm[iParam].Integral() << std::endl;
	    
	  }

	}

	if(ksize > 0){

	  for(int iParam = 0; iParam < fNSignal; ++iParam){
	    
	        double new_sel_eff = SelEff[iParam];
	    
	    //	    double new_sel_eff = fInterpolHist[iParam].Integral()/415307.97; // new Lumi
	    
	    
	    if(iParam == 0){
	      
	      for(int l = 1; l <= fNBins; ++l)
	       	hist_sum.SetBinContent(l, 0);
	      
	    }
	    
	    hist_help = fInterpolHistNorm[iParam];
	    
	    //	    par[iParam] = 4713.11*89.2;

	    hist_help.Scale(par[iParam]*new_sel_eff);
	    
	    hist_sum = hist_sum + hist_help;
	    
	  }
	  for(int iParam = 0; iParam < fNBkg; ++iParam){
	    
	    if(fBkgMode == "BkgFit"){
	    //std::cout<<"######## fBkgMode is BkgFit" << std::endl;
	      hist_help = fInterpolHistNorm[iParam+fNSignal];
	      hist_help.Scale(par[iParam+fNSignal]);
	    
	    }
	    else{
	    //	std::cout<<"######## fBkgMode is not BkgFit" << std::endl;
	      hist_help = fInterpolHist[iParam+fNSignal];

	    }
	 
	    
	    hist_sum = hist_sum + hist_help;
	    
	  }
	}
	else{ // ksize= 0
	  for(int iParam = 0; iParam < fNSignal; ++iParam){
	    
            if(iParam == 0){
	      
              for(int l = 1; l <= fNBins; ++l)
		hist_sum.SetBinContent(l, 0.0);
	      
	    }

	    hist_help = fParameterHistNorm[iParam];

        hist_help.Scale(par[iParam]*SelEff[iParam]);
	    
	    // std::cout<< "***** par[iParam]= " << par[iParam] << std::endl;
	    // std::cout<< "***** Scale= " << par[iParam]*SelEff[iParam] <<std::endl;

	    hist_sum = hist_sum + hist_help;

       }
       //std::cout<< std::endl;

          for(int iParam = 0; iParam < fNBkg; ++iParam){

           if(fBkgMode == "BkgFit" || fBkgMode == ""){
	      hist_help = fParameterHistNorm[iParam+fNSignal];
	      hist_help.Scale(par[iParam+fNSignal]);
	          
	    }
	    else{

	      hist_help = fInterpolHist[iParam+fNSignal];

	    }
	    
            hist_sum = hist_sum + hist_help;

	  }

	  //	     std::cout << hist_sum.Integral() << std::endl;
	}



	//for(int iParam = 0; iParam < fNSignal+fNBkg+fNNui; ++iParam)
	//std::cout << "hist sum  " << hist_sum.Integral() << "\t" << hist_ensemble.Integral()  << std::endl;


	//FIT
	long double like = 0.0; //like corresponds to -2*ln(L), error poisson sqrt
	
	//Calculation of likelihood
	for (int i = 1; i != fNBins+1; i++){
	  
	  double val   = hist_sum.GetBinContent(i);
	  double data  = hist_ensemble.GetBinContent(i);
	  double sigma = sqrt(data);

	  //std::cout<< "***** val= " << val << std::endl;
	  //std::cout<< "***** data= " << data << std::endl;

	  if(val > 0)
	    like -= BCMath::LogPoisson(data, val); // LogPoisson::Calculate the natural logarithm of a poisson distribution
	  else
	    like += 10e10;

	}	

	//	std::cout << "like sig " << like << std::endl;

	
	std::vector<double> num_bkg;
	num_bkg.clear();
	for(int i = fNSignal; i < fNSignal+fNBkg; ++i)
	  num_bkg.push_back(par[i]);
	
	if(fBkgMode == "BkgFit" || fBkgMode == ""){
	  //Gaussian constraints on background
	  for (int i = 0; i != fNBkg; i++)
	    {
	      
	      //    std::cout << fBkgNorm[i] << std::endl;
	      like -= BCMath::LogGaus(num_bkg[i], fBkgNorm[i], fBkgUnc[i], true); // Calculate the natural logarithm of a gaussian function with mean and sigma. If norm=true (default is false) the result is multiplied by the normalization constant, i.e. divided by sqrt(2*Pi)*sigma.
	      
	      
	    }
	}

	//	std::cout << "like bkg " << like << std::endl;

	for (int i = 0; i != ksize; i++)//start with 0 since not sure whether cos_theta without penalty
	{
	  like -= BCMath::LogGaus(k[i], 0.0, 1.0, true);

	}

	//	std::cout << "like final   " << like << std::endl;

	f = 2.0*like;
	
}


