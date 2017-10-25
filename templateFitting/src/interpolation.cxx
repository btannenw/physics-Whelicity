/*
 interpolation.cxx
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
#include <cstdlib> 

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
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMinuit.h"


#include "interpolation.h"
#include "fit.h"
#include "functions.h"
#include "TemplateInfo.h"
#include "SystematicInfo.h"
#include "ProfilingClass.h"
#include "ValidationClass.h"
#include "StatusLogbook.h"

using namespace std;

double kfix = 10.0;

//Constructor 
//tinterpolation::tinterpolation(TH1D **hist, const int hsize, int mode)
tinterpolation::tinterpolation(std::string TemplateType, int ParamNumber)
{

  _hist.clear();
  _hist_2k.clear();
  _hist_3k.clear();

  if(TemplateType == "Signal"){

    _hist.push_back(fSignalInterpolation[0][ParamNumber].hist[0]);

    fParamNumber = ParamNumber;
    
    for (int iSys = 0; iSys < fSignalInterpolation.size(); ++iSys)
      {

	_hist.push_back(fSignalInterpolation[iSys][ParamNumber].hist[2]); // down
        _hist.push_back(fSignalInterpolation[iSys][ParamNumber].hist[1]); // up

        if(fSignalInterpolation[iSys][ParamNumber].hist.size() > 3){

          _hist_2k.push_back(fSignalInterpolation[iSys][ParamNumber].hist[4]); // down
          _hist_2k.push_back(fSignalInterpolation[iSys][ParamNumber].hist[3]); // up

        }
        if(fSignalInterpolation[iSys][ParamNumber].hist.size() > 5){

          _hist_3k.push_back(fSignalInterpolation[iSys][ParamNumber].hist[6]); // down
          _hist_3k.push_back(fSignalInterpolation[iSys][ParamNumber].hist[5]); // up

        }

      }
  }
  else{

    _hist.push_back(fBkgInterpolation[0][ParamNumber].hist[0]);

    // has to be made more general!!!
    fParamNumber = ParamNumber + 3;

    for (int iSys = 0; iSys < fBkgInterpolation.size(); ++iSys)
      {
	_hist.push_back(fBkgInterpolation[iSys][ParamNumber].hist[2]); // down
        _hist.push_back(fBkgInterpolation[iSys][ParamNumber].hist[1]); // up

        if(fBkgInterpolation[iSys][ParamNumber].hist.size() > 3){

          _hist_2k.push_back(fBkgInterpolation[iSys][ParamNumber].hist[4]); // down
          _hist_2k.push_back(fBkgInterpolation[iSys][ParamNumber].hist[3]); // up

        }
        if(fSignalInterpolation[iSys][ParamNumber].hist.size() > 5){

          _hist_3k.push_back(fBkgInterpolation[iSys][ParamNumber].hist[6]); // down
          _hist_3k.push_back(fBkgInterpolation[iSys][ParamNumber].hist[5]); // up

        }


      }


  }

}

tinterpolation::~tinterpolation()
{

}

//Interpolate: Returns a histogram using k as the interpolation parameter (-1: hist_low, +1: hist_high, 0: hist)
TH1D tinterpolation::interpolate(std::string interpol_method, std::vector<double> &k)
{

  //  WriteParameterStatus("interpolation", "Using the interpolation method "+fInterpolationMethod);
  
  //  fNBins = _hist[0].GetNbinsX();

  TH1D result00 = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

        const int ksize = k.size();
	
	//	int nbins = _hist[0].GetNbinsX()+1;

	for (int i =1; i != fNBins+1; i++)
	{
	
	  	//cout << "interp i = " << i << endl; 
		
		double value   = 0.0;//0 for linear interpolation, 1 for quadratic!!!
		double central = _hist[0].GetBinContent(i);
		
		if(fInterpolationMethod == "PiecewiseLinear"){

		  //Old linear interpolation
		  for(int j = 0; j != ksize; j++){

		    //linear interpolation
		    if (k[j]>0.0)
		      {
			double slope_pos=_hist[2*j+2].GetBinContent(i);
			value += (slope_pos-central)*k[j]; 
		      }
		    else if (k[j]<0.0)
		      {
			double slope_neg=_hist[2*j+1].GetBinContent(i);
			value += (central-slope_neg)*k[j];
		      }
		    else if (k[j]==0.0)
		      {
			value += 0.0; //_hist[0]->GetBinContent(i);
		      }	
		  }
		    value += central;
	
/*	  
		  double central_func = fFitFunc[0][fParamNumber][i-1].Eval(0);

		  for(int j = 0; j != ksize; j++){

		    value += fFitFunc[j][fParamNumber][i-1].Eval(k[j]) - central_func;
		  }
		  value += central_func;
*/
		}
		else if(fInterpolationMethod == "LinearInterp"){

		  //New linear interpolation
		  // ynom + 0.5*delta*(yup - ydown) 
		    for(int j = 0; j != ksize; j++)
		    {
		    double slope_pos=_hist[2*j+2].GetBinContent(i);
		    double slope_neg=_hist[2*j+1].GetBinContent(i);
		    
		    value += 0.5*k[j]*(slope_pos-slope_neg);
		    }
		    value += central;

/*		   double central_func = fFitFunc[0][fParamNumber][i-1].Eval(0);

		  for(int j = 0; j != ksize; j++){

		    value += fFitFunc[j][fParamNumber][i-1].Eval(k[j]) - central_func;
			}
		  value += central_func;
*/
		}
		else if(fInterpolationMethod == "QuadraticInterp"){
/*		  
		  for(int j = 0; j != ksize; j++){
		    
		    double slope_pos=_hist[2*j+2].GetBinContent(i);
		    double slope_neg=_hist[2*j+1].GetBinContent(i);
		    
		    if (k[j]>1.0)
                      value += slope_pos + (k[j]-1.0)*(1.5*slope_pos + 0.5*slope_neg-2.0*central ) - central;
                    else if(k[j]<-1.0)
                      value += slope_neg + (k[j]+1.0)*(-0.5*slope_pos - 1.5*slope_neg + 2.0*central ) - central;
                    else
                      value += k[j]*k[j]*( (slope_neg+slope_pos)/2 - central) + k[j]*(slope_pos-slope_neg)/2;

		  }
		  
		  value += central;
*/		  

		  double central_func = fFitFunc[0][fParamNumber][i-1].Eval(0);

		  for(int j = 0; j != ksize; j++){

		    value += fFitFunc[j][fParamNumber][i-1].Eval(k[j]) - central_func;
		    
		  }

		  value += central_func;


		}
		else if(fInterpolationMethod == "QuadraticFit"){

		  double central_func = fFitFunc[0][fParamNumber][i-1].Eval(0);

		  for(int j = 0; j != ksize; j++){

		    value += fFitFunc[j][fParamNumber][i-1].Eval(k[j]) - central_func;
		    
		  }

		  value += central_func;
		  
		}
		else{

		  WriteErrorStatus("interpolation", "No Interpolation Method defined!");
		  WriteErrorStatus("interpolation", "EXIT");

		  exit(1);

		}

		result00.SetBinContent(i,value);
		//cout << value << endl;

	}


	return result00;
}

// return histograms for final plot with final conditions
TH1D ReturnInterpolationHistogramScaled(int ParamNumber, std::string interpol_method, std::vector<double> &k, double norm, bool BkgNui)
{
  
  TH1D result00 = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
  
  const int ksize = k.size();
  
  for (int i =1; i != fNBins+1; i++){

    double value = 0.0;
    
    //if(fInterpolationMethod == "QuadraticFit"){
      
      // indices are: Systematic, Parameter, Bin number

      double central_func = fFitFunc[0][ParamNumber][i-1].Eval(0);
      
      for(int j = 0; j != ksize; j++){
	
	value += fFitFunc[j][ParamNumber][i-1].Eval(k[j]) - central_func;
	
      }
      
      value += central_func;
      
    //}
    /*else{
      
      WriteErrorStatus("interpolation", "No Interpolation Method defined!");
      WriteErrorStatus("interpolation", "EXIT");

      exit(1);
      
    }*/
    
    
    result00.SetBinContent(i,value);
    
    
  }

  if(!BkgNui) {
    double seleff = result00.Integral()/415307.97;
    double integral = result00.Integral();
    result00.Scale(norm*seleff/integral);
  }
  else {
    if(fBkgMode == "BkgFit"){
    double integral = result00.Integral();
    result00.Scale(norm/integral);
    }	
  } 
  return result00;

 
}


//Interpolation function to plot contributions of single nuisance parameters
TH1D * interpolate_single(std::vector <TH1D*> &hist_vec, TH1D* result00, double k)
{	
	double value;
	
	for (int i =1; i != hist_vec[0]->GetNbinsX()+1; i++)
	{
		
		value =0.0;//0 for linear interpolation, 1 for relative!!!
		double central=hist_vec[0]->GetBinContent(i);
		
			//linear interpolation
			if (k>0.0)
			{
			double slope_pos=hist_vec[2]->GetBinContent(i);
			//double central=_hist[0]->GetBinContent(i);
			value += (slope_pos-central)*k; 
			}
			else if (k<0.0)
			{
			double slope_neg=hist_vec[1]->GetBinContent(i);
			//double central=_hist[0]->GetBinContent(i);
			value += (central-slope_neg)*k;
			}
			else if (k==0.0)
			{
			value += 0.0; //_hist[0]->GetBinContent(i);
			}	
	

		//Linear interpolation
		value += central;
		result00->SetBinContent(i,value);
	}
	return result00;
}


//Sign function
double sign(double j)
{
	if(j != 0) return (j/fabs(j));
	else return 0;
}

double max(double a, double b)//with ternary operator "?:"
{
	return (a<b)?b:a;     // or: return comp(a,b)?b:a; for the comp version
}

/*
//Function for interpolation with uncertainty according to rules of error propagation
TGraphAsymmErrors * tinterpolation::interpolate_histuncertainty(double* matrix, const int msize,std::vector<TH1D*> &hist_sig, std::vector<TH1D*> &hist_bg, const int hsize, std::vector<double> &k, const int ksize, double best_alpha_s, double best_alpha_b1, double best_alpha_b2, double best_alpha_b3, TH1D* mask1, TH1D* mask2, TH1D* mask3, int mode)
{
	TH1D * result00= new TH1D("result00", "result00", 15, -1.0, 1.0);

  	//Recalculation errormatrix
	double errormatrix[msize][msize];
	for (int m =0; m != msize; m++)
	{
	  	for (int n =0; n != msize; n++)
		{
		errormatrix[m][n] = matrix[m*msize+n]; 
		}		  
	}

  
	TGraphAsymmErrors * result = new TGraphAsymmErrors();
	result->Set(hist_sig[0]->GetNbinsX());
	
	//h0 as best fit hist_total
	TH1D * h0 = GetPlotHistogram(hist_sig, hist_bg, hist_sig, hist_bg, hist_sig, hist_bg, hsize, k, ksize, best_alpha_s, best_alpha_b1, best_alpha_b2, best_alpha_b3, mask1, mask2, mask3, mode=1);
	
	//Definition of hs and hb as histograms for signal and background resulting from the fit
	tinterpolation inter_obj_s (hist_sig, hsize);
	tinterpolation inter_obj_b (hist_bg, hsize);

 	//Interpolation separately for signal and background
	TH1D * hs = inter_obj_s.interpolate(result00, k, ksize, "hist_sig");
	TH1D * hb = inter_obj_b.interpolate(result00, k, ksize, "hist_bg");
	
	//loop over all bins
	for (int i =1; i != h0->GetNbinsX()+1; i++)
	{
		result->SetPoint(i,h0->GetBinCenter(i),h0->GetBinContent(i));
		result->SetPointEXlow(i,h0->GetBinWidth(i)/2);
		result->SetPointEXhigh(i,h0->GetBinWidth(i)/2);
		
		//Definition of derivatives
		double deriv [ksize];//alpha_b's included
		const int dsize = sizeof(deriv) / sizeof(deriv[0]);
		
		//Calculation of derivatives (see formula error propagation)
		//Calculation according to number of parameters
		
		//k[0] = alpha_wt and alpha_b's included in k array
		for(int j = 0; j != ksize; j++)
		{
			if (j == 100)
			//if (j == symm_inter_num[0] || j == symm_inter_num[1])
			{
			//Smoother interpolation:
			double hskup = hist_sig[2*j+2]->GetBinContent(i);
			double hbkup = hist_bg[2*j+2]->GetBinContent(i);
			double hskdown = hist_sig[2*j+1]->GetBinContent(i);
			double hbkdown = hist_bg[2*j+1]->GetBinContent(i);
			
			deriv[j] = 0.5*(hskup-hskdown) + 0.5*(hbkup-hbkdown);
			}
			else
			{
			if (k[j]>0.0)
			{
			double hsk = hist_sig[2*j+2]->GetBinContent(i);
			double hbk = hist_bg[2*j+2]->GetBinContent(i);
			double hs0 = hist_sig[0]->GetBinContent(i);
			double hb0 = hist_bg[0]->GetBinContent(i);
			
			//For linear interpolation:
			deriv[j] = (hsk-hs0) + (hbk-hb0);
			}
			else if (k[j]<0.0)
			{
			double hsk = hist_sig[2*j+1]->GetBinContent(i);
			double hbk = hist_bg[2*j+1]->GetBinContent(i);
			double hs0 = hist_sig[0]->GetBinContent(i);
			double hb0 = hist_bg[0]->GetBinContent(i);
			
			//For linear interpolation:
			deriv[j] = (hs0-hsk) + (hb0-hbk);
			}
			else if (k[j]==0.0)
			{
			deriv[j] = 0.0; //_hist[0]->GetBinContent(i);
			}
			}
		}
		
		double error = 0.0;
		
		//Error propagation, sum over indices as two foor loops
		for (int m =0; m != dsize; m++)
		{
		  	for (int n =0; n != dsize; n++)
			{
			error += deriv[m]*deriv[n]*errormatrix[m][n]; 
			}		  
		}
		//cout << error << endl;
		//Link bin and corresponding error
		result->SetPointEYlow(i,sqrt(error));
		result->SetPointEYhigh(i,sqrt(error));
		//cout << "error a " << sqrt(error) << endl;
		//cout << "error b " << result->GetErrorYhigh(i) << endl;
		
		//cout << i << " error " << sqrt(error) << endl;
		
	
	}
	
	return result;
}


//New function for interpolation with uncertainty according to rules of error propagation (only for background histogram)

TGraphAsymmErrors * tinterpolation::interpolate_bguncertainty(double* matrix, const int msize, std::vector<TH1D*> &hist_sig, std::vector<TH1D*> &hist_bg, const int hsize, std::vector<double> &k, const int ksize, double best_alpha_s, double best_alpha_b1, double best_alpha_b2, double best_alpha_b3, TH1D* mask1, TH1D* mask2, TH1D* mask3, int mode)
{
	TH1D * result00= new TH1D("result00", "result00", 15, -1.0, 1.0);

	//Recalculation errormatrix
	double errormatrix[msize][msize];
	for (int m =0; m != msize; m++)
	{
	  	for (int n =0; n != msize; n++)
		{
		errormatrix[m][n] = matrix[m*msize+n]; 
		}		  
	}
	
	//Definitions for background histogram
	//best_alpha_s = 0;
	
	TGraphAsymmErrors * result = new TGraphAsymmErrors();
	result->Set(hist_bg[0]->GetNbinsX());
	
	//h0 as best fit hist_total
	//TH1D * h0 = GetPlotHistogram(hist_sig, hist_bg, hsize, k, ksize, best_alpha_s, best_alpha_b, mode=1);
	
	//Definition of hs and hb as histograms for signal and background resulting from the fit
	tinterpolation inter_obj_s (hist_sig, hsize);
	tinterpolation inter_obj_b (hist_bg, hsize);

 	//Interpolation for background
	TH1D * hs = inter_obj_s.interpolate(result00, k, ksize, "hist_sig");
	TH1D * hb = inter_obj_b.interpolate(result00, k, ksize, "hist_bg");
	
	//loop over all bins
	for (int i =1; i != hb->GetNbinsX()+1; i++)
	{
		result->SetPoint(i,hb->GetBinCenter(i),hb->GetBinContent(i));
		result->SetPointEXlow(i,hb->GetBinWidth(i)/2);
		result->SetPointEXhigh(i,hb->GetBinWidth(i)/2);
		
		//Definition of derivatives
		double deriv [ksize];//without k_fsr and k_isr
		const int dsize = sizeof(deriv) / sizeof(deriv[0]);
		
		//Calculation of derivatives (see formula error propagation
		//Calculation according to number of parameters
		//H_S not used for background!!!
		
		//k[0] = alpha_wt and alpha_b's included in k array
		for(int j = 0; j != ksize; j++)
		{
			if (j == 100)
			//if (j == symm_inter_num[0] || j == symm_inter_num[1])
			{
			//Smoother interpolation:
			double hskup = hist_sig[2*j+2]->GetBinContent(i);
			double hbkup = hist_bg[2*j+2]->GetBinContent(i);
			double hskdown = hist_sig[2*j+1]->GetBinContent(i);
			double hbkdown = hist_bg[2*j+1]->GetBinContent(i);
			
			deriv[j] = 0.5*(hbkup-hbkdown);
			}
			else
			{
			if (k[j]>0.0)
			{
			double hbk = hist_bg[2*j+2]->GetBinContent(i);
			double hb0 = hist_bg[0]->GetBinContent(i);
			
			//For linear interpolation:
			deriv[j] = (hbk-hb0);
			}
			else if (k[j]<0.0)
			{
			double hbk = hist_bg[2*j+1]->GetBinContent(i);
			double hb0 = hist_bg[0]->GetBinContent(i);
			
			//For linear interpolation
			deriv[j] = (hb0-hbk);
			}
			else if (k[j]==0.0)
			{
			deriv[j] = 0.0; //_hist[0]->GetBinContent(i);
			}
			}
		}
		
		double error = 0.0;
		
		//Error propagation, sum over indices as two foor loops
		for (int m =0; m != dsize; m++)
		{
		  	for (int n =0; n != dsize; n++)
			{
			error += deriv[m]*deriv[n]*errormatrix[m][n]; 
			}		  
		}
		//cout << error << endl;
		//Link bin and corresponding error
		result->SetPointEYlow(i,sqrt(error));
		result->SetPointEYhigh(i,sqrt(error));
	}
	
	return result;
}
*/



//Get histogram which can directly be used for various plots
//TH1D * GetPlotHistogram(TH1D **hist_sig, TH1D **hist_bg, int hsize, double *k, int ksize, double best_alpha_s, double best_alpha_b1, double best_alpha_b2, double best_alpha_b3, TH1D* mask1, TH1D* mask2, TH1D* mask3,  int mode)
/*
TH1D * GetPlotHistogram(std::vector<TH1D*> &hist_F0_vec_sig,std::vector<TH1D*> &hist_F0_vec_bg, std::vector<TH1D*> &hist_FL_vec_sig,std::vector<TH1D*> &hist_FL_vec_bg, std::vector<TH1D*> &hist_FR_vec_sig,std::vector<TH1D*> &hist_FR_vec_bg, int hsize, std::vector<double> &k, int ksize, double best_alpha_s, double best_alpha_b1, double best_alpha_b2, double best_alpha_b3, TH1D* mask1, TH1D* mask2, TH1D* mask3, int mode)
{
	TH1D * result00= new TH1D("result00", "result00", 15, -1.0, 1.0);

	//Interpolation objects
	//tinterpolation inter_obj_s (hist_sig, hsize, mode=1);
	//tinterpolation inter_obj_b (hist_bg, hsize, mode=1);
	tinterpolation inter_obj_F0_s (hist_F0_vec_sig, hsize);
	tinterpolation inter_obj_F0_b (hist_F0_vec_bg, hsize);
	tinterpolation inter_obj_FL_s (hist_FL_vec_sig, hsize);
	tinterpolation inter_obj_FL_b (hist_FL_vec_bg, hsize);
	tinterpolation inter_obj_FR_s (hist_FR_vec_sig, hsize);
	tinterpolation inter_obj_FR_b (hist_FR_vec_bg, hsize);

	//Interpolation separately for signal and background
	//TH1D * hist_int_sig = inter_obj_s.interpolate(k, ksize, symm_inter_num, "hist_sig");
	//TH1D * hist_int_bg = inter_obj_b.interpolate(k, ksize, symm_inter_num, "hist_bg");
	TH1D * hist_F0_sig = inter_obj_F0_s.interpolate(result00, k, ksize, "hist_F0_sig");
	TH1D * hist_F0_bg = inter_obj_F0_b.interpolate(result00, k, ksize, "hist_F0_bg");
	TH1D * hist_FL_sig = inter_obj_FL_s.interpolate(result00, k, ksize, "hist_FL_sig");
	TH1D * hist_FL_bg = inter_obj_FL_b.interpolate(result00, k, ksize, "hist_FL_bg");
	TH1D * hist_FR_sig = inter_obj_FR_s.interpolate(result00, k, ksize, "hist_FR_sig");
	TH1D * hist_FR_bg = inter_obj_FR_b.interpolate(result00, k, ksize, "hist_FR_bg");
			
	//Combination of signal and background histograms
	TH1D * hist_F0_sum = (TH1D*) hist_F0_sig->Clone("F0_sig");
	hist_F0_sum->Scale(1.0);
	hist_F0_sum->Add(hist_F0_bg,1.0);
	TH1D * hist_FL_sum = (TH1D*) hist_FL_sig->Clone("FL_sig");
	hist_FL_sum->Scale(1.0);
	hist_FL_sum->Add(hist_FL_bg,1.0);
	TH1D * hist_FR_sum = (TH1D*) hist_FR_sig->Clone("FR_sig");
	hist_FR_sum->Scale(1.0);
	hist_FR_sum->Add(hist_FR_bg,1.0);
	
	//Combine F0, FL, FR:
	//hist_sum = hist_sum_F0 + hist_sum_FL + hist_sum_FR with command "add"
	TH1D * hist_sum = (TH1D*) hist_F0_sum->Clone("F0_sum");
	hist_sum->Add(hist_FL_sum,1.0);
	hist_sum->Add(hist_FR_sum,1.0);
	
	
	return hist_sum;

}
*/

TH1D TemplateFitErrors ( double* matrix, std::vector<TH1D> Hist, std::vector<double> EffVec, int NSignal)
{

  int nParam  = Hist.size();
  int nSignal = NSignal;
  int nBkg    = nParam-nSignal;
  
  int MatrixSize = nParam;
  
  //Recalculation errormatrix
  double errormatrix[MatrixSize][MatrixSize];
  for (int m =0; m != MatrixSize; m++)
    {
      for (int n =0; n != MatrixSize; n++)
	{
	  errormatrix[m][n] = matrix[m*MatrixSize+n];
	}
    }
  
  std::vector<TH1D> NormHist;

  for(int iHist = 0; iHist < nParam; ++iHist){

    Hist[iHist].Scale(1.0/Hist[iHist].Integral());
    
    NormHist.push_back(Hist[iHist]);

  }

  int nBins = Hist[0].GetNbinsX();
  
  TH1D UncHist = *(TH1D*) Hist[0].Clone("");

  //Loop over all bins
  for (int i =1; i != nBins+1; i++){


    std::vector<double> BinEntries;
    BinEntries.clear();
 
    for (int k = 0; k != nSignal; k++) {

      BinEntries.push_back(NormHist[k].GetBinContent(i)*EffVec[k]);

    }
    for (int k = 0; k != nBkg; k++) {

      BinEntries.push_back(NormHist[k + nSignal].GetBinContent(i));

    }


    double error = 0.0;

    //Error propagation, sum over indices as two foor loops
    for (int m = 0; m != nParam; m++){
      
      for (int n = 0; n != nParam; n++){

	error += BinEntries[m]*BinEntries[n]*errormatrix[m][n];

	//	std::cout << m << "\t" << n << "\t" << BinEntries[m]*BinEntries[n]*errormatrix[m][n] << std::endl;

      }

    }
  
    UncHist.SetBinContent(i, sqrt(error));
    
  }

  return UncHist;

}



TGraphAsymmErrors * TemplateFitErrors ( double* matrix, const int msize, TH1D hist_sum0, std::vector<TH1D> &hist_norm_vec0, std::vector<double> &eff_vec, TGraphAsymmErrors* result, int Nsignal, int Nbkg) 
{
	TH1D * hist_sum = &hist_sum0;
	
	std::vector<TH1D*> hist_norm_vec;

	for (int i = 0; i != hist_norm_vec0.size(); i++) {
		
		hist_norm_vec.push_back( &hist_norm_vec0[i] );
	}
	

  	//Recalculation errormatrix
	double errormatrix[msize][msize];
	for (int m =0; m != msize; m++)
	{
	  	for (int n =0; n != msize; n++)
		{
		errormatrix[m][n] = matrix[m*msize+n]; 
		}		  
	}

	//Define TGraphAsymmErrors  
	result->Set(hist_sum->GetNbinsX());

	//Loop over all bins
	for (int i =1; i != hist_sum->GetNbinsX()+1; i++)
	{
		result->SetPoint(i,hist_sum->GetBinCenter(i),hist_sum->GetBinContent(i));
		result->SetPointEXlow(i,hist_sum->GetBinWidth(i)/2);
		result->SetPointEXhigh(i,hist_sum->GetBinWidth(i)/2);
		
		//Definition of derivatives
		double deriv[6];//3 signal and background templates
		const int dsize = sizeof(deriv) / sizeof(deriv[0]);
		
		//Calculation of derivatives (see formula error propagation)
		//Calculation according to number of parameters

		for (int k = 0; k != Nsignal; k++) {

			deriv[k] = hist_norm_vec[k]->GetBinContent(i)*eff_vec[k];
		
		}
		for (int k = 0; k != Nbkg; k++) {

			deriv[k+Nsignal] = hist_norm_vec[k + Nsignal]->GetBinContent(i);
		
		}

		/*deriv[0] = hist_norm_vec[0]->GetBinContent(i)*eff_vec[0];
		deriv[1] = hist_norm_vec[1]->GetBinContent(i)*eff_vec[1];
		deriv[2] = hist_norm_vec[2]->GetBinContent(i)*eff_vec[2];
		deriv[3] = hist_norm_vec[3]->GetBinContent(i);
		deriv[4] = hist_norm_vec[4]->GetBinContent(i);
		deriv[5] = hist_norm_vec[5]->GetBinContent(i);*/
		
		//interpolation part NOT included !!!
		
		/*cout << deriv[0] << endl;
		cout << deriv[1] << endl;
		cout << deriv[2] << endl;
		cout << deriv[3] << endl;
		cout << deriv[4] << endl;
		cout << deriv[5] << endl;*/
		
		double error = 0.0;
		
		//Error propagation, sum over indices as two foor loops
		for (int m =0; m != hist_norm_vec0.size(); m++)
		{
		  	for (int n =0; n != dsize; n++)
			{
			error += deriv[m]*deriv[n]*errormatrix[m][n]; 
			}		  
		}

		//cout << "error " << error << endl;

		//Link bin and corresponding error
		result->SetPointEYlow(i,sqrt(error));
		result->SetPointEYhigh(i,sqrt(error));
	}

	return result;

}



