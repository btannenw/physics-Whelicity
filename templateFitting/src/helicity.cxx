/*
 helicity.cxx
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

#include "helicity.h"

 using namespace std;  

///____________________CALCULATE HELICITY FRACTIONS___________________________

double eval_fi (TH1D hist_nom, double alpha, double alpha_sum )
{
	//double hist_int = hist_nom->GetSumOfWeights();
  //	double hist_int = hist_nom.Integral();

	//cout << "hist_int " << hist_int << endl;

	double Fi = alpha/alpha_sum;
	
	return Fi;
}


double fi_error ( double* matrix, const int msize, double alpha, double alpha_sum, int index)
{
	//Recalculation errormatrix
	double errormatrix[msize][msize];
	for (int m =0; m != msize; m++)
	{
	  	for (int n =0; n != msize; n++)
		{
		errormatrix[m][n] = matrix[m*msize+n]; 
		}		  
	}


	/* 	for (int m =0; m != 3; m++)
	  {
	    for (int n =0; n != 3; n++)
	      {
		
 		cout << m << "\t" << n << "\t" << errormatrix[m][n] <<  endl;
	      }
	      } */

	double deriv[3]; //three helicity states

	//Calculate derivatives depending on alpha_fi
	if (index == 0) //F0
	{
	deriv[0] = (alpha_sum - alpha)/(alpha_sum*alpha_sum);
	deriv[1] = -alpha/(alpha_sum*alpha_sum);
	deriv[2] = -alpha/(alpha_sum*alpha_sum);
	}
	else if (index == 1) //FL
	{
	deriv[0] = -alpha/(alpha_sum*alpha_sum);
	deriv[1] = (alpha_sum - alpha)/(alpha_sum*alpha_sum);
	deriv[2] = -alpha/(alpha_sum*alpha_sum);
	}
	else if (index == 2) //FR
	{
	deriv[0] = -alpha/(alpha_sum*alpha_sum);
	deriv[1] = -alpha/(alpha_sum*alpha_sum);
	deriv[2] = (alpha_sum - alpha)/(alpha_sum*alpha_sum);
	}
	else 
	{
	cout << index << " is a wrong index - not possible to calculate derivatives !" << endl;
	return 0;
	}

	double error = 0.0;

	for (int m =0; m != 3; m++)
	{
  		for (int n =0; n != 3; n++)
		{
		/*if (m != n) error += deriv[m]*deriv[n]*errormatrix[m][n]; 
		else error += deriv[m]*deriv[m]*k_err[m]*k_err[m];*/
		error += deriv[m]*deriv[n]*errormatrix[m][n];
		}		  
	}
	double sigma = sqrt(error); 

	//cout << "Helicity error squared : " << error << endl;

	return sigma;
}




