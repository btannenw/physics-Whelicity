/*
 validation.cxx
*/

#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <ios>

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
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TRandom.h"
#include <TRandom3.h>

#include "ValidationClass.h"
#include "helicity.h"
#include "plots.h"
#include "interpolation.h"
#include "ProfilingClass.h"
#include "StatusLogbook.h"

#include "fit.h"
#include "fit_validation.h"
#include "functions.h"

using namespace std;

//TH1D hist_ensemble;
//TH1D hist_data;
std::vector<tinterpolation> inter_obj;
std::vector<double> fBkgNorm;
std::vector<double> fBkgUnc;
std::vector<TH1D> fParameterHist;
std::vector<double> SelEff;

std::string fBkgMode;

//Constructor 
tvalidation::tvalidation()
{
  
  fHistNomVec.clear();
  fHistNormVec.clear();
  fEffVec.clear();
  fBkgExpVec.clear();
  fBkgUncVec.clear();

  fNNui          = 0;
  fNScaleFactors = 0;
  
  fRandom = new TRandom3(0);
  gaussRandom = new TRandom3(0);

  fit   = new TF1();

  graph = new TGraph(50);

  fitter_prof = new TMinuit(50);
  fitter_prof->SetPrintLevel(-1);

  //Set interpolation object for fit (at first only linear)                                                                                                    
  fitter_prof->SetFCN(eval_chisqrt_vali);

  out_nui_nom.clear();

  gErrorIgnoreLevel = kError;


}

void tvalidation::clean()
{

   out_F0             = 0.0; out_FL           = 0.0;  out_FR             = 0.0;
   out_N0             = 0.0; out_NL           = 0.0;  out_NR             = 0.0;
   out_Nges           = 0.0; out_Nges_nom     = 0.0;
   out_Wjets          = 0.0; out_QCD          = 0.0; out_RemBkg          = 0.0;
   out_F0_nom         = 0.0; out_FL_nom       = 0.0; out_FR_nom          = 0.0;
   out_N0_nom         = 0.0; out_NL_nom       = 0.0; out_NR_nom          = 0.0;
   out_Wjets_nom      = 0.0; out_QCD_nom      = 0.0; out_RemBkg_nom      = 0.0;
   out_F0_err         = 0.0; out_FL_err       = 0.0; out_FR_err          = 0.0;
   out_N0_err         = 0.0; out_NL_err       = 0.0; out_NR_err          = 0.0;
   out_Wjets_err      = 0.0; out_QCD_err      = 0.0; out_RemBkg_err      = 0.0;
   out_N0_prof_err    = 0.0; out_NL_prof_err  = 0.0; out_NR_prof_err     = 0.0;
   out_Wjets_prof_err = 0.0; out_QCD_prof_err = 0.0; out_RemBkg_prof_err = 0.0;

   out_nui.clear();
   // out_nui_nom.clear();
   out_nui_err.clear();
   out_nui_prof_err.clear();

   out_MinuitStatus    = -1;
   out_NumNuisancePar  =  0;



}


void tvalidation::FillInterpolationVector()
{

  for(int iParam = 0; iParam < fNSignal; ++iParam){

    tinterpolation inter_obj_help("Signal", iParam);

    if(iParam < fNSignal){

      fParameterHist.push_back(fSignalInterpolation[0][iParam].hist[0]);

      fNBins = fSignalInterpolation[0][iParam].hist[0].GetNbinsX();

      SelEff.push_back(fSignalInterpolation[0][iParam].SelEff);

    }

    inter_obj.push_back(inter_obj_help);

  }

  for(int iParam = 0; iParam < fNBkg; ++iParam){

    tinterpolation inter_obj_help("Bkg", iParam);
    fParameterHist.push_back(fBkgInterpolation[0][iParam].hist[0]);

    fBkgNorm.push_back(fBkgInterpolation[0][iParam].BkgExp);
    fBkgUnc.push_back(fBkgInterpolation[0][iParam].BkgUnc);

    inter_obj.push_back(inter_obj_help);

  }

}

void tvalidation::FillTemplateVector(std::vector<TemplateInfo::MySignalTemplate> SigTempl, std::vector<TemplateInfo::MyBkgTemplate> BkgTempl)
{

  for(int iParam = 0; iParam < fNSignal; ++iParam){

    if(iParam < fNSignal){
      fParameterHist.push_back(SigTempl[iParam].hist);
      SelEff.push_back(SigTempl[iParam].SelEff);

    }

  }
  for(int iParam = 0; iParam < fNBkg; ++iParam){

    fParameterHist.push_back(BkgTempl[iParam].hist);
    fBkgNorm.push_back(BkgTempl[iParam].Exp);
    fBkgUnc.push_back(BkgTempl[iParam].Unc);

  }


  /*  for(int iParam = 0; iParam < fNSignal+fNBkg; ++iParam){

    cout << iParam << "\t" << fParameterHist[iParam].GetNbinsY() << "\t" << fParameterHist[iParam].GetYaxis()->GetXmin() << "\t" << fParameterHist[iParam].GetYaxis()->GetXmax()<< endl;

    cout << iParam << "\t" << fParameterHist[iParam].GetNbinsX() << "\t" << fParameterHist[iParam].GetXaxis()->GetXmin() << "\t" << fParameterHist[iParam].GetXaxis()->GetXmax()<< endl;
    
    cout << iParam << "\t" << fParameterHist[iParam].GetNbinsZ() << "\t" << fParameterHist[iParam].GetZaxis()->GetXmin() << "\t" << fParameterHist[iParam].GetZaxis()->GetXmax()<< endl;



    }*/

}

void tvalidation::SetFitParameters(std::string BkgMode)
{

  fBkgMode = BkgMode;
  
  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){

    fHistNomVec.push_back(fSignalParameter[iSig].hist);
    fHistNormVec.push_back(normalise(fSignalParameter[iSig].hist));
    fEffVec.push_back(fSignalParameter[iSig].SelEff);

  }
  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){

    fHistNomVec.push_back(fBkgParameter[iBkg].hist);
    fHistNormVec.push_back(normalise(fBkgParameter[iBkg].hist));
    fBkgExpVec.push_back(fBkgParameter[iBkg].Exp);
    fBkgUncVec.push_back(fBkgParameter[iBkg].Unc);

  }

  if(BkgMode == "BkgFit"){
    fNScaleFactors = fHistNomVec.size();

    fNNui          = fNuisanceParameter.size();
  }
  else if (BkgMode == "BkgNui"){

    fNScaleFactors = fSignalParameter.size();
    fNNui          = fNuisanceParameter.size();

  }
    

}


//Get a histogram containing pseudo-data based on certain input values for F_i
TH1D tvalidation::GetPseudodataFit(double F0, double FL, double FR, std::vector<double> vec_nui)
{		
  std::vector<double> Fractions;
  Fractions.push_back(F0);
  Fractions.push_back(FL);
  Fractions.push_back(FR);

  double lower_edge  = fHistNormVec[0].GetBinLowEdge(1);
  double bin_width   = fHistNormVec[0].GetBinWidth(1);
  double number_bins = fHistNormVec[0].GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;

  //  std::cout << "lower = " << lower_edge << "   upper = " << upper_edge << std::endl;

  //Create histogram with pseudodata
  TH1D hist_pseudo0   = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
  TH1D hist_help      = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
  TH1D hist_test      = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);


  for(int iParam = 0; iParam < fNSignal+fNBkg; ++iParam){
    
    for(int iBin = 1; iBin <= number_bins; ++iBin){
            
      double bin_sum = 0.0;
      
      if(iParam < fNSignal)
	//bin_sum = fFitFunc[0][iParam][iBin-1].Eval(0)*Fractions[iParam];
	bin_sum = fFitFunc_QuadFit[0][iParam][iBin-1].Eval(0)*Fractions[iParam];
      else
	//bin_sum = fFitFunc[0][iParam][iBin-1].Eval(0);
	bin_sum = fFitFunc_QuadFit[0][iParam][iBin-1].Eval(0);
      
      //double central = fFitFunc[0][iParam][iBin-1].Eval(0);
      double central = fFitFunc_QuadFit[0][iParam][iBin-1].Eval(0);
      

      for(int iNui = 0; iNui < fNNui; ++iNui){
		
	if(iParam < fNSignal){
	  
	  //bin_sum += (fFitFunc[iNui][iParam][iBin-1].Eval(vec_nui[iNui]) - central)*Fractions[iParam];
	  bin_sum += (fFitFunc_QuadFit[iNui][iParam][iBin-1].Eval(vec_nui[iNui]) - central)*Fractions[iParam];	

	}
	else	  
	  //bin_sum += fFitFunc[iNui][iParam][iBin-1].Eval(vec_nui[iNui]) - central;
	  bin_sum += fFitFunc_QuadFit[iNui][iParam][iBin-1].Eval(vec_nui[iNui]) - central;

      }
      
      //hist_test.SetBinContent(iBin, fFitFunc[0][iParam][iBin-1].Eval(0));
      hist_test.SetBinContent(iBin, fFitFunc_QuadFit[0][iParam][iBin-1].Eval(0));

      hist_help.SetBinContent(iBin, bin_sum);
      

    }
    
    //    std::cout << iParam << "\t" << hist_help.Integral() << std::endl;

    fNominalScaleParameters.push_back(hist_test.Integral());
    
    hist_pseudo0 = hist_pseudo0 + hist_help;
    
    if(iParam >= fNSignal)
      fBkgNorm[iParam-fNSignal] = hist_help.Integral();
    
  }
  
  std::stringstream oss;
  oss << hist_pseudo0.Integral();
  
  WriteParameterStatus("ValidationClass", "Integral pseudo data (fit): "+oss.str());
  
  return hist_pseudo0;

}

//Get a histogram containing pseudo-data based on certain input values for F_i
TH1D tvalidation::GetPseudodata(double F0, double FL, double FR)
{		
  std::vector<double> Fractions;
  Fractions.push_back(F0);
  Fractions.push_back(FL);
  Fractions.push_back(FR);

  double lower_edge  = fHistNormVec[0].GetBinLowEdge(1);
  double bin_width   = fHistNormVec[0].GetBinWidth(1);
  double number_bins = fHistNormVec[0].GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;

  std::cout << "lower = " << lower_edge << "   upper = " << upper_edge << std::endl;


  //Create histogram with pseudodata
  TH1D hist_pseudo0   = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

  TH1D hist_help      = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){
    
    hist_help = fHistNormVec[iSig];
    hist_help.Scale(fXsec*fLumi*Fractions[iSig]*fEffVec[iSig]);

    //    std::cout << iSig << "\t" <<  hist_help.Integral() << std::endl;
  
    hist_pseudo0 = hist_pseudo0+hist_help;

    fNominalScaleParameters.push_back(hist_help.Integral());

  }
  
  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
    
    hist_help = fHistNormVec[iBkg+fSignalParameter.size()];
    hist_help.Scale(fHistNomVec[iBkg+fSignalParameter.size()].Integral());

    fNominalScaleParameters.push_back(hist_help.Integral());

    //    std::cout << iBkg << "\t" <<  hist_help.Integral() << std::endl;

    fBkgNorm[iBkg] = hist_help.Integral();
      
    hist_pseudo0 = hist_pseudo0+hist_help;

  }

  std::stringstream oss;
  oss << hist_pseudo0.Integral();

  WriteParameterStatus("ValidationClass", "Integral pseudo data (templates): "+oss.str());
  
  return hist_pseudo0;

}


//Get variation of pseudo-data
TH1D tvalidation::RandomPseudodata(TH1D Histo)
{

  double lower_edge  = Histo.GetBinLowEdge(1);
  double bin_width   = Histo.GetBinWidth(1);
  double number_bins = Histo.GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;
  
  hist_ensemble = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

  //Loop over bins and fill them
  for(int ibin = 1; ibin <= number_bins; ibin++)
    {
      double nexp = Histo.GetBinContent(ibin);
      
      double nobs = fRandom->PoissonD(nexp);
      
      //Set the bin content
      hist_ensemble.SetBinContent(ibin, nobs);
    }

  std::stringstream oss;
  oss << hist_ensemble.Integral();

  //  WriteParameterStatus("ValidationClass", "Integral ensemble: "+oss.str());
  
  //Return the ensemble histogram
  return hist_ensemble;

}

///_______________________________________________________________________
///_____________________________VALIDATION________________________________
///_______________________________________________________________________


void tvalidation::CallValidation(std::string validation_mode, int nui_var, std::string input_channel)
{

  if(validation_mode == "CaliCurvesFraction"){

    WriteInfoStatus("Validation", "Run validation chain with option: "+validation_mode);

    std::vector<double> Nui_Vec;

    for(int iNui = 0; iNui < fNNui; ++iNui)
      Nui_Vec.push_back(0.0);

    double diff = 0.1;

    ///___________For-loop for different helicity fractions___________
    //for(int i = 0; i < 7; ++i){
    for(int i = 3; i < 4; ++i){ //for tests with F_0 = 0.7, etc...

      // Calculate fractions F0, FL, FR
      // SM configuration is now realized for counter_fractions = 3       
      double counter_fractions = i;
      double diff_counter = diff*counter_fractions;
  
      // makes calibration curves for F0 between 0.4 and 1.0, FL between 0.45 and 0.15 and FR between +- 0.15
      double F0_true = 0.4  + diff_counter;
      double FL_true = 0.45 - diff_counter*0.5;
      double FR_true = 0.15 - diff_counter*0.5;

      std::string outputfile_mode = "Validation";

      Validation(i, 0, Nui_Vec, diff, F0_true, FL_true, FR_true, input_channel, outputfile_mode);
      
    }

  }
  else if(validation_mode == "CaliCurvesNuisanceParam"){

    WriteInfoStatus("Validation", "Run validation chain with option: "+validation_mode);

    std::vector<double> Nui_Vec;

    for(int i = 0; i < 7; ++i){
    //int i = 0;

      out_nui_nom.clear();
      Nui_Vec.clear();

      for(int iNui = 0; iNui < fNNui; ++iNui)
	Nui_Vec.push_back(0.0);

      Nui_Vec[nui_var-1] = -1.5 + (double)i*0.5; 

      std::cout << "Vary nui " << nui_var-1 << "\t" << i << "\t" <<  Nui_Vec[nui_var-1] << std::endl;

      for(int iNui = 0; iNui < fNNui; ++iNui)
	out_nui_nom.push_back(Nui_Vec[iNui]);

      double F0_true = 0.7;
      double FL_true = 0.3;
      double FR_true = 0.0;

      std::string outputfile_mode = "Validation";

      // for counter fractions == 3 we get the SM configuration
      Validation(3, nui_var, Nui_Vec, 0.1, F0_true, FL_true, FR_true, input_channel, outputfile_mode);

    }
  }
  else if(validation_mode == "FluctuateNuisanceParam"){
    
    std::vector<double> Nui_Vec;

    //for(int i = 0; i < 7; ++i){
    int i = 0;

      out_nui_nom.clear();
      Nui_Vec.clear();

      for(int iNui = 0; iNui < fNNui; ++iNui)
	Nui_Vec.push_back(0.0);

      //Nui_Vec[nui_var-1] = -1.5 + (double)i*0.5; 
      Nui_Vec[nui_var-1] = 0.0; //gaussRandom->Gaus(0.0, 1.0);

      std::cout << "Vary nui " << nui_var-1 << "\t" << i << "\t" <<  Nui_Vec[nui_var-1] << std::endl;

      for(int iNui = 0; iNui < fNNui; ++iNui)
	out_nui_nom.push_back(Nui_Vec[iNui]);

      double F0_true = 0.7;
      double FL_true = 0.3;
      double FR_true = 0.0;

      std::string outputfile_mode = "ValidationFluctuate";

      // for counter fractions == 3 we get the SM configuration
      Validation(3, nui_var, Nui_Vec, 0.1, F0_true, FL_true, FR_true, input_channel, outputfile_mode);


  }
  else if(validation_mode == "Datafit"){

    WriteInfoStatus("Validation", "Run validation chain with option: "+validation_mode);

    std::vector<double> Nui_Vec;

    for(int iNui = 0; iNui < fNNui; ++iNui)
      Nui_Vec.push_back(0.0);

    Datafit(Nui_Vec, input_channel, "Full");

  }
  else if(validation_mode == "Datafit_Single"){

    WriteInfoStatus("Validation", "Run validation chain with option: "+validation_mode);

    std::vector<double> Nui_Vec;

    for(int iNui = 0; iNui < fNNui; ++iNui)
      Nui_Vec.push_back(0.0);

    Datafit(Nui_Vec, input_channel, "Single");

  }
  else{

    WriteErrorStatus("Validation", "No such option for the validation mode!");
    WriteErrorStatus("Validation", "Note: options are CaliCurvesFraction or CaliCurvesNuisanceParam");
    WriteErrorStatus("Validation", "EXIT");

    exit(-1);

  }

}

//Parameter distributions based on pseudo-data: nui_var is the index of the parameter that is varied, Nui_Vec contains all values for the nui param that are used to generate the pseudo data
void tvalidation::Validation(int counter_fractions, int nui_var, std::vector<double> Nui_Vec, double diff, double F0_true, double FL_true, double FR_true, std::string input_channel, std::string outputfile_mode)
{

      double nui_value = 0.0;

      for(int bla = 0; bla < Nui_Vec.size(); ++bla)
	std::cout << bla << "\t" << Nui_Vec[bla] << std::endl;

      if(nui_var <= Nui_Vec.size() && Nui_Vec.size() != 0)
	nui_value = Nui_Vec[nui_var-1];

      std::cout << nui_var << "\t" << nui_value << std::endl;

      WriteInfoStatus("ValidationClass", "Get pseudo data histogram");
      
      cout << "Interpolation method   " << fInterpolationMethod.c_str() << endl;

      TH1D hist_pseudo;

      //Generate peudo data according to input files
      if( Nui_Vec.size() > 0)
	if (outputfile_mode == "Validation") {
	hist_pseudo = GetPseudodataFit(F0_true, FL_true, FR_true, Nui_Vec);
	}
      else
	hist_pseudo = GetPseudodata(F0_true, FL_true, FR_true);

      WriteInfoStatus("ValidationClass", "Obtained pseudo data histogram");

      int k_counter    = 0;
      int conv_counter = 0;
      
      std::stringstream F0_label, nr_file, nui_num,output_nr,nui_var_str,nui_value_str;
      F0_label      << F0_true;
      nui_num       << fNNui;
      output_nr     << fNOutputSample;
      nui_var_str   << nui_var;
      if (fabs(nui_value) < 0.001) nui_value_str << 0;
      else nui_value_str << fabs(nui_value);
      
      std::string output_file;

      if(outputfile_mode == "Validation"){

      	if(nui_value < 0)
	  output_file = ("Output_"+input_channel+"_F0_"+F0_label.str()+"_NrNui"+nui_num.str()+"_NuiVar"+nui_var_str.str()+"_m"+nui_value_str.str()+"_"+output_nr.str()+".root");
      	else
	  output_file = ("Output_"+input_channel+"_F0_"+F0_label.str()+"_NrNui"+nui_num.str()+"_NuiVar"+nui_var_str.str()+"_p"+nui_value_str.str()+"_"+output_nr.str()+".root");

      }
      if(outputfile_mode == "ValidationFluctuate"){

	  output_file = ("OutputFluctuate_"+input_channel+"_F0_"+F0_label.str()+"_NrNui"+nui_num.str()+"_NuiVar"+nui_var_str.str()+"_"+output_nr.str()+".root");

      }

      else if (outputfile_mode == "Datafit") {

      	if(nui_value < 0)
	  output_file = ("DatafitPseudo_Output_"+input_channel+"_NrNui"+nui_num.str()+"_NuiVar"+nui_var_str.str()+"_m"+nui_value_str.str()+"_"+output_nr.str()+".root");
      	else
	  output_file = ("DatafitPseudo_Output_"+input_channel+"_NrNui"+nui_num.str()+"_NuiVar"+nui_var_str.str()+"_p"+nui_value_str.str()+"_"+output_nr.str()+".root");
      }

      // create file for new ensemble test
      
      WriteInfoStatus("ValidationClass", "Create new file: "+output_file);
      
      TFile *fNewFile = new TFile(output_file.c_str(), "RECREATE");
      
      WriteInfoStatus("ValidationClass", "Create new output tree");
      
      TTree *fNewTree = InitializeOutputTree();
      
      ///_____________________________FIT_______________________________
      
      WriteInfoStatus("ValidationClass", "Initialize TMinuit");

      //Setting up the fitting procedure; initialize TMinuit
      TMinuit *fitter = new TMinuit(50);
      
      //      hist_ensemble = 0;	

      int counter_good = 0;
      
      while(counter_good < fNPseudoExp)
	{

	  // clear the output tree
	  clean();

	  std::stringstream oss;
	  oss << counter_good;

	  WriteInfoStatus("ValidationClass", "Get new ensemble for pseudo experiments "+oss.str());

	  if (outputfile_mode == "ValidationFluctuate") {
	    //Vary/fluctuate nuisance parameter in each pseudo-experiment
            Nui_Vec[nui_var-1] = gaussRandom->Gaus(0.0, 1.0);
	    hist_pseudo = GetPseudodataFit(F0_true, FL_true, FR_true, Nui_Vec);
            //hist_ensemble = hist_pseudo;
	    hist_ensemble = RandomPseudodata(hist_pseudo);
	    cout << "Fluctuate nuisance parameter : " << Nui_Vec[nui_var-1] << endl;
	  }
	  else {
	  hist_ensemble = RandomPseudodata(hist_pseudo);
          //cout << "alles gut!!! " << endl;
	  }

	  //	  hist_ensemble = hist_pseudo;

	  //Setting up the fitting procedure; initialize TMinuit
	  if(fitter) fitter->Delete();
	  fitter = new TMinuit(50);
      	  fitter->SetPrintLevel(-1);

	  int ierflg_err=0, ierflg = 0, ierflg_hesse = 0, par_num;
	  double start_val, start_step=0.001, A=0.0, B=0.0;  
	  
	  double arglist[50];

	  //Set interpolation object for fit (at first only linear)
	  fitter->SetFCN(eval_chisqrt_vali);	
	  
	  //Generate new start values
	  std::vector <double> start_val_gauss;
	  
	  double Ntotal = fXsec*fLumi;
	  
	  // has to be made more general!!!
	  //	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*F0_true, 20000.0 ) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*F0_true, 20000.0) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FL_true, 20000.0) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FR_true, 20000.0) );
	  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg)
	    start_val_gauss.push_back( gaussRandom->Gaus(fBkgExpVec[iBkg], fBkgUncVec[iBkg]) );
	  
	  //Start values
	  double start_val_array[] = {Ntotal*F0_true};
	  // double start_val_array[] = {Ntotal*F0_true, Ntotal*FL_true, Ntotal*FR_true, fBkgExpVec[0], fBkgExpVec[1], fBkgExpVec[2]};
	  
	  //Define 6 scale parameters alpha for signal and background
	  for (int k = 0; k != fNScaleFactors; k++)
	    {
	      
	      std::string s;
	      std::stringstream out;
	      out << k;
	      s = out.str();

	      fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, A, B, ierflg);
	      
	      /*	      if(k < 2)
		fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, 50000, 550000, ierflg);
	      else if(k == 2)
		fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, -50000, 50000, ierflg);
	      else
	      fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, 0.0, 7500.0, ierflg); */


	    }
	  
	  //Define nuisance parameters k as free pararameters
	  for (int k = 0; k != fNNui; k++){

	    //probably nui-param for all three states + three backgrounds 	
	    //-> 6 nui-alphas + all nuisance parameters k 
	    //if(i < 6) fitter->mnparm(par_num=i, "alpha_" + used_nui_parameters[i], start_val=0.0, start_step, A, B, ierflg); 
	    
	    std::string s;
	    std::stringstream out;
	    out << k;
	    s = out.str();
	    
	    start_val = gaussRandom->Gaus(0.0, 1.0);
	    //start_val = 0.0;

	    //	    fitter->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, A, B, ierflg);
	    // debug A.Knue
	    fitter->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, -5.0, 5.0, ierflg);
	    
	    
	  }
	  
	  const int par_size = fNScaleFactors+fNNui;//6 scale parameters
	  
	  arglist[0] = 1.0;
	  fitter->mnexcm("SET ERR", arglist, 1, ierflg_err);

	  //Starting the fit
	  arglist[0] = 10000.0;
	  arglist[1] = 0.01;
	  fitter->mnexcm("MIGRAD", arglist, 2, ierflg);

	  arglist[0] = 1000.0;
	  arglist[1] = 0;
	  //not needed for tests in goe, but consider errors larger than 1 then
	  fitter->mnexcm("HESSE", arglist, 2, ierflg_hesse);//use HESSE/MINOS alternatively
	  
	  //Print Results; access to fit parameter
	  double fval,edm,errdef;
	  //fval: the best function value found so far; edm: the estimated vertical distance remaining to minimum; errdef: the value of UP defining parametEr uncertainties
	  int nvpar,nparx,icstat;//nvpar: the number of currently variable parameters, nparx: the highest (external) parameter number defined by user, icstat: a status integer indicating how good is the covariance
	  
	  fitter->mnstat(fval,edm,errdef,nvpar,nparx,icstat);
	  //fitter->mnprin(3,fval);
	  //3 refers to: values, errors, step sizes, first derivs.

	  //Get fit parameters
	  std::vector<double> best_par; 
	  std::vector<double> par_err;
	  double par_value, par_err_value;
	  best_par.clear();
	  par_err.clear();
	  	  
	  for (int k = 0; k != par_size; k++)
	    { 
	      fitter->GetParameter(k,par_value,par_err_value);
	      best_par.push_back(par_value);
	      par_err.push_back(par_err_value);
	      //      std::cout << par_err_value << std::endl;

	      std::stringstream oss_k,oss_par,oss_err;
	      oss_k << k;
	      oss_par << best_par[k];
	      oss_err << par_err[k];

	      WriteInfoStatus("ValidationClass", "Best parameter "+oss_k.str()+", par: "+oss_par.str()+", err: "+oss_err.str());

	    }
	  
	  //Get error matrix / covariance matrix
	  const int emsize = par_size;//k values and scale factors
	  double errormatrix[emsize][emsize];
	  fitter->mnemat(&errormatrix[0][0],emsize);

	  out_MinuitStatus   = ierflg;
	  out_NumNuisancePar = fNNui;

	  //add new error calculation: profiling plot
	  std::vector<double> error_profile_vec;
	  for (int k = 0; k != par_size; k++)
	    error_profile_vec.push_back( error_profile_PE(hist_pseudo, best_par[k], k, best_par, par_err, counter_good, 1, F0_true, FL_true, FR_true, 0) ) ;
	  
	    //	    error_profile_vec.push_back(1.0);
	    //	    error_profile_vec.push_back( error_profile_PE(hist_pseudo, best_par[k], k, best_par, par_err, i, counter, F0_true, FL_true, FR_true) ) ;
	      
	  //Sum of all fitted signal parameters alpha_(F0,FL,FR)
	  double alpha_sum = 0.0;
	  
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig)
	    alpha_sum += best_par[iSig];
	  
	  //Calculate helicity fractions ...
	  
	  if(fNSignal > 2){
	    out_F0 = eval_fi (fHistNomVec[0], best_par[0], alpha_sum);
	    out_FL = eval_fi (fHistNomVec[1], best_par[1], alpha_sum);
	    out_FR = eval_fi (fHistNomVec[2], best_par[2], alpha_sum);
	    
	    //... and related errors:
	    out_F0_err = fi_error (&errormatrix[0][0], emsize, best_par[0], alpha_sum, 0);
	    out_FL_err = fi_error (&errormatrix[0][0], emsize, best_par[1], alpha_sum, 1);
	    out_FR_err = fi_error (&errormatrix[0][0], emsize, best_par[2], alpha_sum, 2);
	    //last int represents index according to: 0 corresponds to F0, 1 to FL, 2 to FR!
	    
	  }
	  
	  // store nominal (true) values;
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){
	    
	    std::string signal = fSignalParameter[iSig].ParamName;
	    
	    if(signal == "N0")
	      out_F0_nom     = F0_true;
	    else if(signal == "NL")
	      out_FL_nom     = FL_true;
	    else if(signal == "NR")
	      out_FR_nom     = FR_true;
	  }
	  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
	    
	    std::string bkg = fBkgParameter[iBkg].ParamName;
	    
	    if(bkg == "Wjets")
	      out_Wjets_nom  = fBkgExpVec[0];
	    else if(bkg == "QCD")
	      out_QCD_nom    = fBkgExpVec[1];
	    else if(bkg == "RemBkg")
	      out_RemBkg_nom = fBkgExpVec[2];
	    
	  }
	  

	  out_N0         = best_par[0];
	  out_NL         = best_par[1];
	  out_NR         = best_par[2];
	  out_Wjets      = best_par[3];
	  out_QCD        = best_par[4];
	  out_RemBkg     = best_par[5];
	  out_Nges       = alpha_sum;
	  out_N0_err     = par_err[0];
	  out_NL_err     = par_err[1];
	  out_NR_err     = par_err[2];
	  out_Wjets_err  = par_err[3];
	  out_QCD_err    = par_err[4];
	  out_RemBkg_err = par_err[5]; 
	  
	  for (int k = 0; k != fNNui; k++){
	    
	    // has to be filled with the correct values later
	    
	    out_nui_nom[k] = Nui_Vec[k];
	    out_nui[k]     = best_par[k+fNScaleFactors];
	    out_nui_err[k] =  par_err[k+fNScaleFactors];
	    
	  }
	  
	  out_N0_nom     = fXsec*fLumi*F0_true;
	  out_NL_nom     = fXsec*fLumi*FL_true;
	  out_NR_nom     = fXsec*fLumi*FR_true;
	  out_Wjets_nom  = fNominalScaleParameters[3];
	  out_QCD_nom    = fNominalScaleParameters[4];
	  out_RemBkg_nom = fNominalScaleParameters[5];
	  
	  out_Nges_nom = out_N0_nom+out_NL_nom+out_NR_nom;
	  
	  // has to be made more general!!!!
	  	  
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){
	    
	    std::string signal = fSignalParameter[iSig].ParamName;
	    
	      if(signal == "N0")
		out_N0_prof_err     = error_profile_vec[0];
	      else if(signal == "NL")
		out_NL_prof_err     = error_profile_vec[1];
	      else if(signal == "NR")
		out_NR_prof_err     = error_profile_vec[2];
	      
	    }
	    
	    for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
	      
	      std::string bkg = fBkgParameter[iBkg].ParamName;
	      
	      if(bkg == "Wjets")
		out_Wjets_prof_err  = error_profile_vec[3];
	      else if(bkg == "QCD")
		out_QCD_prof_err    = error_profile_vec[4];
	      else if(bkg == "RemBkg")
		out_RemBkg_prof_err = error_profile_vec[5];   
	      
	    }
	    

	    for (int iNui = 0; iNui != fNNui; iNui++){
	    
	      out_nui_prof_err[iNui] = error_profile_vec[iNui+fNScaleFactors];
	    
	    }
	    
	    if(best_par[fNScaleFactors] > 0) k_counter++;
	 

	    if(out_MinuitStatus == 0){
	      
	      fNewTree -> Fill();
	      counter_good++;
	      
	    }
	//counter_good++;

      //if(best_par[6] > 0) k_counter++;
      
      ///_________________Interpolation histograms______________
	  
	  //plot interpolation histograms for one nuisance parameter
	  //use just one nuisance parameter otherwise too many output histograms!!!
	  
	  /*	    TH1D * hist_interp_F0 = new TH1D("hist_interp_F0", "hist_interp_F0", 15, -1.0, 1.0);
	    TH1D * hist_interp_FL = new TH1D("hist_interp_FL", "hist_interp_FL", 15, -1.0, 1.0);
	    TH1D * hist_interp_FR = new TH1D("hist_interp_FR", "hist_interp_FR", 15, -1.0, 1.0);		
	    
	    hist_interp_F0 = interpolate_single(hist_vec_F0_scaled, hist_interp_F0, best_par[fNScaleFactors]);
	    hist_interp_FL = interpolate_single(hist_vec_FL_scaled, hist_interp_FL, best_par[fNScaleFactors]);
	    hist_interp_FR = interpolate_single(hist_vec_FR_scaled, hist_interp_FR, best_par[fNScaleFactors]);
	    
	    char k[100] = "PE = %d  /  k = %.4f #pm %.4f";// use %f for double, %d for int
	    char latex_k[sizeof k + 100];
	    sprintf(latex_k, k, i, best_par[fNScaleFactors], par_err[fNScaleFactors] );
	    
	    hist_interp_F0->SetTitle(latex_k);	
	    hist_interp_FL->SetTitle(latex_k);
	    hist_interp_FR->SetTitle(latex_k);
	    
	    char title_F0[100] = "F0_%d";// use %f for double, %d for int
	    char latex_title_F0[sizeof title_F0 + 100];
	    sprintf(latex_title_F0, title_F0, i );
	    char title_FL[100] = "FL_%d";// use %f for double, %d for int
	    char latex_title_FL[sizeof title_FL + 100];
	    sprintf(latex_title_FL, title_FL, i );
	    char title_FR[100] = "FR_%d";// use %f for double, %d for int
	    char latex_title_FR[sizeof title_FR + 100];
	    sprintf(latex_title_FR, title_FR, i );
	    
	    hist_interp_F0->Write(latex_title_F0);
	    hist_interp_FL->Write(latex_title_FL);
	    hist_interp_FR->Write(latex_title_FR);
	    
	    hist_interp_F0->Delete();
	    hist_interp_FL->Delete();
	    hist_interp_FR->Delete();
	  
	     */
      

	}

      //cout << "Close tree" << endl;
      
      fNewTree -> AutoSave();
      
      fNewFile -> Close();
      
      fitter   -> Delete();
      
    
      //    }  

}




///_________________________Start Data Fit______________________________

//Parameter distributions based on pseudo-data: nui_var is the index of the parameter that is varied, Nui_Vec contains all values for the nui param that are used to generate the pseudo data
void tvalidation::Datafit(std::vector<double> Nui_Vec, std::string input_channel, std::string datafit_mode)
{
  
     //Set start values to SM values
      double F0_true = 0.698; 
      double FL_true = 0.301; 
      double FR_true = 0.00041;

      std::string output_file;

      output_file = ("Datafit_Output_"+input_channel+".root");


      // create file for datafit
      
      WriteInfoStatus("ValidationClass", "Create new file: "+output_file);
      
      TFile *fNewFile = new TFile(output_file.c_str(), "RECREATE");
      
      WriteInfoStatus("ValidationClass", "Create new output tree");
      
      TTree *fNewTree = InitializeOutputTree();
      


	//Define output data files
	std::string output = "txt_output/datafit_results_" + input_channel + ".txt";

	const char *filename_output;
	filename_output=output.c_str();

	ofstream ofile(filename_output);

	ofile.setf(ios::fixed); //in order to count leading zeros in set precision 


      ///_____________________________FIT_______________________________
      
      WriteInfoStatus("ValidationClass", "Initialize TMinuit");

      //Setting up the fitting procedure; initialize TMinuit
      TMinuit *fitter = new TMinuit(50);
      
      //int counter_good = 0;
      
      //while(counter_good < fNPseudoExp)
//	{

	  // clear the output tree
	  clean();

	  //Setting up the fitting procedure; initialize TMinuit
	  if(fitter) fitter->Delete();
	  fitter = new TMinuit(50);
      	  fitter->SetPrintLevel(-1);

	  int ierflg_err=0, ierflg = 0, ierflg_hesse = 0, par_num;
	  double start_val, start_step=0.001, A=0.0, B=0.0;  
	  
	  double arglist[50];

	  //Set interpolation object for fit (at first only linear)
	  fitter->SetFCN(eval_chisqrt);	
	  
	  //Generate new start values
	  std::vector <double> start_val_gauss;
	  std::vector <double> start_val_vec;
	  
	  double Ntotal = fXsec*fLumi;
	  
	  start_val_vec.push_back(Ntotal*F0_true);
	  start_val_vec.push_back(Ntotal*FL_true);
	  start_val_vec.push_back(Ntotal*FR_true);

	  if(fBkgMode == "BkgFit"){
	
	    for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg)
	      start_val_vec.push_back(fBkgExpVec[iBkg] );

	  }

	  //Define 6 scale parameters alpha for signal and background
	  for (int k = 0; k != fNScaleFactors; k++)
	    {
	      
	      std::string s;
	      std::stringstream out;
	      out << k;
	      s = out.str();

	      fitter->mnparm(par_num=k, "alpha_" + s, start_val_vec[k], start_step, A, B, ierflg);
	      
	      /*	      if(k < 2)
		fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, 50000, 550000, ierflg);
	      else if(k == 2)
		fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, -50000, 50000, ierflg);
	      else
	      fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, 0.0, 7500.0, ierflg); */


	    }
	  
	  //Define nuisance parameters k as free pararameters
	  for (int k = 0; k != fNNui; k++){

	    //probably nui-param for all three states + three backgrounds 	
	    //-> 6 nui-alphas + all nuisance parameters k 
	    //if(i < 6) fitter->mnparm(par_num=i, "alpha_" + used_nui_parameters[i], start_val=0.0, start_step, A, B, ierflg); 
	    
	    std::string s;
	    std::stringstream out;
	    out << k;
	    s = out.str();
	    
	    //start_val = gaussRandom->Gaus(0.0, 1.0);
	    start_val = 0.0;	    

	    fitter->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, A, B, ierflg);
	    // debug A.Knue
	    //	    fitter->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, -4.0, 4.0, ierflg);
	    
	    
	  }
	  
	  const int par_size = fNScaleFactors+fNNui;//6 scale parameters

	  std::stringstream oss;
	  oss << par_size;

	  WriteParameterStatus("ValidationClass", "Number of parameters in data fit: "+oss.str());
	  
	  arglist[0] = 1.0;
	  fitter->mnexcm("SET ERR", arglist, 1, ierflg_err);

	  //Starting the fit
	  arglist[0] = 10000.0;
	  arglist[1] = 0.01;
	  fitter->mnexcm("MIGRAD", arglist, 2, ierflg);

	  arglist[0] = 1000.0;
	  arglist[1] = 0;
	  //not needed for tests in goe, but consider errors larger than 1 then
	  fitter->mnexcm("HESSE", arglist, 2, ierflg_hesse);//use HESSE/MINOS alternatively
	  
	  //Print Results; access to fit parameter
	  double fval,edm,errdef;
	  //fval: the best function value found so far; edm: the estimated vertical distance remaining to minimum; errdef: the value of UP defining parametEr uncertainties
	  int nvpar,nparx,icstat;//nvpar: the number of currently variable parameters, nparx: the highest (external) parameter number defined by user, icstat: a status integer indicating how good is the covariance
	  
	  fitter->mnstat(fval,edm,errdef,nvpar,nparx,icstat);
	  //fitter->mnprin(3,fval);
	  //3 refers to: values, errors, step sizes, first derivs.

	  //Get fit parameters
	  std::vector<double> best_par; 
	  std::vector<double> par_err;
	  double par_value, par_err_value;
	  best_par.clear();
	  par_err.clear();
	  	  
	  for (int k = 0; k != par_size; k++)
	    { 
	      fitter->GetParameter(k,par_value,par_err_value);
	      best_par.push_back(par_value);
	      par_err.push_back(par_err_value);
	      //      std::cout << par_err_value << std::endl;

	      std::stringstream oss_k,oss_par,oss_err;
	      oss_k << k;
	      oss_par << best_par[k];
	      oss_err << par_err[k];

	      WriteInfoStatus("ValidationClass", "Best parameter "+oss_k.str()+", par: "+oss_par.str()+", err: "+oss_err.str());

	    }
	  
	  //Get error matrix / covariance matrix
	  const int emsize = par_size;//k values and scale factors
	  double errormatrix[emsize][emsize];
	  fitter->mnemat(&errormatrix[0][0],emsize);

	  out_MinuitStatus   = ierflg;
	  out_NumNuisancePar = fNNui;

	int counter_good = 0;

	  //add new error calculation: profiling plot
	  std::vector<double> error_profile_vec;

	//error_profile_vec.push_back(par_err[0]);
	//error_profile_vec.push_back(par_err[1]);
	//error_profile_vec.push_back(par_err[2]);

	  for (int k = 0; k != par_size; k++) {
	    error_profile_vec.push_back( error_profile_PE(hist_data, best_par[k], k, best_par, par_err, counter_good, 1, F0_true, FL_true, FR_true, 1) ) ;
	  
	    //error_profile_vec.push_back(1.0);
	    }
	      
	

	  //Sum of all fitted signal parameters alpha_(F0,FL,FR)
	  double alpha_sum = 0.0;
	  
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig)
	    alpha_sum += best_par[iSig];
	  
	  //Calculate helicity fractions ...
	  
	  if(fNSignal > 2){
	    out_F0 = eval_fi (fHistNomVec[0], best_par[0], alpha_sum);
	    out_FL = eval_fi (fHistNomVec[1], best_par[1], alpha_sum);
	    out_FR = eval_fi (fHistNomVec[2], best_par[2], alpha_sum);
	    
	    //... and related errors:
	    out_F0_err = fi_error (&errormatrix[0][0], emsize, best_par[0], alpha_sum, 0);
	    out_FL_err = fi_error (&errormatrix[0][0], emsize, best_par[1], alpha_sum, 1);
	    out_FR_err = fi_error (&errormatrix[0][0], emsize, best_par[2], alpha_sum, 2);
	    //last int represents index according to: 0 corresponds to F0, 1 to FL, 2 to FR!
	    
	  }
	  


	  // store nominal (true) values;
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){
	    
	    std::string signal = fSignalParameter[iSig].ParamName;
	    
	    if(signal == "N0")
	      out_F0_nom     = F0_true;
	    else if(signal == "NL")
	      out_FL_nom     = FL_true;
	    else if(signal == "NR")
	      out_FR_nom     = FR_true;
	  }
	  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
	    
	    std::string bkg = fBkgParameter[iBkg].ParamName;
	    
	    if(bkg == "Wjets")
	      out_Wjets_nom  = fBkgExpVec[0];
	    else if(bkg == "QCD")
	      out_QCD_nom    = fBkgExpVec[1];
	    else if(bkg == "RemBkg")
	      out_RemBkg_nom = fBkgExpVec[2];
	    
	  }
	  

	  out_N0         = best_par[0];
	  out_NL         = best_par[1];
	  out_NR         = best_par[2];
	  out_Nges       = alpha_sum;
          out_N0_err     = par_err[0];
          out_NL_err     = par_err[1];
          out_NR_err     = par_err[2];

	  if(fBkgMode == "BkgFit"){
	    out_Wjets      = best_par[3];
	    out_QCD        = best_par[4];
	    out_RemBkg     = best_par[5];
	    out_Wjets_err  = par_err[3];
	    out_QCD_err    = par_err[4];
	    out_RemBkg_err = par_err[5];

	  }
	  
	  for (int k = 0; k != fNNui; k++){
	    
	    // has to be filled with the correct values later
	    
	    out_nui_nom[k] = Nui_Vec[k];
	    out_nui[k]     = best_par[k+fNScaleFactors];
	    out_nui_err[k] =  par_err[k+fNScaleFactors];
	    
	  }
	  
	  
	  out_N0_nom     = fXsec*fLumi*F0_true;
	  out_NL_nom     = fXsec*fLumi*FL_true;
	  out_NR_nom     = fXsec*fLumi*FR_true;

	 // out_Wjets_nom  = fNominalScaleParameters[3];
	 // out_QCD_nom    = fNominalScaleParameters[4];
	 // out_RemBkg_nom = fNominalScaleParameters[5];
	  
	  out_Nges_nom = out_N0_nom+out_NL_nom+out_NR_nom;
	  
	  // has to be made more general!!!!
	  	  
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){
	    
	    std::string signal = fSignalParameter[iSig].ParamName;
	    
	      if(signal == "N0")
		out_N0_prof_err     = error_profile_vec[0];
	      else if(signal == "NL")
		out_NL_prof_err     = error_profile_vec[1];
	      else if(signal == "NR")
		out_NR_prof_err     = error_profile_vec[2];
	      
	    }
	    
	  
	  if(fBkgMode == "BkgFit"){

	    for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
	      
	      std::string bkg = fBkgParameter[iBkg].ParamName;
	      
	      if(bkg == "Wjets")
		out_Wjets_prof_err  = error_profile_vec[3];
	      else if(bkg == "QCD")
		out_QCD_prof_err    = error_profile_vec[4];
	      else if(bkg == "RemBkg")
		out_RemBkg_prof_err = error_profile_vec[5];   
	      
	    }
	  
	  }

	  for (int iNui = 0; iNui != fNNui; iNui++){
	    
	    out_nui_prof_err[iNui] = error_profile_vec[iNui+fNScaleFactors];
	    
	  }

	     
	       //if(best_par[fNScaleFactors] > 0) k_counter++;
	       
	       
	       if(out_MinuitStatus == 0){
	       
	       fNewTree -> Fill();
	       //counter_good++;
	       
	       }
	  
	  //Data output in datafile ...:



	///__________________________RESULTS IN FILE________________________



	ofile << endl;
	ofile << "*********__________RESULTS OF THE TMINUIT-FIT_________**********" << endl;
	ofile << endl;
	ofile << "*************************Fit Parameters*************************" << endl;
	ofile << "****************************************************************" << endl;
	ofile << setw(25) << left << "fval : " << setw(14) << right << setprecision(3) << fval << endl;
	//Signal
	for (int i = 0; i != fSignalParameter.size(); i++)
	{	
		ofile << setw(3) << left << i+1 << setw(22) << left << " alpha (" + fSignalParameter[i].ParamName + ") :" << setw(14) << right << setprecision(3) << best_par[i] << "  with error: " << setw(11) << right << setprecision(3) << par_err[i] << " and " << setw(11) << right << setprecision(3) << error_profile_vec[i] << endl;
	}

	if(fBkgMode == "BkgFit"){

	  //Background
	  for (int i = 0; i != fBkgParameter.size(); i++)
	    {	
	      ofile << setw(3) << left << i+1+fSignalParameter.size()<< setw(22) << left << " alpha (" + fBkgParameter[i].ParamName + ") :" << setw(14) << right << setprecision(3) << best_par[i+fSignalParameter.size()] << "  with error: " << setw(11) << right << setprecision(3) << par_err[i+fSignalParameter.size()] << " and " << setw(11) << right << setprecision(3) << error_profile_vec[i+fSignalParameter.size()] << endl;
	    }
	}

	//Nuisance parameters
	for (int i = 0; i != fNNui; i++)
	{	
		std::stringstream istr;
		istr << i;
		//if(i < 6) ofile << setw(3) << left << i+1 << setw(16) << left <<  " alpha_" + used_parameters[i] + " : " << setw(14) << right << setprecision(10) << best_k[i] << "  with error: " << setw(14) << right << setprecision(10) << k_err[i] << endl;
		ofile << setw(3) << left << i+1+fNScaleFactors << setw(22) << left << (" k (" + fNuisanceParameter[i] + ") :" ).c_str()  << setw(14) << right << setprecision(5) << best_par[i+fNScaleFactors] << "  with error: " << setw(11) << right << setprecision(5) << par_err[i+fNScaleFactors] << " and " << setw(11) << right << setprecision(3) << error_profile_vec[i+fNScaleFactors] << endl;
	}

	ofile << "****************************************************************" << endl;
	ofile << endl;



	char line[200];
	double vline[25];
	double eline[25];

	ofile << "*********************** Correlation Matrix: ********************" << endl;
	ofile << endl;
	ofile << "Par. No.  Global           ";
	for (int id = 1; id <= gMinuit->fNu; ++id)  ofile << setw(14) << left << gMinuit->fNexofi[id-1];
	ofile << endl;
	for (int i = 1; i <= gMinuit->fNu; ++i) 
	{
		int ix  = gMinuit->fNexofi[i-1];
		int   ndi = i*(i + 1) / 2;
    		for (int j = 1; j <= gMinuit->fNu; ++j) 
		{
      			int m = TMath::Max(i,j);
      			int n = TMath::Min(i,j);
      			int ndex = m*(m-1) / 2 + n;
      			int  ndj  = j*(j + 1) / 2;
      			vline[j-1] = gMinuit->fVhmat[ndex-1] / TMath::Sqrt(TMath::Abs(gMinuit->fVhmat[ndi-1]*gMinuit->fVhmat[ndj-1]));
    		}
    		//fprintf(ofile, line," %1d   %14.6g ",ix,gMinuit->fGlobcc[i-1]);
    		ofile << setw(5) << left << ix << setw(12) << right << setprecision(6) << gMinuit->fGlobcc[i-1];
		//ofile << line;
    		for (int it = 1; it <= gMinuit->fNu; ++it) 
		{
			//fprintf(ofile, line, "  %12.5g",vline[it-1]);
			ofile << "  " << setw(12) << right << setprecision(5) << right << vline[it-1];
	 	}
		ofile << endl;
	}
	
	ofile << "*****************************************************************" << endl;
	ofile << endl;

	ofile << "**************HELICITY RESULTS***************" << endl;
	ofile << setw(3) << left << setw(10) << left <<  " F0 : " << setw(10) << right << setprecision(5) << out_F0 << "  with error: " << setw(10) << right << setprecision(5) << out_F0_err << endl;
	ofile << setw(3) << left << setw(10) << left <<  " FL : " << setw(10) << right << setprecision(5) << out_FL << "  with error: " << setw(10) << right << setprecision(5) << out_FL_err << endl;
	ofile << setw(3) << left << setw(10) << left <<  " FR : " << setw(10) << right << setprecision(5) << out_FR << "  with error: " << setw(10) << right << setprecision(5) << out_FR_err << endl;
	ofile << "*********************************************" << endl;
	ofile << endl;

///_____________________CALCULATE CROSS-SECTION___________________________


	double xsec, xsec_err;
	double BR = 0.543, BR_err = 0.0037;
	double lumi_corr = 4655.74, lumi_corr_err = 0.0;//pb^-1

	xsec = out_Nges/(lumi_corr*BR);

	double N_err = sqrt( out_N0_err*out_N0_err + out_NL_err*out_NL_err + out_NR_err*out_NR_err);
	//error caluclation based on error in Nges (lumi and BR too small)
	xsec_err = xsec*sqrt(   N_err*N_err/(out_Nges*out_Nges) + BR_err*BR_err/(BR*BR) +  lumi_corr_err*lumi_corr_err/(lumi_corr*lumi_corr) );

	cout << "Nges : " << out_Nges << endl;

	ofile << "************CROSS-SECTION RESULT***************" << endl;
	ofile << setw(3) << left << setw(10) << left <<  " XSEC (pb) : " << setw(10) << right << setprecision(5) << xsec << "  with error: " << setw(10) << right << setprecision(5) << xsec_err << endl;






	//Perform pseudo-experiments
	
	if (datafit_mode == "Full") {
	  
 	  std::vector<double> Nui_Vec;
	  
      	  Nui_Vec.clear();
	  
          for(int iNui = 0; iNui < fNNui; ++iNui) {     
	    Nui_Vec.push_back(best_par[iNui+fNScaleFactors]);
	  }
	  
	  std::string outputfile_mode = "Datafit";
	  
	//out_F0 = 0.08513;
	//out_FL = 0.05038;
	//out_FR = 0.03935;
	  
       	  // for counter fractions == 3 we get the SM configuration
      	  Validation(3, 0, Nui_Vec, 0.1, out_F0, out_FL, out_FR, input_channel, outputfile_mode);
	  
	}
	
	
	
	///_____________________Correlation based on Profiling___________________
	
	if (datafit_mode == "Single") {
	  
	  std::vector<double> vec_corr;
	  
	  std::vector<double> Nui_Vec;
	
	  Nui_Vec.clear();
	  
          for(int iNui = 0; iNui < fNNui; ++iNui) {     
	    Nui_Vec.push_back(best_par[iNui+fNScaleFactors]);
	  }
	  
	  vec_corr = Validation_corr(3, 0, Nui_Vec, 0.1, out_F0, out_FL, out_FR, input_channel);
	  
	  //vec_corr.push_back(-0.776);
	  //vec_corr.push_back(-0.625);
	  //vec_corr.push_back(0.854);
	  
	  
	//part from fillhistograms ...

	  std::vector<double> vec_error;
	  vec_error.clear();
	  vec_error.push_back(error_profile_vec[0]);
	  vec_error.push_back(error_profile_vec[1]);
	  vec_error.push_back(error_profile_vec[2]);
	  
	  
	  const int emsize_corr=3;
	  double errormatrix_corr[emsize_corr][emsize_corr];
	  
	  for(int i = 0; i < 3; ++i){
	    errormatrix_corr[i][i] = vec_error[i]*vec_error[i];
	    for (int j = i+1; j<3; ++j) {
	      errormatrix_corr[i][j] = vec_error[i]*vec_error[j]*vec_corr[i+j-1];
	      errormatrix_corr[j][i] = vec_error[i]*vec_error[j]*vec_corr[i+j-1];
	      //cout << errormatrix[i][j] << endl;
	    }
	  }
	  
	  
	  
	  // 	h_cov_N0NL->Fill(errormatrix[0][1]);
	  // 	h_cov_N0NR->Fill(errormatrix[0][2]);
	  // 	h_cov_NLNR->Fill(errormatrix[1][2]);
	  
	  //cout << best_par[0] << " " << best_par[1] << " " << best_par[2] << " " << endl;
	  //cout << alpha_sum << endl;
	  
	  //... determine errors on Fi then :
	  double sigma_F0_corr = fi_error (&errormatrix_corr[0][0], emsize_corr, best_par[0], alpha_sum, 0);
	  double sigma_FL_corr = fi_error (&errormatrix_corr[0][0], emsize_corr, best_par[1], alpha_sum, 1);
	  double sigma_FR_corr = fi_error (&errormatrix_corr[0][0], emsize_corr, best_par[2], alpha_sum, 2);
	  
	  ofile << endl;
	  ofile << "**************CORRELATION****************" << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " N0NL : " << setw(10) << right << setprecision(5) << vec_corr[0] << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " N0NR : " << setw(10) << right << setprecision(5) << vec_corr[1] << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " NLNR : " << setw(10) << right << setprecision(5) << vec_corr[2] << endl;
	  ofile << endl;
	  ofile << "**************NEW CORRELATION****************" << endl;
	  ofile << "**************HELICITY RESULTS***************" << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " F0 : " << setw(10) << right << setprecision(5) << out_F0 << "  with error: " << setw(10) << right << setprecision(5) << sigma_F0_corr << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " FL : " << setw(10) << right << setprecision(5) << out_FL << "  with error: " << setw(10) << right << setprecision(5) << sigma_FL_corr << endl;
	  ofile << setw(3) << left << setw(10) << left <<  " FR : " << setw(10) << right << setprecision(5) << out_FR << "  with error: " << setw(10) << right << setprecision(5) << sigma_FR_corr << endl;
	  ofile << "*********************************************" << endl;
	  ofile << endl;
	  
	  
	}
	
	ofile.close();
	
///_________________________Plot based on data !!!_____________________:
	
	
	std::vector<double> k_final;
	
	cout << endl;
	cout << "Hist_output" << endl;
	cout << endl;

	std::vector<TH1D> fParameterHistNorm;
	fParameterHistNorm.clear();

	//Define histograms and set bins
	TH1D hist_sum  = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
	TH1D hist_helper;
	TH1D hist_bkg = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);	
	//TH1D hist_ratio = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

	for(int iParam = 0; iParam < fNSignal + fNBkg; ++iParam)
	{
            fParameterHistNorm.push_back(normalise_fast(fParameterHist[iParam]));
	    //cout << "hist bins " << fParameterHistNorm[iParam].GetNbinsX() << endl;
	}

	for(int k = fNScaleFactors; k < fNScaleFactors+fNNui; ++k)
	  k_final.push_back(best_par[k]);
	
	//Consider interpolation since bkg fitted via nuisance parameters
	//with additional correction factor because error sometimes smaller than 100% 


	if (fNNui > 0) {
		hist_sum = ReturnInterpolationHistogramScaled(0, fInterpolationMethod, k_final, best_par[0], false);//out_Nges_nom
	}
	else {
		hist_helper = fParameterHistNorm[0];
            	hist_helper.Scale(best_par[0]*SelEff[0]);
	    
	  	hist_sum = hist_sum + hist_helper;
	}


	//cout << hist_sum.GetBinContent(1) << endl;

	for (int i = 1; i != fSignalParameter.size(); i++) {

		if (fNNui > 0) {
		hist_sum = hist_sum + ReturnInterpolationHistogramScaled(i, fInterpolationMethod, k_final, best_par[i], false);
		}
		else {
		hist_helper = fParameterHistNorm[i];
            	hist_helper.Scale(best_par[i]*SelEff[i]);
	    
	  	hist_sum = hist_sum + hist_helper;

		std::cout << i << "\t" << hist_sum.Integral() << std::endl;

		}

	//cout << hist_sum.GetBinContent(1) << endl;
	}

	//Bkg histogram
	//TH1D hist_bkg = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
	if (fBkgParameter.size() != 0){

	if (fNNui > 0) {
	    	hist_bkg = ReturnInterpolationHistogramScaled(fNSignal, fInterpolationMethod, k_final, best_par[fNSignal], true);
	//cout << hist_bkg.GetBinContent(1) << endl;
		//hist_bkg.Scale(best_par[0+fSignalParameter.size()]);
	  //		hist_bkg = hist_bkg*best_par[0+fSignalParameter.size()];
	}
	else {
	        hist_helper = fParameterHistNorm[fNSignal];
	        hist_helper.Scale(best_par[fNSignal]);
	 
		hist_bkg = hist_bkg + hist_helper;

		std::cout << hist_bkg.Integral() << std::endl;

	}

	for (int i = 1; i != fBkgParameter.size(); i++) {

		if (fNNui > 0) {		
		hist_bkg = hist_bkg + ReturnInterpolationHistogramScaled(i+fNSignal, fInterpolationMethod, k_final, best_par[i+fNSignal], true);
		}
		else {
	        hist_helper = fParameterHistNorm[i+fNSignal];
	      	hist_helper.Scale(best_par[i+fNSignal]);
	  	hist_bkg = hist_bkg + hist_helper;

		std::cout << i << "\t" << hist_bkg.Integral() << std::endl;

		}
	}


	}

	hist_sum = hist_sum + hist_bkg;
	//cout << hist_sum.GetBinContent(1) << endl;

	//SM histogram

	//double sum_input = hist_F0_nom->Integral() + hist_FL_nom->Integral() + hist_FR_nom->Integral();

	//cout << "Sum input " << sum_input << endl;

	std::cout << "Alles aufaddieren...." << std::endl;

	TH1D hist_sum_SM;
	hist_sum_SM = fHistNormVec[0];
	//hist_sum_SM.Scale(fXsec*fLumi*0.698*fEffVec[0]);


	hist_sum_SM = hist_sum_SM*fXsec*fLumi*0.698*fEffVec[0];
	cout << hist_sum_SM.Integral() <<  "\t" << fXsec*fLumi*0.698*fEffVec[0]  << "\t" << fEffVec[0] << endl;
	hist_sum_SM = hist_sum_SM + fHistNormVec[1]*fXsec*fLumi*0.301*fEffVec[1];
	
	cout << fHistNormVec[1].Integral() <<  "\t" << fXsec*fLumi*0.301*fEffVec[1] << "\t" << fEffVec[1] << endl;
	//cout << hist_sum_SM.GetBinContent(1) << " " << fEffVec[1] << endl;
	hist_sum_SM = hist_sum_SM + fHistNormVec[2]*fXsec*fLumi*0.00041*fEffVec[2];

	cout << fHistNormVec[2].Integral() <<  "\t" << fXsec*fLumi*0.00041*fEffVec[2] << "\t" << fEffVec[2] << endl;

	//hist_sum_SM->Add(fHistNormVec[1],fXsec*fLumi*0.301*fEffVec[1]);
	//hist_sum_SM->Add(fHistNormVec[2],fXsec*fLumi*0.00041*fEffVec[2]);
	//cout << hist_sum_SM.GetBinContent(1) << " " << fEffVec[2] << endl;
	

	for (int i = 0; i != fBkgParameter.size(); i++) {
		hist_sum_SM = hist_sum_SM + fHistNormVec[i+fSignalParameter.size()]*fHistNomVec[i+fSignalParameter.size()].Integral();

		cout << (fHistNormVec[i+fSignalParameter.size()]*fHistNomVec[i+fSignalParameter.size()].Integral()).Integral() << endl;

		//cout << "int : " <<  fHistNormVec[i+fSignalParameter.size()].Integral() << endl;
		//hist_sum_SM->Add( fHistNormVec[i+fSignalParameter.size()],fHistNormVec[i+fSignalParameter.size()].Integral() );
	}

	//std::cout << hist_sum_SM.Integral() << std::endl;

	//TH1D hist_ratio = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);
	//TH1D hist_ratio_err = TH1D("", "", Nbins, -1.0, 1.0);


	//Redefine histos with new bin widths
	
	TH1D hist_sum_adj  = *(TH1D*) hist_data.Clone("");
	TH1D hist_bkg_adj = *(TH1D*) hist_data.Clone("");	
	TH1D hist_ratio = *(TH1D*) hist_data.Clone("");
	for  (int i = 1; i != fNBins+1; i++) {
		hist_sum_adj.SetBinContent(i, hist_sum.GetBinContent(i) );
		hist_bkg_adj.SetBinContent(i, hist_bkg.GetBinContent(i) );
	} 

	//Calculate histogram errors
	//Think about this function later
	TGraphAsymmErrors * hist_unc = new TGraphAsymmErrors();
	hist_unc = TemplateFitErrors(&errormatrix[0][0], emsize, hist_sum_adj, fHistNormVec, fEffVec, hist_unc, fSignalParameter.size(), fBkgParameter.size() ); 

	for (int i = 1; i != fNBins+1; i++) {

		double fit = hist_sum_adj.GetBinContent(i);
		double data = hist_data.GetBinContent(i);
		double fit_error = hist_unc->GetErrorY(i);
		double data_error = hist_data.GetBinError(i);
		double ratio=1.0, ratio_error=0.0;

		if (fit != 0) ratio = data/fit;
		if (fit != 0 || data != 0) ratio_error = ratio*sqrt( data_error*data_error/(data*data) + fit_error*fit_error/(fit*fit) );

		hist_ratio.SetBinContent(i,ratio);
		hist_ratio.SetBinError(i,ratio_error);
	}


	std::string output_folder = "txt_output";

	//Plot fit results
	//Need pointer for the plot ...
	PlotDataFit(hist_sum_adj, hist_bkg_adj, hist_sum_SM, hist_data, hist_unc, hist_ratio, input_channel, output_folder);


	//hist_sum->Delete();
	//hist_sum_SM->Delete();
	//hist_bkg->Delete();
	//hist_unc->Delete();



      //if(best_par[6] > 0) k_counter++;
      
     
      fNewTree -> AutoSave();
      
      fNewFile -> Close();
      
      fitter   -> Delete();
      

      //    }  

}









///________________________________________________________________________
///_______________________Profiling to calculate errors____________________
///________________________________________________________________________

double tvalidation::error_profile_PE (TH1D hist_pseudo, double par_int, int number, std::vector<double> &par, std::vector<double> &par_err, int num_pe, int counter, double F0_true, double FL_true, double FR_true, int switcher)
{

///_____________________________Definitions______________________________

	//number of likelihood calculations
	int nmax = 50;

	//Choose parameters for likelihood plots
	//first k value + k_stepsize
	double par_step, par_stepsize;

	if(fBkgMode == "BkgFit"){

	  if (number < 3) {par_step = par_int-20000; par_stepsize = 800;}
	  else if (number > 2 && number < 6) {par_step = par_int-2000; par_stepsize = 80;}
	  else if (number > 5) {par_step = par_int-3.0; par_stepsize = 0.12;}

	}
	else if(fBkgMode == "BkgNui"){

          if (number < 3) {par_step = par_int-20000; par_stepsize = 800;}
	  else if (number > 2) {par_step = par_int-3.0; par_stepsize = 0.12;}

        }

	

///_________________________Calculate likelihood value___________________


	double par_step_ini = par_step;	

	std::vector<double> func_min;
	std::vector<int> check_min;
		
	//double k_var[ksize];
	std::vector<double> par_var;
	for (int j = 0; j != par.size(); j++)
	{
		par_var.push_back(par[j]);
	}
	//Change k-factor which should be plotted
	par_var[number] = par_step;

	int conv_counter=0;	
	
	//Loop for varying of the signal parameter
	for (int n = 0; n != nmax; n++)
	{
		par_var[number] += par_stepsize;

				
		//Redo the fit with fixed parameter "number"

///__________________________FITTING PROCEDURE_____________________________
		
		//Setting up the fitting procedure; initialize TMinuit
		TMinuit* fitter_prof = new TMinuit(50);
		fitter_prof->SetPrintLevel(-1);
	
		int ierflg_err=0, ierflg = 0, ierflg_hesse = 0, par_num;
		double start_val, start_step=0.01, A=0.0, B=0.0;  
		
		double arglist[50];
		
		//Set interpolation object for fit (at first only linear)
		if (switcher == 0) fitter_prof->SetFCN(eval_chisqrt_vali);	
		else fitter_prof->SetFCN(eval_chisqrt);	
		
		//Generate new start values
		std::vector <double> start_val_gauss;

		start_val_gauss.push_back( gaussRandom->Gaus(fXsec*fLumi*F0_true, 20000.0 ) );
		start_val_gauss.push_back( gaussRandom->Gaus(fXsec*fLumi*FL_true, 20000.0) );
		start_val_gauss.push_back( gaussRandom->Gaus(fXsec*fLumi*FR_true, 20000.0) );

		if(fBkgMode == "BkgFit" && fBkgExpVec.size() > 0){

		  start_val_gauss.push_back( gaussRandom->Gaus( fBkgExpVec[0], fBkgUncVec[0]) );
		  start_val_gauss.push_back( gaussRandom->Gaus( fBkgExpVec[1], fBkgUncVec[1]) );
		  start_val_gauss.push_back( gaussRandom->Gaus( fBkgExpVec[2], fBkgUncVec[2]) );
		
		}

		//Start values
		//double start_val_array[] = {_xsec*_lumi*F0_true, _xsec*_lumi*FL_true, _xsec*_lumi*FR_true, fBkgExpVec[0], fBkgExpVec[1], fBkgExpVec[2]};

		std::stringstream oss, oss2;
		oss  << fNScaleFactors;
		oss2 << fNNui;

		//		WriteParameterStatus("ValidationClass", "ProfileError: Number of scale factors:  " + oss.str()  );
		//		WriteParameterStatus("ValidationClass", "ProfileError: Number of nui parameters: " + oss2.str() );

		
		//Define 6 scale parameters alpha for signal and background
		for (int k = 0; k != fNScaleFactors; k++)
		{
			std::string s;
			std::stringstream out;
			out << k;
			s = out.str();

			if (k == number) {
			  //	cout << "here number scale : " << k << endl; 
			fitter_prof->mnparm(par_num=k, "alpha_" + s, par_var[number], 0.0, A, B, ierflg); 
			}
			else {
			fitter_prof->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, A, B, ierflg); 
			}
		}

		//Define nuisance parameters k as free pararameters
		for (int k = 0; k != fNNui; k++)
		{
			//probably nui-param for all three states + three backgrounds 	
			//-> 6 nui-alphas + all nuisance parameters k 
			//if(i < 6) fitter->mnparm(par_num=i, "alpha_" + used_nui_parameters[i], start_val=0.0, start_step, A, B, ierflg); 
			
			std::string s;
			std::stringstream out;
			out << k;
			s = out.str();
			
			start_val = gaussRandom->Gaus(0.0, 1.0);
			//start_val = 0.0;

			if (k == number-fNScaleFactors) {
			  //	cout << "here number : " << i << endl; 
			fitter_prof->mnparm(par_num=k+fNScaleFactors, "k_" + s, par_var[number], 0.0, A, B, ierflg); 
			}
			else {
		 	fitter_prof->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, A, B, ierflg);
			} 
		}

		const int par_size = fNNui+fNScaleFactors;//6 scale parameters

		//		cout << "par size here  " << par_size << endl;

		arglist[0] = 1.0;
		fitter_prof->mnexcm("SET ERR", arglist, 1, ierflg_err);
	

		//Starting the fit
		arglist[0] = 10000.0;
		arglist[1] = 0.01;
		fitter_prof->mnexcm("MIGRAD", arglist, 2, ierflg);
		

		arglist[0] = 1000.0;
		arglist[1] = 0;
		//not needed for tests in goe, but consider errors larger than 1 then
		fitter_prof->mnexcm("HESSE", arglist, 2, ierflg_hesse);//use HESSE/MINOS alternatively
		
		//Print Results; access to fit parameter
		double fval,edm,errdef;
		//fval: the best function value found so far; edm: the estimated vertical distance remaining to minimum; errdef: the value of UP defining parametEr uncertainties
		int nvpar,nparx,icstat;//nvpar: the number of currently variable parameters, nparx: the highest (external) parameter number defined by user, icstat: a status integer indicating how good is the covariance

		fitter_prof->mnstat(fval,edm,errdef,nvpar,nparx,icstat);
		//fitter_prof->mnprin(3,fval);
		//3 refers to: values, errors, step sizes, first derivs.
		

		//if (ierflg_mig == 0) func_min.push_back(fval);
		//else conv_counter++;
		func_min.push_back(fval);
		check_min.push_back(ierflg);
		if (ierflg != 0)  conv_counter++;

		//Get fit parameters
		std::vector<double> best_par; 
		std::vector<double> par_err;
		double par_value, par_err_value;
		best_par.clear();
		par_err.clear();
		//double par_value, par_err_value;
		
		for (int k = 0; k != par_size; k++)
		{ 
			fitter_prof->GetParameter(k,par_value,par_err_value);
			best_par.push_back(par_value);
			par_err.push_back(par_err_value);
		}
			

		//cout << "Par " << 0 << " : "  << setprecision(6) << best_par[0] << " Par " << 7 << " : "  << setprecision(6) << best_par[7] << " fval : "  << setprecision(10) << fval << endl; 

		fitter_prof->Delete();

	}
	
	
	///_______________________ERROR CALCULATIONS_____________________________
	
	//Get smallest likelihood value to shift distribution
	double min = *min_element( func_min.begin(), func_min.end() );
	
	//Define TGraph for likelihood values
	TGraph *graph = new TGraph(3);
	TGraph *graph_test = new TGraph(3);
	
	par_step = par_step_ini;
	
	/*for (int n = 0; n != nmax; n++)
	  {
	  //cout << i << endl;
	  par_step += par_stepsize;
	  //if (check_min[n] == 0) {
	  graph->SetPoint( n+1,par_step,func_min[n]-min);	
	  //}
	  }*/
	
	//graph->RemovePoint(0);
	int checker = 0;
	int ngraph = 1;
	par_step += par_stepsize;
	
	int check_number = 0;

	if(fBkgMode == "BkgFit") check_number = 5;
	else if(fBkgMode == "BkgNui") check_number = 2;

	//Add criteria to remove values which differ from expected curve value probably due to a problems with convergence
	if (number > check_number) {
	  
	  for (int n = 1; n != nmax; n++)
	    {
	      par_step += par_stepsize;

	      //here additional stronger criterion in comparison to else case
	      if (check_min[n] == 0 && func_min[n]-min < 5) {
		
		if (checker == 0) {
		  if ( fabs(func_min[n] - func_min[n-1]) > 1.5) {
		    //graph->RemovePoint(n);
		    checker = 1;	
		  }
		  else {
		    graph->SetPoint( ngraph,par_step,func_min[n]-min);
		    ngraph++;
		  }
		}
		else if (checker == 1){
		  checker = 0;
		  if ( fabs(func_min[n] - func_min[n-2]) > 1.5) {
		    //graph->RemovePoint(n);
		  }
		  else {

		    graph->SetPoint( ngraph,par_step,func_min[n]-min);

		    ngraph++;
		  
		  }
		}
		
	      }
	      
	    }
	  
	}
	else {
	  
	  for (int n = 1; n != nmax; n++)
	    {
	      par_step += par_stepsize;
	      //par_step += par_stepsize;
	      
	      if (check_min[n] == 0) {
		
		if (checker == 0) {
		  if ( fabs(func_min[n] - func_min[n-1]) > 1.5) {
		    //graph->RemovePoint(n);
		    checker = 1;	
		  }
		  else {
		    graph->SetPoint( ngraph,par_step,func_min[n]-min);
		    ngraph++;
		}
		}
		else if (checker == 1){
		  checker = 0;
		  if ( fabs(func_min[n] - func_min[n-2]) > 1.5) {
		    //graph->RemovePoint(n);
		  }
		  else {
		    graph->SetPoint( ngraph,par_step,func_min[n]-min);
		    ngraph++;
		  }
		}
		
	      }
	      
	    }
	  
	}
	
	graph->RemovePoint(0);
	
	//Ordinate	
	//cout << number << " : " << func_min[0]-min << endl;
	//cout << par_int << endl;
	TF1 *fit;
	TF1 *fit_test;

	//	std::cout << fBkgMode.c_str() << std::endl;

	//Fit to curve, approximated with quadratic polynomial
	//Define fit function here
	if(fBkgMode == "BkgFit"){
	  
	  if (number < 3)                    fit      = new TF1("fit","[0]*x*x + [1]*x + [2]", par_int-20000,par_int+20000);
	  else if (number > 2 && number < 6) fit      = new TF1("fit","[0]*x*x + [1]*x + [2]", par_int-500,par_int+500);
	  else if (number > 5)               fit      = new TF1("fit","[0]*x*x + [1]*x + [2]", par_int-2.0,par_int+2.0);
	  
	  if (number < 3)                    fit_test = new TF1("fit_test","[0]*x*x + [1]*x + [2]", par_int-2000,par_int+2000);
	  else if (number > 2 && number < 6) fit_test = new TF1("fit_test","[0]*x*x + [1]*x + [2]", par_int-500,par_int+500);
	  else if (number > 5)               fit_test = new TF1("fit_test","[0]*x*x + [1]*x + [2]", par_int-2.0,par_int+2.0);
	  
	}
	else if (fBkgMode == "BkgNui"){

	  if (number < 3)      fit = new TF1("fit","[0]*x*x + [1]*x + [2]", par_int-20000,par_int+20000);
	  else if (number > 2) fit = new TF1("fit","[0]*x*x + [1]*x + [2]", par_int-2.0,par_int+2.0);

          if (number < 3)      fit_test = new TF1("fit_test","[0]*x*x + [1]*x + [2]", par_int-2000,par_int+2000);
	  else if (number > 2) fit_test = new TF1("fit_test","[0]*x*x + [1]*x + [2]", par_int-2.0,par_int+2.0);


	}

	//fit->SetParameter(0,1);
	//fit->SetParameter(2,100000); 
	fit -> SetLineColor(kCyan);
	fit -> SetLineWidth(2);

	graph -> Fit("fit","Q");


	//fit->GetParameter(1), fit->GetParError(1)
	double a = fit -> GetParameter(0);
	double b = fit -> GetParameter(1);
	double c = fit -> GetParameter(2);
	//double a = fit_test->GetParameter(0);
	//double b = fit_test->GetParameter(1);
	//double c = fit_test->GetParameter(2);

	//Determine error on par with help of quadratic fit
	double sigma_par = sqrt( (1.0-c)/a + b*b/(4*a*a) );
	
	//	std::cout << par_int << "\t" << sigma_par << "\t" << a << "\t" << b << "\t" << c << std::endl;


	///_______________________________PLOTS__________________________________
	
	//...usually not used if pseudo-experiments are performed:
	//otherwise remove /* and */
	///MAYBE OTHER OPTION NEEDED TO CALL THIS PART!!!	

	//Additional plot for profiling the likelihood if needed
	/*
	  TCanvas *c1 = new TCanvas("c1","c1", 700, 700);
	  c1->SetFillColor(0);
	  c1->SetGrid();
	  
	  //Define all necessary labels/titles
	  //char slabel_x[100] = "#scale[1.20]{Parameter (%d)}";
	  char slabel_x[100] = "#scale[1.20]{Parameter N}";
	  char label_x[sizeof slabel_x + 100];
	  sprintf(label_x, slabel_x, number);
	  
	  char stitle[100] = "plots_scans/scan_par%d_counter%d_pe%d.png";//%d for int, %f for double
	  char title[sizeof stitle + 100];
	  sprintf(title, stitle, number, counter, num_pe);
	  
	  char hist_stitle[100] = "Scan Par=%d Counter=%d PE=%d ";//%d for int, %f for double
	  char hist_title[sizeof stitle + 100];
	  sprintf(hist_title, hist_stitle, number, counter, num_pe);
	  
	  //graph->SetTitle(hist_title);
	  graph->SetTitle("");
	  graph->GetXaxis()->SetTitle(label_x);
	  graph->GetYaxis()->SetTitle("#scale[1.2]{-2ln(L)}");
	  graph->GetXaxis()->SetTitleOffset(1.2);
	  graph->GetYaxis()->SetTitleOffset(1.2);
	  graph->GetXaxis()->SetLabelSize(0.035);
	  graph->GetYaxis()->SetLabelSize(0.035);
	  graph->SetMinimum(-0.15);//(-20285.0);
	  //graph->SetMaximum(16.0);//-20281.7);
	  //graph->GetXaxis()->SetLimits(min,max);
	  
	  //graph->SetStats(kFALSE);
	  graph->SetFillStyle(0);
	  graph->SetLineColor(kBlue);
	  graph->SetMarkerStyle(20);
	  graph->SetMarkerColor(kBlue);
	  graph->Draw("AP");
	  
	  //not needed now with setlimits (see above)
	  //graph->GetXaxis()->SetRangeUser(-2.0,2.0);
	  //graph->Draw("AP");
	  
	  fit->SetLineColor(kCyan);
	  fit->SetLineWidth(2);
	  fit->Draw("same");
	  
	  TLegend *legend = new TLegend(0.70,0.77,0.89,0.89);
	  //legend->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
	  legend->SetTextFont(42);
	  legend->AddEntry(graph,"#scale[1.5]{Likelihood}","p");
	  legend->AddEntry(graph,"#scale[1.5]{value}",NULL);	
	  legend->AddEntry(fit,"#scale[1.5]{Fit}","l");	
	  //legend->AddEntry(hb_c,"Reconstructed as a b-jet","f");
	  legend->SetFillColor(0);
	  legend->SetBorderSize(0);
	  legend->Draw("same");
	  
	  
	  char config[50] = "Parameter Number %d";// use %f for double, %d for int
	  char subtitle[sizeof config + 100];
	  sprintf(subtitle, config, number);
	  
	  TLatex l; //l.SetTextAlign(12); 
	  l.SetTextSize(0.045); 
	  l.SetNDC();
	  l.SetTextFont(72);
	  //l.SetTextColor();
	  //l.DrawLatex(0.14,0.83,"ATLAS");
	  l.SetTextFont(42);
	  //l.DrawLatex(0.29,0.83,"Work in progress");
	  l.SetTextSize(0.041); 
	  //l.DrawLatex(0.14,0.77,"#alpha_{S}-Scan");
	  l.DrawLatex(0.18,0.85,"Parameter Scan");
	  l.SetTextSize(0.041); 
	  if(number == 0) l.DrawLatex(0.18, 0.79, "N_{0}");
	  if(number == 1) l.DrawLatex(0.18, 0.79, "N_{L}");
	  if(number == 2) l.DrawLatex(0.18, 0.79, "N_{R}");
	  if(number == 3) l.DrawLatex(0.18, 0.79, "N_{Wjets}");
	  if(number == 4) l.DrawLatex(0.18, 0.79, "N_{QCD}");
	  if(number == 5) l.DrawLatex(0.18, 0.79, "N_{RemBkg}");
	  //if(number == 6) l.DrawLatex(0.14, 0.76, "EER");
	  if(number == 6) l.DrawLatex(0.18, 0.79, "JES");
	  
	  c1->Print(title);	
	//*/
	


	delete fit;
	delete graph;

	//return value
	return sigma_par;
} 

//For correlation tests
std::vector<double> tvalidation::Validation_corr(int counter_fractions, int nui_var, std::vector<double> Nui_Vec, double diff, double F0_true, double FL_true, double FR_true, std::string input_channel)
{

      double nui_value = 0.0;

      for(int bla = 0; bla < Nui_Vec.size(); ++bla)
	std::cout << bla << "\t" << Nui_Vec[bla] << std::endl;

      if(nui_var <= Nui_Vec.size() && Nui_Vec.size() != 0)
	nui_value = Nui_Vec[nui_var-1];
      std::cout << nui_var << "\t" << nui_value << std::endl;

      WriteInfoStatus("ValidationClass", "Corr: Get pseudo data histogram");
      
      TH1D hist_pseudo;

      if( Nui_Vec.size() > 0)
	hist_pseudo = GetPseudodataFit(F0_true, FL_true, FR_true, Nui_Vec);
      else
	hist_pseudo = GetPseudodata(F0_true, FL_true, FR_true);

      WriteInfoStatus("ValidationClass", "Corr: Obtained pseudo data histogram");

      int k_counter    = 0;
      int conv_counter = 0;
      
      std::stringstream F0_label, nr_file, nui_num,output_nr,nui_var_str,nui_value_str;
      F0_label      << F0_true;
      nui_num       << fNNui;
      output_nr     << fNOutputSample;
      nui_var_str   << nui_var;
      if (fabs(nui_value) < 0.001) nui_value_str << 0;
      else nui_value_str << fabs(nui_value);
      
      TH2D *h2_N0NL = new TH2D("h2_N0NL", "", 100, fXsec*fLumi*F0_true-110000, fXsec*fLumi*F0_true+110000, 100, fXsec*fLumi*FL_true-110000, fXsec*fLumi*FL_true+110000);
      TH2D *h2_N0NR = new TH2D("h2_N0NR", "", 100, fXsec*fLumi*F0_true-110000, fXsec*fLumi*F0_true+110000, 100, fXsec*fLumi*FR_true-110000, fXsec*fLumi*FR_true+110000);
      TH2D *h2_NLNR = new TH2D("h2_NLNR", "", 100, fXsec*fLumi*FL_true-80000,  fXsec*fLumi*FL_true+80000,  100, fXsec*fLumi*FR_true-80000,  fXsec*fLumi*FR_true+80000);
      
      ///_____________________________FIT_______________________________
      
      WriteInfoStatus("ValidationClass", "Corr: Initialize TMinuit");

      //Setting up the fitting procedure; initialize TMinuit
      TMinuit *fitter_corr = new TMinuit(50);
      
      //      hist_ensemble = 0;	

      int counter_good = 0;
      int NumberPseudoExp = 1000;      //usually 1000
      //int NumberPseudoExp = 1;

      while(counter_good < NumberPseudoExp)
	{
	  std::stringstream oss;
	  oss << counter_good;

	  //WriteInfoStatus("ValidationClass", "Get new ensemble for pseudo experiments "+oss.str());

	  hist_ensemble = RandomPseudodata(hist_pseudo);
	  
	  if (counter_good % 10 == 0) cout << "Number of pseudo-experiments: " << counter_good << endl;

	  //	  hist_ensemble = hist_pseudo;

	  //Setting up the fitting procedure; initialize TMinuit
	  if(fitter_corr) fitter_corr->Delete();
	  fitter_corr = new TMinuit(50);
      	  fitter_corr->SetPrintLevel(-1);

	  int ierflg_err=0, ierflg = 0, ierflg_hesse = 0, par_num;
	  double start_val, start_step=0.001, A=0.0, B=0.0;  
	  
	  double arglist[50];

	  //Set interpolation object for fit (at first only linear)
	  fitter_corr->SetFCN(eval_chisqrt_vali);	
	  
	  //Generate new start values
	  std::vector <double> start_val_gauss;
	  
	  double Ntotal = fXsec*fLumi;
	  
	  // has to be made more general!!!
	  //	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*F0_true, 20000.0 ) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*F0_true, 20000.0 ) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FL_true, 20000.0) );
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FR_true, 20000.0) );
	  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg)
	    start_val_gauss.push_back( gaussRandom->Gaus(fBkgExpVec[iBkg], fBkgUncVec[iBkg]) );
	  
	  //Start values
	  double start_val_array[] = {Ntotal*F0_true};
	  // double start_val_array[] = {Ntotal*F0_true, Ntotal*FL_true, Ntotal*FR_true, fBkgExpVec[0], fBkgExpVec[1], fBkgExpVec[2]};
	  
	  //Define 6 scale parameters alpha for signal and background
	  for (int k = 0; k != fNScaleFactors; k++)
	    {
	      
	      std::string s;
	      std::stringstream out;
	      out << k;
	      s = out.str();

	      fitter_corr->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, A, B, ierflg);
	      
	    }
	  
	  //Define nuisance parameters k as free pararameters
	  for (int k = 0; k != fNNui; k++){
	    
	    std::string s;
	    std::stringstream out;
	    out << k;
	    s = out.str();
	    
	    start_val = gaussRandom->Gaus(0.0, 1.0);
	    
	    //	    fitter->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, A, B, ierflg);
	    // debug A.Knue
	    fitter_corr->mnparm(par_num=k+fNScaleFactors, "k_" + s, start_val, start_step, -5.0, 5.0, ierflg);
	    
	    
	  }
	  
	  const int par_size = fNScaleFactors+fNNui;//6 scale parameters
	  
	  arglist[0] = 1.0;
	  fitter_corr->mnexcm("SET ERR", arglist, 1, ierflg_err);

	  //Starting the fit
	  arglist[0] = 10000.0;
	  arglist[1] = 0.01;
	  fitter_corr->mnexcm("MIGRAD", arglist, 2, ierflg);

	  arglist[0] = 1000.0;
	  arglist[1] = 0;
	  //not needed for tests in goe, but consider errors larger than 1 then
	  fitter_corr->mnexcm("HESSE", arglist, 2, ierflg_hesse);//use HESSE/MINOS alternatively
	  
	  //Print Results; access to fit parameter
	  double fval,edm,errdef;
	  //fval: the best function value found so far; edm: the estimated vertical distance remaining to minimum; errdef: the value of UP defining parametEr uncertainties
	  int nvpar,nparx,icstat;//nvpar: the number of currently variable parameters, nparx: the highest (external) parameter number defined by user, icstat: a status integer indicating how good is the covariance
	  
	  fitter_corr->mnstat(fval,edm,errdef,nvpar,nparx,icstat);
	  //fitter->mnprin(3,fval);
	  //3 refers to: values, errors, step sizes, first derivs.

	  //Get fit parameters
	  std::vector<double> best_par; 
	  std::vector<double> par_err;
	  double par_value, par_err_value;
	  best_par.clear();
	  par_err.clear();
	  	  
	  for (int k = 0; k != par_size; k++)
	    { 
	      fitter_corr->GetParameter(k,par_value,par_err_value);
	      best_par.push_back(par_value);
	      par_err.push_back(par_err_value);
	      //      std::cout << par_err_value << std::endl;

	    }
	  
	  //Get error matrix / covariance matrix
	  const int emsize = par_size;//k values and scale factors
	  double errormatrix[emsize][emsize];
	  fitter_corr->mnemat(&errormatrix[0][0],emsize);

	  int MinuitStatus = ierflg;

	  //add new error calculation: profiling plot - not needed here!
	  std::vector<double> error_profile_vec;
	  //for (int k = 0; k != par_size; k++)
	  //  error_profile_vec.push_back( error_profile_PE(hist_pseudo, best_par[k], k, best_par, par_err, counter_good, 1, F0_true, FL_true, FR_true, 0) ) ;
	 

	  double value_N0 = best_par[0];
	  double value_NL = best_par[1];
	  double value_NR = best_par[2];

	  h2_N0NL->Fill(value_N0, value_NL);
	  h2_N0NR->Fill(value_N0, value_NR);
	  h2_NLNR->Fill(value_NL, value_NR);

	  if(MinuitStatus == 0){
	      
	      //fNewTree -> Fill();
	      counter_good++;
	      
	  }

	//end while
	}

      TFile *output = new TFile("Correlations_JES.root", "RECREATE");
      
      h2_N0NL->Write("N0NL");
      h2_N0NR->Write("N0NR");
      h2_NLNR->Write("NLNR");

      output -> Close();

  	std::vector<double> vec_corr;
  	vec_corr.push_back(h2_N0NL->GetCorrelationFactor() );
  	vec_corr.push_back(h2_N0NR->GetCorrelationFactor() );
  	vec_corr.push_back(h2_NLNR->GetCorrelationFactor() );

	fitter_corr -> Delete();
	
	return vec_corr;
}


TTree *tvalidation::InitializeOutputTree(){

  TTree *OutputTree = new TTree("EnsembleTree", "EnsembleTree");

  OutputTree -> Branch("F0",              &out_F0);
  OutputTree -> Branch("FL",              &out_FL);
  OutputTree -> Branch("FR",              &out_FR);
  OutputTree -> Branch("N0",              &out_N0);
  OutputTree -> Branch("NL",              &out_NL);
  OutputTree -> Branch("NR",              &out_NR);
  OutputTree -> Branch("Nges",            &out_Nges);
  OutputTree -> Branch("Wjets",           &out_Wjets);
  OutputTree -> Branch("QCD",             &out_QCD);
  OutputTree -> Branch("RemBkg",          &out_RemBkg);
  OutputTree -> Branch("F0_nom",          &out_F0_nom);
  OutputTree -> Branch("FL_nom",          &out_FL_nom);
  OutputTree -> Branch("FR_nom",          &out_FR_nom);
  OutputTree -> Branch("N0_nom",          &out_N0_nom);
  OutputTree -> Branch("NL_nom",          &out_NL_nom);
  OutputTree -> Branch("NR_nom",          &out_NR_nom);
  OutputTree -> Branch("Nges_nom",        &out_Nges_nom);
  OutputTree -> Branch("Wjets_nom",       &out_Wjets_nom);
  OutputTree -> Branch("QCD_nom",         &out_QCD_nom);
  OutputTree -> Branch("RemBkg_nom",      &out_RemBkg_nom);
  OutputTree -> Branch("F0_err",          &out_F0_err);
  OutputTree -> Branch("FL_err",          &out_FL_err);
  OutputTree -> Branch("FR_err",          &out_FR_err);
  OutputTree -> Branch("N0_err",          &out_N0_err);
  OutputTree -> Branch("NL_err",          &out_NL_err);
  OutputTree -> Branch("NR_err",          &out_NR_err);
  OutputTree -> Branch("Wjets_err",       &out_Wjets_err);
  OutputTree -> Branch("QCD_err",         &out_QCD_err);
  OutputTree -> Branch("RemBkg_err",      &out_RemBkg_err);
  OutputTree -> Branch("N0_prof_err",     &out_N0_prof_err);
  OutputTree -> Branch("NL_prof_err",     &out_NL_prof_err);
  OutputTree -> Branch("NR_prof_err",     &out_NR_prof_err);
  OutputTree -> Branch("Wjets_prof_err",  &out_Wjets_prof_err);
  OutputTree -> Branch("QCD_prof_err",    &out_QCD_prof_err);
  OutputTree -> Branch("RemBkg_prof_err", &out_RemBkg_prof_err);

  out_nui.clear();
  out_nui_nom.clear();
  out_nui_err.clear();
  out_nui_prof_err.clear();

  int fNui = fNNui;

  out_nui          = std::vector<double>(fNui);
  // out_nui_nom      = Nui_Vec;
  out_nui_nom      = std::vector<double>(fNui);
  out_nui_err      = std::vector<double>(fNui); 
  out_nui_prof_err = std::vector<double>(fNui);

  for(int k = 0; k < fNui; ++k){

    std::stringstream out;
    out << k;
    std::string nr = out.str();

    std::string nui          = "nui_" + nr;
    std::string nui_nom      = "nui_nom_" + nr;
    std::string nui_err      = "nui_err_" + nr;
    std::string nui_prof_err = "nui_prof_err_" + nr;

    OutputTree -> Branch(nui.c_str(),          &(out_nui[k]));
    OutputTree -> Branch(nui_nom.c_str(),      &(out_nui_nom[k]));
    OutputTree -> Branch(nui_err.c_str(),      &(out_nui_err[k]));
    OutputTree -> Branch(nui_prof_err.c_str(), &(out_nui_prof_err[k]));

  }

  OutputTree -> Branch("MinuitStatus",   &out_MinuitStatus);
  OutputTree -> Branch("NumNuisancePar", &out_NumNuisancePar);

  return OutputTree;

}
  
//Destructor
tvalidation::~tvalidation()
{

  delete fRandom;
  delete gaussRandom;
  delete fit;
  delete graph;
  delete fitter_prof; 

}






