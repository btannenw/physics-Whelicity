//---------------------------------------------------------------------
//
// Class used to evaluate external systematics and produce directly
// the final plots for the analysis,
// partly derived from ValidationClass
//
// 28.08.2012 by A.Knue
//---------------------------------------------------------------------


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

#include "AtlasStyle.h"

#include "ExternalSystematics.h"
#include "ValidationClass.h"
#include "helicity.h"
#include "plots.h"
#include "interpolation.h"
#include "ProfilingClass.h"
#include "StatusLogbook.h"
#include "PlotInterpolationCurves.h"

#include "fit.h"
#include "fit_validation.h" // contains the minimization function
#include "functions.h"

#include <fstream>

using namespace std;

std::vector<double> fBkgNorm2;
std::vector<double> fBkgUnc2;
std::vector<TH1D> fParameterHist2;
std::vector<TH1D> fParameterHistTS;
std::vector<double> SelEff2;

//Constructor 
ExternalSystematics::ExternalSystematics(std::string fittingMode)
{

  fFittingMode = fittingMode;
  
  fHistNomVec.clear();
  fHistNormVec.clear();
  fEffVec.clear();
  fBkgExpVec.clear();
  fBkgUncVec.clear();

  fNNui          = 0;
  fNScaleFactors = 0;
  
  //fRandom = new TRandom3(0);
  fRandom = new TRandom3(2101); // (24 Nov.2015: usig same random seed to make PD)
  gaussRandom = new TRandom3(0);

  fit   = new TF1();

  graph = new TGraph(50);

  fitter_prof = new TMinuit(50); // TMinuit(Int_t maxpar): maxpar is the maximum number of parameters used with this TMinuit object.
  //fitter_prof->SetPrintLevel(-1); // printlevel = -1  quiet (also suppresse all warnings)
  fitter_prof->SetPrintLevel(1); // printlevel = -1  quiet (also suppresse all warnings)

  //Set interpolation object for fit (at first only linear)                                                                                                    
  fitter_prof->SetFCN(eval_chisqrt_vali);

  fFractionNames.push_back("F0");
  fFractionNames.push_back("FL");
  if(fFittingMode == "3D")
    fFractionNames.push_back("FR");

  gErrorIgnoreLevel = kError;

}

void ExternalSystematics::clean()
{

}


void ExternalSystematics::FillTemplateVector(std::vector<TemplateInfo::MySignalTemplate> SigTempl, std::vector<TemplateInfo::MyBkgTemplate> BkgTempl)
{

  fParameterNames.clear();
  fParameterValues.clear();
  fParameterValues_nom.clear();

  for(int iParam = 0; iParam < fNSignal; ++iParam){

    if(iParam < fNSignal){
      fParameterHist.push_back(SigTempl[iParam].hist);
      fParameterHistTS.push_back(SigTempl[iParam].histTS);
      SelEff.push_back(SigTempl[iParam].SelEff);

      fParameterNames.push_back(SigTempl[iParam].ParamName);
      fParameterValues_nom.push_back(1.0);
      // fFractionValues_nom.push_back(1.0);

    }

  }
  for(int iParam = 0; iParam < fNBkg; ++iParam){

    fParameterHist.push_back(BkgTempl[iParam].hist);
    fParameterHistTS.push_back(BkgTempl[iParam].histTS);
    fBkgNorm.push_back(BkgTempl[iParam].Exp);
    fBkgUnc.push_back(BkgTempl[iParam].Unc);

    fParameterNames.push_back(BkgTempl[iParam].ParamName);
    fParameterValues_nom.push_back(BkgTempl[iParam].Exp);
    //  fParameterValues.push_back(-1.0);

  }


}

void ExternalSystematics::SetFitParameters(std::string BkgMode)
{

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


  fNScaleFactors = fHistNomVec.size();
    

}

//Get variation of pseudo-data
TH1D ExternalSystematics::RandomPseudodata(TH1D Histo)
{

  double lower_edge  = Histo.GetBinLowEdge(1);
  double bin_width   = Histo.GetBinWidth(1);
  double number_bins = Histo.GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;
  
  TH1D hist_ensemble = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

  //Loop over bins and fill them
  for(int ibin = 1; ibin <= number_bins; ibin++)
    {
      double nexp = Histo.GetBinContent(ibin);
      
      double nobs = fRandom->PoissonD(nexp);
      
      //      std::cout << nexp << "\t" << nobs << std::endl;

      //Set the bin content
      hist_ensemble.SetBinContent(ibin, nobs);
    }

  std::stringstream oss;
  oss << hist_ensemble.Integral();

    //WriteParameterStatus("ExternalSystematics", "Integral ensemble: "+oss.str());
  
  //Return the ensemble histogram
  return hist_ensemble;

}

//Get fluctuated templates
void ExternalSystematics::ProduceRandomTemplates()
{

  // fluctuate each template seperately (old)
  // test: fluctuate all templates at the same time

  int nParam = fParameterHistTS.size(); 

  for(int iParam = 0; iParam < nParam; ++iParam){
    
    TH1D help = fParameterHistTS[iParam];
  
    //    int nEntries = help.GetEntries();
  
    //    std::cout << iParam << "\t" << help.Integral() << "\t" << nEntries << std::endl;
  
    //    double norm  = help.Integral(); 
  
    //    help.Scale(nEntries/help.Integral());
  
    help = RandomPseudodata(help);
  
    //    std::cout << iParam << "\t" << help.Integral() << std::endl;
  
    fParameterHist[iParam] = help;
  
  }
  
}

///_______________________________________________________________________
///_____________________________VALIDATION________________________________
///_______________________________________________________________________
void ExternalSystematics::CallEvaluation(std::string OutputTxtFile, std::string mode)
{

  createHTML();

  //  InitializeHistograms(fSystematicPD.size());

  std::vector<std::string> labels;
  std::vector<TH1D> PD;

  int nSyst = fSystematicPD.size(); // = 7  // for DataFit fSystematicPD.size() = 1

  std::cout << "##### fSystematicPD.size= " << nSyst <<std::endl;

  if(mode != "Datafit"){

    for(int iSyst = 0; iSyst < nSyst; ++iSyst){

      Validation(fSystematicPD[iSyst].histoPD, iSyst, mode);
      
      PD.push_back(fSystematicPD[iSyst].histoPD);
      
      labels.push_back(fSystematicPD[iSyst].SystematicType);
      
    }
    
  }

  PlotInterpolationCurves *fPlot = new PlotInterpolationCurves();

  fPlot -> SetSubPlotLabels(fLeptonLabel, fJetBinLabel, fBTagLabel);

// new : added by MJK
  fPlot -> SetHtmlLable(_htmls[fHTMLLabel]);

  if(mode == "Syst" || mode == "PDF" || mode == "mass" ){


    if(nSyst == 3){
      fPlot -> DrawFinalPlot(PD, labels, fOutputFolder+"/HistComparison_"+fSystematicLabel+"_"+fChannel+"_PD", true, false);
      fPlot -> DrawFinalPlot(PD, labels, fOutputFolder+"/HistComparison_"+fSystematicLabel+"_"+fChannel+"_PD_normalised", true, true);
    }
    else if(nSyst == 2){
      fPlot -> DrawFinalPlotSingle(PD, labels, fOutputFolder+"/HistComparison_"+fSystematicLabel+"_"+fChannel+"_PD", true, false);
      fPlot -> DrawFinalPlotSingle(PD, labels, fOutputFolder+"/HistComparison_"+fSystematicLabel+"_"+fChannel+"_PD_normalised", true, true);
    }
    else
      std::cout << "ERROR::ExternalSystematics: We do not have a proper Plotting code for this number of variations!!!" << std::endl;
    
    delete fPlot;

    InitializeHistograms(fSystematicPD.size());

    // 25 Feb. 2015 bt: split evaluation by PDF 
    if(mode == "Syst"){
      //EvaluateSystematic(fEvaluationMode, OutputTxtFile);
      EvaluateSystematic(fEvaluationMode, OutputTxtFile, labels); // 7 Dec. 2015 new by mjk: passing labels for txt file
    }
    else if (mode == "PDF"){
      EvaluateSystematic("PDF", OutputTxtFile, labels);
    }
    else if (mode == "mass"){
      EvaluateSystematic("mass", OutputTxtFile, labels);
    }

    PlotSystematicDistributions(fHistoF0,     labels, "F0");
    PlotSystematicDistributions(fHistoFL,     labels, "FL");
    if(fFittingMode == "3D")
      PlotSystematicDistributions(fHistoFR,     labels, "FR");
    
    int nParam = fParameterNames.size();
    
    for(int iParam = 0; iParam < nParam; ++iParam){
      
      PlotSystematicDistributions(fHistoAll[iParam], labels, fParameterNames[iParam].c_str());
      
    }
    
  }
  else if(mode == "Cali"){
    
    std::cout << "InitializeAnalysis" << std::endl;
    
    InitializeGraphs();  
    
    
    //    for(int i = 0; i < fSystList.size(); ++i){
    
    for(int i = 0; i < fParameterNames.size(); ++i){
      
      for(int k = 0; k < fSystList.size(); ++k){
	AnalyseParameter(i, k, fOutputFolder);
      }
      
      ProduceAllPlots(i);
      
    }
    
    for(int i = 0; i < fSystList.size(); ++i)
      AnalyseFractions(i, "Cali");
    
    ProduceAllFractionPlots();
  }
  
  else if(mode == "Datafit"){
    
    fDataHist = fSystematicPD[0].histoPD; // mjk: it is the data histo in the template file
    
    Validation(fSystematicPD[0].histoPD, 0, mode);
    
    AnalyseDatafit();
    
  }
  
}

void ExternalSystematics::ProduceAllPlots(int iParam)
{

  int nParam = fParameterNames.size();

  std::string xtitle   = "";
  std::string ytitle   = "";
  std::string ytitle2  = "";
  std::string filename = "";
  double xmin = 0.0; // used in (plots.cxx) as: ==> TF1 *fit = new TF1("fit", "[0]*x + [1]",xmin,xmax);
  double xmax = 0.0;
  TGraphErrors fNewGraph,fNewGraphDiff,fNewGraphPull,fNewGraphRMS;
  std::string F0Label           = "";
  std::string DistLabel         = "";
  std::string SigmaLabel        = "";
  std::string PullLabel         = "";
  std::string DistFilename      = "";
  std::string SigmaFilename     = "";
  std::string PullFilename      = "";
  std::string PullDistFilename  = "";
  std::string RMSLabel          = "";
  std::string RMSFilename       = "";
  std::string RMSDistFilename   = "";


  std::string ParamName = fParameterNames[iParam];

  if(ParamName == "N0"){
    
    xtitle        = "N_{0} - N_{0,nom}";  filename = "Calicurve_N0_"+fChannel; 

    if(fFittingMode == "3D"){

      xmin = fXsec*fLumi*(-0.365); xmax = fXsec*fLumi*0.365;

    }
    else{

      xmin = fXsec*fLumi*(-0.425); xmax = fXsec*fLumi*0.425;
      
    }

    ytitle2       = "N_{0} - N_{0,nom}";
    ytitle        = "N_{0}";
    F0Label       = "N_{0,input}=";
    DistFilename  = "Distr_N0_F0_";
    SigmaFilename = "Distr_N0_err_F0_";
    PullFilename  = "Pull_Dependence_N0_"+fChannel;
    PullDistFilename  = "Pull_Dist_N0_"+fChannel;
    DistLabel     = "N_{0}";
    SigmaLabel    = "#sigma(N_{0})";
    PullLabel     = "Pull distribution for N_{0}";
    RMSFilename     = "RMS_Dependence_N0_"+fChannel;
    RMSDistFilename = "RMS_Dist_N0_"+fChannel;
    RMSLabel        = "RMS distribution for N_{0}";

        
  }
  else if(ParamName == "NL"){
    
    xtitle        = "N_{L} - N_{L,nom}";  
    
    if(fFittingMode == "3D"){

      filename = "Calicurve_NL_"+fChannel; xmin = fXsec*fLumi*(-0.20); xmax = fXsec*fLumi*0.2;

    }
    else if(fFittingMode == "2D"){

      filename = "Calicurve_NL_"+fChannel; xmin = fXsec*fLumi*(-0.25); xmax = fXsec*fLumi*0.25;

    }    
    
    ytitle2       = "N_{L} - N_{L,nom}";
    ytitle        = "N_{L}";
    F0Label       = "N_{L,input}=";
    DistFilename  = "Distr_NL_FL_";
    SigmaFilename = "Distr_NL_err_FL_";
    PullFilename  = "Pull_Dependence_NL_"+fChannel;
    PullDistFilename  = "Pull_Dist_NL_"+fChannel;
    DistLabel     = "N_{L}";
    SigmaLabel    = "#sigma(N_{L})";
    PullLabel     = "Pull distribution for N_{L}";
    RMSFilename     = "RMS_Dependence_NL_"+fChannel;
    RMSDistFilename = "RMS_Dist_NL_"+fChannel;
    RMSLabel        = "RMS distribution for N_{L}";
 
  }
  else if(ParamName == "NR"){
    
    xtitle        = "N_{R} - N_{R,nom}";  filename = "Calicurve_NR_"+fChannel; xmin = fXsec*fLumi*(-0.21); xmax = fXsec*fLumi*0.21;
    ytitle2       = "N_{R} - N_{R,nom}";
    ytitle        = "N_{R}";
    F0Label       = "N_{R,input}=";
    DistFilename  = "Distr_NR_FR_";
    SigmaFilename = "Distr_NR_err_FR_";
    PullFilename  = "Pull_Dependence_NR_"+fChannel;
    PullDistFilename  = "Pull_Dist_NR_"+fChannel;
    DistLabel     = "N_{R}";
    SigmaLabel    = "#sigma(N_{R})";
    PullLabel     = "Pull distribution for N_{R}";
    RMSFilename     = "RMS_Dependence_NR_"+fChannel;
    RMSDistFilename = "RMS_Dist_NR_"+fChannel;
    RMSLabel        = "RMS distribution for N_{R}";

    
  }
  else if(ParamName == "Wjets"){
      
      std::stringstream oss;
      oss << fBkgExpVec[iParam-fFractionNames.size()];

      xtitle        = "F_{0}";  filename = "Calicurve_Wjets_"+fChannel; xmin = 0.35; xmax = 1.05;
      ytitle        = "N_{W+jets}-N_{W+jets,nom}";
      F0Label       = "N_{W+jets,input}="+oss.str();
      DistFilename  = "Distr_Wjets_F0_";
      SigmaFilename = "Distr_Wjets_err_F0_";
      PullFilename  = "Pull_Dependence_Wjets_"+fChannel;
      PullDistFilename  = "Pull_Dist_Wjets_"+fChannel;
      DistLabel     = "N_{Wjets}";
      SigmaLabel    = "#sigma(N_{W+jets})";
      PullLabel     = "Pull distribution for N_{W+jets}";
      RMSFilename     = "RMS_Dependence_Wjets_"+fChannel;
      RMSDistFilename = "RMS_Dist_Wjets_"+fChannel;
      RMSLabel        = "RMS distribution for Wjets";

      if(fFittingMode == "3D"){

	xmin = 0.35; xmax = 1.075;

      }
      else{

	xmin = 0.425; xmax = 0.975;

      }

  }
  else if(ParamName == "QCD"){

      std::stringstream oss;
      oss << fBkgExpVec[iParam-fFractionNames.size()];

      xtitle        = "F_{0}";  filename = "Calicurve_QCD_"+fChannel; xmin = 0.35; xmax = 1.05;
      ytitle        = "N_{QCD} - N_{QCD,nom}";
      F0Label       = "N_{QCD,input}="+oss.str();
      DistFilename  = "Distr_QCD_F0_";
      SigmaFilename = "Distr_QCD_err_F0_";
      PullFilename  = "Pull_Dependence_QCD_"+fChannel;
      PullDistFilename  = "Pull_Dist_QCD_"+fChannel;
      DistLabel     = "N_{QCD}";
      SigmaLabel    = "#sigma(N_{QCD})";
      PullLabel     = "Pull distribution for N_{QCD}";

      RMSFilename     = "RMS_Dependence_QCD_"+fChannel;
      RMSDistFilename = "RMS_Dist_QCD_"+fChannel;
      RMSLabel        = "RMS distribution for QCD";


      if(fFittingMode == "3D"){

	xmin = 0.35;  xmax = 1.075;

      }
      else{

	xmin = 0.425; xmax = 0.975;

      }
      
  }
  else if(ParamName == "RemBkg"){

      std::stringstream oss;
      oss << fBkgExpVec[iParam-fFractionNames.size()];

      xtitle        = "F_{0}";  filename = "Calicurve_RemBkg_"+fChannel; xmin = 0.35; xmax = 1.05;
      ytitle        = "N_{RemBkg} - N_{RemBkg,nom}";
      F0Label       = "N_{RemBkg,input}="+oss.str();
      DistFilename  = "Distr_RemBkg_F0_";
      SigmaFilename = "Distr_RemBkg_err_F0_";
      PullFilename  = "Pull_Dependence_RemBkg_"+fChannel;
      PullDistFilename  = "Pull_Dist_RemBkg_"+fChannel;
      DistLabel     = "N_{RemBkg}";
      SigmaLabel    = "#sigma(N_{RemBkg})";
      PullLabel     = "Pull distribution for N_{RemBkg}";

      RMSFilename     = "RMS_Dependence_RemBkg_"+fChannel;
      RMSDistFilename = "RMS_Dist_RemBkg_"+fChannel;
      RMSLabel        = "RMS distribution for RemBkg";

      if(fFittingMode == "3D"){

	xmin = 0.35; xmax = 1.075;

      }
      else{

	xmin = 0.425; xmax = 0.975;

      }
      
  }
  else if(ParamName == "Wjets_4incl_1incl_el"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
    xtitle        = "F_{0}";  filename = "Calicurve_Wjets_el_"+fChannel; xmin = 0.35; xmax = 1.05;
    ytitle        = "#Delta N_{W+jets,e+jets}";
    F0Label       = "N_{W+jets,e+jets,input}="+oss.str();
    DistFilename  = "Distr_Wjets_el_F0_";
    SigmaFilename = "Distr_Wjets_err_el_F0_";
    PullFilename  = "Pull_Dependence_Wjets_el_"+fChannel;
    PullDistFilename  = "Pull_Dist_Wjets_el_"+fChannel;
    DistLabel     = "N_{Wjets,e+jets}";
    SigmaLabel    = "#sigma(N_{W+jets,e+jets})";
    PullLabel     = "Pull distribution for N_{W+jets,e+jets}";

    RMSFilename     = "RMS_Dependence_Wjets_e_"+fChannel;
    RMSDistFilename = "RMS_Dist_Wjets_e_"+fChannel;
    RMSLabel        = "RMS distribution for Wjets e";


    if(fFittingMode == "3D"){

      xmin = 0.35; xmax = 1.075;

    }
    else{

      xmin = 0.425; xmax = 0.975;

    }

    
  }
  else if(ParamName == "QCD_4incl_1incl_el"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
    xtitle        = "F_{0}";  filename = "Calicurve_QCD_el_"+fChannel; xmin = 0.35; xmax = 1.05;
    ytitle        = "#Delta N_{QCD,e+jets}";
    F0Label       = "N_{QCD,e+jets,input}="+oss.str();
    DistFilename  = "Distr_QCD_el_F0_";
    SigmaFilename = "Distr_QCD_err_el_F0_";
    PullFilename  = "Pull_Dependence_QCD_el_"+fChannel;
    PullDistFilename  = "Pull_Dist_QCD_el_"+fChannel;
    DistLabel     = "N_{QCD,e+jets}";
    SigmaLabel    = "#sigma(N_{QCD,e+jets})";
    PullLabel     = "Pull distribution for N_{QCD,e+jets}";
    RMSFilename     = "RMS_Dependence_QCD_e_"+fChannel;
    RMSDistFilename = "RMS_Dist_QCD_e_"+fChannel;
    RMSLabel        = "RMS distribution for QCD e";

    if(fFittingMode == "3D"){

      xmin = 0.35; xmax = 1.075;

    }
    else{

      xmin = 0.425; xmax = 0.975;

    }
    
  }
  else if(ParamName == "RemBkg_4incl_1incl_el"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
      xtitle        = "F_{0}";  filename = "Calicurve_RemBkg_el_"+fChannel; xmin = 0.35; xmax = 1.05;
      ytitle        = "#Delta N_{RemBkg,e+jets}";
      F0Label       = "N_{RemBkg,e+jets,input}="+oss.str();
      DistFilename  = "Distr_RemBkg_el_F0_";
      SigmaFilename = "Distr_RemBkg_err_el_F0_";
      PullFilename  = "Pull_Dependence_RemBkg_el_"+fChannel;
      PullDistFilename  = "Pull_Dist_RemBkg_el_"+fChannel;
      DistLabel     = "N_{RemBkg,e+jets}";
      SigmaLabel    = "#sigma(N_{RemBkg,e+jets})";
      PullLabel     = "Pull distribution for N_{RemBkg,e+jets}";

      RMSFilename     = "RMS_Dependence_RB_e_"+fChannel;
      RMSDistFilename = "RMS_Dist_RB_e_"+fChannel;
      RMSLabel        = "RMS distribution for RB e";


      if(fFittingMode == "3D"){

        xmin = 0.35; xmax = 1.075;

      }
      else{

        xmin = 0.425; xmax = 0.975;

      }
      
  }
  else if(ParamName == "Wjets_4incl_1incl_mu"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
    xtitle        = "F_{0}";  filename = "Calicurve_Wjets_mu_"+fChannel; xmin = 0.35; xmax = 1.05;
    ytitle        = "#Delta N_{W+jets,#mu+jets}";
    F0Label       = "N_{W+jets,#mu+jets,input}="+oss.str();
    DistFilename  = "Distr_Wjets_mu_F0_";
    SigmaFilename = "Distr_Wjets_err_mu_F0_";
    PullFilename  = "Pull_Dependence_Wjets_mu_"+fChannel;
    PullDistFilename  = "Pull_Dist_Wjets_mu_"+fChannel;
    DistLabel     = "N_{Wjets,#mu+jets}";
    SigmaLabel    = "#sigma(N_{W+jets,#mu+jets})";
    PullLabel     = "Pull distribution for N_{W+jets,#mu+jets}";

    RMSFilename     = "RMS_Dependence_Wjets_mu_"+fChannel;
    RMSDistFilename = "RMS_Dist_Wjets_mu_"+fChannel;
    RMSLabel        = "RMS distribution for Wjets mu";


    if(fFittingMode == "3D"){

      xmin = 0.35; xmax = 1.075;

    }
    else{

      xmin = 0.425; xmax = 0.975;

    }
    
  }
  else if(ParamName == "QCD_4incl_1incl_mu"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
    xtitle        = "F_{0}";  filename = "Calicurve_QCD_mu_"+fChannel; xmin = 0.35; xmax = 1.05;
    ytitle        = "#Delta N_{QCD,#mu+jets}";
    F0Label       = "N_{QCD,#mu+jets,input}="+oss.str();
    DistFilename  = "Distr_QCD_mu_F0_";
    SigmaFilename = "Distr_QCD_err_mu_F0_";
    PullFilename  = "Pull_Dependence_QCD_mu_"+fChannel;
    PullDistFilename  = "Pull_Dist_QCD_mu_"+fChannel;
    DistLabel     = "N_{QCD,#mu+jets}";
    SigmaLabel    = "#sigma(N_{QCD,#mu+jets})";
    PullLabel     = "Pull distribution for N_{QCD,#mu+jets}";
    
    RMSFilename     = "RMS_Dependence_QCD_m_"+fChannel;
    RMSDistFilename = "RMS_Dist_QCD_m_"+fChannel;
    RMSLabel        = "RMS distribution for QCD m";


    if(fFittingMode == "3D"){

      xmin = 0.35; xmax = 1.075;

    }
    else{

      xmin = 0.425; xmax = 0.975;

    }
    
  }
  else if(ParamName == "RemBkg_4incl_1incl_mu"){
    
    std::stringstream oss;
    oss << fBkgExpVec[iParam-fFractionNames.size()];
    
      xtitle        = "F_{0}";  filename = "Calicurve_RemBkg_mu_"+fChannel; xmin = 0.35; xmax = 1.05;
      ytitle        = "#Delta N_{RemBkg,#mu+jets}";
      F0Label       = "N_{RemBkg,#mu+jets,input}="+oss.str();
      DistFilename  = "Distr_RemBkg_mu_F0_";
      SigmaFilename = "Distr_RemBkg_err_mu_F0_";
      PullFilename  = "Pull_Dependence_RemBkg_mu_"+fChannel;
      PullDistFilename  = "Pull_Dist_RemBkg_mu_"+fChannel;
      DistLabel     = "N_{RemBkg,#mu+jets}";
      SigmaLabel    = "#sigma(N_{RemBkg,#mu+jets})";
      PullLabel     = "Pull distribution for N_{RemBkg,#mu+jets}";

      RMSFilename     = "RMS_Dependence_RB_mu_"+fChannel;
      RMSDistFilename = "RMS_Dist_RB_mu_"+fChannel;
      RMSLabel        = "RMS distribution for RB mu";


      if(fFittingMode == "3D"){

        xmin = 0.35; xmax = 1.075;

      }
      else{

        xmin = 0.425; xmax = 0.975;

      }
            
  }

  fNewGraph = fGraphs[iParam][0]; fNewGraphDiff = fGraphs[iParam][1]; fNewGraphPull = fGraphs[iParam][2]; fNewGraphRMS = fGraphs[iParam][3];

  int nFrac = fFractionNames.size();

  std::cout << "Calicurve for : " << ParamName.c_str() << std::endl;
  
  if(iParam < nFrac)
    PlotCalDistribution_N(&fNewGraph,  &fNewGraphDiff, xtitle, ytitle2, filename,  xmin, xmax, fChannel, fOutputFolder, "");
  else
    PlotCalDistribution_N(&fNewGraph,  &fNewGraphDiff, xtitle, ytitle, filename,  xmin, xmax, fChannel, fOutputFolder, F0Label);

  std::cout << "\t" << std::endl;
  
  std::cout << "Pullcurve for : " << ParamName.c_str() << std::endl;

  PlotCalPullDistribution(&fNewGraphPull,    xtitle, ytitle, PullFilename, xmin, xmax, -0.099, 0.099, fChannel, "", fOutputFolder);

  std::cout << "\t" << std::endl;

  std::cout << "RMS-curve for : " << ParamName.c_str() << std::endl;

  //void PlotCalRMSDistribution (TGraphErrors* graph, std::string xlabel, std::string ylabel, std::string xtitle, double xmin, double xmax, double ymin, double ymax, std::string input_channel, std::string smode, std::string output_folder)
  PlotCalRMSDistribution(&fNewGraphRMS,      xtitle, ytitle, RMSFilename,  xmin, xmax, 0.5,      1.1, fChannel, "", fOutputFolder);

  std::cout << "\t" << std::endl;


  for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){
    
    std::stringstream oss;
    
    if(iParam == 1)
      oss << fFraction[iSyst][1]*fXsec*fLumi; // value for FL;
    else if(iParam == 2 && fFittingMode == "3D")
      oss << fFraction[iSyst][2]*fXsec*fLumi; // value for FR;
    else
      oss << fFraction[iSyst][0]*fXsec*fLumi; // value for F0;
    
    
    if(iParam < nFrac){
      PlotDistribution_N(&fHistoAnalysis[iSyst][0], 0.4, fXsec, fLumi, DistLabel,  DistFilename+"_"+oss.str()+"_"+fChannel,      fChannel, F0Label+oss.str(), 0, fOutputFolder, "");
      PlotDistribution_N(&fHistoAnalysis[iSyst][1], 0.4, fXsec, fLumi, SigmaLabel, SigmaFilename+"_"+oss.str()+"_"+fChannel,     fChannel, F0Label+oss.str(), 0, fOutputFolder, "");
      PlotDistribution_N(&fHistoAnalysis[iSyst][2], 0.4, fXsec, fLumi, PullLabel,  PullDistFilename+"_"+oss.str()+"_"+fChannel,  fChannel, F0Label+oss.str(), 0, fOutputFolder, "Pull");
          
    }
    else{
      PlotDistribution_N(&fHistoAnalysis[iSyst][0], 0.4, fXsec, fLumi, DistLabel,  DistFilename+"_"+oss.str()+"_"+fChannel,      fChannel, F0Label, 0, fOutputFolder, "");
      PlotDistribution_N(&fHistoAnalysis[iSyst][1], 0.4, fXsec, fLumi, SigmaLabel, SigmaFilename+"_"+oss.str()+"_"+fChannel,     fChannel, F0Label, 0, fOutputFolder, "");
      PlotDistribution_N(&fHistoAnalysis[iSyst][2], 0.4, fXsec, fLumi, PullLabel,  PullDistFilename+"_"+oss.str()+"_"+fChannel,  fChannel, F0Label, 0, fOutputFolder, "Pull");
          
    }
    
    
  }
  
  fHistoAnalysis.clear();
}

void ExternalSystematics::ProduceAllFractionPlots()
{

    std::cout << "F0" << std::endl;

  // F0 
  if(fFittingMode == "3D")
    PlotCalDistribution_N(&fFractionGraphs[0][0],   &fFractionGraphs[0][1], "F_{0}", "F_{0}", "Calicurve_F0_"+fChannel,  0.35, 1.05, fChannel, fOutputFolder, "");
  else
    PlotCalDistribution_N(&fFractionGraphs[0][0],   &fFractionGraphs[0][1], "F_{0}", "F_{0}", "Calicurve_F0_"+fChannel,  0.425, 0.975, fChannel, fOutputFolder, "");

  
  //  PlotCalPullDistribution(&fFractionGraphs[0][2],                         "F_{0}",  "(F_{0})", "Pull_Dependence_F0_"+fChannel, 0.35, 1.05, -0.1, 0.1, fChannel, "", fOutputFolder);
  
  for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){

    std::stringstream oss;
    oss << fFraction[iSyst][0]; // value for F0;
 
    
    PlotDistribution_N(&fFractionAnalysis[iSyst][0][0], 0.4, fXsec, fLumi, "F_{0}", "Distr_F0_F0_"+oss.str()+"_"+fChannel, fChannel, "F_{0,input}="+oss.str(), 0, fOutputFolder, "Pull\
");
    
    PlotDistribution_N(&fFractionAnalysis[iSyst][0][1], 0.4, fXsec, fLumi, "#sigma(F_{0})", "Distr_F0_err_F0_"+oss.str()+"_"+fChannel,  fChannel, "F_{0,input}="+oss.str(), 0, fOutputFolder, "Pull\
");
    
    PlotDistribution_N(&fFractionAnalysis[iSyst][0][2], 0.4, fXsec, fLumi, "Pull distribution for F_{0}", "Pull_F0_F0_"+oss.str()+"_"+fChannel, fChannel, "F_{0,input}="+oss.str(), 0, fOutputFolder, "Pull\
");
    

  }

    std::cout << "FL" << std::endl;

  // FL

  if(fFittingMode == "3D"){

    PlotCalDistribution_N(&fFractionGraphs[1][0],   &fFractionGraphs[1][1], "F_{L}", "F_{L}", "Calicurve_FL_"+fChannel,  0.10, 0.5, fChannel, fOutputFolder, "");
    //    PlotCalPullDistribution(&fFractionGraphs[1][2],                         "F_{L}",  "(F_{L})", "Pull_Dependence_FL_"+fChannel, 0.10, 0.5, -0.1, 0.1, fChannel, "", fOutputFolder);
 
  }
  else if(fFittingMode == "2D"){
    
    PlotCalDistribution_N(&fFractionGraphs[1][0],   &fFractionGraphs[1][1], "F_{L}", "F_{L}", "Calicurve_FL_"+fChannel,          0.025, 0.575, fChannel, fOutputFolder, "");
    //  PlotCalPullDistribution(&fFractionGraphs[1][2],                         "F_{L}",  "(F_{L})", "Pull_Dependence_FL_"+fChannel, -0.05, 0.65, -0.1, 0.1, fChannel, "", fOutputFolder);
    
  }


  for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){

    std::stringstream oss;
    oss << fFraction[iSyst][1]; // value for FL;
   
    PlotDistribution_N(&fFractionAnalysis[iSyst][1][0], 0.4, fXsec, fLumi, "F_{L}", "Distr_FL_FL_"+oss.str()+"_"+fChannel,                      fChannel, "F_{L,input}="+oss.str(), 0, fOutputFolder, "Pull\
");
    PlotDistribution_N(&fFractionAnalysis[iSyst][1][1], 0.4, fXsec, fLumi, "#sigma(F_{L})", "Distr_FL_err_FL_"+oss.str()+"_"+fChannel,          fChannel, "F_{L,input}="+oss.str(), 0, fOutputFolder, "Pull\
");
    PlotDistribution_N(&fFractionAnalysis[iSyst][1][2], 0.4, fXsec, fLumi, "Pull distribution for F_{L}", "Pull_FL_FL_"+oss.str()+"_"+fChannel, fChannel, "F_{L,input}="+oss.str(), 0, fOutputFolder, "Pull");
    
  }

  if(fFittingMode == "3D"){
  
        std::cout << "FR" << std::endl;

    // FR
    PlotCalDistribution_N(&fFractionGraphs[2][0],   &fFractionGraphs[2][1], "F_{R}", "F_{R}", "Calicurve_FR_"+fChannel,  -0.2, 0.2, fChannel, fOutputFolder, "");
    //    PlotCalPullDistribution(&fFractionGraphs[2][2],                         "F_{R}",  "(F_{R})", "Pull_Dependence_FR_"+fChannel, -0.2, 0.2, -0.1, 0.1, fChannel, "", fOutputFolder);
    
    for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){
      
      std::stringstream oss;
      oss << fFraction[iSyst][2]; // value for FR;
      
      PlotDistribution_N(&fFractionAnalysis[iSyst][2][0], 0.4, fXsec, fLumi, "F_{R}", "Distr_FR_FR_"+oss.str()+"_"+fChannel,                      fChannel, "F_{R,input}="+oss.str(), 0, fOutputFolder, "Pull");
      
      //    std::cout << "Mean error " << fFractionAnalysis[iSyst][2][1].GetMean() << std::endl;
      
      PlotDistribution_N(&fFractionAnalysis[iSyst][2][1], 0.4, fXsec, fLumi, "#sigma(F_{R})", "Distr_FRerr_FR_"+oss.str()+"_"+fChannel,           fChannel, "F_{R,input}="+oss.str(), 0, fOutputFolder, "Pull");
      
      PlotDistribution_N(&fFractionAnalysis[iSyst][2][2], 0.4, fXsec, fLumi, "Pull distribution for F_{R}", "Pull_FR_FR_"+oss.str()+"_"+fChannel, fChannel, "F_{R,input}="+oss.str(), 0, fOutputFolder, "Pull");
      
    }
    
  }

}

//Parameter distributions based on pseudo-data: nui_var is the index of the parameter that is varied, Nui_Vec contains all values for the nui param that are used to generate the pseudo data
void ExternalSystematics::Validation(TH1D histPD, int nSyst, std::string mode)
{

      //Set start values to SM values
  //      double F0_true = 0.698;
  //    double FL_true = 0.301;
  //    double FR_true = 0.00041;
    
      double F0_true = fFraction[nSyst][0];
      double FL_true = fFraction[nSyst][1];
      double FR_true = 0.0;
      if(fFittingMode == "3D")
	           FR_true = fFraction[nSyst][2];
      
      // if(mode == "Datafit"){

	     //   F0_true = 0.698;
	     //   FL_true = 0.301;
      //    FR_true = 0.00041;
      // }

      int k_counter    = 0;
      int conv_counter = 0;
      
      ///_____________________________FIT_______________________________
      
      WriteInfoStatus("ExternalSystematics", "Initialize TMinuit");
      
      //Setting up the fitting procedure; initialize TMinuit
      TMinuit *fitter = new TMinuit(50); // TMinuit(Int_t maxpar): maxpar is the maximum number of parameters used with this TMinuit object
      
      //      hist_ensemble = 0;	
      
      int counter_good = 0;

      if(mode=="Datafit")
      {
        std::cout<<"$$$$$ (1):nSyst"<<nSyst <<" !!!"<<endl;
        std::cout<<"$$$$$ (2):fSystematicPD.size()"<<fSystematicPD.size() <<" !!!"<<endl;
        std::cout<<"$$$$$ (3):fSystematicPD[0].SystematicType"<<fSystematicPD[0].SystematicType <<" !!!"<<endl;
      }
      // if(mode=="DatafitPD")
      // {
      //   std::cout<<"$$$$$ (1):nSyst"<<nSyst <<" !!!"<<endl;
      //   std::cout<<"$$$$$ (2):fSystematicPD.size()"<<fSystematicPD.size() <<" !!!"<<endl;
      //   std::cout<<"$$$$$ (3):fSystematicPD[0].SystematicType"<<fSystematicPD[0].SystematicType <<" !!!"<<endl;
      // }
        
      std::string outputFile;

      if(mode != "DatafitPD")
        outputFile = fOutputFolder+"/Syst_"+fSystematicPD[nSyst].SystematicType+"_"+fChannel+".root";
      

      else
	       outputFile = fOutputFolder+"/Syst_DatafitPD_"+fChannel+".root";
      
      //std::cout<<"$$$$$ (1):outputFile= "<<outputFile<<" !!!"<<endl;
      fSystList.push_back(outputFile);
      //std::cout<<"$$$$$ (2):fSystList.size()= "<<fSystList.size()<<" !!!"<<endl;

      //std::cout<<"%%%%%%%%%%%%%% 0000 %%%%%%%%%%%%%%%"<<std::endl;
      //if(fChannel=="Nominal")
      //cout<<"!!!!!!!!!!!! fChannel="<<fChannel<<" !!!!!!!!!!!!!"<<endl;
      
      if(gSystem->AccessPathName(outputFile.c_str())){
	//std::cout<<"%%%%%%%%%%%%%% 11111 %%%%%%%%%%%%%%%"<<std::endl;
	std::cout<<"Creating file "<<outputFile<<" !!!"<<endl;
	TFile *fNewFile = new TFile(outputFile.c_str(), "RECREATE");
	fNewFile->cd();
	}
      else{
	std::cout<<"Output file "<<outputFile<<" already exists!!! creating temp file"<<endl;
	TFile *fNewFile = new TFile("temp.root", "RECREATE");
	fNewFile->cd();
	}
    std::cout<<"$$$ initializeOutputTree, mode="<<mode<<std::endl;
      TTree *fNewTree = InitializeOutputTree();

      // std::cout << fParameterValues.size() << std::endl;

      // std::cout << "Integral of Pseudo Data = " << histPD.Integral() << std::endl;

      if(fSystematicLabel == "TemplateStat")
	       fParameterHist2 = fParameterHist;

      int NumberPseudoExp = fNPseudoExp;

      if(mode == "Datafit") NumberPseudoExp = 1;

      //TFile* f_ensemblesIn;
      //TFile* f_ensemblesOut;

      // std::string make_use_PD = "make"; // use , make
      // std::string ensemblesDir = "/afs/cern.ch/work/m/mkareem/FitPackage/Outputs/15bins/ensembles";

      // if(make_use_PD=="make"){
      //   f_ensemblesOut = new TFile((ensemblesDir+ "/hist_ensembles_"+fSystematicPD[nSyst].SystematicType+".root").c_str(),"RECREATE");
      //   f_ensemblesOut->cd();  
      // }

      // else
      //   f_ensemblesIn = new TFile((ensemblesDir+ "/hist_ensembles_"+fSystematicPD[nSyst].SystematicType+".root").c_str(), "READ");

      while(counter_good < NumberPseudoExp)
	{ 
	  
	  // clear the output tree
	  clean();
	  
	  std::stringstream oss;
	  oss << counter_good;

	  if(counter_good%100 == 0 && mode != "Datafit")
	    WriteInfoStatus("ExternalSystematics", "Get new ensemble for pseudo experiments "+oss.str());
	  //std::cout<<"%%%%%%%%%%%%%% 22222 %%%%%%%%%%%%%%%"<<std::endl;

	  // std::cout << "Hier 1a" << std::endl;

	  if(fSystematicLabel == "TemplateStat"){

	    hist_ensemble = histPD; // defined as extern variable in fit_validation.h
	    ProduceRandomTemplates();
	    
	  }
	  else{
      if(mode != "Datafit"){

/*        if(make_use_PD == "make"){
            hist_ensemble = RandomPseudodata(histPD); // PoissonD fluctuated histPD histogram
            hist_ensemble.Write( ("hist_ensemble_"+ std::to_string(counter_good)).c_str() );
        }
          else
            hist_ensemble = *(TH1D*)f_ensemblesIn->Get( ("hist_ensemble_"+ std::to_string(counter_good)).c_str() );
            */
          hist_ensemble = RandomPseudodata(histPD); // PoissonD fluctuated histPD histogram
      }
        
      
      else hist_ensemble = histPD;
    }

    //hist_ensemble.Write( ("hist_ensemble_"+ std::to_string(counter_good)).c_str() );
    
	  //if(mode == "Datafit")
	  //  hist_ensemble = histPD; 

	  // std::cout << "Hier 1e" << std::endl;

	  //Setting up the fitting procedure; initialize TMinuit
	  if(fitter) fitter->Delete();
	  fitter = new TMinuit(50);
	  //fitter->SetPrintLevel(-1); // quiet (also suppresse all warnings)
	  fitter->SetPrintLevel(1); // quiet (also suppresse all warnings)
	  
	  // std::cout << "Hier 1f" << std::endl;

	  int ierflg_err=0, ierflg = 0, ierflg_hesse = 0, par_num;
	  double start_val, start_step=0.001, A=0.0, B=0.0;  
	  
	  double arglist[50];
	  
	  // std::cout << "Hier 1g" << std::endl;

	  //Set interpolation object for fit (at first only linear)
	  fitter->SetFCN(eval_chisqrt_vali);	//To set the address of the minimization function
	  
	  // std::cout << "Hier 1h" << std::endl;

	  //Generate new start values
	  std::vector <double> start_val_gauss;
	  
	  double Ntotal = fXsec*fLumi;

    //std::cout<<"Ntotal= " << Ntotal << std::endl;

	  //	  std::cout << Ntotal << "\t" << F0_true << std::endl;
	  
	  // has to be made more general!!!
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*F0_true, 50000.0)); // Gaus(Double_t mean, Double_t sigma)
	  start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FL_true, 50000.0));
	  if(fFittingMode == "3D")
	    start_val_gauss.push_back( gaussRandom->Gaus(Ntotal*FR_true, 50000.0));

    //-------------------------
	  
	  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){
	    start_val_gauss.push_back( gaussRandom->Gaus(fBkgNorm[iBkg], 1.0*fBkgUnc[iBkg]) );
	  
	    //	    std::cout << fBkgNorm[iBkg] << "\t" << fBkgUnc[iBkg] << std::endl;

	  }

	  //Define 6 scale parameters alpha for signal and background
	  for (int k = 0; k != fNScaleFactors; k++){
	    
	    //	    // std::cout << k << std::endl;

	    std::string s;
	    std::stringstream out;
	    out << k;
	    s = out.str();


	    // MJK TEST 18 Feb. 2016 (setting limits for A and B for bkg norm)
	    /*
      if(k<fSignalParameter.size())
	      fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, A, B, ierflg); //mjk: https://root.cern.ch/root/html/TMinuit.html#TMinuit:mnparm
	    else {
	      //limiting bkg norm to be within its uncertainty
	      //fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, fBkgNorm[k-fBkgParameter.size()]-fBkgUnc[k-fBkgParameter.size()], fBkgNorm[k-fBkgParameter.size()]+fBkgUnc[k-fBkgParameter.size()], ierflg); 
	      fitter->mnparm(par_num=k, "alpha_" + s, fBkgNorm[k-fSignalParameter.size()] , start_step, fBkgNorm[k-fSignalParameter.size()]-fBkgUnc[k-fSignalParameter.size()], fBkgNorm[k-fSignalParameter.size()]+fBkgUnc[k-fSignalParameter.size()], ierflg); 
	      //std::cout<<"$$$$ --- fBkgNorm[" << k << "]= "<<fBkgNorm[k-fSignalParameter.size()] << "\tfBkgUnc= "<< fBkgUnc[k-fSignalParameter.size()] << std::endl;
	    }
      */
	    
	    // ORIGINAL --- Feb 18, 2016
	    //	    // std::cout << "Start value " << k << "\t" << start_val_gauss[k] << std::endl;
	    
	    fitter->mnparm(par_num=k, "alpha_" + s, start_val_gauss[k], start_step, A, B, ierflg); //mjk: https://root.cern.ch/root/html/TMinuit.html#TMinuit:mnparm
	    
	  }

    //-----------> FIXME (added)
    //fitter->SetErrorDef(0.5); // Set error Definition: 1 for Chi square , 0.5 for negative log likelihood


	  const int par_size = fNScaleFactors;
	  
	  arglist[0] = 1.0;
    
	  fitter->mnexcm("SET ERR", arglist, 1, ierflg_err);

    //Starting the fit
    arglist[0] = 10000.0; // number of iterations 

	  arglist[1] = 0.01;
	  fitter->mnexcm("MIGRAD", arglist, 2, ierflg);  // mjk: MIGRAD: Causes minimization of the function by the method of Migrad, the most efficient and complete single method, recommended for general functions    
	  
	  arglist[0] = 1000.0;
	  arglist[1] = 0;
	  
    //not needed for tests in goe, but consider errors larger than 1 then
    // mjk: Instructs Minuit to calculate, by finite differences, the Hessian or error matrix
  
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
	  	  
	  // std::cout << "Hier 1k" << std::endl;

	  for (int k = 0; k != par_size; k++)
	    { 

	      fitter->GetParameter(k,par_value,par_err_value);
	      best_par.push_back(par_value);
	      par_err.push_back(par_err_value);
	      //      std::cout << par_err_value << std::endl;
	      
	      std::stringstream oss_k,oss_par,oss_err;
	      oss_k   << k;
	      oss_par << best_par[k];
	      oss_err << par_err[k];
	      
	      //	      WriteInfoStatus("ExternalSystematics", "Best parameter "+oss_k.str()+", par: "+oss_par.str()+", err: "+oss_err.str());
	      
	    }
	  
	  //Get error matrix / covariance matrix
	  const int emsize = par_size;//k values and scale factors
	  double errormatrix[emsize][emsize];

	  fCorrelationmatrix[emsize][emsize];

	  fitter->mnemat(&errormatrix[0][0],emsize);
	  
	  //Sum of all fitted signal parameters alpha_(F0,FL,FR)
	  double N_sum = 0.0;
	  
	  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig)
	    N_sum += best_par[iSig];

	  //std::cout <<"ierflg(3)= " <<ierflg << std::endl;
	  
	  if(ierflg == 0){

	    int nParam = fParameterNames.size();

	    for(int iParam = 0; iParam < nParam; ++iParam){

	      fParameterValues[iParam]     = best_par[iParam];
	      fParameterValues_err[iParam] = par_err[iParam];

        //std::cout<<"fParameterValues[iParam]= " << fParameterValues[iParam] << std::endl;
	         
	      	 //     if(iParam < 3)
	      	//	std::cout << "###iparam: " <<iParam << "\t" << best_par[iParam] << "\t" << best_par[iParam]*SelEff[iParam] << "\t"<< SelEff[iParam] <<std::endl;

	      double unc_factor;
	      double sum = 0.0;
	      if(fFittingMode == "3D")
		      sum = best_par[0]+best_par[1]+best_par[2];
	      else{

		      sum = best_par[0]+best_par[1];

		      double termN0 = best_par[0]/sum/sum*par_err[1];
		      double termNL = best_par[1]/sum/sum*par_err[0];
		
		      unc_factor = (termN0 + termNL);
		
	      }

	  if(iParam == 0){
		  //fParameterValues_nom[iParam] = F0_true*fLumi*fXsec;
      fParameterValues_nom[iParam] = F0_true*fLumi*fXsec*fKfactor*fF0eff;    // ----------- TEST (13 Nov 2015 by MJK)
      //std::cout << "### fParameterValues_nom = " << fParameterValues_nom[iParam] << std::endl;
		  fFractionValues[iParam]      = best_par[0]/sum;
		  if(fFittingMode == "3D")
		    fFractionValues_err[iParam]  = fi_error(&errormatrix[0][0], emsize, best_par[0], sum, 0);
		  else
		    fFractionValues_err[iParam]  = unc_factor;
		  fFractionValues_nom[iParam]  = F0_true;

		fVarCorrN0NL = errormatrix[1][0]/par_err[0]/par_err[1];
		
		if(mode == "Datafit")
		  std::cout << "F0 = " << fFractionValues[iParam] << "\t" << fFractionValues_err[iParam] << std::endl;

	      }
	      if(iParam == 1){
		//fParameterValues_nom[iParam] = FL_true*fLumi*fXsec;
    fParameterValues_nom[iParam] = FL_true*fLumi*fXsec*fKfactor*fFLeff;    // ----------- TEST (13 Nov 2015)
		fFractionValues[iParam]      = best_par[1]/sum;
		if(fFittingMode == "3D")
		  fFractionValues_err[iParam]  = fi_error(&errormatrix[0][0], emsize, best_par[1], sum, 1);
		else
		  fFractionValues_err[iParam]  = unc_factor;
		fFractionValues_nom[iParam]  = FL_true;

		if(mode == "Datafit")
		  std::cout << "FL = " << fFractionValues[iParam] << "\t" << fFractionValues_err[iParam] << std::endl;

	      }
    if(iParam == 2 && fFittingMode == "3D"){
		//fParameterValues_nom[iParam] = FR_true*fLumi*fXsec;
    fParameterValues_nom[iParam] = FR_true*fLumi*fXsec*fKfactor*fFReff;    // ----------- TEST (13 Nov 2015)
		fFractionValues[iParam]      = best_par[2]/sum;
                fFractionValues_err[iParam]  = fi_error(&errormatrix[0][0], emsize, best_par[2], sum, 2);
                fFractionValues_nom[iParam]  = FR_true;

		fVarCorrN0NR = errormatrix[2][0]/par_err[0]/par_err[2];
		fVarCorrNRNL = errormatrix[2][1]/par_err[2]/par_err[1];

		if(mode == "Datafit")
		  std::cout << "FR = " << fFractionValues[iParam] << "\t" << fFractionValues_err[iParam] << std::endl;

	      }

	      //	      std::cout << best_par[iParam] << "\t" << par_err[iParam] << std::endl;



	      //	      std::cout << "corr " << errormatrix[1][0]/par_err[0]/par_err[1] << "\t" << errormatrix[2][0]/par_err[0]/par_err[2] << "\t" << errormatrix[1][2]/par_err[2]/par_err[1] << std::endl;
	           
	      //	      std::cout << "corr diag" << errormatrix[0][0]/par_err[0]/par_err[0] << std::endl;
	      //    std::cout << "corr diag" << errormatrix[1][1]/par_err[1]/par_err[1] << std::endl;
	      //  std::cout << "corr diag" << errormatrix[2][2]/par_err[2]/par_err[2] << std::endl;

	      //	      std::cout << fParameterValues[iParam] << "\t" << best_par[iParam] << std::endl;
	      
	    } // end loop over iParam

	    fNewTree -> Fill();

	    counter_good++;

	  
	    if(mode == "Datafit"){


	      fFitErrorHist = TemplateFitErrors(&errormatrix[0][0], fParameterHist, SelEff, fSignalParameter.size());

	      for(int iParam = 0; iParam < nParam; ++iParam){

		for(int jParam = 0; jParam < nParam; ++jParam){

		  //		  std::cout << "2" << std::endl;

		  fCorrelationmatrix[iParam][jParam] = errormatrix[iParam][jParam]/par_err[iParam]/par_err[jParam];

		}
		
	      }

	    }
	    
	  }
	  
	  // std::cout << "Hier 1m" << std::endl;

	} // end while(counter_good < NumberPseudoExp)
  
     // fNewFile->cd();

      // std::cout << "Hier 1n" << std::endl;

      fNewTree -> Write();
      //      fNewFile -> Close();

      fitter   -> Delete();


}

//void ExternalSystematics::EvaluateSystematic(std::string mode, std::string OutputTxtFile)
void ExternalSystematics::EvaluateSystematic(std::string mode, std::string OutputTxtFile, std::vector<std::string> Label)
{


  for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){

    TFile *fFile = new TFile(fSystList[iSyst].c_str(), "READ");

    TTree *tree = (TTree*) fFile -> Get("EnsembleTree");
    
    int nParam = fParameterNames.size();
    
    std::vector<double> ParamValues(nParam);

    for(int iParam = 0; iParam < nParam; ++iParam){
      
      tree -> SetBranchAddress(fParameterNames[iParam].c_str(), &(ParamValues[iParam]));
      std::cout<<"!!!! Branch: "<<fParameterNames[iParam].c_str()<<"\t Value: "<<&(ParamValues[iParam])<<std::endl;
	  
    }
    
    //TH1D protoN0 = TH1D("protoN0",     "", 75,  180000.0,  490000.0);
    //TH1D protoNL = TH1D("protoNL",     "", 75,   40000.0,  330000.0);
    //TH1D protoNR = TH1D("protoNR",     "", 75, -125000.0,  165000.0);
    
    TH2D *N0NL   = new TH2D("", "", 100, 180000.0, 490000.0, 100,  40000.0, 490000.0);
    TH2D *N0NR   = new TH2D("", "", 100, 180000.0, 490000.0, 100, -65000.0,  65000.0);
    TH2D *NRNL   = new TH2D("", "", 100, -65000.0,  65000.0, 100,  40000.0, 230000.0);

    int nentries = tree -> GetEntries();

    for(int i = 0; i < nentries; ++i){
      
      tree -> GetEntry(i);

      double N0 = ParamValues[0];
      double NL = ParamValues[1];
      double NR = ParamValues[2];

      N0NL -> Fill(N0, NL);

      if(fFittingMode == "3D"){
	N0NR -> Fill(N0, NR);
	NRNL -> Fill(NR, NL);
      }

      double N_sum = 0.0;
      if(fFittingMode == "3D")
	N_sum = N0+NL+NR;
      else
	N_sum = N0+NL;

      for(int iParam = 0; iParam < fParameterNames.size(); ++iParam){

	fHistoAll[iParam][iSyst].Fill(ParamValues[iParam]);

      }

      fHistoF0[iSyst].Fill(N0/N_sum);                                                                                                                 
      fHistoFL[iSyst].Fill(NL/N_sum);
      if(fFittingMode == "3D")
	     fHistoFR[iSyst].Fill(NR/N_sum);

    }

    double corrN0NL = N0NL -> GetCovariance();
    double corrN0NR = N0NR -> GetCovariance();
    double corrNRNL = NRNL -> GetCovariance();

    delete tree;
    delete fFile;

  }

  std::cout<<"###@@@ fHistoF0[0].Integral()= " << fHistoF0[0].Integral() << std::endl;
  std::cout<<"###@@@ fHistoFL[0].Integral()= " << fHistoFL[0].Integral() << std::endl;
  std::cout<<"###@@@ fHistoFR[0].Integral()= " << fHistoFR[0].Integral() << std::endl;
  fHistoF0[0].Fit("gaus","Q"); 
  fHistoFL[0].Fit("gaus","Q");
  fHistoFR[0].Fit("gaus","Q");

  std::cout<< "#####@@@ fSystList.size(): " <<fSystList.size() <<std::endl;
  // Feb 25 BT: Store means for all variations only when running PDF systematics
  if( mode=="PDF" || mode=="mass") {

    for(int i = 0; i < fSystList.size(); ++i){
      fHistoF0[i].Fit("gaus","Q");
      fHistoFL[i].Fit("gaus","Q");
      fHistoFR[i].Fit("gaus","Q");
      
      meanF0.push_back(fHistoF0[i].GetFunction("gaus")->GetParameter(1) );
      meanFL.push_back(fHistoFL[i].GetFunction("gaus")->GetParameter(1) );
      if(fFittingMode == "3D")
	       meanFR.push_back(fHistoFR[i].GetFunction("gaus")->GetParameter(1) );
    }
  } // end PDF fill
  else{ // follow normal nominal - variation scheme for all other systematics
    
    for(int i = 1; i < fSystList.size(); ++i){

      fHistoF0[i].Fit("gaus","Q");
      fHistoFL[i].Fit("gaus","Q");
      fHistoFR[i].Fit("gaus","Q");
      
      // diffF0.push_back(fabs(fHistoF0[0].GetMean() - fHistoF0[i].GetMean()));
      // diffFL.push_back(fabs(fHistoFL[0].GetMean() - fHistoFL[i].GetMean()));
      // if(fFittingMode == "3D")
      //   diffFR.push_back(fabs(fHistoFR[0].GetMean() - fHistoFR[i].GetMean()));
      
      
      //--new: mjk: using a gaussian fit to read out the mean val. instead of GetMean().
      std::cout<< "#####@@@ EvaluateSystematic: " << " fHistoF0[0]= "<<fHistoF0[0].GetFunction("gaus")->GetParameter(1) <<std::endl;
      //std::cout<< "#####@@@ EvaluateSystematic: " << " fHistoF0[0]= "<<fHistoF0[0].GetFunction("gaus")->GetParameter(1) << "\t fHistoF0[i]= "<<fHistoF0[i].GetFunction("gaus")->GetParameter(1) << "\t diff= "<< fHistoF0[i].GetFunction("gaus")->GetParameter(1) - fHistoF0[0].GetFunction("gaus")->GetParameter(1) <<std::endl;
      //std::cout<< "##### EvaluateSystematic: " << " fHistoFL[0]= "<<fHistoFL[0].GetFunction("gaus")->GetParameter(1) << "\t fHistoFL[i]= "<<fHistoFL[i].GetFunction("gaus")->GetParameter(1) << "\t diff= "<< fHistoFL[i].GetFunction("gaus")->GetParameter(1) - fHistoFL[0].GetFunction("gaus")->GetParameter(1) <<std::endl;
      // std::cout<< "##### EvaluateSystematic: " << " fHistoFR[0]= "<<fHistoFR[0].GetFunction("gaus")->GetParameter(1) << "\t fHistoFR[i]= "<<fHistoFR[i].GetFunction("gaus")->GetParameter(1) << "\t diff= "<< fHistoFR[i].GetFunction("gaus")->GetParameter(1) - fHistoFR[0].GetFunction("gaus")->GetParameter(1) <<std::endl;
      
      diffF0.push_back(fHistoF0[i].GetFunction("gaus")->GetParameter(1) - fHistoF0[0].GetFunction("gaus")->GetParameter(1));
      diffFL.push_back(fHistoFL[i].GetFunction("gaus")->GetParameter(1) - fHistoFL[0].GetFunction("gaus")->GetParameter(1));
      if(fFittingMode == "3D")
	     diffFR.push_back(fHistoFR[i].GetFunction("gaus")->GetParameter(1) - fHistoFR[0].GetFunction("gaus")->GetParameter(1));
    }
  }
  


  double FinalF0 = 0.0;
  double FinalFL = 0.0;
  double FinalFR = 0.0;

  if(mode != "Dependence" && mode != "FractionWidth"){

    for(int i = 0; i < diffF0.size(); ++i){
      
      if(diffF0[i] > FinalF0)
	FinalF0 = diffF0[i];
      
      if(diffFL[i] > FinalFL)
        FinalFL = diffFL[i];
      
      if(fFittingMode == "3D"){
	if(diffFR[i] > FinalFR)
	  FinalFR = diffFR[i];
      }
      
    } // largest difference


    if(mode == "HalfDiff"){
      
      FinalF0 = FinalF0*0.5;
      FinalFL = FinalFL*0.5;

      if(fFittingMode == "3D")
	FinalFR = FinalFR*0.5;
      
    }
      

  }
  else if(mode == "FractionWidth"){
    
    for(int i = 0; i < fSystList.size(); ++i){
      
      FinalF0 += fHistoF0[i].GetRMS()*fHistoF0[i].GetRMS();
      FinalFL += fHistoFL[i].GetRMS()*fHistoFL[i].GetRMS();
    
      if(fFittingMode == "3D")
	FinalFR += fHistoFR[i].GetRMS()*fHistoFR[i].GetRMS();
      
    }
    
  }

  std::ofstream file;

  file.open(OutputTxtFile.c_str(),fstream::app);

  //file << fSystematicLabel << "\t" << sqrt(FinalF0) << "\t" << sqrt(FinalFL) << "\t" << sqrt(FinalFR) << std::endl;
  
  // Feb 22, BT
  if(mode != "PDF" && mode != "mass") { // normal output of difference between variation and nominal
    // keep up and down with sign
    for(int i = 0; i < fSystList.size()-1; ++i){
      //std::cout<< "@@@@@ EvaluateSystematic: " <<diffF0[i] << "\t" << diffFL[i]<< "\t" << diffFR[i]  <<std::endl;
      file << Label[i+1] << "\t& " << diffF0[i] << "\t& " << diffFL[i] << "\t& " << diffFR[i] << " \\ q"<< std::endl;
    }
  }
  else{ // Print values of each member in vector (mean values, not differences!)
    for(int i = 0; i < fSystList.size(); ++i){
      file << Label[i] << "\t& " << meanF0[i] << "\t& " << meanFL[i] << "\t& " << meanFR[i] << " \\ q"<< std::endl; 
    }
  }

  file.close();

}

std::pair<double,double> ExternalSystematics::CalculateParameterRanges(int Param, bool err, bool fraction)
{

  double xmin = 1000000.0;
  double xmax =  -10000.0;

  for(int iSyst = 0; iSyst < fSystList.size(); ++iSyst){

    TFile *fFile = new TFile(fSystList[iSyst].c_str(), "READ");

    TTree *tree  = (TTree*) fFile -> Get("EnsembleTree");

    int nParam   = fParameterNames.size();

    double Value;

    std::cout<<"2 !!!! Branch: "<<fParameterNames[iSyst].c_str()<<"\t Value: "<<&Value<<std::endl;

    if(!fraction){

      if(err)
	tree -> SetBranchAddress((fParameterNames[Param]+"_err").c_str(), &Value);
      else
	tree -> SetBranchAddress(fParameterNames[Param].c_str(),          &Value);

    }
    else{
      
      if(err)
        tree -> SetBranchAddress((fFractionNames[Param]+"_err").c_str(), &Value);
      else{
        //std::cout<<"&&& === fFractionNames["<<Param << "]= " <<fFractionNames[Param] << std::endl;
        tree -> SetBranchAddress(fFractionNames[Param].c_str(), &Value);
      }
        
      
    }

    int nentries = tree -> GetEntries();

    for(int i = 0; i < nentries; ++i){

      tree -> GetEntry(i);

      if(xmin > Value) xmin = Value;
      if(xmax < Value) xmax = Value;
      
    }
    
    fFile -> Close();

  }

  double range  = xmax-xmin;
  double middle = range/2.0;
 
  xmin = xmin-0.25*middle;
  xmax = xmax+0.25*middle;

  std::pair<double,double> edges;
  edges.first  = xmin;
  edges.second = xmax;
 
  return edges;
  
}

void ExternalSystematics::InitializeHistogramsAnalysis(int nrHisto, int iParam)
{

  //  std::cout << "Number Histo " << nrHisto << std::endl;

  std::vector<double> min,max;

  std::pair<double, double> result = CalculateParameterRanges(iParam, false, false);

  TH1D help_val     = TH1D("", "", 75, result.first, result.second);

  std::pair<double, double> result_unc = CalculateParameterRanges(iParam, true, false);

  TH1D help_val_unc = TH1D("", "", 75, result_unc.first, result_unc.second);

  TH1D help_pull    = TH1D("", "", 71, -3.55, 3.55);

  std::vector<TH1D> hist_summary;
  hist_summary.push_back(help_val);
  hist_summary.push_back(help_val_unc);
  hist_summary.push_back(help_pull);

  fHistoAnalysis.push_back(hist_summary);

}

void ExternalSystematics::InitializeHistogramsFractions()
{

  //  std::cout << "Number Histo " << nrHisto << std::endl;

  std::pair<double, double> resultF0 = CalculateParameterRanges(0, false, true);

  TH1D help_val_F0  = TH1D("", "", 75, resultF0.first-0.2, resultF0.second+0.2);
  

  std::pair<double, double> resultFL = CalculateParameterRanges(1, false, true);

  TH1D help_val_FL  = TH1D("", "", 75, resultFL.first-0.2, resultFL.second+0.2);

  //  std::pair<double, double> resultFR = CalculateParameterRanges(2, false, true);

  TH1D help_val_FR;
  if(fFittingMode == "3D"){
  
    std::pair<double, double> resultFR = CalculateParameterRanges(2, false, true);
    help_val_FR = TH1D("", "", 75, resultFR.first-0.2, resultFR.second+0.2);
  
  }

  //  std::cout << result.first << "\t" << result.second << std::endl;

  //  std::pair<double, double> result_unc = CalculateParameterRanges(iParam, true, true);

  //  std::cout << result_unc.first << "\t" << result_unc.second << std::endl;

  TH1D help_val_unc_F0 = TH1D("", "", 80, 0.025,  0.065);
  TH1D help_val_unc_FL = TH1D("", "", 80, 0.0075, 0.0375);
  TH1D help_val_unc_FR = TH1D("", "", 80, 0.0165,  0.0445);

  if(fFittingMode == "2D"){
    
    if(fParameterNames.size() > 6){

      help_val_unc_F0 = TH1D("", "", 80, 0.012, 0.023);
      help_val_unc_FL = TH1D("", "", 80, 0.012, 0.023);

    }
    else{

      if(fChannel == "el"){

	help_val_unc_F0 = TH1D("", "", 80, 0.014, 0.038);
	help_val_unc_FL = TH1D("", "", 80, 0.014, 0.038);

      }
      else{

	help_val_unc_F0 = TH1D("", "", 80, 0.014, 0.0285);
        help_val_unc_FL = TH1D("", "", 80, 0.014, 0.0285);

      }

    }

  }
  else{

    if(fParameterNames.size() > 6){

      help_val_unc_F0 = TH1D("", "", 80, 0.018, 0.042);
      help_val_unc_FL = TH1D("", "", 80, 0.006, 0.024);
      help_val_unc_FR = TH1D("", "", 80, 0.005, 0.032);

    }

  }

  TH1D help_pull    = TH1D("", "", 71, -3.55, 3.55);

  std::vector<TH1D> hist_summary_F0;
  hist_summary_F0.push_back(help_val_F0);
  hist_summary_F0.push_back(help_val_unc_F0);
  hist_summary_F0.push_back(help_pull);

  std::vector<TH1D> hist_summary_FL;
  hist_summary_FL.push_back(help_val_FL);
  hist_summary_FL.push_back(help_val_unc_FL);
  hist_summary_FL.push_back(help_pull);

  std::vector<TH1D> hist_summary_FR;
  hist_summary_FR.push_back(help_val_FR);
  hist_summary_FR.push_back(help_val_unc_FR);
  hist_summary_FR.push_back(help_pull);

  std::vector<std::vector<TH1D> > help;
  help.push_back(hist_summary_F0);
  help.push_back(hist_summary_FL);
 
  if(fFittingMode == "3D")
    help.push_back(hist_summary_FR);

  fFractionAnalysis.push_back(help);

}


void ExternalSystematics::InitializeHistograms(int nrHisto){
    
  //define prototype histograms
  TH1D protoF0     = TH1D("protoF0",     "", 75,      0.5,      0.9);
  TH1D protoFL     = TH1D("protoFL",     "", 75,      0.1,      0.5);
  TH1D protoFR     = TH1D("protoFR",     "", 75,     -0.25,      0.25);

  //TH1D protoF0     = TH1D("protoF0",     "", 40,      0.6,      0.8);
  //TH1D protoFL     = TH1D("protoFL",     "", 40,      0.2,      0.4);
  //TH1D protoFR     = TH1D("protoFR",     "", 40,     -0.05,      0.05);
  
  std::vector<double> min,max;

  std::vector<TH1D> Protos;

  for(int iParam = 0; iParam < fParameterNames.size(); ++iParam){

    std::pair<double, double> result = CalculateParameterRanges(iParam, false, false);

    //    std::cout << result.first << "\t" << result.second << std::endl;
    
    TH1D help = TH1D("", "", 75, result.first, result.second);
 
    Protos.push_back(help);

  }

  for(int i = 0; i < nrHisto; ++i){

    fHistoF0.push_back(protoF0);
    fHistoFL.push_back(protoFL);
    fHistoFR.push_back(protoFR);

  }

  for(int iParam = 0; iParam < fParameterNames.size(); ++iParam){

    std::vector<TH1D> helpList;

    for(int i = 0; i < nrHisto; ++i){
    
      helpList.push_back(Protos[iParam]);

    }

    fHistoAll.push_back(helpList);

  }
  
}

void ExternalSystematics::PlotSystematicDistributions(std::vector<TH1D> hist, std::vector<std::string> Label, std::string Variable)
{

  //SetAtlasStyle();

  std::vector<TH1D*> histo;

  for (int i = 0; i != hist.size(); i++) {
    TH1D * hist_dummy = &hist[i];
    histo.push_back(hist_dummy);
  }

  TCanvas *c0 = new TCanvas("", "", 820, 820);
  c0->SetLeftMargin(0.10);
  c0->SetRightMargin(-0.05);

  double nmax = histo[0] -> GetBinContent(histo[0] -> GetMaximumBin())*1.75;

  int nbins   = histo[0] -> GetNbinsX();

  double lower_edge  = histo[0] -> GetBinLowEdge(1);
  double bin_width   = histo[0] -> GetBinWidth(1);
  double number_bins = histo[0] -> GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;

  for(int iHist = 0; iHist < histo.size(); ++iHist){

    histo[iHist] -> SetFillColor(kWhite);
    histo[iHist] -> SetLineWidth(2);

    if(iHist==2) histo[iHist] -> SetLineColor(kCyan+2);
    if(iHist==1) histo[iHist] -> SetLineColor(kOrange+8);

    // if(iHist!=4) histo[iHist] -> SetLineColor(iHist+1);
    //    else histo[iHist] -> SetLineColor(kMagenta);
    histo[iHist] -> SetMarkerSize(0);
    TF1 *fit;
    if(iHist == 0){
      histo[iHist] -> Draw();
      histo[iHist] ->Fit("gaus","Q");
      fit = histo[iHist]->GetFunction("gaus");
      fit->SetLineColor(histo[iHist]->GetLineColor());
    }
      
    else{
      histo[iHist] -> Draw("SAME");
      histo[iHist] ->Fit("gaus","Q");
      fit = histo[iHist]->GetFunction("gaus");
      fit->SetLineColor(histo[iHist]->GetLineColor());
    }
      

  }
  histo[0] -> GetXaxis() -> SetTitle(Variable.c_str());
  histo[0] -> GetYaxis() -> SetTitle("Ensembles");

  histo[0] -> GetXaxis()->SetTitleSize(0.035);
  histo[0] -> GetYaxis()->SetTitleSize(0.035);

  histo[0] -> GetXaxis()->SetLabelSize(0.032);
  histo[0] -> GetYaxis()->SetLabelSize(0.032);

  histo[0] -> GetXaxis() -> SetTitleOffset(1.3);
  histo[0] -> GetYaxis() -> SetTitleOffset(1.3);

  std::stringstream diff;
  //diff  << setprecision(3) << histo[0] -> GetMean();
  //--new: using a gaussian fit to read out the mean val. instead of GetMean().
  diff  << setprecision(3) << histo[0]->GetFunction("gaus")->GetParameter(1);

  Label[0] = Label[0]+" (Mean = "+diff.str()+")";


  for(int i = 1; i < Label.size(); ++i){

    std::stringstream diff;
    //diff  << setprecision(3) << (histo[i] -> GetMean() - histo[0] -> GetMean());

    diff  << setprecision(3) << (histo[i]->GetFunction("gaus")->GetParameter(1) - histo[0]->GetFunction("gaus")->GetParameter(1));

    if(i == 1)
      Label[i] = Label[i]+"    (Diff = "+diff.str()+")";
    else
      Label[i] = Label[i]+" (Diff = "+diff.str()+")";
 
  }


  //define the legend...
  //TLegend *fLegend = new TLegend(0.50, 0.79, 0.80, 0.925);
  TLegend *fLegend = new TLegend(0.50, 0.79, 0.80, 0.9);
  for(unsigned int k = 0; k < histo.size(); ++k)
    fLegend -> AddEntry(histo[k],  Label[k].c_str(), "l");
  fLegend   -> SetFillColor(0);
  fLegend   -> SetLineColor(0);
  fLegend   -> SetBorderSize(0);
  fLegend   -> SetTextFont(42);
  //fLegend   -> SetTextSize(0.0275);
  fLegend   -> SetTextSize(0.02);
  fLegend   -> Draw();

  histo[0]  -> SetMaximum(nmax);
  histo[0]  -> SetMinimum(0.01);

  TLatex l2;
  l2.SetTextAlign(9);

	//TLatex l; //l.SetTextAlign(12); 
	l2.SetTextSize(0.032); 
	l2.SetNDC();
	l2.SetTextFont(72);
	//l.SetTextColor();
	l2.DrawLatex(0.14,0.88,"ATLAS");
	l2.SetTextFont(42);
	l2.DrawLatex(0.25,0.88,"Work in progress");

  l2.SetTextFont(42);
  l2.SetTextSize(0.032);
  l2.SetNDC();
  if(fChannel == "mu")
    l2.DrawLatex(0.14, 0.835, "#mu+jets channel");
  else if(fChannel == "el")
    l2.DrawLatex(0.14, 0.835, "e+jets channel");
  else
    l2.DrawLatex(0.14, 0.835, "Combined channel");

  /*  int NrSubPlots = fLeptonLabel.size();
  
  double SubPlotWidth     = 1.0/double(NrSubPlots);
  double SubPlotHalfWidth = SubPlotWidth/2.0; 
  
  for(int iPlot = 0; iPlot < NrSubPlots; ++iPlot){

    double SubPlotPosition = (1+2*iPlot)*SubPlotHalfWidth;

    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextFont(42);
    l3.SetTextSize(0.032);
    l3.SetNDC();
    l3.DrawLatex(SubPlotPosition, 0.880, (fLeptonLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.860, (fJetBinLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.840, (fBTagLabel[iPlot]).c_str());
    
  }

  */

  std::stringstream ensembles;
  ensembles << fNPseudoExp;

  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextFont(42);
  l3.SetTextSize(0.032);
  l3.SetNDC();
  l3.DrawLatex(0.14, 0.79, (ensembles.str()+" Ensembles").c_str());

  c0 -> Draw();

  std::string outputFile     = fOutputFolder+"/Syst_"+fSystematicLabel+"_"+fChannel+"_"+Variable+".eps";
  std::string outputFile_png = "Syst_"+fSystematicLabel+"_"+fChannel+"_"+Variable+".png";

  c0 -> Print((outputFile).c_str());
  c0 -> Print((fOutputFolder+"/"+outputFile_png).c_str());

  std::ofstream* htmlfile = new std::ofstream();

  htmlfile -> open((_htmls[fHTMLLabel]).c_str(),fstream::app);
  *htmlfile << "<tr>"
	    << "<td>"<< Variable << "</td>";

  *htmlfile << "<td>";
  *htmlfile << "<p><a href=\""   << outputFile_png.c_str() << "\">";
  *htmlfile << "<img src=\""     << outputFile_png.c_str() << "\">";
  *htmlfile << "</td>"           << std::endl;

  htmlfile->close();


}

void ExternalSystematics::createHTML(){

  // create the page
  std::ofstream page;
  std::string   pname = fOutputFolder+"/index_Syst_"+fSystematicLabel+"_"+fChannel+".html";

  fHTMLLabel = "index_Syst_"+fSystematicLabel+"_"+fChannel+".html";

  struct stat buf;
  //  if(stat(pname.c_str(), &buf) == -1){

    page.open(pname.c_str());

    page << "<html><head><title> Output from systematic evaluation </title></head>" << std::endl;
    page << "<body>" << std::endl;
    page << "<h1> Comparison of nominal and systematic ensemble tests </h1>"    << std::endl;
    page << "<a href="  << fSystematicLabel << "> >=4 jets </a>" << std::endl;
    
    page << "<br> <br>" << std::endl;
    
    page << "<table border = 1> <tr>"
         << "<th> name </th>"
         << "<th> variable </th>"
         << "</tr>" << std::endl;

    _htmls[fHTMLLabel] = pname;
    page.close();

    //  }

  return;

}

void ExternalSystematics::InitializeGraphs()
{

  int nSyst = fSystList.size();

  TGraphErrors help = TGraphErrors(nSyst);

  std::vector<TGraphErrors> helpVec;
  helpVec.push_back(help);
  helpVec.push_back(help);
  helpVec.push_back(help);
  helpVec.push_back(help);

  for(int i = 0; i < fParameterValues.size(); ++i){

    fGraphs.push_back(helpVec);
    
  }

  fFractionGraphs.push_back(helpVec);
  fFractionGraphs.push_back(helpVec);

  if(fFittingMode == "3D")
    fFractionGraphs.push_back(helpVec);

}


TTree *ExternalSystematics::InitializeOutputTree(){

  TTree *OutputTree = new TTree("EnsembleTree", "EnsembleTree");

  int nParam = fParameterNames.size();

  fFractionValues.clear();
  fFractionValues_err.clear();
  fFractionValues_nom.clear();

  fParameterValues.clear();
  fParameterValues_err.clear();
  // fParameterValues_nom.clear();

  std::vector<double> ValuesHelp(nParam);
  std::vector<double> ValuesHelp_err(nParam);
  std::vector<double> ValuesHelp_nom(nParam);

  std::vector<double> FractionValuesHelp(3,0);
  std::vector<double> FractionValuesHelp_err(3,0);
  std::vector<double> FractionValuesHelp_nom(3,0);



  std::cout<<"..(1) nParam= "<<nParam<<std::endl;
  // std::cout<<"..(1) fFractionValues.size()= "<<fFractionValues.size()<<std::endl;
  // std::cout<<"..(1) fFractionValues_err.size()= "<<fFractionValues_err.size()<<std::endl;
  // std::cout<<"..(1) fFractionValues_nom.size()= "<<fFractionValues_nom.size()<<std::endl;
  // std::cout<<"..(1) fParameterValues.size()= "<<fParameterValues.size()<<std::endl;
  // std::cout<<"..(1) fParameterValues_err.size()= "<<fParameterValues_err.size()<<std::endl;

  fFractionValues      = FractionValuesHelp;
  fFractionValues_err  = FractionValuesHelp_err;
  fFractionValues_nom  = FractionValuesHelp_nom;

  fParameterValues     = ValuesHelp;
  fParameterValues_err = ValuesHelp_err;
  //  fParameterValues_nom = ValuesHelp_nom;
    
  for(int iParam = 0; iParam < nParam; ++iParam){                               
    std::cout<<"3!!!! Branch: "<<fParameterNames[iParam].c_str()<<"\t Value: "<<fParameterValues_nom[iParam]<<std::endl;
    OutputTree -> Branch(fParameterNames[iParam].c_str(),          &(fParameterValues[iParam]));
    OutputTree -> Branch((fParameterNames[iParam]+"_err").c_str(), &(fParameterValues_err[iParam]));
    OutputTree -> Branch((fParameterNames[iParam]+"_nom").c_str(), &(fParameterValues_nom[iParam]));

  }
  
  OutputTree -> Branch("F0",          &(fFractionValues[0]));
  OutputTree -> Branch("F0_err",      &(fFractionValues_err[0]));
  OutputTree -> Branch("F0_nom",      &(fFractionValues_nom[0]));
  OutputTree -> Branch("FL",          &(fFractionValues[1]));
  OutputTree -> Branch("FL_err",      &(fFractionValues_err[1]));
  OutputTree -> Branch("FL_nom",      &(fFractionValues_nom[1]));

  OutputTree -> Branch("CorrN0NL",    &fVarCorrN0NL);
  OutputTree -> Branch("CorrN0NR",    &fVarCorrN0NR);
  OutputTree -> Branch("CorrNRNL",    &fVarCorrNRNL);


  if(fFittingMode == "3D"){

    OutputTree -> Branch("FR",          &(fFractionValues[2]));
    OutputTree -> Branch("FR_err",      &(fFractionValues_err[2]));
    OutputTree -> Branch("FR_nom",      &(fFractionValues_nom[2]));
    
  }

  return OutputTree;

}

void ExternalSystematics::AnalyseDatafit()
{


  //std::cout<<"### fSystList.size()= " << fSystList.size() << std::endl;
  //std::cout<<"### fSystList.[0]= " << fSystList[0] << std::endl;
  TFile *fFile = new TFile(fSystList[0].c_str(), "READ");
  fFile -> cd();

  TTree *tree    = (TTree*) fFile -> Get("EnsembleTree");
  
  TTree *tree_cp = tree -> CloneTree();

  Double_t F0err_corr,FLerr_corr,FRerr_corr;

  TBranch* bF0 = tree_cp -> Branch("F0err_corr", &F0err_corr, "F0err_corr/D");
  TBranch* bFL = tree_cp -> Branch("FLerr_corr", &FLerr_corr, "FLerr_corr/D");
  TBranch* bFR = tree_cp -> Branch("FRerr_corr", &FRerr_corr, "FRerr_corr/D");

  std::vector<double> SignalBranchValues;
  std::vector<double> SignalBranchErrors;
  std::vector<double> BkgBranchValues;
  std::vector<double> BkgBranchErrors;
  std::vector<double> FractionBranchValues;
  std::vector<double> FractionBranchErrors;
  std::vector<double> IntegralFit;

  std::vector<std::string> SignalBranchNames;
  std::vector<std::string> BkgBranchNames;

  int nSignal = fFractionNames.size();
  int nParam  = fParameterNames.size();
  int nBkg    = nParam - nSignal;

  // define branch names...

  // std::cout<<"### fFractionNames.size()= " << fFractionNames.size() << std::endl;
  // std::cout<<"### fParameterNames.size()= " << fParameterNames.size() << std::endl;
  // std::cout<<"### fParameterNames.[0]= " << fParameterNames[0] << std::endl;
  
  for(int i = 0; i < nSignal; ++i){

    SignalBranchValues.push_back(1.0);
    SignalBranchErrors.push_back(1.0);
    SignalBranchNames.push_back(fParameterNames[i]);
    
    FractionBranchValues.push_back(1.0);
    FractionBranchErrors.push_back(1.0);

  }
  
  for(int i = nSignal; i < nParam; ++i){
    
    BkgBranchValues.push_back(1.0);
    BkgBranchNames.push_back(fParameterNames[i]);
    BkgBranchErrors.push_back(1.0);
    
  }
  
  // define branches...
  
  for(int i = 0; i < nSignal; ++i){
    
    tree -> SetBranchAddress(SignalBranchNames[i].c_str(),          &(SignalBranchValues[i]));
    tree -> SetBranchAddress((SignalBranchNames[i]+"_err").c_str(), &(SignalBranchErrors[i]));

    tree -> SetBranchAddress(fFractionNames[i].c_str(),             &(FractionBranchValues[i]));
    tree -> SetBranchAddress((fFractionNames[i]+"_err").c_str(),    &(FractionBranchErrors[i]));
    

  }
  
  for(int i = 0; i < nBkg; ++i){

    tree -> SetBranchAddress(BkgBranchNames[i].c_str(),          &(BkgBranchValues[i]));
    tree -> SetBranchAddress((BkgBranchNames[i]+"_err").c_str(), &(BkgBranchErrors[i]));

  }
  
  tree -> GetEntry(0);
  
  IntegralFit.clear();

  IntegralFit.push_back(SignalBranchValues[0]*SelEff[0]);
  IntegralFit.push_back(SignalBranchValues[1]*SelEff[1]);
  IntegralFit.push_back(SignalBranchValues[2]*SelEff[2]);

  std::cout << SignalBranchValues[0] << "\t" << FractionBranchValues[0] << "\t" << SelEff[0] << std::endl;


  std::cout << SignalBranchValues[0]*FractionBranchValues[0]*SelEff[0] << std::endl;
  std::cout << SignalBranchValues[1]*FractionBranchValues[1]*SelEff[1] << std::endl;
  std::cout << SignalBranchValues[2]*FractionBranchValues[2]*SelEff[2] << std::endl;
  std::cout << BkgBranchValues[0] << std::endl;
  std::cout << BkgBranchValues[1] << std::endl;
  std::cout << BkgBranchValues[2] << std::endl;

  // make now pseudo data with fit values above

  TH1D TestPseudodata = PseudoDataBestFit(IntegralFit, BkgBranchValues);

  //  std::cout << "Integral Data:       " << fDataHist.Integral()      << std::endl;
  std::cout << "Integral PseudoData: " << TestPseudodata.Integral() << std::endl;

  std::vector<double> Frac;  
  if(fFittingMode == "3D"){
    Frac.push_back(0.698); Frac.push_back(0.301); Frac.push_back(0.00041);
    fFraction.push_back(Frac);
  }
  else{

    Frac.push_back(0.700); Frac.push_back(0.300);
    fFraction.push_back(Frac);
    
  }

  std::vector<double> Bkg0; Bkg0.push_back(0.0); Bkg0.push_back(0.0); Bkg0.push_back(0.0);

  std::vector<double> Frac0; Frac0.push_back(0.0);  Frac0.push_back(0.0);  
  
  if(fFittingMode == "3D")
    Frac0.push_back(0.0);

  TH1D SMPseudodata     = PseudoData(Frac,  fBkgNorm);
  TH1D BkgPseudodata    = PseudoData(Frac0, BkgBranchValues);
  TH1D SignalPseudodata = PseudoData(Frac,  Bkg0);

  std::cout << "Integral SM Data:    " << SMPseudodata.Integral()    << std::endl;
  std::cout << "Integral Bkg PD:     " << BkgPseudodata.Integral()   << std::endl;
  std::cout << "Integral Data:       " << fDataHist.Integral()       << std::endl;
  std::cout << "Integral Test PD:    " << TestPseudodata.Integral()    << std::endl;
  std::cout << "Errors:              " << fFitErrorHist.Integral()    << std::endl;


  TGraphErrors data = MakeTGraphErrors(fDataHist);

  PlotDataFitSimple(TestPseudodata, BkgPseudodata, SMPseudodata, data, fChannel, fOutputFolder);

  TFile *fDataFile = new TFile((fOutputFolder+"/DatafitHistograms_"+fChannel+".root").c_str(), "RECREATE");

  BkgPseudodata.Write("SumBkgHist");
  fDataHist.Write("DataHist");
  TestPseudodata.Write("FitPseudoData");
  SMPseudodata.Write("SMPseudoData");
  SignalPseudodata.Write("SumSignalHist");
  fFitErrorHist.Write("StatErrorHist");

  TTree *OutputTree = new TTree("CorrelationTree", "CorrelationTree");

  //  int nParam = fParameterNames.size();

  std::vector<double> correlation(nParam);
  
  int numberParam;

  OutputTree -> Branch("nParam", numberParam);

  for(int i = 0; i < nParam; ++i){
    
    OutputTree -> Branch(Form("corrN%i", i), &(correlation[i]));
    
  }
  
  for(int i = 0; i < nParam; ++i){
    
    for(int j = 0; j < nParam; ++j){
      
      std::cout << i << "\t" << j << "\t" << fCorrelationmatrix[i][j] << std::endl;

      correlation[j] = fCorrelationmatrix[i][j];

    }
    
    OutputTree -> Fill();
    
  }

  OutputTree -> Write();

  fDataFile -> Close();
  
  PlotDataFit(TestPseudodata, BkgPseudodata, SMPseudodata, data, fFitErrorHist, fChannel, fOutputFolder);
  
  /*
  Validation(TestPseudodata, 0, "DatafitPD"); // crash point in 4channel fit !

  AnalyseFractions(0, "Datafit");

  Plot2Dcorrelation_N(fN0NL, FractionBranchValues[0], "N_{0}", "N_{L}", fOutputFolder+"/Correlation_N0_NL_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fN0NR, FractionBranchValues[0], "N_{0}", "N_{R}", fOutputFolder+"/Correlation_N0_NR_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fNRNL, FractionBranchValues[0], "N_{R}", "N_{L}", fOutputFolder+"/Correlation_NR_NL_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fN0N0, FractionBranchValues[0], "N_{0}", "N_{0}", fOutputFolder+"/Correlation_N0_N0_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fNLNL, FractionBranchValues[0], "N_{L}", "N_{L}", fOutputFolder+"/Correlation_NL_NL_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fNRNR, FractionBranchValues[0], "N_{R}", "N_{R}", fOutputFolder+"/Correlation_NR_NR_"+fChannel+".eps", fChannel);

  Plot2Dcorrelation_N(fF0FL, FractionBranchValues[0], "F_{0}", "F_{L}", fOutputFolder+"/Correlation_F0_FL_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fF0FR, FractionBranchValues[0], "F_{0}", "F_{R}", fOutputFolder+"/Correlation_F0_FR_"+fChannel+".eps", fChannel);
  Plot2Dcorrelation_N(fFRFL, FractionBranchValues[0], "F_{R}", "F_{L}", fOutputFolder+"/Correlation_FR_FL_"+fChannel+".eps", fChannel);


  double Nges = 0.0;

  if(fFittingMode == "3D")
    Nges = SignalBranchValues[0]+SignalBranchValues[1]+SignalBranchValues[2];
  else
    Nges = SignalBranchValues[0]+SignalBranchValues[1];

  // calculate uncertainties with improved correlations
  double f0Error = fi_error (&fErrormatrix[0][0],      3, SignalBranchValues[0], Nges, 0);
  double fLError = fi_error (&fErrormatrix[0][0],      3, SignalBranchValues[1], Nges, 1);
  double fRError = 0.0;

  if(fFittingMode == "3D")
    fRError = fabs(fi_error (&fErrormatrix[0][0], 3, SignalBranchValues[2], Nges, 2));
  else{

    f0Error = FractionBranchErrors[0];
    fLError = FractionBranchErrors[1];
  }

  

  // sigma F0
  PlotErrorDataComp(&fFractionAnalysis[0][0][1], 0.698, f0Error, "Statistical uncertainty F_{0}", fChannel, "", fOutputFolder+"/ExpectedUncertainty_F0_"+fChannel+".eps", 1);
  // sigma FL
  PlotErrorDataComp(&fFractionAnalysis[0][1][1], 0.698, fLError, "Statistical uncertainty F_{L}", fChannel, "", fOutputFolder+"/ExpectedUncertainty_FL_"+fChannel+".eps", 1);

  F0err_corr = f0Error;
  FLerr_corr = fLError;
  FRerr_corr = fRError;
 
  tree_cp -> Fill();

  TFile *newFile = new TFile((fOutputFolder+"/DatafitCorrectedUnc_"+fChannel+".root").c_str(), "RECREATE");

  tree_cp -> Write();

  newFile -> Close();

  // sigma FR
  if(fFittingMode == "3D")
    PlotErrorDataComp(&fFractionAnalysis[0][2][1], 0.0, fRError, "Statistical uncertainty F_{R}", fChannel, "", fOutputFolder+"/ExpectedUncertainty_FR_"+fChannel+".eps", 1);
  
  std::ofstream fOut;

  fOut.open((fOutputFolder+"/FitResult_"+fChannel+".tex").c_str(), std::ios::out | std::ios::app );

  //  fOut << "\\begin{table}[!htb]"                << std::endl;
  //  fOut << "\\centering"                         << std::endl;
  fOut << "\\begin{tabular}{l|rrrr|c}"            << std::endl;
  fOut << "\\hline"                             << std::endl;
  fOut << "\\hline"                             << std::endl;
  fOut << "Fraction/Param.   & Sim. & $\\sigma(\\textrm{Sim.})$ & Fit & $\\sigma(\\textrm{Fit})$ & Diff. [$\\sigma_{\\textrm{Sim.}}$] \\\\" << std::endl;
  fOut << "\\hline"                             << std::endl;
  fOut << "F$_{\\textrm{0}}$ &     &  & " << std::setprecision(3) << FractionBranchValues[0] << " & " << std::setprecision(3) << f0Error << " & " << "\\\\" << std::endl;
  fOut << "F$_{\\textrm{L}}$ &     &  & " << std::setprecision(3) << FractionBranchValues[1] << " & " << std::setprecision(3) << fLError << " & " << "\\\\" << std::endl;
  if(fFittingMode == "3D")
    fOut << "F$_{\\textrm{R}}$ &   &  & " << std::setprecision(3) << FractionBranchValues[2] << " & " << std::setprecision(3) << fRError << " & " << "\\\\" << std::endl;
  fOut << "\\hline \\hline"                             << std::endl;
  fOut << "N$_{\\textrm{0}}$ &     &  & " << std::setprecision(7) << SignalBranchValues[0] << " & " << std::setprecision(7) << SignalBranchErrors[0] << " & \\\\" << std::endl;
  fOut << "N$_{\\textrm{L}}$ &     &  & " << std::setprecision(7) << SignalBranchValues[1] << " & " << std::setprecision(7) << SignalBranchErrors[1] << " & \\\\" << std::endl;
  if(fFittingMode == "3D")
    fOut << "N$_{\\textrm{R}}$ &     &  & " << std::setprecision(7) << SignalBranchValues[2] << " & " << std::setprecision(7) << SignalBranchErrors[2] << " & \\\\ " << std::endl;
  fOut << "\\hline"                             << std::endl;
  for(int i = nSignal; i < nParam; ++i) 
    fOut << "N("+fParameterNames[i]+")"  << " & " << std::setprecision(4) << fBkgNorm[i-nSignal] << " & " << std::setprecision(4) << fBkgUnc[i-nSignal] <<  " & " << std::setprecision(4) << BkgBranchValues[i-nSignal] << " & " << std::setprecision(4) << BkgBranchErrors[i-nSignal] << " & " << fabs((fBkgNorm[i-nSignal]-BkgBranchValues[i-nSignal])/fBkgUnc[i-nSignal])  << "\\\\" << std::endl;
 
  fOut << "\\hline"                             << std::endl;
  fOut << "\\hline"                             << std::endl;
  fOut << "\\end{tabular}"                      << std::endl;
  //  fOut << "\\end{table}"                        << std::endl;

  // Calculate fitted xsec...

  double xsec_full = Nges/fLumi/0.543024; // BR from 2012 PDG values...

  fOut << "\t" << std::endl;
  fOut << "\t" << std::endl;

  for(int i = 0; i < nSignal; ++i)
    fOut << "% fitted number of signal events: " << SignalBranchValues[i] << "\t" << SignalBranchErrors[i] << std::endl;

  fOut << "% fitted ttbarXsec with BR = 0.543024 from 2012 PDG and fLumi = " << fLumi << " : " << xsec_full << " pb" << std::endl;

  fOut.close();
*/
}

TGraphErrors ExternalSystematics::MakeTGraphErrors(TH1D hist)
{

  int nBins = hist.GetNbinsX();

  TGraphErrors graph = TGraphErrors(nBins);

  for(int iBin = 1; iBin <= nBins; ++iBin){

    double xval = hist.GetBinCenter(iBin);
    double yval = hist.GetBinContent(iBin);
    double yerr = hist.GetBinError(iBin);
    double xerr = 0.0;

    graph.SetPoint(iBin-1, xval, yval);
    
    graph.SetPointError(iBin-1, xerr, yerr);
    
  }

  return graph;

}

TH1D ExternalSystematics::PseudoDataBestFit(std::vector<double> FractionTimesNorm, std::vector<double> BkgValues)
{

  TH1D HistSum;

  std::vector<TH1D> HistVec = fParameterHist;

  int nSignal = FractionTimesNorm.size();
  int nParam  = fParameterHist.size();
  int nBkg    = nParam-nSignal;

  for(int i = 0; i < nSignal; ++i){

    std::cout << "@@PseudoDataBestFit: "<<SelEff[i] << "\t" << FractionTimesNorm[i] << "\t" << HistVec[i].Integral() << std::endl;

    HistVec[i].Scale(FractionTimesNorm[i]/HistVec[i].Integral());

  }

  for(int i = nSignal; i < nParam; ++i){

    //    std::cout << BkgValues[i-nSignal] << "\t" << HistVec[i].Integral() << std::endl;                                                                              

    HistVec[i].Scale(BkgValues[i-nSignal]/HistVec[i].Integral());

  }

  for(int i = 0; i < nParam; ++i){

    if(i == 0)
      HistSum = HistVec[i];

    else
      HistSum = HistSum + HistVec[i];

  }

  return HistSum;

}

TH1D ExternalSystematics::PseudoData(std::vector<double> Fraction, std::vector<double> BkgValues)
{
  
  TH1D HistSum;

  std::vector<TH1D> HistVec = fParameterHist;

  int nSignal = Fraction.size();
  int nParam  = fParameterHist.size();
  int nBkg    = nParam-nSignal;

  for(int i = 0; i < nSignal; ++i){

    std::cout << Fraction[i] << "\t" << SelEff[i] << "\t" << Fraction[i]*SelEff[i]*fXsec*fLumi << "\t" << HistVec[i].Integral() << std::endl;

    HistVec[i].Scale(Fraction[i]*SelEff[i]*fXsec*fLumi/HistVec[i].Integral());

  }

  for(int i = nSignal; i < nParam; ++i){

    //    std::cout << BkgValues[i-nSignal] << "\t" << HistVec[i].Integral() << std::endl;

    HistVec[i].Scale(BkgValues[i-nSignal]/HistVec[i].Integral());

  }

  for(int i = 0; i < nParam; ++i){

    if(i == 0)
      HistSum = HistVec[i];

    else
      HistSum = HistSum + HistVec[i];
    
  }
    
  return HistSum;
  
}

void ExternalSystematics::AnalyseFractions(int iSyst, std::string mode)
{
  
  InitializeHistogramsFractions();

  TFile *fFile;
  
  if(mode != "Datafit")
    fFile = new TFile(fSystList[iSyst].c_str(), "READ");
  else
    fFile = new TFile(fSystList[1].c_str(), "READ");

  fFile -> cd();

  TTree *tree = (TTree*) fFile -> Get("EnsembleTree");

  double valN0, uncN0, nomF0;
  double valNL, uncNL, nomFL;
  double valNR, uncNR, nomFR;

  std::string Name0 = fParameterNames[0];
  std::string Name1 = fParameterNames[1];
  std::string Name2 = fParameterNames[2];

  tree -> SetBranchAddress(Name0.c_str(), &valN0);
  tree -> SetBranchAddress((Name0+"_err").c_str(), &uncN0);
  tree -> SetBranchAddress("F0_nom",      &nomF0);

  tree -> SetBranchAddress(Name1.c_str(), &valNL);
  tree -> SetBranchAddress((Name1+"_err").c_str(), &uncNL);
  tree -> SetBranchAddress("FL_nom",      &nomFL);

  if(fFittingMode == "3D"){
    
    tree -> SetBranchAddress(Name2.c_str(), &valNR);
    tree -> SetBranchAddress("FR_nom",      &nomFR);

  }

  std::pair<double, double> N0_edges = CalculateParameterRanges(0, false, false);
  std::pair<double, double> NL_edges = CalculateParameterRanges(1, false, false);
  std::pair<double, double> NR_edges = CalculateParameterRanges(2, false, false);

  //  std::cout << "Edges N0: " << N0_edges.first << "\t" << N0_edges.second << std::endl;
  //  std::cout << "Edges NL: " << NL_edges.first << "\t" << NL_edges.second << std::endl;
  //  std::cout << "Edges NR: " << NR_edges.first << "\t" << NR_edges.second << std::endl;

  //  NL_edges.second = 200000.0;
  //  NR_edges.second =  70000.0;

  fN0NL   = new TH2D("", "", 200,  N0_edges.first, N0_edges.second, 200, NL_edges.first, NL_edges.second);
  fN0NR   = new TH2D("", "", 200,  N0_edges.first, N0_edges.second, 200, NR_edges.first, NR_edges.second);
  fNRNL   = new TH2D("", "", 200,  NR_edges.first, NR_edges.second, 200, NL_edges.first, NL_edges.second);
  fN0N0   = new TH2D("", "", 200,  N0_edges.first, N0_edges.second, 200, N0_edges.first, N0_edges.second);
  fNLNL   = new TH2D("", "", 200,  NL_edges.first, NL_edges.second, 200, NL_edges.first, NL_edges.second);
  fNRNR   = new TH2D("", "", 200,  NR_edges.first, NR_edges.second, 200, NR_edges.first, NR_edges.second);

  fF0FL   = new TH2D("", "", 100,  0.45, 0.90, 100,  0.20, 0.45);
  fF0FR   = new TH2D("", "", 100,  0.45, 0.90, 100, -0.20, 0.20);
  fFRFL   = new TH2D("", "", 100, -0.20, 0.20, 100,  0.20, 0.45);

  TH1D N0hist = TH1D("", "", 200,  N0_edges.first, N0_edges.second);
  TH1D NLhist = TH1D("", "", 200,  NL_edges.first, NL_edges.second);
  TH1D NRhist = TH1D("", "", 200,  NR_edges.first, NR_edges.second);

  int entries = tree -> GetEntries();

  // first: get correlation between Ns:
 
  for(int iEvent = 0; iEvent < entries; ++iEvent){
    
    tree -> GetEntry(iEvent);
    
    //    std::cout << "Values " << valN0 << "\t" << valNL << "\t" << valNR << std::endl;

    N0hist.Fill(valN0);
    NLhist.Fill(valNL);
    NRhist.Fill(valNR);

    double Nges = valN0+valNL+valNR;

    fN0NL -> Fill(valN0, valNL);
    fN0N0 -> Fill(valN0, valN0);
    fNLNL -> Fill(valNL, valNL);

    fF0FL-> Fill(valN0/Nges, valNL/Nges);
    fF0FR-> Fill(valN0/Nges, valNR/Nges);
    fFRFL-> Fill(valNR/Nges, valNL/Nges);

    if(fFittingMode == "3D"){
      
      fN0NR -> Fill(valN0, valNR);
      fNRNL -> Fill(valNR, valNL);
      fNRNR -> Fill(valNR, valNR);
      
    }
    
  }
  
  //  fErrormatrix[3][3];
  
  double corrN0NL = fN0NL -> GetCovariance();
  double corrN0NR = fN0NR -> GetCovariance();
  double corrNRNL = fNRNL -> GetCovariance();
  /*
  std::cout << "Covariance " << fN0NL -> GetCorrelationFactor() << "\t" << fN0NR -> GetCorrelationFactor() << "\t" << fNRNL -> GetCorrelationFactor() << std::endl;

  std::cout << fN0N0 -> GetCorrelationFactor() << std::endl;
  std::cout << fNLNL -> GetCorrelationFactor() << std::endl;
  std::cout << fNRNR -> GetCorrelationFactor() << std::endl; */

  //  std::cout << "Covariance " << corrN0NL << "\t" << corrN0NR << "\t" << corrNRNL << std::endl;

  if(fFittingMode == "3D"){

    fErrormatrix[0][0] = fN0N0 -> GetCovariance();
    fErrormatrix[1][1] = fNLNL -> GetCovariance();
    fErrormatrix[2][2] = fNRNR -> GetCovariance();
    fErrormatrix[0][1] = corrN0NL;
    fErrormatrix[0][2] = corrN0NR;
    fErrormatrix[1][0] = corrN0NL;
    fErrormatrix[2][0] = corrN0NR;
    fErrormatrix[1][2] = corrNRNL;
    fErrormatrix[2][1] = corrNRNL;

    //  std::cout << "Cov N0N0: " << fN0N0 -> GetCovariance() << std::endl;
    //  std::cout << "Cov NLNL: " << fNLNL -> GetCovariance() << std::endl;
    /// std::cout << "Cov NRNR: " << fNRNR -> GetCovariance() << std::endl;

    
  }

  double nominal;

  for(int iEvent = 0; iEvent < entries; ++iEvent){
    
    tree -> GetEntry(iEvent);
    
    double Nges = 0.0;

    if(fFittingMode == "3D")
      Nges = valN0+valNL+valNR;
    else
      Nges = valN0+valNL;


    double f0Error = 0.0;
    double fLError = 0.0;
    double fRError = 0.0;

    if(fFittingMode == "3D"){

      f0Error = fi_error (&fErrormatrix[0][0], 3, valN0, Nges, 0);
      fLError = fi_error (&fErrormatrix[0][0], 3, valNL, Nges, 1);
      fRError = fabs(fi_error (&fErrormatrix[0][0], 3, valNR, Nges, 2));
      
    }
    else{

      //      std::cout << valN0 << "\t" << valNL << "\t" << Nges << "\t" << uncN0 << "\t" << uncNL << std::endl;

      double termN0 = valN0/Nges/Nges*uncNL;
      double termNL = valNL/Nges/Nges*uncN0;

      double unc_factor = sqrt(termN0*termN0 + termNL*termNL + 2.0*termN0*termNL);

      //      std::cout << unc_factor << "\t" << termN0 << "\t" << termNL << std::endl;

      f0Error = unc_factor;
      fLError = unc_factor;

    }
    
    fFractionAnalysis[iSyst][0][0].Fill(valN0/Nges);
    fFractionAnalysis[iSyst][1][0].Fill(valNL/Nges);
    
    if(fFittingMode == "3D")
      fFractionAnalysis[iSyst][2][0].Fill(valNR/Nges);
    
    fFractionAnalysis[iSyst][0][1].Fill(f0Error);
    fFractionAnalysis[iSyst][1][1].Fill(fLError);

    if(fFittingMode == "3D")
      fFractionAnalysis[iSyst][2][1].Fill(fRError);
    
    fFractionAnalysis[iSyst][0][2].Fill((valN0/Nges - nomF0)/f0Error);
    fFractionAnalysis[iSyst][1][2].Fill((valNL/Nges - nomFL)/fLError);
    
    if(fFittingMode == "3D")
      fFractionAnalysis[iSyst][2][2].Fill((valNR/Nges - nomFR)/fRError);
    
  }

  int nFrac = fFractionNames.size();

  if(mode != "Datafit"){

    double helpN0 = N0hist.GetMean();
    double helpNL = NLhist.GetMean();
    double helpNR = NRhist.GetMean();
    
    double helpNges = helpN0+helpNL+helpNR;

    double helpF0    = helpN0/helpNges;
    double helpFL    = helpNL/helpNges;
    double helpFR    = helpNR/helpNges;
    
    double helpF0err = fi_error(&fErrormatrix[0][0], 3, helpN0, helpNges, 0);
    double helpFLerr = fi_error(&fErrormatrix[0][0], 3, helpNL, helpNges, 1);
    double helpFRerr = fi_error(&fErrormatrix[0][0], 3, helpNR, helpNges, 2);
    
    if(fFittingMode != "3D"){

      double termN0 = helpN0/helpNges/helpNges*NLhist.GetMeanError();
      double termNL = helpNL/helpNges/helpNges*N0hist.GetMeanError();

      double unc_factor = sqrt(termN0*termN0 + termNL*termNL  + 2.0*termN0*termNL);

      //      std::cout << unc_factor << "\t" << termN0 << "\t" << termNL << std::endl;

      helpF0err = unc_factor;
      helpFLerr = unc_factor;

    }

    if(fFittingMode == "3D"){

       //================= mjk: Change DeltaF/F to DetaF
      fFractionGraphs[0][0].SetPoint(iSyst,      fFraction[iSyst][0], helpF0);
      fFractionGraphs[0][0].SetPointError(iSyst,                 0.0, helpF0err/sqrt(entries));
      

      fFractionGraphs[0][1].SetPoint(iSyst,      fFraction[iSyst][0], (fFraction[iSyst][0] - helpF0));
      fFractionGraphs[0][1].SetPointError(iSyst, 0.0, helpF0err/sqrt(entries));
      
      
      //std::cout << "###FL = " << fFraction[iSyst][1] << "\t" << helpFL << "\t" << helpFLerr << "\t" << entries << std::endl;

      fFractionGraphs[1][0].SetPoint(iSyst,      fFraction[iSyst][1], helpFL);
      fFractionGraphs[1][0].SetPointError(iSyst,                 0.0, helpFLerr/sqrt(entries));

      fFractionGraphs[1][1].SetPoint(iSyst,      fFraction[iSyst][1], (fFraction[iSyst][1] - helpFL));
      fFractionGraphs[1][1].SetPointError(iSyst, 0.0, helpFLerr/sqrt(entries));

      if(fFraction[iSyst][2] != 0){
	
	fFractionGraphs[2][0].SetPoint(iSyst,      fFraction[iSyst][2], helpFR);
	fFractionGraphs[2][0].SetPointError(iSyst,                 0.0, helpFRerr/sqrt(entries));
	
 
 //  fFractionGraphs[2][1].SetPoint(iSyst,      fFraction[iSyst][2], (fFraction[iSyst][2] - helpFR)/fFraction[iSyst][2]);
	// fFractionGraphs[2][1].SetPointError(iSyst, 0.0, helpFRerr/sqrt(entries)/fFraction[iSyst][2]);

  fFractionGraphs[2][1].SetPoint(iSyst,      fFraction[iSyst][2], (fFraction[iSyst][2] - helpFR));
  fFractionGraphs[2][1].SetPointError(iSyst, 0.0, helpFRerr/sqrt(entries));

  //=================
  
      }
      else{ // fFraction[iSyst][2] == 0
	
	fFractionGraphs[2][0].SetPoint(iSyst,      fFraction[iSyst][2], helpFR);
	fFractionGraphs[2][0].SetPointError(iSyst, 0.0, helpFRerr/sqrt(entries));
	
	//================== mjk
  // fFractionGraphs[2][1].SetPoint(iSyst,      fFraction[iSyst][2], (fFraction[iSyst][2] - helpFR)/fFraction[iSyst][2]);
	// fFractionGraphs[2][1].SetPointError(iSyst, 0.0, 0.0);
  fFractionGraphs[2][1].SetPoint(iSyst,      fFraction[iSyst][2], (fFraction[iSyst][2] - helpFR));
  //fFractionGraphs[2][1].SetPointError(iSyst, 0.0, helpFRerr/sqrt(entries)); //FIXME
  fFractionGraphs[2][1].SetPointError(iSyst, 0.0, 0); //test

  //==================
	
      }
      
    }
    else{ // if not 3D

      fFractionGraphs[0][0].SetPoint(iSyst,      fFraction[iSyst][0], helpF0);
      fFractionGraphs[0][0].SetPointError(iSyst,                 0.0, helpF0err);

      fFractionGraphs[0][1].SetPoint(iSyst,      fFraction[iSyst][0], (fFraction[iSyst][0] - helpF0));
      fFractionGraphs[0][1].SetPointError(iSyst, 0.0, helpF0err/fFraction[iSyst][0]);

      //    std::cout << "FL = " << fFraction[iSyst][1] << "\t" << helpFL << "\t" << helpFLerr << "\t" << entries << std::endl;

      fFractionGraphs[1][0].SetPoint(iSyst,      fFraction[iSyst][1], helpFL);
      fFractionGraphs[1][0].SetPointError(iSyst,                 0.0, helpFLerr);

      fFractionGraphs[1][1].SetPoint(iSyst,      fFraction[iSyst][1], (fFraction[iSyst][1] - helpFL));
      fFractionGraphs[1][1].SetPointError(iSyst, 0.0, helpFLerr);
      
    }
    
  }
  
  
}

void ExternalSystematics::AnalyseParameter(int iParam, int iSyst, std::string OutputFolder)
{

  //  std::cout << "Analyse parameter " << iParam << std::endl;

  InitializeHistogramsAnalysis(fSystList.size(), iParam);

  //  std::cout << fSystList[iSyst].c_str() << std::endl;
  
  TFile *fFile = new TFile(fSystList[iSyst].c_str(), "READ");
  fFile -> cd();

  TTree *tree = (TTree*) fFile -> Get("EnsembleTree");
  
  double val, unc, nom;
  
  std::string Name = fParameterNames[iParam];
  
  tree -> SetBranchAddress(Name.c_str(),          &val);
  tree -> SetBranchAddress((Name+"_err").c_str(), &unc);
  tree -> SetBranchAddress((Name+"_nom").c_str(), &nom);
  
  int entries = tree -> GetEntries();
  
  double nominal;
  double nominal1;

  for(int iEvent = 0; iEvent < entries; ++iEvent){
    
    tree -> GetEntry(iEvent);
    
    //    std::cout << nom << "\t" << val << "\t" << iEvent << std::endl;

    fHistoAnalysis[iSyst][0].Fill(val);
    fHistoAnalysis[iSyst][1].Fill(unc);
    fHistoAnalysis[iSyst][2].Fill((val-nom)/unc);
    
    if(iParam > fFractionNames.size()-1)
      nominal = fFraction[iSyst][0];
    else
      nominal = nom;

    nominal1 = nom;

    //    std::cout << fHistoAnalysis[iSyst][0].GetMean() << std::endl;

  }

  if((iParam > 2 && fFittingMode == "3D") || (iParam > 1 && fFittingMode == "2D")){
  
    // fill calibration curves
    fGraphs[iParam][0].SetPoint(iSyst,  nominal, fHistoAnalysis[iSyst][0].GetMean()-nominal1);
    fGraphs[iParam][0].SetPointError(iSyst, 0.0, fHistoAnalysis[iSyst][0].GetMeanError());
    
    fGraphs[iParam][1].SetPoint(iSyst,     nominal, (nominal1 - fHistoAnalysis[iSyst][0].GetMean())/nominal1);
    fGraphs[iParam][1].SetPointError(iSyst,    0.0, fHistoAnalysis[iSyst][0].GetMeanError()/nominal1);
    
    // fill pull curves
    fGraphs[iParam][2].SetPoint(iSyst,  nominal, fHistoAnalysis[iSyst][2].GetMean());
    fGraphs[iParam][2].SetPointError(iSyst, 0.0, fHistoAnalysis[iSyst][2].GetMeanError());
    
    // fill RMS curves
    fGraphs[iParam][3].SetPoint(iSyst,  nominal, fHistoAnalysis[iSyst][2].GetRMS());
    fGraphs[iParam][3].SetPointError(iSyst, 0.0, fHistoAnalysis[iSyst][2].GetRMSError());


  }


  else{

    double delta = 0;

    if(iParam == 0) delta = fFraction[3][0]*fLumi*fXsec;
    if(iParam == 1) delta = fFraction[3][1]*fLumi*fXsec;
    if(iParam == 2) delta = fFraction[3][2]*fLumi*fXsec;

    // fill calibration curves
    fGraphs[iParam][0].SetPoint(iSyst,  nominal-delta, fHistoAnalysis[iSyst][0].GetMean() - delta);
    fGraphs[iParam][0].SetPointError(iSyst, 0.0, fHistoAnalysis[iSyst][0].GetMeanError());

    // residual for calicurves
    fGraphs[iParam][1].SetPoint(iSyst,  nominal-delta, (nominal1 - fHistoAnalysis[iSyst][0].GetMean())/nominal1);
    fGraphs[iParam][1].SetPointError(iSyst,    0.0, fHistoAnalysis[iSyst][0].GetMeanError()/nominal1);

    // fill pull curves
    fGraphs[iParam][2].SetPoint(iSyst,  nominal-delta, fHistoAnalysis[iSyst][2].GetMean());
    fGraphs[iParam][2].SetPointError(iSyst, 0.0,       fHistoAnalysis[iSyst][2].GetMeanError());

    // fill RMS curves
    fGraphs[iParam][3].SetPoint(iSyst,  nominal-delta, fHistoAnalysis[iSyst][2].GetRMS());
    fGraphs[iParam][3].SetPointError(iSyst, 0.0,       fHistoAnalysis[iSyst][2].GetRMSError());

  }
  
  
  fFile -> Close();
  
}

//Destructor
ExternalSystematics::~ExternalSystematics()
{

  delete fRandom;
  std::cout<<"fRandom deleted" << std::endl;
  delete gaussRandom;
  std::cout<<"gaussRandom deleted" << std::endl;
  delete fit;
  std::cout<<"fit deleted" << std::endl;
  delete graph;
  std::cout<<"graph deleted" << std::endl;
  delete fitter_prof; 
  std::cout<<"fitter_prof deleted" << std::endl;

}






