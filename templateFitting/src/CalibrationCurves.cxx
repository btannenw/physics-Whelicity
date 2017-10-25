#include "CalibrationCurves.h"

#include "plots.h"

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

CalibrationCurves::CalibrationCurves()
{

  


}

CalibrationCurves::~CalibrationCurves()
{




}

void CalibrationCurves::FillGraphs(std::string tree_name, std::string var_name, std::string input_channel, bool xval, bool is_fraction)
{

  int n = fInputFolderList.size();

  std::cout << n << std::endl;

  fCalibCurve                  = new TGraphErrors(n+1);
  fPullCurve                   = new TGraphErrors(n+1);
  fRMSCurve                    = new TGraphErrors(n+1);
  fPullCombCurve               = new TGraphErrors(n+1);
  fPullProfileCurve            = new TGraphErrors(n+1);
  fRMSProfileCurve             = new TGraphErrors(n+1);
  fPullCombProfileCurve        = new TGraphErrors(n+1);
  fPullCurveGauss              = new TGraphErrors(n+1);
  fRMSCurveGauss               = new TGraphErrors(n+1);
  fPullCombCurveGauss          = new TGraphErrors(n+1);
  fPullProfileCurveGauss       = new TGraphErrors(n+1);
  fRMSProfileCurveGauss        = new TGraphErrors(n+1);
  fPullCombProfileCurveGauss   = new TGraphErrors(n+1);

  fDiffCalibCurve              = new TGraphErrors(n+1);

	//std::cout << "here 0 " << std::endl;

  for(int iFile = 0; iFile < n; ++iFile)
    ReadInputValues(tree_name, iFile, xval, var_name);
  
	//std::cout << "here 1 " << std::endl;
  
  std::string xlabel      = "";
  std::string ylabel      = "";
  std::string xtitle      = "";
  std::string xtitle_pull = "";
  std::string xtitle_rms  = "";


  if(var_name == "F0"){
    ylabel = "F_{0}";      xtitle = input_channel + "_calibration_F0";     xtitle_pull = input_channel + "_pull_F0";     xtitle_rms = input_channel + "_rms_F0";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "FL"){
    ylabel = "F_{L}";      xtitle = input_channel + "_calibration_FL";     xtitle_pull = input_channel + "_pull_FL";     xtitle_rms = input_channel + "_rms_FL";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "FR"){
    ylabel = "F_{R}";      xtitle = input_channel + "_calibration_FR";     xtitle_pull = input_channel + "_pull_FR";     xtitle_rms = input_channel + "_rms_FR";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "N0"){
    ylabel = "N_{0}";      xtitle = input_channel + "_calibration_N0";     xtitle_pull = input_channel + "_pull_N0";     xtitle_rms = input_channel + "_rms_N0";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "NL"){
    ylabel = "N_{L}";      xtitle = input_channel + "_calibration_NL";     xtitle_pull = input_channel + "_pull_NL";     xtitle_rms = input_channel + "_rms_NL";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "NR"){
    ylabel = "N_{R}";      xtitle = input_channel + "_calibration_NR";     xtitle_pull = input_channel + "_pull_NR";     xtitle_rms = input_channel + "_rms_NR";

    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "Wjets"){
    ylabel = "N(W+jets)";  xtitle = input_channel + "_calibration_Wjets";  xtitle_pull = input_channel + "_pull_Wjets";  xtitle_rms = input_channel + "_rms_Wjets";

    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "QCD"){
    ylabel = "N(QCD)";    xtitle = input_channel + "_calibration_QCD";    xtitle_pull = input_channel + "_pull_QCD";    xtitle_rms = input_channel + "_rms_QCD";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "RemBkg"){ 
    ylabel = "N(RemBkg)"; xtitle = input_channel + "_calibration_RemBkg"; xtitle_pull = input_channel + "_pull_RemBkg"; xtitle_rms = input_channel + "_rms_RemBkg";
  
    if(!is_fraction) xlabel = "Nui";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui1"){
    ylabel = "Nui1";       xtitle = input_channel + "_calibration_Nui1";   xtitle_pull = input_channel + "_pull_Nui1";   xtitle_rms = input_channel + "_rms_Nui1";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui2"){
    ylabel = "Nui2";       xtitle = input_channel + "_calibration_Nui2";   xtitle_pull = input_channel + "_pull_Nui2";   xtitle_rms = input_channel + "_rms_Nui2";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }
  else if(var_name == "Nui3"){
    ylabel = "Nui3";       xtitle = input_channel + "_calibration_Nui3";   xtitle_pull = input_channel + "_pull_Nui3";   xtitle_rms = input_channel + "_rms_Nui3";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui4"){
    ylabel = "Nui4";       xtitle = input_channel + "_calibration_Nui4";   xtitle_pull = input_channel + "_pull_Nui4";   xtitle_rms = input_channel + "_rms_Nui4";
    
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui5"){
    ylabel = "Nui5";       xtitle = input_channel + "_calibration_Nui5";   xtitle_pull = input_channel + "_pull_Nui5";   xtitle_rms = input_channel + "_rms_Nui5";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui6"){
    ylabel = "Nui6";       xtitle = input_channel + "_calibration_Nui6";   xtitle_pull = input_channel + "_pull_Nui6";   xtitle_rms = input_channel + "_rms_Nui6";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui7"){
    ylabel = "Nui7";       xtitle = input_channel + "_calibration_Nui7";   xtitle_pull = input_channel + "_pull_Nui7";   xtitle_rms = input_channel + "_rms_Nui7";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;

  }
  else if(var_name == "Nui8"){
    ylabel = "Nui8";       xtitle = input_channel + "_calibration_Nui8";   xtitle_pull = input_channel + "_pull_Nui8";   xtitle_rms = input_channel + "_rms_Nui8";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }
  else if(var_name == "Nui9"){
    ylabel = "Nui9";       xtitle = input_channel + "_calibration_Nui9";   xtitle_pull = input_channel + "_pull_Nui9";   xtitle_rms = input_channel + "_rms_Nui9";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }
  else if(var_name == "Nui10"){
    ylabel = "Nui10";       xtitle = input_channel + "_calibration_Nui10";   xtitle_pull = input_channel + "_pull_Nui10";   xtitle_rms = input_channel + "_rms_Nui10";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }
  else if(var_name == "Nui11"){
    ylabel = "Nui11";       xtitle = input_channel + "_calibration_Nui11";   xtitle_pull = input_channel + "_pull_Nui11";   xtitle_rms = input_channel + "_rms_Nui11";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }
  else if(var_name == "Nui12"){
    ylabel = "Nui12";       xtitle = input_channel + "_calibration_Nui12";   xtitle_pull = input_channel + "_pull_Nui12";   xtitle_rms = input_channel + "_rms_Nui12";
  
    if(is_fraction) xlabel = "F_{0}";
    else xlabel = ylabel;
    
  }



  double max_y = TMath::MaxElement(n, fCalibCurveVecY);
  double min_y = TMath::MinElement(n, fCalibCurveVecY);

  double max_pull_y = TMath::MaxElement(n, fPullCurve -> GetY());
  double min_pull_y = TMath::MinElement(n, fPullCurve -> GetY());

  double max_rms_y = TMath::MaxElement(n, fRMSCurve -> GetY());
  double min_rms_y = TMath::MinElement(n, fRMSCurve -> GetY());

  double max_x = TMath::MaxElement(n, fCalibCurveVecX);
  double min_x = TMath::MinElement(n, fCalibCurveVecX);

  double diff_x = (max_x - min_x)*0.1;
  double diff_y = (max_y - min_y)*0.1;
 
  max_x += diff_x;
  max_y += diff_y*4;
  min_x -= diff_x;
  min_y -= diff_y*4;

  //  if(min_y < 0) min_y = min_y*1.25;
  // else min_y = min_y*0.75;

  //  if(min_x < 0) min_x = min_x*1.25;
  // else min_x = min_x*0.75;

  if(min_pull_y < 0) min_pull_y = min_pull_y*1.5;
  else min_pull_y = min_pull_y*0.5;

  if(max_pull_y < 0) max_pull_y = max_pull_y*0.5;
  else max_pull_y = max_pull_y*1.5;
  
  if(min_rms_y < 0) min_rms_y = min_rms_y*1.5;
  else min_rms_y = min_rms_y*0.5;

  if(max_rms_y < 0) max_rms_y = max_rms_y*0.5;
  else max_rms_y = max_rms_y*1.5;

  std::cout << min_x << "\t" << max_x << std::endl;
  std::cout << min_y << "\t" << max_y << std::endl;

///________________________PLOT CALIBRATION__________________________

  if (var_name == "N0" || var_name == "NL" || var_name == "NR") {
  	PlotCalDistribution(fCalibCurve, fDiffCalibCurve, xlabel, ylabel, xtitle, min_x, max_x, min_y, max_y, input_channel, 1, fOutputFolder);
  }
  else if (var_name == "F0" || var_name == "FL" || var_name == "FR") {
	PlotCalDistribution(fCalibCurve, fDiffCalibCurve, xlabel, ylabel, xtitle, min_x, max_x, min_y, max_y, input_channel, 0, fOutputFolder);
  }
  else {
	PlotCalDistribution(fCalibCurve, fDiffCalibCurve, xlabel, ylabel, xtitle, min_x, max_x, min_y, max_y, input_channel, 2, fOutputFolder);
  }

  //  std::cout << max_pull_y << "\t" << min_pull_y << std::endl;


///_____________________________PLOT PULL___________________________

  PlotCalPullDistribution (fPullCurve, xlabel, ylabel, xtitle_pull, min_x, max_x, -0.12, 0.12, input_channel, "Mean - w/o profiling", fOutputFolder);
  
  //PlotCalPullDistribution (fRMSCurve, xlabel, xtitle_rms, min_x, max_x, min_rms_y, max_rms_y, input_channel, "RMS", fOutputFolder);
  PlotCalPullDistribution (fRMSCurve, xlabel, ylabel, xtitle_rms, min_x, max_x, 0.5, 1.1, input_channel, "Sigma - w/o profiling", fOutputFolder);  

  PlotCalPullDistribution (fPullCombCurve, xlabel, ylabel, xtitle_pull+"_comb", min_x, max_x, -1.2, 1.2, input_channel, "Mean + Sigma - w/o profiling", fOutputFolder);

  //PlotCalPullDistribution (fPullProfileCurve, xlabel, ylabel, xtitle_pull+"_profile", min_x, max_x, min_pull_y, max_pull_y, input_channel, "Mean Profiling", fOutputFolder);
  PlotCalPullDistribution (fPullProfileCurve, xlabel, ylabel, xtitle_pull+"_profile", min_x, max_x, -0.12, 0.12, input_channel, "Mean", fOutputFolder);

  PlotCalPullDistribution (fRMSProfileCurve, xlabel, ylabel, xtitle_rms+"_profile", min_x, max_x, 0.5, 1.1, input_channel, "Sigma", fOutputFolder);
  
  PlotCalPullDistribution (fPullCombProfileCurve, xlabel, ylabel, xtitle_pull+"_comb_profile", min_x, max_x, -1.2, 1.2, input_channel, "Mean + Sigma", fOutputFolder);

///GAUSS
/*
  PlotCalPullDistribution (fPullCurveGauss, xlabel, ylabel, xtitle_pull+"_gauss", min_x, max_x, -0.15, 0.15, input_channel, "Mean (Gauss)", fOutputFolder);
  
  //PlotCalPullDistribution (fRMSCurve, xlabel, xtitle_rms, min_x, max_x, min_rms_y, max_rms_y, input_channel, "RMS", fOutputFolder);
  PlotCalPullDistribution (fRMSCurveGauss, xlabel, ylabel, xtitle_rms+"_gauss", min_x, max_x, 0.6, 1.2, input_channel, "RMS (Gauss)", fOutputFolder);  

  PlotCalPullDistribution (fPullCombCurveGauss, xlabel, ylabel, xtitle_pull+"_comb_gauss", min_x, max_x, -1.2, 1.2, input_channel, "Mean + RMS Error (Gauss)", fOutputFolder);

  PlotCalPullDistribution (fPullProfileCurveGauss, xlabel, ylabel, xtitle_pull+"_profile_gauss", min_x, max_x, -0.15, 0.15, input_channel, "Mean - Profiling (Gauss)", fOutputFolder);

  PlotCalPullDistribution (fRMSProfileCurveGauss, xlabel, ylabel, xtitle_rms+"_profile_gauss", min_x, max_x, 0.6, 1.2, input_channel, "RMS - Profiling (Gauss)", fOutputFolder);
  
  PlotCalPullDistribution (fPullCombProfileCurveGauss, xlabel, ylabel, xtitle_pull+"_comb_profile_gauss", min_x, max_x, -1.2, 1.2, input_channel, "Mean + RMS Error - Profiling (Gauss)", fOutputFolder);
*/
}








void CalibrationCurves::ReadInputValues(std::string tree_name, int counter, bool xval, std::string var_name)
{
  
  std::string input_path = fInputFolderList[counter]+"/CalibrationTree.root";

  //  std::cout << input_path.c_str() << "\t" << tree_name.c_str() << std::endl;

	//std::cout << "here func 0 " << std::endl;

  TFile* file = new TFile(input_path.c_str(), "open");

	//std::cout << "here func 0b " << std::endl;

  TTree *tree = (TTree*) file -> Get(tree_name.c_str());

	//std::cout << "here func 1 " << std::endl;

  Double_t        Nominal;
  Double_t        Mean;
  Double_t        MeanErr;
  Double_t        GaussMean;
  Double_t        GaussMeanErr;
  Double_t        PullMean;
  Double_t        PullMeanErr;
  Double_t        PullRMS;
  Double_t        PullRMSErr;
  Double_t        GaussPullMean;
  Double_t        GaussPullMeanErr;
  Double_t        GaussPullRMS;
  Double_t        GaussPullRMSErr;
  Double_t        MeanProf;
  Double_t        MeanProfErr;
  Double_t        PullMeanProf;
  Double_t        PullMeanProfErr;
  Double_t        PullRMSProf;
  Double_t        PullRMSProfErr;
  Double_t        GaussPullMeanProf;
  Double_t        GaussPullMeanProfErr;
  Double_t        GaussPullRMSProf;
  Double_t        GaussPullRMSProfErr;

  tree -> SetBranchAddress("Nominal",              &Nominal);
  tree -> SetBranchAddress("Mean",                 &Mean);
  tree -> SetBranchAddress("MeanErr",              &MeanErr);
  tree -> SetBranchAddress("GaussMean",            &Mean);
  tree -> SetBranchAddress("GaussMeanErr",         &MeanErr);
  tree -> SetBranchAddress("PullMean",             &PullMean);
  tree -> SetBranchAddress("PullMeanErr",          &PullMeanErr);
  tree -> SetBranchAddress("PullRMS",              &PullRMS);
  tree -> SetBranchAddress("PullRMSErr",           &PullRMSErr);
  tree -> SetBranchAddress("GaussPullMean",        &GaussPullMean);
  tree -> SetBranchAddress("GaussPullMeanErr",     &GaussPullMeanErr);
  tree -> SetBranchAddress("GaussPullRMS",         &GaussPullRMS);
  tree -> SetBranchAddress("GaussPullRMSErr",      &GaussPullRMSErr);
  tree -> SetBranchAddress("MeanProf",             &MeanProf);
  tree -> SetBranchAddress("MeanProfErr",          &MeanProfErr);
  tree -> SetBranchAddress("PullMeanProf",         &PullMeanProf);
  tree -> SetBranchAddress("PullMeanProfErr",      &PullMeanProfErr);
  tree -> SetBranchAddress("PullRMSProf",          &PullRMSProf);
  tree -> SetBranchAddress("PullRMSProfErr",       &PullRMSProfErr);
  tree -> SetBranchAddress("GaussPullMeanProf",    &GaussPullMeanProf);
  tree -> SetBranchAddress("GaussPullMeanProfErr", &GaussPullMeanProfErr);
  tree -> SetBranchAddress("GaussPullRMSProf",     &GaussPullRMSProf);
  tree -> SetBranchAddress("GaussPullRMSProfErr",  &GaussPullRMSProfErr);

  // Get the Event
  tree -> GetEntry(0);



  //  std::cout << counter << "\t" << Nominal << "\t" << Mean << std::endl;

  double   x = 0.0;

  if(xval) x = fXValueList[counter];
  else     x = Nominal;
  
  double add_mean = PullMean;
  double add_rms = PullRMS;
  double add_mean_prof = PullMeanProf;
  double add_rms_prof = PullRMSProf;

  //  std::cout << input_path.c_str() << std::endl;
  //  std::cout << counter << "\t" << x << "\t" << Mean << std::endl;

  fCalibCurve                 -> SetPoint(counter+1, x, Mean);
  fCalibCurve                 -> SetPointError(counter+1, 0, fabs(MeanErr));
  fPullCurve                  -> SetPoint(counter+1, x, PullMean);
  fPullCurve                  -> SetPointError(counter+1, 0, PullMeanErr);
  fRMSCurve                   -> SetPoint(counter+1, x, PullRMS);
  fRMSCurve                   -> SetPointError(counter+1, 0, PullRMSErr);
  fPullCombCurve              -> SetPoint(counter+1, x, PullMean);
  fPullCombCurve              -> SetPointError(counter+1, 0, PullRMS);
  fPullProfileCurve           -> SetPoint(counter+1, x, PullMeanProf);
  fPullProfileCurve           -> SetPointError(counter+1, 0, PullMeanProfErr);
  fRMSProfileCurve            -> SetPoint(counter+1, x, PullRMSProf);
  fRMSProfileCurve            -> SetPointError(counter+1, 0, PullRMSProfErr);
  fPullCombProfileCurve       -> SetPoint(counter+1, x, PullMeanProf);
  fPullCombProfileCurve       -> SetPointError(counter+1, 0, PullRMSProf);
  fPullCurveGauss             -> SetPoint(counter+1, x, GaussPullMean);
  fPullCurveGauss             -> SetPointError(counter+1, 0, GaussPullMeanErr);
  fRMSCurveGauss              -> SetPoint(counter+1, x, GaussPullRMS);
  fRMSCurveGauss              -> SetPointError(counter+1, 0, GaussPullRMSErr);
  fPullCombCurveGauss         -> SetPoint(counter+1, x, GaussPullMean);
  fPullCombCurveGauss         -> SetPointError(counter+1, 0, GaussPullRMS);
  fPullProfileCurveGauss      -> SetPoint(counter+1, x, GaussPullMeanProf);
  fPullProfileCurveGauss      -> SetPointError(counter+1, 0, GaussPullMeanProfErr);
  fRMSProfileCurveGauss       -> SetPoint(counter+1, x, GaussPullRMSProf);
  fRMSProfileCurveGauss       -> SetPointError(counter+1, 0, GaussPullRMSProfErr);
  fPullCombProfileCurveGauss  -> SetPoint(counter+1, x, GaussPullMeanProf);
  fPullCombProfileCurveGauss  -> SetPointError(counter+1, 0, GaussPullRMSProf);



  fCalibCurveVecX[counter]              = x;
  fCalibCurveVecY[counter]              = Mean;
  fPullCurveVecX[counter]               = x;
  fPullCurveVecY[counter]               = PullMean;
  //fPullCombCurveVecX[counter]           = x;
  //fPullCombCurveVecY[counter]           = PullMean;
  fRMSCurveVecX[counter]                = x;
  fRMSCurveVecY[counter]                = PullRMS;
  fPullProfileCurveVecX[counter]        = x;
  fPullProfileCurveVecY[counter]        = PullMeanProfErr;
  fRMSProfileCurveVecX[counter]         = x;
  fRMSProfileCurveVecY[counter]         = PullRMSProfErr;
  //fPullCombProfileCurveVecX[counter]    = x;
  //fPullCombProfileCurveVecY[counter]    = PullMeanProfErr;



  if (var_name == "F0" || var_name == "FL" || var_name == "FR") {
    fDiffCalibCurve   -> SetPoint(counter+1, x, (Mean-Nominal));
  }
  else if (Nominal == 0) {
    fDiffCalibCurve   -> SetPoint(counter+1, x, (Mean-Nominal));
  }
  else {
    fDiffCalibCurve   -> SetPoint(counter+1, x, (Mean-Nominal)/Nominal);
  }

  if (var_name == "F0" || var_name == "FL" || var_name == "FR") {
    fDiffCalibCurve   -> SetPointError(counter+1, 0, fabs(MeanErr));
  }
  else if (Nominal == 0) {
    fDiffCalibCurve   -> SetPointError(counter+1, 0, fabs(MeanErr));
  }
  else {
    fDiffCalibCurve   -> SetPointError(counter+1, 0, fabs(MeanErr)/Nominal);
  }

  //  std::cout << (Mean-Nominal)/Nominal << "\t" << MeanErr/Nominal << std::endl;

  if(counter == fInputFolderList.size()){

    fCalibCurve                 -> RemovePoint(0);
    fPullCurve                  -> RemovePoint(0);
    fRMSCurve                   -> RemovePoint(0);
    fPullCombCurve              -> RemovePoint(0);
    fPullProfileCurve           -> RemovePoint(0); 
    fRMSProfileCurve            -> RemovePoint(0);
    fPullCombProfileCurve       -> RemovePoint(0);
    fPullCurveGauss             -> RemovePoint(0);
    fRMSCurveGauss              -> RemovePoint(0);
    fPullCombCurveGauss         -> RemovePoint(0);
    fPullProfileCurveGauss      -> RemovePoint(0); 
    fRMSProfileCurveGauss       -> RemovePoint(0);
    fPullCombProfileCurveGauss  -> RemovePoint(0);
  }
    

  delete tree;

  delete file;

}
