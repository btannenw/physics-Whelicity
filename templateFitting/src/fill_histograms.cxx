
#include "fill_histograms.h"
#include "OutputTreeReader.h"
#include "helicity.h"
#include "plots.h"

#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include <sys/stat.h>
#include <dirent.h>

#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TError.h"

using namespace std;


histograms::histograms()
{

  fFileNames.clear();

}

histograms::~histograms()
{

}

///_______________________FILL HISTOGRAMS__________________________
void histograms::FillHistograms(std::string input_channel, int num_of_nui_parameters)
{

  gErrorIgnoreLevel = kError;

  this -> FindRanges(num_of_nui_parameters);
  
  //Define vectors needed for correlation and errormatrix
  std::vector<double> vec_F0, vec_FL, vec_FR; 
  std::vector<double> vec_N0, vec_NL, vec_NR;
  std::vector<double> vec_N0_prof_err, vec_NL_prof_err, vec_NR_prof_err;

  double xsec = 96.0;    //pb  //better 89.45 ... check later
  double lumi = 4655.74; //i_pb

  double F0_true, FL_true, FR_true;

  for(int param = 0; param < fParameter.size(); ++param){
    
    std::string label = fParameter[param].label;
    
    //    std::cout << fParameter[param].Low_prof_unc << "\t" << fParameter[param].Up_prof_unc << std::endl;
    
    if (fParameter[param].Low > 0.0) {
      fParameter[param].hist_param     = new TH1D(label.c_str(),                "", 100, fParameter[param].Low*0.8,          fParameter[param].Up*1.2);
    }
    else {
      fParameter[param].hist_param     = new TH1D(label.c_str(),                "", 100, fParameter[param].Low*1.2,          fParameter[param].Up*1.2);
    }
    fParameter[param].hist_unc       = new TH1D((label+"_unc").c_str(),       "", 100, fParameter[param].Low_unc*0.8,      fParameter[param].Up_unc*1.2);
    fParameter[param].hist_prof_unc  = new TH1D((label+"+unc_prof").c_str(),  "", 100, fParameter[param].Low_prof_unc*0.8, fParameter[param].Up_prof_unc*1.2);
    fParameter[param].hist_pull      = new TH1D((label+"_pull").c_str(),      "", 100, -7.0, 7.0);
    fParameter[param].hist_pull_prof = new TH1D((label+"_pull_prof").c_str(), "", 100, -7.0, 7.0);
    
  }

  for(int param1 = 0; param1 < fParameter.size(); ++param1){
    
    for(int param2 = 0; param2 < fParameter.size(); ++param2){
      
      fParameter[param1].hist2D.push_back(TH2D("", "", 100, fParameter[param1].Low*0.8, fParameter[param1].Up*1.2, 100, fParameter[param2].Low*0.8, fParameter[param2].Up*1.2));

    }
    
  }
  

  std::cout << "Fill Histograms" << std::endl;

  for(int k = 0; k < fTotalEntries; ++k){

    fInputTree -> fChain -> GetEntry(k);

    for(int param = 0; param < fParameter.size(); ++param){
      
      double central  = 0.0;
      double unc      = 0.0;
      double prof_unc = 0.0;
      double nominal  = 0.0;
      
      if(fParameter[param].name == "F0"){
	
	central = fInputTree -> F0; unc = fInputTree -> F0_err; nominal = fInputTree -> F0_nom;
	vec_F0.push_back(central);
	fParameter[param].val_vec.push_back(central);
	
	if(param == 0)
	  F0_true = fInputTree -> F0_nom;

      }
      else if(fParameter[param].name == "FL"){

	central = fInputTree -> FL; unc = fInputTree -> FL_err; nominal = fInputTree -> FL_nom;
	vec_FL.push_back(central);
	fParameter[param].val_vec.push_back(central);

	if(param == 1)
          FL_true = fInputTree -> FL_nom;

      }
      else if(fParameter[param].name == "FR"){

	central = fInputTree -> FR; unc = fInputTree -> FR_err; nominal = fInputTree -> FR_nom;
	vec_FR.push_back(central);
	fParameter[param].val_vec.push_back(central);

	if(param == 2)
          FR_true = fInputTree -> FR_nom;

      }
      else if(fParameter[param].name == "N0"){

	central = fInputTree -> N0; unc = fInputTree -> N0_err; prof_unc = fInputTree -> N0_prof_err; nominal = fInputTree -> N0_nom;
	vec_N0.push_back(central);
	vec_N0_prof_err.push_back(prof_unc);

	fParameter[param].val_vec.push_back(central);

      }
      else if(fParameter[param].name == "NL"){

	central = fInputTree -> NL; unc = fInputTree -> NL_err; prof_unc = fInputTree -> NL_prof_err; nominal = fInputTree -> NL_nom;
  	vec_NL.push_back(central);
	vec_NL_prof_err.push_back(prof_unc);

	fParameter[param].val_vec.push_back(central);

      }
      else if(fParameter[param].name == "NR"){

	central = fInputTree -> NR; unc = fInputTree -> NR_err; prof_unc = fInputTree -> NR_prof_err; nominal = fInputTree -> NR_nom;
	vec_NR.push_back(central);
	vec_NR_prof_err.push_back(prof_unc);

	fParameter[param].val_vec.push_back(central);

      }
      else if(fParameter[param].name == "Nges"){

	central = fInputTree -> Nges; nominal = fInputTree->N0_nom + fInputTree->NL_nom + fInputTree->NR_nom;

	fParameter[param].val_vec.push_back(central);

      }
      else if(fParameter[param].name == "Wjets"){

        central = fInputTree -> Wjets; unc = fInputTree -> Wjets_err; prof_unc = fInputTree -> Wjets_prof_err; nominal = fInputTree -> Wjets_nom;

	fParameter[param].val_vec.push_back(central);

      }
      else if(fParameter[param].name == "QCD"){

	central = fInputTree -> QCD; unc = fInputTree -> QCD_err; prof_unc = fInputTree -> QCD_prof_err; nominal = fInputTree -> QCD_nom;

	fParameter[param].val_vec.push_back(central);
	
      }
      else if(fParameter[param].name == "RemBkg"){

        central = fInputTree -> RemBkg; unc = fInputTree -> RemBkg_err; prof_unc = fInputTree -> RemBkg_prof_err; nominal = fInputTree -> RemBkg_nom;

	fParameter[param].val_vec.push_back(central);

      }
      else {
	
	//	std::cout << "Find reason for nui problems   " << i << "\t" << param << "\t" << fParameter[param].name << std::endl;
	
	if (fParameter[param].name == "nui_0") {
	  central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0; 
	  fParameter[param].val_vec.push_back(central);
	}
	if (fParameter[param].name == "nui_1") {
	  central = fInputTree -> nui_1; unc = fInputTree -> nui_err_1; prof_unc = fInputTree -> nui_prof_err_1; nominal = fInputTree -> nui_nom_1; 
	  fParameter[param].val_vec.push_back(central);
	}
	if (fParameter[param].name == "nui_2") {
	  central = fInputTree -> nui_2; unc = fInputTree -> nui_err_2; prof_unc = fInputTree -> nui_prof_err_2; nominal = fInputTree -> nui_nom_2; 
	  fParameter[param].val_vec.push_back(central);
	}
	if (fParameter[param].name == "nui_3") {
	  central = fInputTree -> nui_3; unc = fInputTree -> nui_err_3; prof_unc = fInputTree -> nui_prof_err_3; nominal = fInputTree -> nui_nom_3; 
	  fParameter[param].val_vec.push_back(central);
	}
	if (fParameter[param].name == "nui_4") {
	  central = fInputTree -> nui_4; unc = fInputTree -> nui_err_4; prof_unc = fInputTree -> nui_prof_err_4; nominal = fInputTree -> nui_nom_4; 
	  fParameter[param].val_vec.push_back(central);
	}
	if (fParameter[param].name == "nui_5") {
	  central = fInputTree -> nui_5; unc = fInputTree -> nui_err_5; prof_unc = fInputTree -> nui_prof_err_5; nominal = fInputTree -> nui_nom_5; 
	  fParameter[param].val_vec.push_back(central);
	}		
	if (fParameter[param].name == "nui_6") {
	  central = fInputTree -> nui_6; unc = fInputTree -> nui_err_6; prof_unc = fInputTree -> nui_prof_err_6; nominal = fInputTree -> nui_nom_6; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_7") {
	  central = fInputTree -> nui_7; unc = fInputTree -> nui_err_7; prof_unc = fInputTree -> nui_prof_err_7; nominal = fInputTree -> nui_nom_7; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_8") {
	  central = fInputTree -> nui_8; unc = fInputTree -> nui_err_8; prof_unc = fInputTree -> nui_prof_err_8; nominal = fInputTree -> nui_nom_8; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_9") {
	  central = fInputTree -> nui_9; unc = fInputTree -> nui_err_9; prof_unc = fInputTree -> nui_prof_err_9; nominal = fInputTree -> nui_nom_9; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_10") {
	  central = fInputTree -> nui_10; unc = fInputTree -> nui_err_10; prof_unc = fInputTree -> nui_prof_err_10; nominal = fInputTree -> nui_nom_10; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_11") {
	  central = fInputTree -> nui_11; unc = fInputTree -> nui_err_11; prof_unc = fInputTree -> nui_prof_err_11; nominal = fInputTree -> nui_nom_11; 
	  fParameter[param].val_vec.push_back(central);
	}	
	if (fParameter[param].name == "nui_12") {
	  central = fInputTree -> nui_12; unc = fInputTree -> nui_err_12; prof_unc = fInputTree -> nui_prof_err_12; nominal = fInputTree -> nui_nom_12; 
	  fParameter[param].val_vec.push_back(central);
	}	



	//      	}
	
      }
      
      
      if (fParameter[param].name == "Nges") {
	fParameter[param].hist_param     -> Fill(central);
      }
      else if(fParameter[param].name == "F0" || fParameter[param].name == "FL" || fParameter[param].name == "FR" ) {

	fParameter[param].hist_param     -> Fill(central);      
        fParameter[param].hist_unc       -> Fill(unc);
        fParameter[param].hist_pull      -> Fill((central - nominal)/unc);

      }
      else {
	fParameter[param].hist_param     -> Fill(central);      
        fParameter[param].hist_unc       -> Fill(unc);
        fParameter[param].hist_prof_unc  -> Fill(prof_unc); 
        fParameter[param].hist_pull      -> Fill((central - nominal)/unc);
        fParameter[param].hist_pull_prof -> Fill((central - nominal)/prof_unc);

      }
    }
  }


  // fill 2D distributions

  for(int param1 = 0; param1 < fParameter.size(); ++param1){
    
    for(int param2 = 0; param2 < fParameter.size(); ++param2){
      
      for(int k = 0; k < fParameter[param1].val_vec.size(); ++k){
	
	double x = fParameter[param1].val_vec[k];
	double y = fParameter[param2].val_vec[k];
	
	fParameter[param1].hist2D[param2].Fill(x,y); 
	
      }
      
    }
    
  }
  
  
  std::cout << "Print out Plots" << std::endl;


  for(int param = 0; param < fParameter.size(); ++param){

    std::stringstream oss;
    oss << fParameter[param].nom;
    
    std::stringstream oss2;
    if (F0_true < 0.1) oss2 << 0.0;
    else oss2 << F0_true;
    
    std::string name  = fParameter[param].name;
    std::string label = fParameter[param].label;
    

    double      mean;
    double      mean_err;
    double      gauss_mean;
    double      gauss_mean_err;
    double      pull_mean;
    double      pull_mean_err;
    double      pull_rms;
    double      pull_rms_err;
    double      gauss_pull_mean;
    double      gauss_pull_mean_err;
    double      gauss_pull_rms;
    double      gauss_pull_rms_err;
    double      mean_prof;
    double      mean_prof_err;
    double      gauss_pull_mean_prof;
    double      gauss_pull_mean_prof_err;
    double      gauss_pull_rms_prof;
    double      gauss_pull_rms_prof_err;

    //Gauss fits to all histograms:
    //so far in MakePlots
    std::vector<double> gauss_fits_param = MakeFits(fParameter[param].hist_param, label, fParameter[param].Low, fParameter[param].Up);
    std::vector<double> gauss_fits_pull = MakeFits(fParameter[param].hist_pull, "Pull distribution for "+label,-2.0, 2.0);


    fParameter[param].mean               = fParameter[param].hist_param     -> GetMean();
    fParameter[param].mean_err           = fParameter[param].hist_param     -> GetMeanError();
    fParameter[param].pull_mean          = fParameter[param].hist_pull      -> GetMean();
    fParameter[param].pull_mean_err      = fParameter[param].hist_pull      -> GetMeanError();
    fParameter[param].pull_rms           = fParameter[param].hist_pull      -> GetRMS();
    fParameter[param].pull_rms_err       = fParameter[param].hist_pull      -> GetRMSError();
    // fParameter[param].mean_prof          = fParameter[param].hist_ -> Get();
    // fParameter[param].mean_prof_err      = fParameter[param].hist_prof_unc  -> Get();
    
    fParameter[param].gauss_mean               = gauss_fits_param[0];
    fParameter[param].gauss_mean_err           = gauss_fits_param[1];
    fParameter[param].gauss_pull_mean          = gauss_fits_pull[0];
    fParameter[param].gauss_pull_mean_err      = gauss_fits_pull[1];
    fParameter[param].gauss_pull_rms           = gauss_fits_pull[2];
    fParameter[param].gauss_pull_rms_err       = gauss_fits_pull[3];
    
    
    
    if(fParameter[param].name != "F0" && fParameter[param].name != "FL" && fParameter[param].name != "FR" ) {
      
      std::vector<double> gauss_fits_pull_prof = MakeFits(fParameter[param].hist_pull_prof, "Pull distribution for "+label,-2.0, 2.0);
      
      fParameter[param].pull_mean_prof     = fParameter[param].hist_pull_prof -> GetMean();
      fParameter[param].pull_mean_prof_err = fParameter[param].hist_pull_prof -> GetMeanError();
      fParameter[param].pull_rms_prof      = fParameter[param].hist_pull_prof -> GetRMS();
      fParameter[param].pull_rms_prof_err  = fParameter[param].hist_pull_prof -> GetRMSError();
      
      fParameter[param].gauss_pull_mean_prof     = gauss_fits_pull_prof[0];
      fParameter[param].gauss_pull_mean_prof_err = gauss_fits_pull_prof[1];
      fParameter[param].gauss_pull_rms_prof      = gauss_fits_pull_prof[2];
      fParameter[param].gauss_pull_rms_prof_err  = gauss_fits_pull_prof[3];
    }
    
    fParameter[param].hist_param;
    fParameter[param].hist_unc;
    fParameter[param].hist_prof_unc;
    fParameter[param].hist_pull;
    fParameter[param].hist_pull_prof;
    
    
    ///________________________Datafit values if available________________
    
    std::vector<double> F_err_vec;    
    
    if (fValidationMode == "Datafit"){	
      
      std::string filename_data = "Datafit_Output/Datafit_Output_" + input_channel + ".root";
      TFile *file_data = GetInputFile(filename_data);
      
      F_err_vec = GetErrorVector(file_data);
      
      // 	cout << "-------------------------------------" << endl;
      // 	cout << "-------------------------------------" << endl;
      // 	cout << "-------------------------------------" << endl;
      // 	cout << F_err_vec[0] << " " << F_err_vec[1] << " " << F_err_vec[2] << endl;
      // 	cout << "-------------------------------------" << endl;
      // 	cout << "-------------------------------------" << endl;
      // 	cout << "-------------------------------------" << endl;
    }
    
    
    ///____________________________Start making plots_______________________
    
    if (fParameter[param].name == "Nges") {
      std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str();
      
      fFileNames.push_back(filename+".png");
      
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".png", input_channel, 1, 0);
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".eps", input_channel, 1, 0);
    }
    else if (fParameter[param].name == "F0" || fParameter[param].name == "FL" || fParameter[param].name == "FR" ) {
      
      // parameter distribution
      std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".png", input_channel, 0, 0);
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".eps", input_channel, 0, 0);
      
      fFileNames.push_back(filename+".png");

      // error distribution
      filename = fOutputFolder+"/plot_"+input_channel+"_parerror_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 7, 1);
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 7, 1);

      fFileNames.push_back(filename+".png");

	//error distributions including datafit error
      if (fValidationMode == "Datafit"){	
	
	filename = fOutputFolder+"/plot_"+input_channel+"_parerror_data_"+name+"_"+oss2.str();
	
	if (name == "F0") {      
	  
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[0], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0);
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[0], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0);
	}
	else if (name == "FL") {      
	  
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[1], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0);
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[1], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0);
	}
	else if (name == "FR") {      
	  
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[2], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0);
	  MakeErrorDataPlots(fParameter[param].hist_unc, F0_true, F_err_vec[2], oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0);
	}
	
	fFileNames.push_back(filename+".png");
      }
      
      // profile error distribution
      /*filename = fOutputFolder+"/plot_"+input_channel+"_proferror_"+name+"_"+oss2.str();
	
      //    if(fParameter[param].Up_prof_unc > 35000) fParameter[param].Up_prof_unc = 35000;

      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+") (prof.)", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0);
      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+") (prof.)", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0);

      fFileNames.push_back(filename+".png");
	*/

      // pull distribution minuit error
      filename = fOutputFolder+"/plot_"+input_channel+"_pull_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // pull distribution profiling error
      /*filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_"+name+"_"+oss2.str();
	
	MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label+" (prof.)", -2.0, 2.0, filename+".png", input_channel, 0);
	MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label+" (prof.)", -2.0, 2.0, filename+".eps", input_channel, 0);
	
	fFileNames.push_back(filename+".png");
      */
    }
    else if  (fParameter[param].name == "N0" || fParameter[param].name == "NL" || fParameter[param].name == "NR" )  {
      
      // parameter distribution
      std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".png", input_channel, 1, 0);
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".eps", input_channel, 1, 0);
      
      fFileNames.push_back(filename+".png");
      
      // error distribution
      filename = fOutputFolder+"/plot_"+input_channel+"_parerror_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 1, 1);
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 1, 1);
      
      fFileNames.push_back(filename+".png");
      
      // profile error distribution
      filename = fOutputFolder+"/plot_"+input_channel+"_proferror_"+name+"_"+oss2.str();
      
      //      if(fParameter[param].Up_prof_unc > 35000) fParameter[param].Up_prof_unc = 35000;
      
      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 1, 1);
      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 1, 1);
      
      fFileNames.push_back(filename+".png");
      
      // pull distribution minuit error
      filename = fOutputFolder+"/plot_"+input_channel+"_pull_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // pull distribution profiling error   
      filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      
      filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_gauss_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 10);
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 10);
      
      fFileNames.push_back(filename+".png");
      
    }
    else {
      
      // parameter distribution
      std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // error distribution
      filename = fOutputFolder+"/plot_"+input_channel+"_parerror_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // profile error distribution
      filename = fOutputFolder+"/plot_"+input_channel+"_proferror_"+name+"_"+oss2.str();
      
      //    if(fParameter[param].Up_prof_unc > 35000) fParameter[param].Up_prof_unc = 35000;
      
      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // pull distribution minuit error
      filename = fOutputFolder+"/plot_"+input_channel+"_pull_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
      // pull distribution profiling error   
      filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_"+name+"_"+oss2.str();
      
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".png", input_channel, 0, 1);
      MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename+".eps", input_channel, 0, 1);
      
      fFileNames.push_back(filename+".png");
      
    }
    
  }
  
  //Handle covariance plots
  TH2D *h2_F0FL  = new TH2D("h2_F0FL", "", 100, F0_true-0.19,F0_true+0.2, 100, FL_true-0.2,FL_true+0.2);
  TH2D *h2_F0FR  = new TH2D("h2_F0FR", "", 100, F0_true-0.19,F0_true+0.2, 100, FR_true-0.2,FR_true+0.2);
  TH2D *h2_FLFR  = new TH2D("h2_FLFR", "", 100, FL_true-0.19,FL_true+0.2, 100, FR_true-0.2,FR_true+0.2);
  
  TH2D *h2_N0NL  = new TH2D("h2_N0NL", "", 100, xsec*lumi*F0_true-110000,xsec*lumi*F0_true+110000, 100, xsec*lumi*FL_true-110000,xsec*lumi*FL_true+110000);
  TH2D *h2_N0NR  = new TH2D("h2_N0NR", "", 100, xsec*lumi*F0_true-110000,xsec*lumi*F0_true+110000, 100, xsec*lumi*FR_true-110000,xsec*lumi*FR_true+110000);
  TH2D *h2_NLNR  = new TH2D("h2_NLNR", "", 100, xsec*lumi*FL_true-80000,xsec*lumi*FL_true+80000, 100, xsec*lumi*FR_true-80000,xsec*lumi*FR_true+80000);
  
  for(int k = 0; k < fTotalEntries; ++k){
    
    h2_F0FL->Fill(vec_F0[k], vec_FL[k]);	
    h2_F0FR->Fill(vec_F0[k], vec_FR[k]);	
    h2_FLFR->Fill(vec_FL[k], vec_FR[k]);
    
    h2_N0NL->Fill(vec_N0[k], vec_NL[k]);
    h2_N0NR->Fill(vec_N0[k], vec_NR[k]);
    h2_NLNR->Fill(vec_NL[k], vec_NR[k]);
  }
  
  std::string label_x, label_y, title;
  
  std::stringstream oss2;
  if (F0_true < 0.1) oss2 << 0.0;
  else oss2 << F0_true;
  
  std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_F0FL_"+oss2.str();
  Plot2Dcorrelation_F(h2_F0FL, F0_true, FL_true, label_x = "F_{0}", label_y = "F_{L}", filename+".png", input_channel);
  Plot2Dcorrelation_F(h2_F0FL, F0_true, FL_true, label_x = "F_{0}", label_y = "F_{L}", filename+".eps", input_channel);
  
  //  fFileNames.push_back(filename+".png");

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_F0FR_"+oss2.str();
  Plot2Dcorrelation_F(h2_F0FR, F0_true, FR_true, label_x = "F_{0}", label_y = "F_{R}", filename+".png", input_channel);
  Plot2Dcorrelation_F(h2_F0FR, F0_true, FR_true, label_x = "F_{0}", label_y = "F_{R}", filename+".eps", input_channel);

  // fFileNames.push_back(filename+".png");

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_FLFR_"+oss2.str();
  Plot2Dcorrelation_F(h2_FLFR, FL_true, FR_true, label_x = "F_{L}", label_y = "F_{R}", filename+".png", input_channel);
  Plot2Dcorrelation_F(h2_FLFR, FL_true, FR_true, label_x = "F_{L}", label_y = "F_{R}", filename+".eps", input_channel);

  // fFileNames.push_back(filename+".png");

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_N0NL_"+oss2.str();
  Plot2Dcorrelation_N(h2_N0NL, F0_true, label_x = "N_{0}", label_y = "N_{L}", filename+".png", input_channel);
  Plot2Dcorrelation_N(h2_N0NL, F0_true, label_x = "N_{0}", label_y = "N_{L}", filename+".eps", input_channel);

  // fFileNames.push_back(filename+".png");

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_N0NR_"+oss2.str();
  Plot2Dcorrelation_N(h2_N0NR, F0_true, label_x = "N_{0}", label_y = "N_{R}", filename+".png", input_channel);
  Plot2Dcorrelation_N(h2_N0NR, F0_true, label_x = "N_{0}", label_y = "N_{R}", filename+".eps", input_channel);

  // fFileNames.push_back(filename+".png");

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_NLNR_"+oss2.str();

  Plot2Dcorrelation_N(h2_NLNR, F0_true, label_x = "N_{L}", label_y = "N_{R}", filename+".png", input_channel);
  Plot2Dcorrelation_N(h2_NLNR, F0_true, label_x = "N_{L}", label_y = "N_{R}", filename+".eps", input_channel);
  //
  // fFileNames.push_back(filename+".png");

  //Calculate correlation
  std::vector<double> vec_corr;
  vec_corr.push_back(h2_N0NL->GetCorrelationFactor() );
  vec_corr.push_back(h2_N0NR->GetCorrelationFactor() );
  vec_corr.push_back(h2_NLNR->GetCorrelationFactor() );

  //Calculate new errormatrix based on profiling likelihood
  //Calculation for each PE

  //Histograms for errors on different helicity states
  TH1D *h_prof_err_F0 = new TH1D("h_prof_err_F0","",100,0.03,+0.075);
  TH1D *h_prof_err_FL = new TH1D("h_prof_err_FL","",100,0.02,+0.055);
  TH1D *h_prof_err_FR = new TH1D("h_prof_err_FR","",100,0.01,+0.045);
  //Helicity fractions
  TH1D *h_prof_pull_F0 = new TH1D("h_prof_pull_F0","",100,-7,+7);
  TH1D *h_prof_pull_FL = new TH1D("h_prof_pull_FL","",100,-7,+7);
  TH1D *h_prof_pull_FR = new TH1D("h_prof_pull_FR","",100,-7,+7);

  const int emsize=3;
  double errormatrix[emsize][emsize];
  std::vector<double> vec_error;

  for(int k = 0; k < fTotalEntries; ++k){

	vec_error.clear();
	vec_error.push_back(vec_N0_prof_err[k]);
	vec_error.push_back(vec_NL_prof_err[k]);
	vec_error.push_back(vec_NR_prof_err[k]);

	for(int i = 0; i < 3; ++i){
		errormatrix[i][i] = vec_error[i]*vec_error[i];
		for (int j = i+1; j<3; ++j) {
			errormatrix[i][j] = vec_error[i]*vec_error[j]*vec_corr[i+j-1];
			errormatrix[j][i] = vec_error[i]*vec_error[j]*vec_corr[i+j-1];
		}
	}
// 	h_cov_N0NL->Fill(errormatrix[0][1]);
// 	h_cov_N0NR->Fill(errormatrix[0][2]);
// 	h_cov_NLNR->Fill(errormatrix[1][2]);

	double sum = vec_N0[k]+vec_NL[k]+vec_NR[k];
	//... determine errors on Fi then :
	double sigma_F0 = fi_error (&errormatrix[0][0], emsize, vec_N0[k], sum, 0);
	double sigma_FL = fi_error (&errormatrix[0][0], emsize, vec_NL[k], sum, 1);
	double sigma_FR = fi_error (&errormatrix[0][0], emsize, vec_NR[k], sum, 2);

	h_prof_err_F0->Fill(sigma_F0);
	h_prof_err_FL->Fill(sigma_FL);
	h_prof_err_FR->Fill(sigma_FR);

	//Corresponding pull values
	h_prof_pull_F0->Fill( (vec_F0[k]-F0_true)/sigma_F0);
	h_prof_pull_FL->Fill( (vec_FL[k]-FL_true)/sigma_FL);
	h_prof_pull_FR->Fill( (vec_FR[k]-FR_true)/sigma_FR);
  }

  std::stringstream oss;
  oss << F0_true;

  //Plot profile error plots for Fi
  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_F0_"+oss2.str();
  MakePlots(h_prof_err_F0, F0_true, oss.str(), "#sigma(F0)", 0.0, 0.08, filename+".png", input_channel, 7, 1);
  MakePlots(h_prof_err_F0, F0_true, oss.str(), "#sigma(F0)", 0.0, 0.08, filename+".eps", input_channel, 7, 1);

  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_FL_"+oss2.str();
  MakePlots(h_prof_err_FL, F0_true, oss.str(), "#sigma(FL)", 0.0, 0.08, filename+".png", input_channel, 7, 1);
  MakePlots(h_prof_err_FL, F0_true, oss.str(), "#sigma(FL)", 0.0, 0.08, filename+".eps", input_channel, 7, 1);

  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_FR_"+oss2.str();
  MakePlots(h_prof_err_FR, F0_true, oss.str(), "#sigma(FR)", 0.0, 0.08, filename+".png", input_channel, 7, 1);
  MakePlots(h_prof_err_FR, F0_true, oss.str(), "#sigma(FR)", 0.0, 0.08, filename+".eps", input_channel, 7, 1);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_F0_"+oss2.str();
  MakePlots(h_prof_pull_F0, F0_true, oss.str(), "Pull distribution for F0", -2.0, 2.0, filename+".png", input_channel, 0, 1);
  MakePlots(h_prof_pull_F0, F0_true, oss.str(), "Pull distribution for F0", -2.0, 2.0, filename+".eps", input_channel, 0, 1);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_FL_"+oss2.str();
  MakePlots(h_prof_pull_FL, F0_true, oss.str(), "Pull distribution for FL", -2.0, 2.0, filename+".png", input_channel, 0, 1);
  MakePlots(h_prof_pull_FL, F0_true, oss.str(), "Pull distribution for FL", -2.0, 2.0, filename+".eps", input_channel, 0, 1);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_FR_"+oss2.str();
  MakePlots(h_prof_pull_FR, F0_true, oss.str(), "Pull distribution for FR", -2.0, 2.0, filename+".png", input_channel, 0, 1);
  MakePlots(h_prof_pull_FR, F0_true, oss.str(), "Pull distribution for FR", -2.0, 2.0, filename+".eps", input_channel, 0, 1);


  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_gauss_F0_"+oss2.str();
  MakePlots(h_prof_pull_F0, F0_true, oss.str(), "Pull distribution for F0", -2.0, 2.0, filename+".png", input_channel, 0, 10);
  MakePlots(h_prof_pull_F0, F0_true, oss.str(), "Pull distribution for F0", -2.0, 2.0, filename+".eps", input_channel, 0, 10);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_gauss_FL_"+oss2.str();
  MakePlots(h_prof_pull_FL, F0_true, oss.str(), "Pull distribution for FL", -2.0, 2.0, filename+".png", input_channel, 0, 10);
  MakePlots(h_prof_pull_FL, F0_true, oss.str(), "Pull distribution for FL", -2.0, 2.0, filename+".eps", input_channel, 0, 10);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_gauss_FR_"+oss2.str();
  MakePlots(h_prof_pull_FR, F0_true, oss.str(), "Pull distribution for FR", -2.0, 2.0, filename+".png", input_channel, 0, 10);
  MakePlots(h_prof_pull_FR, F0_true, oss.str(), "Pull distribution for FR", -2.0, 2.0, filename+".eps", input_channel, 0, 10);


    std::vector<double> F_err_vec;    

 /*   if (fValidationMode == "Datafit"){	

    	std::string filename_data = "Datafit_Output/Datafit_Output_" + input_channel + ".root";
    	TFile *file_data = GetInputFile(filename_data);

  	F_err_vec = GetErrorVector(file_data);
    } */


	//error distributions including datafit error
	if (fValidationMode == "Datafit"){	

		std::string filename_data = "Datafit_Output/Datafit_Output_" + input_channel + ".root";
    		TFile *file_data = GetInputFile(filename_data);

	  	//F_err_vec = GetErrorVector(file_data);
		F_err_vec.clear();
		F_err_vec.push_back(0.08513);
		F_err_vec.push_back(0.05038);
		F_err_vec.push_back(0.03935);
		//F_err_vec.push_back(0.07848);
		//F_err_vec.push_back(0.04517);
		//F_err_vec.push_back(0.03844);

		filename = fOutputFolder+"/plot_"+input_channel+"_proferror_data_F0_"+oss2.str();

		MakeErrorDataPlots(h_prof_err_F0, F0_true, F_err_vec[0], oss.str(), "#sigma(F0)", 0.0, 0.08, filename+".png", input_channel, 0);
      		MakeErrorDataPlots(h_prof_err_F0, F0_true, F_err_vec[0], oss.str(), "#sigma(F0)", 0.0, 0.08, filename+".eps", input_channel, 0);

		filename = fOutputFolder+"/plot_"+input_channel+"_proferror_data_FL_"+oss2.str();

		MakeErrorDataPlots(h_prof_err_FL, F0_true, F_err_vec[1], oss.str(), "#sigma(FL)",  0.0, 0.08, filename+".png", input_channel, 0);
      		MakeErrorDataPlots(h_prof_err_FL, F0_true, F_err_vec[1], oss.str(), "#sigma(FL)",  0.0, 0.08, filename+".eps", input_channel, 0);

		filename = fOutputFolder+"/plot_"+input_channel+"_proferror_data_FR_"+oss2.str();

		MakeErrorDataPlots(h_prof_err_FR, F0_true, F_err_vec[2], oss.str(), "#sigma(FR)",  0.0, 0.08, filename+".png", input_channel, 0);
      		MakeErrorDataPlots(h_prof_err_FR, F0_true, F_err_vec[2], oss.str(), "#sigma(FR)",  0.0, 0.08, filename+".eps", input_channel, 0);
 
	}





  //For filling calibration tree
  fParameter[0].pull_mean_prof     = h_prof_pull_F0 -> GetMean();
  fParameter[0].pull_mean_prof_err = h_prof_pull_F0 -> GetMeanError();
  fParameter[0].pull_rms_prof      = h_prof_pull_F0 -> GetRMS();
  fParameter[0].pull_rms_prof_err  = h_prof_pull_F0 -> GetRMSError();
  fParameter[1].pull_mean_prof     = h_prof_pull_FL -> GetMean();
  fParameter[1].pull_mean_prof_err = h_prof_pull_FL -> GetMeanError();
  fParameter[1].pull_rms_prof      = h_prof_pull_FL -> GetRMS();
  fParameter[1].pull_rms_prof_err  = h_prof_pull_FL -> GetRMSError();
  fParameter[2].pull_mean_prof     = h_prof_pull_FR -> GetMean();
  fParameter[2].pull_mean_prof_err = h_prof_pull_FR -> GetMeanError();
  fParameter[2].pull_rms_prof      = h_prof_pull_FR -> GetRMS();
  fParameter[2].pull_rms_prof_err  = h_prof_pull_FR -> GetRMSError();

  std::vector<double> gauss_fits_pull_prof;
  gauss_fits_pull_prof.clear();
  gauss_fits_pull_prof = MakeFits(h_prof_pull_F0, "Pull distribution for F0",-2.0, 2.0);

  fParameter[0].gauss_pull_mean_prof     = gauss_fits_pull_prof[0];
  fParameter[0].gauss_pull_mean_prof_err = gauss_fits_pull_prof[1];
  fParameter[0].gauss_pull_rms_prof      = gauss_fits_pull_prof[2];
  fParameter[0].gauss_pull_rms_prof_err  = gauss_fits_pull_prof[3];

  gauss_fits_pull_prof.clear();
  gauss_fits_pull_prof = MakeFits(h_prof_pull_FL, "Pull distribution for FL",-2.0, 2.0);

  fParameter[1].gauss_pull_mean_prof     = gauss_fits_pull_prof[0];
  fParameter[1].gauss_pull_mean_prof_err = gauss_fits_pull_prof[1];
  fParameter[1].gauss_pull_rms_prof      = gauss_fits_pull_prof[2];
  fParameter[1].gauss_pull_rms_prof_err  = gauss_fits_pull_prof[3];

  gauss_fits_pull_prof.clear();
  gauss_fits_pull_prof = MakeFits(h_prof_pull_FR, "Pull distribution for FR",-2.0, 2.0);

  fParameter[2].gauss_pull_mean_prof     = gauss_fits_pull_prof[0];
  fParameter[2].gauss_pull_mean_prof_err = gauss_fits_pull_prof[1];
  fParameter[2].gauss_pull_rms_prof      = gauss_fits_pull_prof[2];
  fParameter[2].gauss_pull_rms_prof_err  = gauss_fits_pull_prof[3];
  
  gauss_fits_pull_prof.clear();


  // plot all 2D histograms
  
  for(int param1 = 0; param1 < fParameter.size(); ++param1){

    for(int param2 = 0; param2 < fParameter.size(); ++param2){

      if(param1 != param2){
	if(param1 == 10){

	  std::string namex = fParameter[param1].label;
	  std::string namey = fParameter[param2].label;

	  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_"+namex+"_"+namey+"_"+oss2.str();

	  TH2D *histo = (TH2D*) &fParameter[param1].hist2D[param2];

	  Plot2Dcorrelation_N(histo, F0_true, namex, namey, filename+".png", input_channel);
	  Plot2Dcorrelation_N(histo, F0_true, namex, namey, filename+".eps", input_channel);


	  fFileNames.push_back(filename+".png");

	}
      }

    }

  }

}


void histograms::MakeOutputTrees(int num_nui)
{

  TFile *outputFile = new TFile((fOutputFolder+"/CalibrationTree.root").c_str(), "RECREATE");

  TTree *F0_tree = InitializeOutputTree("F0_tree", 0);
  
  F0_tree -> Fill();

  F0_tree -> Write();

  TTree *FL_tree = InitializeOutputTree("FL_tree", 1);

  FL_tree -> Fill();

  FL_tree -> Write();

  TTree *FR_tree = InitializeOutputTree("FR_tree", 2);

  FR_tree -> Fill();

  FR_tree -> Write();

  TTree *N0_tree = InitializeOutputTree("N0_tree", 3);

  N0_tree -> Fill();

  N0_tree -> Write();

  TTree *NL_tree = InitializeOutputTree("NL_tree", 4);

  NL_tree -> Fill();

  NL_tree -> Write();

  TTree *NR_tree = InitializeOutputTree("NR_tree", 5);

  NR_tree -> Fill();

  NR_tree -> Write();

  TTree *Wjets_tree = InitializeOutputTree("Wjets_tree", 6);

  Wjets_tree -> Fill();

  Wjets_tree -> Write();

  TTree *QCD_tree = InitializeOutputTree("QCD_tree", 7);

  QCD_tree -> Fill();

  QCD_tree -> Write();

  TTree *RemBkg_tree = InitializeOutputTree("RemBkg_tree", 8);

  RemBkg_tree -> Fill();

  RemBkg_tree -> Write();

  if(num_nui > 0){

    TTree *nui1_tree = InitializeOutputTree("nui1_tree", 10);
    
    nui1_tree -> Fill();
    
    nui1_tree -> Write();

  }

  if(num_nui > 1){

    TTree *nui2_tree = InitializeOutputTree("nui2_tree", 11);
    
    nui2_tree -> Fill();
    
    nui2_tree -> Write();

  }

  if(num_nui > 2){

    TTree *nui3_tree = InitializeOutputTree("nui3_tree", 12);
    
    nui3_tree -> Fill();
    
    nui3_tree -> Write();
    
  }

  if(num_nui > 3){

    TTree *nui4_tree = InitializeOutputTree("nui4_tree", 13);
    
    nui4_tree -> Fill();
    
    nui4_tree -> Write();

  }

  if(num_nui > 4){

    TTree *nui5_tree = InitializeOutputTree("nui5_tree", 14);

    nui5_tree -> Fill();
    
    nui5_tree -> Write();
    
  }

  if(num_nui > 5){

    TTree *nui6_tree = InitializeOutputTree("nui6_tree", 15);
    
    nui6_tree -> Fill();
    
    nui6_tree -> Write();
    
  }

  if(num_nui > 6){
    
    TTree *nui7_tree = InitializeOutputTree("nui7_tree", 16);
    
    nui7_tree -> Fill();
    
    nui7_tree -> Write();

  }

  if(num_nui > 7){

    TTree *nui8_tree = InitializeOutputTree("nui8_tree", 17);
    
    nui8_tree -> Fill();
    
    nui8_tree -> Write();
    
  }

  if(num_nui > 8){

    TTree *nui9_tree = InitializeOutputTree("nui9_tree", 18);
    
    nui9_tree -> Fill();
    
    nui9_tree -> Write();
    
  }

  if(num_nui > 9){

    TTree *nui10_tree = InitializeOutputTree("nui10_tree", 19);
    
    nui10_tree -> Fill();
    
    nui10_tree -> Write();
    
  }

  if(num_nui > 10){

    TTree *nui11_tree = InitializeOutputTree("nui11_tree", 20);
    
    nui11_tree -> Fill();
    
    nui11_tree -> Write();
    
  }

  if(num_nui > 11){

    TTree *nui12_tree = InitializeOutputTree("nui12_tree", 21);
    
    nui12_tree -> Fill();
    
    nui12_tree -> Write();
    
  }

  if(num_nui > 12){

    TTree *nui13_tree = InitializeOutputTree("nui13_tree", 22);
    
    nui13_tree -> Fill();
    
    nui13_tree -> Write();
    
  }


  outputFile -> Close();

}


TTree *histograms::InitializeOutputTree(std::string ParameterTree, int param)
{

  TTree *OutputTree = new TTree(ParameterTree.c_str(), ParameterTree.c_str());

  OutputTree -> Branch("Nominal",              &fParameter[param].nom);    
  OutputTree -> Branch("Mean",                 &fParameter[param].mean);
  OutputTree -> Branch("MeanErr",              &fParameter[param].mean_err);
  OutputTree -> Branch("GaussMean",            &fParameter[param].gauss_mean);
  OutputTree -> Branch("GaussMeanErr",         &fParameter[param].gauss_mean_err);
  OutputTree -> Branch("PullMean",             &fParameter[param].pull_mean);
  OutputTree -> Branch("PullMeanErr",          &fParameter[param].pull_mean_err);
  OutputTree -> Branch("PullRMS",              &fParameter[param].pull_rms);
  OutputTree -> Branch("PullRMSErr",           &fParameter[param].pull_rms_err);
  OutputTree -> Branch("GaussPullMean",        &fParameter[param].gauss_pull_mean);
  OutputTree -> Branch("GaussPullMeanErr",     &fParameter[param].gauss_pull_mean_err);
  OutputTree -> Branch("GaussPullRMS",         &fParameter[param].gauss_pull_rms);
  OutputTree -> Branch("GaussPullRMSErr",      &fParameter[param].gauss_pull_rms_err);
  OutputTree -> Branch("MeanProf",             &fParameter[param].mean_prof);
  OutputTree -> Branch("MeanProfErr",          &fParameter[param].mean_prof_err);
  OutputTree -> Branch("PullMeanProf",         &fParameter[param].pull_mean_prof);
  OutputTree -> Branch("PullMeanProfErr",      &fParameter[param].pull_mean_prof_err);
  OutputTree -> Branch("PullRMSProf",          &fParameter[param].pull_rms_prof);
  OutputTree -> Branch("PullRMSProfErr",       &fParameter[param].pull_rms_prof_err);
  OutputTree -> Branch("GaussPullMeanProf",    &fParameter[param].gauss_pull_mean_prof);
  OutputTree -> Branch("GaussPullMeanProfErr", &fParameter[param].gauss_pull_mean_prof_err);
  OutputTree -> Branch("GaussPullRMSProf",     &fParameter[param].gauss_pull_rms_prof);
  OutputTree -> Branch("GaussPullRMSProfErr",  &fParameter[param].gauss_pull_rms_prof_err);

  return OutputTree;

}


///_______________________________MAKE FITS___________________________
std::vector<double> histograms::MakeFits(TH1D *hist, std::string label, double low, double up)
{
  std::vector<double> gauss_fit;

  TF1 *fit = new TF1(label.c_str(), "gaus", low, up);
  hist -> Fit(label.c_str(), "RQ"); // "R" means range

  gauss_fit.push_back( fit->GetParameter(1) );
  gauss_fit.push_back( fit->GetParError(1)  );
  gauss_fit.push_back( fit->GetParameter(2) );
  gauss_fit.push_back( fit->GetParError(2) );

  delete fit;

  return gauss_fit;

}


///________________________________MAKE PLOTS___________________________
void histograms::MakePlots(TH1D *hist, double F0_true, std::string nominal, std::string label, double low, double up, std::string filename, std::string input_channel, int mode, int mode_gauss)
{
  TF1 *fit = new TF1(label.c_str(), "gaus", -4.0, +4.0);

  if(mode_gauss == 10) {

    int max_bin = hist -> GetMaximumBin();
    double max  = hist -> GetBinContent(max_bin);
 
    //double nom = hist->Integral();
    double nom = hist->GetSumOfWeights();
	//cout << "hist " << nom << endl;
    //hist->Scale(0.4/max);
    hist->Scale(1/(nom*hist->GetBinWidth(1)) );
    const double PI = std::atan(1.0)*4;
	//cout << "pi " << PI << endl;
    fit->SetParameter(0, 1/sqrt(2*PI));
    fit->SetParameter(1, 0.0);
    fit->SetParameter(2, 1.0);
    fit->FixParameter(1, 0.0);
    fit->FixParameter(2, 1.0);

//    hist -> Fit(label.c_str(), "RBQ"); // "R" means range
  }
    int max_bin = hist -> GetMaximumBin();
    double max  = hist -> GetBinContent(max_bin);
 

  // hardcoded for a first test!!!                                                                                                                               
  //std::string input_channel = "mu";

  //cout << "low and up position " << up << endl;

    //  PlotDistribution(hist, F0_true, label.c_str(), input_channel + "_pseudo_"+label+"=%.2f", max*1.2, input_channel, nominal, fit, filename, low, up, mode, mode_gauss);

  delete fit;

}



///_________________________MAKE ERROR DATA PLOTS_________________________
void histograms::MakeErrorDataPlots(TH1D *hist, double F0_true, double F_err, std::string nominal, std::string label, double low, double up, std::string filename, std::string input_channel, int mode)
{
  TF1 *fit = new TF1(label.c_str(), "gaus", low, up);

  //hist -> Fit(label.c_str(), "RQ"); // "R" means range

  int max_bin = hist -> GetMaximumBin();
  double max  = hist -> GetBinContent(max_bin);

  // hardcoded for a first test!!!                                                                                                                               
  //std::string input_channel = "mu";

  //cout << "low and up position " << up << endl;

  //  PlotErrorDataComp(hist, F0_true, F_err, label.c_str(), input_channel + "_pseudo_"+label+"=%.2f", max*1.2, input_channel, nominal, fit, filename, mode);
  //void PlotErrorDataComp (TH1D* hist, double F_true, double F_data, std::string xlabel, std::string xtitle, double max, std::string input_channel, std::string sFi_number, int mode)

  delete fit;

}


///__________________________FIND RANGES_______________________________
void histograms::FindRanges(int num_of_nui_parameters)
{

  fTotalEntries = fInputTree -> fChain -> GetEntries();

   //double test_out = fInputTree->test;
   //cout << "test output : " << test_out << endl;
   //cout << "Start loop for entries" << endl;
  
  for(int k = 0; k < fTotalEntries; ++k){

    //cout << "Get entry " << k << endl;
    fInputTree -> fChain -> GetEntry(k);

    for(int param = 0; param < fParameter.size(); ++param){

      double central, unc, prof_unc, nominal;

      if(fParameter[param].name == "F0"){

	central = fInputTree -> F0; unc = fInputTree -> F0_err; nominal = fInputTree -> F0_nom;
      }
      else if(fParameter[param].name == "FL"){

        central = fInputTree -> FL; unc = fInputTree -> FL_err; nominal = fInputTree -> FL_nom;
      }
      else if(fParameter[param].name == "FR"){

        central = fInputTree -> FR; unc = fInputTree -> FR_err; nominal = fInputTree -> FR_nom;
      }
      else if(fParameter[param].name == "N0"){

	central = fInputTree -> N0; unc = fInputTree -> N0_err; prof_unc = fInputTree -> N0_prof_err; nominal = fInputTree -> N0_nom;
      }
      else if(fParameter[param].name == "NL"){

        central = fInputTree -> NL; unc = fInputTree -> NL_err; prof_unc = fInputTree -> NL_prof_err; nominal = fInputTree -> NL_nom;
      }
      else if(fParameter[param].name == "NR"){

        central = fInputTree -> NR; unc = fInputTree -> NR_err; prof_unc = fInputTree -> NR_prof_err; nominal = fInputTree -> NR_nom;
      }
      else if(fParameter[param].name == "Nges"){

        central = fInputTree -> Nges; nominal = fInputTree->NR_nom + fInputTree->NL_nom + fInputTree->N0_nom;
      }
      else if(fParameter[param].name == "Wjets"){

        central = fInputTree -> Wjets; unc = fInputTree -> Wjets_err; prof_unc = fInputTree -> Wjets_prof_err; nominal = fInputTree -> Wjets_nom;
      }
      else if(fParameter[param].name == "QCD"){

        central = fInputTree -> QCD; unc = fInputTree -> QCD_err; prof_unc = fInputTree -> QCD_prof_err; nominal = fInputTree -> QCD_nom;
      }
      else if(fParameter[param].name == "RemBkg"){

        central = fInputTree -> RemBkg; unc = fInputTree -> RemBkg_err; prof_unc = fInputTree -> RemBkg_prof_err; nominal = fInputTree -> RemBkg_nom;
      }
      else {
	
 	//cout << "else condition " << endl;

	//for (int i=0; i != num_of_nui_parameters; i++)
	//{
	/*std::stringstream nui_label;
	nui_label << i;
	std::string nui_name  = "nui_"+nui_label.str() ; 

	if(fParameter[param].name == nui_name){
		
		//cout << "here 0 : " << i << endl;
		//central = nui_value[i];
		//cout << central << endl; 
		//unc = fInputTree -> nui_err[i]; 
		//prof_unc = fInputTree -> nui_prof_err[i]; 
		//nominal = fInputTree -> nui_nom[i];
		//cout << "here 1 " << endl;
        	central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0;

      	}*/

	  if (fParameter[param].name == "nui_0") {
	    central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0; }
	  if (fParameter[param].name == "nui_1") {
	    central = fInputTree -> nui_1; unc = fInputTree -> nui_err_1; prof_unc = fInputTree -> nui_prof_err_1; nominal = fInputTree -> nui_nom_1; }
	  if (fParameter[param].name == "nui_2") {
	    central = fInputTree -> nui_2; unc = fInputTree -> nui_err_2; prof_unc = fInputTree -> nui_prof_err_2; nominal = fInputTree -> nui_nom_2; }
	  if (fParameter[param].name == "nui_3") {
	    central = fInputTree -> nui_3; unc = fInputTree -> nui_err_3; prof_unc = fInputTree -> nui_prof_err_3; nominal = fInputTree -> nui_nom_3; }
	  if (fParameter[param].name == "nui_4") {
	    central = fInputTree -> nui_4; unc = fInputTree -> nui_err_4; prof_unc = fInputTree -> nui_prof_err_4; nominal = fInputTree -> nui_nom_4; }
	  if (fParameter[param].name == "nui_5") {
	    central = fInputTree -> nui_5; unc = fInputTree -> nui_err_5; prof_unc = fInputTree -> nui_prof_err_5; nominal = fInputTree -> nui_nom_5; }
	  if (fParameter[param].name == "nui_6") {
	    central = fInputTree -> nui_6; unc = fInputTree -> nui_err_6; prof_unc = fInputTree -> nui_prof_err_6; nominal = fInputTree -> nui_nom_6; }
	  if (fParameter[param].name == "nui_7") {
	    central = fInputTree -> nui_7; unc = fInputTree -> nui_err_7; prof_unc = fInputTree -> nui_prof_err_7; nominal = fInputTree -> nui_nom_7; }
	  if (fParameter[param].name == "nui_8") {
	    central = fInputTree -> nui_8; unc = fInputTree -> nui_err_8; prof_unc = fInputTree -> nui_prof_err_8; nominal = fInputTree -> nui_nom_8; }
	  if (fParameter[param].name == "nui_9") {
	    central = fInputTree -> nui_9; unc = fInputTree -> nui_err_9; prof_unc = fInputTree -> nui_prof_err_9; nominal = fInputTree -> nui_nom_9; }
	  if (fParameter[param].name == "nui_10") {
	    central = fInputTree -> nui_10; unc = fInputTree -> nui_err_10; prof_unc = fInputTree -> nui_prof_err_10; nominal = fInputTree -> nui_nom_10; }
	  if (fParameter[param].name == "nui_11") {
	    central = fInputTree -> nui_11; unc = fInputTree -> nui_err_11; prof_unc = fInputTree -> nui_prof_err_11; nominal = fInputTree -> nui_nom_11; }
	  if (fParameter[param].name == "nui_12") {
	    central = fInputTree -> nui_12; unc = fInputTree -> nui_err_12; prof_unc = fInputTree -> nui_prof_err_12; nominal = fInputTree -> nui_nom_12; }
      	//}
	
      }
      
      if(central  < fParameter[param].Low)          fParameter[param].Low          = central;
      if(central  > fParameter[param].Up)           fParameter[param].Up           = central;
      if(unc      < fParameter[param].Low_unc)      fParameter[param].Low_unc      = unc;
      if(unc      > fParameter[param].Up_unc)       fParameter[param].Up_unc       = unc;
      if(prof_unc < fParameter[param].Low_prof_unc) fParameter[param].Low_prof_unc = prof_unc;
      if(prof_unc > fParameter[param].Up_prof_unc)  fParameter[param].Up_prof_unc  = prof_unc;

      if(k == 0) fParameter[param].nom = nominal;

    }
  }
}


///__________________________MAKE CHAIN OF FILES______________________
void histograms::MakeChainOfFiles(int num_of_nui_parameters)
{

  int counter_lines = 0;

  // initialize Tree...
  TChain *fHelpChain = new TChain("EnsembleTree");

  // find out which files are in folder  
  DIR *dpdf;
  struct dirent *epdf;

  dpdf = opendir(fInputFolder.c_str());
  if (dpdf != NULL){
    while (epdf = readdir(dpdf)){

      // first two entries are . and ..
      // if(counter_lines > 1){

      std::stringstream oss;
      oss << epdf->d_name;
      
      if(strlen((oss.str()).c_str()) > 3){

	//	std::cout << (fInputFolder+"/"+oss.str()).c_str() << std::endl;
	
	fHelpChain -> Add((fInputFolder+"/"+oss.str()).c_str());

      }

      counter_lines++;

    }
  }

  cout << "Define output tree reader " << endl;

  fInputTree = new OutputTreeReader(fHelpChain, num_of_nui_parameters);

  cout << "Get entries " << endl;

  fInputTree -> fChain -> GetEntries();

}


//Get input file	
TFile * histograms::GetInputFile(std::string filename)
{		
	const char *name;
	name = filename.c_str();

	//Open files
	TFile* file = new TFile(name,"open");
	return file;
}


//Obtain information from tree
std::vector<double> histograms::GetErrorVector(TFile *input_file)
{

	//float F0_err, FL_err, FR_err;
	double F0_err, FL_err, FR_err;
	std::vector<double> err_vec;

	// Get Tree
	TTree *tree = (TTree*)input_file->Get("EnsembleTree");

	// Link to Branch in Tree
	tree->SetBranchAddress("F0_err",&F0_err);
	tree->SetBranchAddress("FL_err",&FL_err);
	tree->SetBranchAddress("FR_err",&FR_err);

	// Get the Event
	tree->GetEntry(0);

	//cout << "F0eff: " << F0eff << endl; 
	//cout << "FLeff: " << FLeff << endl; 
	//cout << "FReff: " << FReff << endl; 

	//Information saved in tree in a vector:
	err_vec.push_back((double) F0_err);
	err_vec.push_back((double) FL_err);
	err_vec.push_back((double) FR_err);

	tree->Delete();

	return err_vec;
}


///______________________________MAKE HTML__________________________________
void histograms::MakeHTML(std::string name)
{

  // create the page
  std::ofstream page;
  std::string   pname = fOutputFolder+"/"+name;

  struct stat buf;
  if(stat(pname.c_str(), &buf) == -1){

    std::cout << "Creating file " << pname.c_str() << " ... " << std::endl;
    
    page.open(pname.c_str());
          
    page << "<html><head><title> Parameter and pull distributions </title></head>" << std::endl;
    page << "<body>" << std::endl;
    page << "<h1> Parameter and pull distributions </h1>" << std::endl;
       
    page << "<table border = 1> <tr>" << "<th> parameter </th>" << "<th> uncertainty  </th>" << "<th> uncertainty (prof.) </th>"  << "<th> pull distribution  </th>" << "<th> pull distribution (prof) </th>"  << "</tr>" << std::endl;

    int nges_val = 45;

    int counter  = 0;

    while(counter < fFileNames.size()){

      if(counter != nges_val){
	
	std::string output1 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()),     strlen((fOutputFolder+"/").c_str()), "");
 
	counter++;
	std::string output2;
	std::string output3;
	std::string output4;
	std::string output5;

	if(counter < fFileNames.size()){
	  output2 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()), strlen((fOutputFolder+"/").c_str()), "");
	  counter++;
	}
	if(counter < fFileNames.size()){
	  output3 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()), strlen((fOutputFolder+"/").c_str()), "");
	  counter++;
	}
	if(counter < fFileNames.size()){
	  output4 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()), strlen((fOutputFolder+"/").c_str()), "");
	  counter++;
	}
	if(counter < fFileNames.size()){
	  output5 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()), strlen((fOutputFolder+"/").c_str()), "");
	  counter++;
	}

	page << "<tr><td>"      << "<p><a href=\"" << output1.c_str() << "\">";
	page << "<img src=\""   << output1.c_str()                    << "\" width=450/>";
	page << "<td>"          << "<p><a href=\"" << output2.c_str() << "\">";
	page << "<img src=\""   << output2.c_str()                    << "\" width=450/>";
	page << "<td>"          << "<p><a href=\"" << output3.c_str() << "\">";
	page << "<img src=\""   << output3.c_str()                    << "\" width=450/>";
	page << "<td>"          << "<p><a href=\"" << output4.c_str() << "\">";
	page << "<img src=\""   << output4.c_str()                    << "\" width=450/>";
	page << "<td>"          << "<p><a href=\"" << output5.c_str() << "\">";
	page << "<img src=\""   << output5.c_str()                    << "\" width=450/>";
	page << "</td>"         << std::endl;


      }
      else{

	std::string output1 = fFileNames[counter].replace(fFileNames[counter].find((fOutputFolder+"/").c_str()),     strlen((fOutputFolder+"/").c_str()), "");
	
	page << "<tr><td>"          << "<p><a href=\"" << output1.c_str() << "\">";
	page << "<img src=\""   << output1.c_str()                        << "\" width=450/>";
	page << "<td>" << "" << "";
	page << "<td>" << "" << "";
	page << "<td>" << "" << "";
	page << "<td>" << "" << "" << std::endl;
	
	counter++;

      }

    }

    //    page << "</tr>"         << std::endl;
    
    page.close();
    
  }
  else{

    std::cout << pname.c_str() << " already exists!" << std::endl;

  }

  return;



}

