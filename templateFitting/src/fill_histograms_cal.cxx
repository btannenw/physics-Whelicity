#include "fill_histograms_cal.h"
#include "OutputTreeReader_cal.h"
#include "helicity.h"
#include "plots.h"


#include <iostream>
#include <string>
#include <sstream>
#include <cmath>


#include <dirent.h>

#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"

using namespace std;


histograms::histograms()
{

}

histograms::~histograms()
{

}

///_______________________FILL HISTOGRAMS__________________________
std::vector<double> histograms::FillHistograms(int counter, std::string input_channel, int num_of_nui_parameters, double xsec, double lumi)
{
  
  std::vector<double> calib_info;

  this -> FindRanges(num_of_nui_parameters);
  
  //Define vectors needed for correlation and errormatrix
  std::vector<double> vec_F0, vec_FL, vec_FR; 
  std::vector<double> vec_N0, vec_NL, vec_NR;
  std::vector<double> vec_N0_prof_err, vec_NL_prof_err, vec_NR_prof_err;



  double F0_true, FL_true, FR_true;

  if (counter == 0) {F0_true = 0.7; FL_true=0.3; FR_true=0.0;}
  else {F0_true = 0.7 - counter*0.1; FL_true=0.3 + counter*0.05 ; FR_true=0.0 + counter*0.05;}

  for(int param = 0; param < fParameter.size(); ++param){

    std::string label = fParameter[param].label;

    std::cout << fParameter[param].Low_prof_unc << "\t" << fParameter[param].Up_prof_unc << std::endl;

    fParameter[param].hist_param     = new TH1D(label.c_str(),                "", 100, fParameter[param].Low*0.8,          fParameter[param].Up*1.2);
    fParameter[param].hist_unc       = new TH1D((label+"_unc").c_str(),       "", 100, fParameter[param].Low_unc*0.8,      fParameter[param].Up_unc*1.2);
    fParameter[param].hist_prof_unc  = new TH1D((label+"+unc_prof").c_str(),  "", 100, fParameter[param].Low_prof_unc*0.8, fParameter[param].Up_prof_unc*1.2);
    fParameter[param].hist_pull      = new TH1D((label+"_pull").c_str(),      "", 100, -8.0, 8.0);
    fParameter[param].hist_pull_prof = new TH1D((label+"_pull_prof").c_str(), "", 100, -8.0, 8.0);

  }

  std::cout << "Fill Histograms" << std::endl;

  for(int k = 0; k < fTotalEntries; ++k){

    fInputTree -> fChain -> GetEntry(k);

    for(int param = 0; param < fParameter.size(); ++param){

      double central, unc, prof_unc, nominal;

      if(fParameter[param].name == "F0"){
	
	central = fInputTree -> F0; unc = fInputTree -> F0_err; nominal = fInputTree -> F0_nom;
	vec_F0.push_back(central);
	//cout << central << endl;
      }
      else if(fParameter[param].name == "FL"){

	central = fInputTree -> FL; unc = fInputTree -> FL_err; nominal = fInputTree -> FL_nom;
	vec_FL.push_back(central);
      }
      else if(fParameter[param].name == "FR"){

	central = fInputTree -> FR; unc = fInputTree -> FR_err; nominal = fInputTree -> FR_nom;
	vec_FR.push_back(central);
      }
      else if(fParameter[param].name == "N0"){

	central = fInputTree -> N0; unc = fInputTree -> N0_err; prof_unc = fInputTree -> N0_prof_err; nominal = fInputTree -> N0_nom;
	vec_N0.push_back(central);
	vec_N0_prof_err.push_back(prof_unc);
      }
      else if(fParameter[param].name == "NL"){

	central = fInputTree -> NL; unc = fInputTree -> NL_err; prof_unc = fInputTree -> NL_prof_err; nominal = fInputTree -> NL_nom;
  	vec_NL.push_back(central);
	vec_NL_prof_err.push_back(prof_unc);
      }
      else if(fParameter[param].name == "NR"){

	central = fInputTree -> NR; unc = fInputTree -> NR_err; prof_unc = fInputTree -> NR_prof_err; nominal = fInputTree -> NR_nom;
	vec_NR.push_back(central);
	vec_NR_prof_err.push_back(prof_unc);
      }
      else if(fParameter[param].name == "Nges"){

	central = fInputTree -> Nges; nominal = fInputTree->N0_nom + fInputTree->NL_nom + fInputTree->NR_nom;
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
	
	for (int i=0; i != num_of_nui_parameters; i++)
	{
	/*std::stringstream nui_label;
	nui_label << i;
	std::string nui_name  = "nui_"+nui_label.str() ; 

	if(fParameter[param].name == nui_name){

		//central = fInputTree -> nui_value[i]; unc = fInputTree -> nui_err[i]; prof_unc = fInputTree -> nui_prof_err[i]; nominal = fInputTree -> nui_nom[i];
        	central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0;

      	}*/

	if (i == 0) {
	central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0; }
	if (i == 1) {
	central = fInputTree -> nui_1; unc = fInputTree -> nui_err_1; prof_unc = fInputTree -> nui_prof_err_1; nominal = fInputTree -> nui_nom_1; }
	if (i == 2) {
	central = fInputTree -> nui_2; unc = fInputTree -> nui_err_2; prof_unc = fInputTree -> nui_prof_err_2; nominal = fInputTree -> nui_nom_2; }
	if (i == 3) {
	central = fInputTree -> nui_3; unc = fInputTree -> nui_err_3; prof_unc = fInputTree -> nui_prof_err_3; nominal = fInputTree -> nui_nom_3; }
	if (i == 4) {
	central = fInputTree -> nui_4; unc = fInputTree -> nui_err_4; prof_unc = fInputTree -> nui_prof_err_4; nominal = fInputTree -> nui_nom_4; }
	if (i == 5) {
	central = fInputTree -> nui_5; unc = fInputTree -> nui_err_5; prof_unc = fInputTree -> nui_prof_err_5; nominal = fInputTree -> nui_nom_5; }		

      	}

      }


      if (fParameter[param].name == "Nges") {
	fParameter[param].hist_param     -> Fill(central);
      }
      else if(fParameter[param].name == "F0" || fParameter[param].name == "FL" || fParameter[param].name == "FR" ) {
	//cout << central << endl;
	fParameter[param].hist_param     -> Fill(central);      
        fParameter[param].hist_unc       -> Fill(unc);
        fParameter[param].hist_pull      -> Fill((central - fParameter[param].nom)/unc);
      }
      else {
	fParameter[param].hist_param     -> Fill(central);      
        fParameter[param].hist_unc       -> Fill(unc);
        fParameter[param].hist_prof_unc  -> Fill(prof_unc); 
        fParameter[param].hist_pull      -> Fill((central - fParameter[param].nom)/unc);
        fParameter[param].hist_pull_prof -> Fill((central - fParameter[param].nom)/prof_unc);
	}
    }
  }

  std::cout << "Print out Plots" << std::endl;

  for(int param = 0; param < fParameter.size(); ++param){

    std::stringstream oss;
    oss << fParameter[param].nom;

    std::stringstream oss2;
    oss2 << F0_true;

    std::string name  = fParameter[param].name;
    std::string label = fParameter[param].label;
/*
    if (fParameter[param].name == "Nges") {
 	std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str()+".eps";
  	cout << oss.str() << endl;
     	MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename, input_channel);
	//MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, 400000, 440000, filename, input_channel);
    }
    else if (fParameter[param].name == "F0" || fParameter[param].name == "FL" || fParameter[param].name == "FR" ) {
	
    // parameter distribution
    std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename, input_channel);

    // error distribution
    filename = fOutputFolder+"/plot_"+input_channel+"_parerror_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename, input_channel);

    // pull distribution minuit error
    filename = fOutputFolder+"/plot_"+input_channel+"_pull_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename, input_channel);


    }
    else {

    // parameter distribution
    std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_"+name+"_"+oss2.str()+".eps";


    MakePlots(fParameter[param].hist_param, F0_true, oss.str(), label, fParameter[param].Low, fParameter[param].Up, filename, input_channel);

    // error distribution
    filename = fOutputFolder+"/plot_"+input_channel+"_parerror_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_unc, fParameter[param].Up_unc, filename, input_channel);

    // profile error distribution
    filename = fOutputFolder+"/plot_"+input_channel+"_proferror_"+name+"_"+oss2.str()+".eps";

    if(fParameter[param].Up_prof_unc > 35000) fParameter[param].Up_prof_unc = 35000;

    MakePlots(fParameter[param].hist_prof_unc, F0_true, oss.str(), "#sigma("+label+")", fParameter[param].Low_prof_unc, fParameter[param].Up_prof_unc, filename, input_channel);

    // pull distribution minuit error
    filename = fOutputFolder+"/plot_"+input_channel+"_pull_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_pull, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename, input_channel);

    // pull distribution profiling error   
    filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_"+name+"_"+oss2.str()+".eps";

    MakePlots(fParameter[param].hist_pull_prof, F0_true, oss.str(), "Pull distribution for "+label, -2.0, 2.0, filename, input_channel);
*/
    

  }

  //Handle covariance plots
  TH2D *h2_F0FL = new TH2D("h2_F0FL", "", 100, F0_true-0.2,F0_true+0.2, 100, FL_true-0.2,FL_true+0.2);
  TH2D *h2_F0FR = new TH2D("h2_F0FR", "", 100, F0_true-0.2,F0_true+0.2, 100, FR_true-0.2,FR_true+0.2);
  TH2D *h2_FLFR = new TH2D("h2_FLFR", "", 100, FL_true-0.2,FL_true+0.2, 100, FR_true-0.2,FR_true+0.2);
		
  TH2D *h2_N0NL = new TH2D("h2_N0NL", "", 100, xsec*lumi*F0_true-110000,xsec*lumi*F0_true+110000, 100, xsec*lumi*FL_true-110000,xsec*lumi*FL_true+110000);
  TH2D *h2_N0NR = new TH2D("h2_N0NR", "", 100, xsec*lumi*F0_true-110000,xsec*lumi*F0_true+110000, 100, xsec*lumi*FR_true-110000,xsec*lumi*FR_true+110000);
  TH2D *h2_NLNR = new TH2D("h2_NLNR", "", 100, xsec*lumi*FL_true-80000,xsec*lumi*FL_true+80000, 100, xsec*lumi*FR_true-80000,xsec*lumi*FR_true+80000);


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
  oss2 << F0_true;
/*
  std::string filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_F0FL_"+oss2.str()+".eps";
  Plot2Dcorrelation_F(h2_F0FL, F0_true, FL_true, label_x = "F_{0}", label_y = "F_{L}", filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_F0FR_"+oss2.str()+".eps";
  Plot2Dcorrelation_F(h2_F0FR, F0_true, FR_true, label_x = "F_{0}", label_y = "F_{R}", filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_FLFR_"+oss2.str()+".eps";
  Plot2Dcorrelation_F(h2_FLFR, FL_true, FR_true, label_x = "F_{L}", label_y = "F_{R}", filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_N0NL_"+oss2.str()+".eps";
  Plot2Dcorrelation_N(h2_N0NL, F0_true, label_x = "N_{0}", label_y = "N_{L}", filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_N0NR_"+oss2.str()+".eps";
  Plot2Dcorrelation_N(h2_N0NR, F0_true, label_x = "N_{0}", label_y = "N_{R}", filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_par_corr_NLNR_"+oss2.str()+".eps";
  Plot2Dcorrelation_N(h2_NLNR, F0_true, label_x = "N_{L}", label_y = "N_{R}", filename, input_channel);
*/
  //Calculate correlation
  std::vector<double> vec_corr;
  vec_corr.push_back(h2_N0NL->GetCorrelationFactor() );
  vec_corr.push_back(h2_N0NR->GetCorrelationFactor() );
  vec_corr.push_back(h2_NLNR->GetCorrelationFactor() );

  //Calculate new errormatrix based on profiling likelihood
  //Calculation for each PE

  //Histograms for errors on different helicity states
  TH1D *h_prof_err_F0 = new TH1D("h_prof_err_F0","",100,0.00,+0.12);
  TH1D *h_prof_err_FL = new TH1D("h_prof_err_FL","",100,0.00,+0.1);
  TH1D *h_prof_err_FR = new TH1D("h_prof_err_FR","",100,0.00,+0.1);
  //Helicity fractions
  TH1D *h_prof_pull_F0 = new TH1D("h_prof_pull_F0","",100,-8,+8);
  TH1D *h_prof_pull_FL = new TH1D("h_prof_pull_FL","",100,-8,+8);
  TH1D *h_prof_pull_FR = new TH1D("h_prof_pull_FR","",100,-8,+8);

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
 	//h_cov_N0NL->Fill(errormatrix[0][1]);
 	//h_cov_N0NR->Fill(errormatrix[0][2]);
 	//h_cov_NLNR->Fill(errormatrix[1][2]);

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
/*  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_F0_"+oss2.str()+".eps";
  MakePlots(h_prof_err_F0, F0_true, oss.str(), "#sigma(F0)", 0.0, 0.08, filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_FL_"+oss2.str()+".eps";
  MakePlots(h_prof_err_FL, F0_true, oss.str(), "#sigma(FL)", 0.0, 0.08, filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_proferror_FR_"+oss2.str()+".eps";
  MakePlots(h_prof_err_FR, F0_true, oss.str(), "#sigma(FR)", 0.0, 0.08, filename, input_channel);

  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_F0_"+oss2.str()+".eps";
  MakePlots(h_prof_pull_F0, F0_true, oss.str(), "Pull distribution for F0", -2.0, 2.0, filename, input_channel);
  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_FL_"+oss2.str()+".eps";
  MakePlots(h_prof_pull_FL, F0_true, oss.str(), "Pull distribution for FL", -2.0, 2.0, filename, input_channel);
  filename = fOutputFolder+"/plot_"+input_channel+"_prof_pull_FR_"+oss2.str()+".eps";
  MakePlots(h_prof_pull_FR, F0_true, oss.str(), "Pull distribution for FR", -2.0, 2.0, filename, input_channel);
*/
  
  ///____________________FILL CALIB INFO______________________________

  //N's
  for(int param = 0; param < fParameter.size(); ++param){

  if( fParameter[param].name == "N0" || fParameter[param].name == "NL" || fParameter[param].name == "NR" ) {
	 
	calib_info.push_back( fParameter[param].hist_param->GetMean() );
	calib_info.push_back( fParameter[param].hist_param->GetMeanError() );
	calib_info.push_back( fParameter[param].hist_pull_prof->GetMean()); 
	calib_info.push_back( fParameter[param].hist_pull_prof->GetMeanError());
	calib_info.push_back( fParameter[param].hist_pull_prof->GetRMS());
	calib_info.push_back( fParameter[param].hist_pull_prof->GetRMSError());
  }
  
  }	

  //F's
  for(int param = 0; param < fParameter.size(); ++param){

  if( fParameter[param].name == "F0") {
	
	calib_info.push_back( fParameter[param].hist_param->GetMean() );
	calib_info.push_back( fParameter[param].hist_param->GetMeanError() );
	calib_info.push_back( h_prof_pull_F0->GetMean() );
	calib_info.push_back( h_prof_pull_F0->GetMeanError() );
	calib_info.push_back( h_prof_pull_F0->GetRMS() );
	calib_info.push_back( h_prof_pull_F0->GetRMSError() );
  }	
  else if ( fParameter[param].name == "FL") {

	calib_info.push_back( fParameter[param].hist_param->GetMean() );
	calib_info.push_back( fParameter[param].hist_param->GetMeanError() );
	calib_info.push_back( h_prof_pull_FL->GetMean() );
	calib_info.push_back( h_prof_pull_FL->GetMeanError() );
	calib_info.push_back( h_prof_pull_FL->GetRMS() );
	calib_info.push_back( h_prof_pull_FL->GetRMSError() );
  }
  else if ( fParameter[param].name == "FR") {

	calib_info.push_back( fParameter[param].hist_param->GetMean() );
	calib_info.push_back( fParameter[param].hist_param->GetMeanError() );
	calib_info.push_back( h_prof_pull_FR->GetMean() );
	calib_info.push_back( h_prof_pull_FR->GetMeanError() );
	calib_info.push_back( h_prof_pull_FR->GetRMS() );
	calib_info.push_back( h_prof_pull_FR->GetRMSError() );
  } 
  }


  for(int param = 0; param < fParameter.size(); ++param){

    //std::string label = fParameter[param].label;

    //std::cout << fParameter[param].Low_prof_unc << "\t" << fParameter[param].Up_prof_unc << std::endl;

    fParameter[param].hist_param->Delete();
    fParameter[param].hist_unc->Delete();
    fParameter[param].hist_prof_unc->Delete();  
    fParameter[param].hist_pull->Delete();    
    fParameter[param].hist_pull_prof->Delete(); 

  }

  h2_F0FL->Delete();
  h2_F0FR->Delete();
  h2_FLFR->Delete();
  h2_N0NL->Delete();
  h2_N0NR->Delete();
  h2_NLNR->Delete();
  h_prof_err_F0->Delete();
  h_prof_err_FL->Delete();
  h_prof_err_FR->Delete();
  h_prof_pull_F0->Delete();
  h_prof_pull_FL->Delete();
  h_prof_pull_FR->Delete();


  return calib_info;

}

///________________________________MAKE PLOTS___________________________
void histograms::MakePlots(TH1D *hist, double F0_true, std::string nominal, std::string label, double low, double up, std::string filename, std::string input_channel, int mode)
{
  TF1 *fit = new TF1(label.c_str(), "gaus", low, up);

  hist -> Fit(label.c_str(), "RQ"); // "R" means range

  int max_bin = hist -> GetMaximumBin();
  double max  = hist -> GetBinContent(max_bin);

  // hardcoded for a first test!!!                                                                                                                               
  //std::string input_channel = "mu";

  //cout << "low and up position " << up << endl;

  PlotDistribution(hist, F0_true, label.c_str(), input_channel + "_pseudo_"+label+"=%.2f", max*1.2, input_channel, nominal, fit, filename, low, up, mode);

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

	for (int i=0; i != num_of_nui_parameters; i++)
	{
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

	if (i == 0) {
	central = fInputTree -> nui_0; unc = fInputTree -> nui_err_0; prof_unc = fInputTree -> nui_prof_err_0; nominal = fInputTree -> nui_nom_0; }
	if (i == 1) {
	central = fInputTree -> nui_1; unc = fInputTree -> nui_err_1; prof_unc = fInputTree -> nui_prof_err_1; nominal = fInputTree -> nui_nom_1; }
	if (i == 2) {
	central = fInputTree -> nui_2; unc = fInputTree -> nui_err_2; prof_unc = fInputTree -> nui_prof_err_2; nominal = fInputTree -> nui_nom_2; }
	if (i == 3) {
	central = fInputTree -> nui_3; unc = fInputTree -> nui_err_3; prof_unc = fInputTree -> nui_prof_err_3; nominal = fInputTree -> nui_nom_3; }
	if (i == 4) {
	central = fInputTree -> nui_4; unc = fInputTree -> nui_err_4; prof_unc = fInputTree -> nui_prof_err_4; nominal = fInputTree -> nui_nom_4; }
	if (i == 5) {
	central = fInputTree -> nui_5; unc = fInputTree -> nui_err_5; prof_unc = fInputTree -> nui_prof_err_5; nominal = fInputTree -> nui_nom_5; }		

      	}

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
      if(counter_lines > 1){

	std::stringstream oss;
	oss << epdf->d_name;

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
