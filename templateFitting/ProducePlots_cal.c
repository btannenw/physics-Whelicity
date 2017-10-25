#include <unistd.h>
#include <dirent.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraph.h"
#include "TMinuit.h"
//#include "OutputTreeReader.h"

#include "fill_histograms_cal.h"
#include "plots.h"

using namespace std;

int main(int argc, char *argv[])
{

  if(argc != 7){

    cout << "\t" << endl;
    cout << "ERROR: Number of input variables is wrong!!!" << endl;
    cout << "\t" << endl;

    return 1;
  }
  else{

    cout << endl;
    cout << "Input Folder:    "               << argv[1] << endl;
    cout << "Input channel:   "               << argv[2] << endl;
    cout << "Configuration:   "               << argv[3] << endl;
    cout << "Number of nuisance parameters: " << argv[4] << endl;
    cout << "Validation mode: "               << argv[5] << endl;
    cout << "Output Folder:   "               << argv[6] << endl;
    cout << endl;
  }

  // find current directory                                                                                                                                          
  char path_cwd[MAXPATHLEN];     // aus Inernet geklaut                                                                                                            
  getcwd(path_cwd, MAXPATHLEN);

  stringstream cwd;
  cwd << path_cwd;

  //std::string fInput            = cwd.str()+"/"+argv[1]+"/";
  std::string input_channel     = argv[2];
  std::string scounter		= argv[3];
  std::string snumber		= argv[4];
  std::string validation_mode   = argv[5];
  std::string fOutput           = cwd.str()+"/"+argv[6]+"/";

  int counter = atoi(scounter.c_str());
  int num_of_nui_parameters = atoi(snumber.c_str());

  double xsec = 89.2;//pb  //better 89.45 ... check later
  double lumi = 4655.73;//i_pb

  std::vector<std::string> fInput;

  for (int k=0; k != counter+1; k++){

	std::stringstream s_count;
        s_count << k;
	
	fInput.push_back(cwd.str()+"/Output_Sample_"+input_channel+"_"+s_count.str()+"/");
  
	cout << fInput[k] << endl;

	// find out if folder exists
	if ( !opendir(fInput[k].c_str()) ){
   
	    std::cout << fInput[k] << " does not exist!!!" << std::endl;
 
	return false;
  	}
  }

  //std::vector<histograms*> fNewHisto(counter);
  histograms *fNewHisto = new histograms();

 //std::vector<histograms::Parameter> nui;
  histograms::Parameter nui[num_of_nui_parameters][counter+1];

  // F0, FL, FR, N0, NL, NR, Wjets, QCD, RemBkg, nui0, nui1
  //as well as additional parameter for number of sugnal evenst Nges

  cout << "Define parameters " << endl;

  histograms::Parameter N0[counter+1],NL[counter+1],NR[counter+1],Wjets[counter+1],QCD[counter+1],RemBkg[counter+1];//,nui0;//, nui1;
  histograms::Parameter F0[counter+1],FL[counter+1],FR[counter+1];
  histograms::Parameter Nges[counter+1];

  //std::vector< std::vector<histograms::Parameter> > ParameterList(counter);
  std::vector<histograms::Parameter> ParameterList;


  std::vector< std::vector<double> > calib_info;

  ///Calibration graphs
	
	//Define calibration graphs
	TGraphErrors *graph_F0   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FL   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FR   = new TGraphErrors(counter+1);

	TGraphErrors *graph_N0   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NL   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NR   = new TGraphErrors(counter+1);
	TGraphErrors *graph_nui1 = new TGraphErrors(counter+1);


	//TGraphErrors *graph_Wjets = new TGraphErrors(counter+1);
	//TGraphErrors *graph_QCD = new TGraphErrors(counter+1);
	//TGraphErrors *graph_RemBkg = new TGraphErrors(counter+1);
	
	//Define graphs for differences
	TGraphErrors *diff_F0 = new TGraphErrors(counter+1);
	TGraphErrors *diff_FL = new TGraphErrors(counter+1);
	TGraphErrors *diff_FR = new TGraphErrors(counter+1);
	
	TGraphErrors *diff_N0 = new TGraphErrors(counter+1);
	TGraphErrors *diff_NL = new TGraphErrors(counter+1);
	TGraphErrors *diff_NR = new TGraphErrors(counter+1);
	
	TGraphErrors *diff_nui1 = new TGraphErrors(counter+1);


	//TGraphErrors *diff_Wjets = new TGraphErrors(counter+1);
	//TGraphErrors *diff_QCD = new TGraphErrors(counter+1);
	//TGraphErrors *diff_RemBkg = new TGraphErrors(counter+1);
	
	//Define calibration pull graphs
	//Mean+RMS values
	TGraphErrors *graph_F0_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FL_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FR_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_N0_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NL_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NR_comb   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NR_nui1   = new TGraphErrors(counter+1);

	//Mean+MeanError
	TGraphErrors *graph_F0_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FL_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FR_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_N0_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NL_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NR_mean   = new TGraphErrors(counter+1);
	TGraphErrors *graph_nui1_mean = new TGraphErrors(counter+1);
	
	//RMS+RMSerror
	TGraphErrors *graph_F0_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FL_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_FR_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_N0_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NL_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_NR_RMS   = new TGraphErrors(counter+1);
	TGraphErrors *graph_nui1_RMS = new TGraphErrors(counter+1);

	double F0_true, FL_true, FR_true;
	
	if (counter == 0) {F0_true = 0.7; FL_true=0.3; FR_true=0.0;}
	else {F0_true = 0.7 - counter*0.1; FL_true=0.3 + counter*0.05 ; FR_true=0.0 + counter*0.05;}	
	
	if (!opendir(fOutput.c_str())){
	  
	  std::cout << "Folder does not exist, make new folder: " << fOutput << std::endl;
	  
	  mkdir(fOutput.c_str(), 0777);
	}
	
	for (int k=0; k != counter+1; k++){
	  
	  fNewHisto = new histograms();
	  fNewHisto -> SetInputFolder(fInput[k]);
	  
	  
	  fNewHisto -> SetOutputFolder(fOutput);
	  
	  cout << "make chain of files " << endl;
	  
	  fNewHisto -> MakeChainOfFiles(num_of_nui_parameters);
	  
	  if (k == 0) {F0_true = 0.7; FL_true=0.3; FR_true=0.0;}
	  else {F0_true = 0.7 - k*0.1; FL_true=0.3 + k*0.05 ; FR_true=0.0 + k*0.05;}	
	  
	  //cout << "Start loop " << k << endl;
	  
	  F0[k].name     = "F0";      F0[k].label     = "F0";
	  FL[k].name     = "FL";      FL[k].label     = "FL";
	  FR[k].name     = "FR";      FR[k].label     = "FR";
	  N0[k].name     = "N0";      N0[k].label     = "N0";
	  NL[k].name     = "NL";      NL[k].label     = "NL";
	  NR[k].name     = "NR";      NR[k].label     = "NR";
	  Nges[k].name   = "Nges";    Nges[k].label   = "Nges";
	  Wjets[k].name  = "Wjets";   Wjets[k].label  = "Wjets";
	  QCD[k].name    = "QCD";     QCD[k].label    = "QCD";
	  RemBkg[k].name = "RemBkg";  RemBkg[k].label = "RemBkg";
	  
	  for (int i=0; i != num_of_nui_parameters; i++)
	    {
	      std::stringstream nui_label;
	      nui_label << i;
	      nui[i][k].name  = "nui_"+nui_label.str() ; 
	      nui[i][k].label = "Nuisance parameter "+nui_label.str() ;
	      //nui1.name   = "nui_1";   nui1.label   = "Nuisance param. 1";
	    }
	  
	  F0[k].Up     =  0.0;      F0[k].Low     = 1.0;
	  FL[k].Up     =  0.0;      FL[k].Low     = 1.0;
	  FR[k].Up     =  0.0;      FR[k].Low     = 1.0;
  	N0[k].Up     =  0.0;      N0[k].Low     = 10e5;
  	NL[k].Up     =  0.0;      NL[k].Low     = 10e5;
  	NR[k].Up     =  0.0;      NR[k].Low     = 10e5;
  	Nges[k].Up   =  0.0;      Nges[k].Low   = 10e5;
  	Wjets[k].Up  =  0.0;      Wjets[k].Low  = 10e5;
  	QCD[k].Up    =  0.0;      QCD[k].Low    = 10e5;
  	RemBkg[k].Up =  0.0;      RemBkg[k].Low = 10e5;
	
	for (int i=0; i != num_of_nui_parameters; i++)
	{
		nui[i][k].Up = -0.001; nui[i][k].Low = 0.001;
	}
	
	F0[k].Up_unc     =  0.0;      F0[k].Low_unc     = 10.0;
	FL[k].Up_unc     =  0.0;      FL[k].Low_unc     = 10.0;
  	FR[k].Up_unc     =  0.0;      FR[k].Low_unc     = 10.0;
  	N0[k].Up_unc     =  0.0;      N0[k].Low_unc     = 10e5;
  	NL[k].Up_unc     =  0.0;      NL[k].Low_unc     = 10e5;
  	NR[k].Up_unc     =  0.0;      NR[k].Low_unc     = 10e5;
  	Nges[k].Up_unc   =  0.0;      Nges[k].Low_unc   = 10e5;
  	Wjets[k].Up_unc  =  0.0;      Wjets[k].Low_unc  = 10e5;
  	QCD[k].Up_unc    =  0.0;      QCD[k].Low_unc    = 10e5;
  	RemBkg[k].Up_unc =  0.0;      RemBkg[k].Low_unc = 10e5;
	
	for (int i=0; i != num_of_nui_parameters; i++)
  	{
		nui[i][k].Up_unc = -0.001; nui[i][k].Low_unc = 0.001;
  	}

  	F0[k].Up_prof_unc     =  0.0;      F0[k].Low_prof_unc     = 10.0;
  	FL[k].Up_prof_unc     =  0.0;      FL[k].Low_prof_unc     = 10.0;
  	FR[k].Up_prof_unc     =  0.0;      FR[k].Low_prof_unc     = 10.0;
  	N0[k].Up_prof_unc     =  0.0;      N0[k].Low_prof_unc     = 10e5;
  	NL[k].Up_prof_unc     =  0.0;      NL[k].Low_prof_unc     = 10e5;
  	NR[k].Up_prof_unc     =  0.0;      NR[k].Low_prof_unc     = 10e5;
  	Nges[k].Up_prof_unc   =  0.0;      Nges[k].Low_prof_unc   = 10e5;
  	Wjets[k].Up_prof_unc  =  0.0;      Wjets[k].Low_prof_unc  = 10e5;
  	QCD[k].Up_prof_unc    =  0.0;      QCD[k].Low_prof_unc    = 10e5;
  	RemBkg[k].Up_prof_unc =  0.0;      RemBkg[k].Low_prof_unc = 10e5;
	
  	for (int i=0; i != num_of_nui_parameters; i++)
  	{
		nui[i][k].Up_prof_unc = -0.001; nui[i][k].Low_prof_unc = 0.001;
  	}

  

  	//cout << "Define Parameter list" << endl;

	ParameterList.clear();

	ParameterList.push_back(F0[k]);
  	ParameterList.push_back(FL[k]);
  	ParameterList.push_back(FR[k]);
  	ParameterList.push_back(N0[k]);
  	ParameterList.push_back(NL[k]);
  	ParameterList.push_back(NR[k]);
  	ParameterList.push_back(Nges[k]);
  	ParameterList.push_back(Wjets[k]);
  	ParameterList.push_back(QCD[k]);
  	ParameterList.push_back(RemBkg[k]);
  
	//cout << "Define Parameter list end" << endl;

  	for (int i=0; i != num_of_nui_parameters; i++)
  	{
		ParameterList.push_back(nui[i][k]);
  	}
  
  	//ParameterList.push_back(nui0);
  	//  ParameterList.push_back(nui1);

  	//cout << "set parameter list " << endl;

  	fNewHisto -> SetParameterList(ParameterList);
  
  	//fNewHisto -> FillHistograms(counter, input_channel, num_of_nui_parameters);


  	calib_info.push_back ( fNewHisto->FillHistograms(k, input_channel, num_of_nui_parameters, xsec, lumi) );

	cout << endl;
 	cout << "Calib info" << endl;
  	cout << endl;
	for (int h = 0; h < 36; h++) {
		cout << calib_info[k][h] << endl;
	}
	cout << endl;

	/*cout << calib_info[k][0] << endl;
 	cout << calib_info[k][6] << endl;
 	cout << calib_info[k][12] << endl;
  	cout << calib_info[k][18] << endl;
	cout << calib_info[k][24] << endl;
 	cout << calib_info[k][30] << endl;
	*/
  

  	//Make calibration plots with information stored in calib_info


	
   	//Fill calibration graphs with help of calib_info:
   ///______________________Fill calibration craphs__________________

   

	graph_F0->SetPoint( k+1,F0_true,calib_info[k][18]);
	graph_FL->SetPoint( k+1,FL_true,calib_info[k][24]);
	graph_FR->SetPoint( k+1,FR_true,calib_info[k][30]);
	
	graph_N0->SetPoint( k+1,xsec*lumi*F0_true,calib_info[k][0]);
	graph_NL->SetPoint( k+1,xsec*lumi*FL_true,calib_info[k][6]);
	graph_NR->SetPoint( k+1,xsec*lumi*FR_true,calib_info[k][12]);
	
	//graph_Wjets->SetPoint( k+1,k+1,h_Wjets->GetMean());
	//graph_QCD->SetPoint( k+1,k+1,h_QCD->GetMean());
	//graph_RemBkg->SetPoint( k+1,k+1,h_RemBkg->GetMean());
	
	graph_F0->SetPointError( k+1,0.0, calib_info[k][19] );
	graph_FL->SetPointError( k+1,0.0, calib_info[k][25] );
	graph_FR->SetPointError( k+1,0.0, calib_info[k][31] );
	
	graph_N0->SetPointError( k+1,0.0, calib_info[k][1] );
	graph_NL->SetPointError( k+1,0.0, calib_info[k][7] );
	graph_NR->SetPointError( k+1,0.0, calib_info[k][13] );
	
	//graph_Wjets->SetPointError( k+1,0.0,h_Wjets->GetMeanError());
	//graph_QCD->SetPointError( k+1,0.0,h_QCD->GetMeanError());
	//graph_RemBkg->SetPointError( k+1,0.0,h_RemBkg->GetMeanError());
	
	diff_F0->SetPoint( k+1, F0_true, (calib_info[k][18]-F0_true)/F0_true );
	diff_FL->SetPoint( k+1, FL_true, (calib_info[k][24]-FL_true)/FL_true );
	diff_FR->SetPoint( k+1, FR_true, (calib_info[k][30]-FR_true)/FR_true );
	//diff_F0->SetPoint( k+1, F0_true, F0_mean/F0_true);
	//diff_FL->SetPoint( k+1, FL_true, FL_mean/FL_true);
	//diff_FR->SetPoint( k+1, FR_true, FR_mean/FR_true);	


	diff_N0->SetPoint( k+1, xsec*lumi*F0_true, (calib_info[k][0] -xsec*lumi*F0_true)/(xsec*lumi*F0_true) );
	diff_NL->SetPoint( k+1, xsec*lumi*FL_true, (calib_info[k][6] -xsec*lumi*FL_true)/(xsec*lumi*FL_true) );
	diff_NR->SetPoint( k+1, xsec*lumi*FR_true, (calib_info[k][12] -xsec*lumi*FR_true)/(xsec*lumi*FR_true) );
	//diff_N0->SetPoint( k+1, xsec*lumi*F0_true, h_N0->GetMean()/xsec*lumi*F0_true);
	//diff_NL->SetPoint( k+1, xsec*lumi*FL_true, h_NL->GetMean() /xsec*lumi*FL_true);
	//diff_NR->SetPoint( k+1, xsec*lumi*FR_true, h_NR->GetMean() /xsec*lumi*FR_true);	

	//diff_Wjets->SetPoint( k+1, k+1, h_Wjets->GetMean()-_bkg_norm_vec[0]);
	//diff_QCD->SetPoint( k+1, k+1, h_QCD->GetMean()-_bkg_norm_vec[1]);
	//diff_RemBkg->SetPoint( k+1, k+1, h_RemBkg->GetMean()-_bkg_norm_vec[2]);
	//diff_Wjets->SetPoint( k+1, k+1, h_Wjets->GetMean()/_bkg_norm_vec[0]);
	//diff_QCD->SetPoint( k+1, k+1, h_QCD->GetMean()/_bkg_norm_vec[1]);
	//diff_RemBkg->SetPoint( k+1, k+1, h_RemBkg->GetMean()/_bkg_norm_vec[2]);
	
	diff_F0->SetPointError( k+1,0.0, calib_info[k][19]/F0_true );
	diff_FL->SetPointError( k+1,0.0, calib_info[k][25]/FL_true );
	diff_FR->SetPointError( k+1,0.0, calib_info[k][31]/FR_true );
	
	diff_N0->SetPointError( k+1,0.0, calib_info[k][1]/(xsec*lumi*F0_true ) );
	diff_NL->SetPointError( k+1,0.0, calib_info[k][7]/(xsec*lumi*FL_true ) );
	diff_NR->SetPointError( k+1,0.0, calib_info[k][13]/(xsec*lumi*FR_true ) );
	
	//cout << "Mean error N0 : " << h_N0->GetMeanError() << endl;
	//cout << "Mean error NL : " << h_NL->GetMeanError() << endl;
	//cout << "Mean error NR : " << h_NR->GetMeanError() << endl;

	//diff_Wjets->SetPointError( k+1,0.0,h_Wjets->GetMeanError());
	//diff_QCD->SetPointError( k+1,0.0,h_QCD->GetMeanError());
	//diff_RemBkg->SetPointError( k+1,0.0,h_RemBkg->GetMeanError());
	
	//Pull calibration curves
	double F0_pull_mean = calib_info[k][20];
	double FL_pull_mean = calib_info[k][26];
	double FR_pull_mean = calib_info[k][32];
	double F0_pull_RMS = calib_info[k][22];
	double FL_pull_RMS = calib_info[k][28];
	double FR_pull_RMS = calib_info[k][34];
	
	double N0_pull_mean = calib_info[k][2];
	double NL_pull_mean = calib_info[k][8];
	double NR_pull_mean = calib_info[k][14];
	double N0_pull_RMS = calib_info[k][4];
	double NL_pull_RMS = calib_info[k][10];
	double NR_pull_RMS = calib_info[k][16];
	
	
	//cout << "Mean pull values: " <<  h_F0_pull->GetMean() << " " << h_FL_pull->GetMean() << " " << h_FR_pull->GetMean() << endl;
	//cout << "RMS pull values: " << h_F0_pull->GetRMS() << " " << h_FL_pull->GetRMS() << " " << h_FR_pull->GetRMS() << endl;
	
	
	//Set points of TGraph
	//Mean+RMS
	graph_F0_comb->SetPoint( k+1,F0_true,F0_pull_mean);
	graph_FL_comb->SetPoint( k+1,FL_true,FL_pull_mean);
	graph_FR_comb->SetPoint( k+1,FR_true,FR_pull_mean);
	graph_N0_comb->SetPoint( k+1,xsec*lumi*F0_true,N0_pull_mean);
	graph_NL_comb->SetPoint( k+1,xsec*lumi*FL_true,NL_pull_mean);
	graph_NR_comb->SetPoint( k+1,xsec*lumi*FR_true,NR_pull_mean);
	//Mean+MeanError
	graph_F0_mean->SetPoint( k+1,F0_true,F0_pull_mean);
	graph_FL_mean->SetPoint( k+1,FL_true,FL_pull_mean);
	graph_FR_mean->SetPoint( k+1,FR_true,FR_pull_mean);
	graph_N0_mean->SetPoint( k+1,xsec*lumi*F0_true,N0_pull_mean);
	graph_NL_mean->SetPoint( k+1,xsec*lumi*FL_true,NL_pull_mean);
	graph_NR_mean->SetPoint( k+1,xsec*lumi*FR_true,NR_pull_mean);
	//RMS+RMSerror
	graph_F0_RMS->SetPoint( k+1,F0_true,F0_pull_RMS);
	graph_FL_RMS->SetPoint( k+1,FL_true,FL_pull_RMS);
	graph_FR_RMS->SetPoint( k+1,FR_true,FR_pull_RMS);
	graph_N0_RMS->SetPoint( k+1,xsec*lumi*F0_true,N0_pull_RMS);
	graph_NL_RMS->SetPoint( k+1,xsec*lumi*FL_true,NL_pull_RMS);
	graph_NR_RMS->SetPoint( k+1,xsec*lumi*FR_true,NR_pull_RMS);
	
	//Set errors of TGraph
	//Mean+RMS
	graph_F0_comb->SetPointError( k+1,0.0,F0_pull_RMS);
	graph_FL_comb->SetPointError( k+1,0.0,FL_pull_RMS);
	graph_FR_comb->SetPointError( k+1,0.0,FR_pull_RMS);
	graph_N0_comb->SetPointError( k+1,0.0,N0_pull_RMS);
	graph_NL_comb->SetPointError( k+1,0.0,NL_pull_RMS);
	graph_NR_comb->SetPointError( k+1,0.0,NR_pull_RMS);
	//Mean+MeanError
	graph_F0_mean->SetPointError( k+1,0.0,calib_info[k][21]);
	graph_FL_mean->SetPointError( k+1,0.0,calib_info[k][27]);
	graph_FR_mean->SetPointError( k+1,0.0,calib_info[k][33]);
	graph_N0_mean->SetPointError( k+1,0.0,calib_info[k][3]);
	graph_NL_mean->SetPointError( k+1,0.0,calib_info[k][9]);
	graph_NR_mean->SetPointError( k+1,0.0,calib_info[k][15]);
	//RMS+RMSerror
	graph_F0_RMS->SetPointError( k+1,0.0,calib_info[k][23]);
	graph_FL_RMS->SetPointError( k+1,0.0,calib_info[k][29]);
	graph_FR_RMS->SetPointError( k+1,0.0,calib_info[k][35]);
	graph_N0_RMS->SetPointError( k+1,0.0,calib_info[k][5]);
	graph_NL_RMS->SetPointError( k+1,0.0,calib_info[k][11]);
	graph_NR_RMS->SetPointError( k+1,0.0,calib_info[k][17]);
	

	delete fNewHisto;

   }

///______________________Calibration plots_____________________

	double min, max, xmin, xmax, ymin, ymax;
	int mode;
	std::string label;
	std::string label_x;
	std::string label_y;
	std::string title;
	std::string string_mode;

	PlotCalDistribution(graph_F0, diff_F0, label = "F_{0}", title = input_channel + "_calibration_F0", min = -0.05, max = 0.9, input_channel, mode=0);
	PlotCalDistribution(graph_FL, diff_FL, label = "F_{L}", title = input_channel + "_calibration_FL", min = 0.25, max = 0.8, input_channel, mode=0);
	PlotCalDistribution(graph_FR, diff_FR, label = "F_{R}", title = input_channel +  "_calibration_FR", min = -0.05, max = 0.5, input_channel, mode=0);	
	//PlotCalDistribution(graph_Wjets, diff_Wjets, label = "Wjets", title = input_channel +  "_calibration_Wjets", min = _bkg_norm_vec[0]-200, max = _bkg_norm_vec[0]+200, input_channel, mode=1);
	//PlotCalDistribution(graph_QCD, diff_QCD, label = "QCD", title = input_channel +  "_calibration_QCD", min = _bkg_norm_vec[1]-200, max = _bkg_norm_vec[1]+200, input_channel, mode=1);
	//PlotCalDistribution(graph_RemBkg, diff_RemBkg, label = "RemBkg", title = input_channel +  "_calibration_RemBkg", min = _bkg_norm_vec[2]-200, max = _bkg_norm_vec[2]+200, input_channel, mode=1);

	PlotCalDistribution_N(graph_N0, diff_N0, label = "N_{0}", title = input_channel + "_calibration_N0", min = -0.05*xsec*lumi, max = xsec*lumi*0.9, input_channel);
	PlotCalDistribution_N(graph_NL, diff_NL, label = "N_{L}", title = input_channel + "_calibration_NL", min = 0.25*xsec*lumi, max = xsec*lumi*0.8, input_channel);
	PlotCalDistribution_N(graph_NR, diff_NR, label = "N_{R}", title = input_channel +  "_calibration_NR", min = -0.05*xsec*lumi, max = xsec*lumi*0.5, input_channel);	


	//Mean+RMS
	PlotCalPullDistribution(graph_F0_comb, label = "F_{0}", title = input_channel + "_calibration_pull_F0_comb", xmin = -0.05, xmax = 0.75, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	PlotCalPullDistribution(graph_FL_comb, label = "F_{L}", title = input_channel +  "_calibration_pull_FL_comb", xmin = 0.25, xmax = 0.7, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	PlotCalPullDistribution(graph_FR_comb, label = "F_{R}", title = input_channel +  "_calibration_pull_FR_comb", xmin = -0.05, xmax = 0.4, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	PlotCalPullDistribution(graph_N0_comb, label = "N_{0}", title = input_channel + "_calibration_pull_N0_comb", xmin = -0.05*xsec*lumi, xmax = 0.75*xsec*lumi, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	PlotCalPullDistribution(graph_NL_comb, label = "N_{L}", title = input_channel +  "_calibration_pull_NL_comb", xmin = 0.25*xsec*lumi, xmax = 0.7*xsec*lumi, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	PlotCalPullDistribution(graph_NR_comb, label = "N_{R}", title = input_channel +  "_calibration_pull_NR_comb", xmin = -0.05*xsec*lumi, xmax = 0.4*xsec*lumi, ymin=-1.5, ymax=1.5, input_channel, string_mode="Mean+RMS");
	//Mean+MeanError
	PlotCalPullDistribution(graph_F0_mean, label = "F_{0}", title = input_channel + "_calibration_pull_F0_mean", xmin = -0.05, xmax = 0.75, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	PlotCalPullDistribution(graph_FL_mean, label = "F_{L}", title = input_channel +  "_calibration_pull_FL_mean", xmin = 0.25, xmax = 0.7, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	PlotCalPullDistribution(graph_FR_mean, label = "F_{R}", title = input_channel +  "_calibration_pull_FR_mean", xmin = -0.05, xmax = 0.4, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	PlotCalPullDistribution(graph_N0_mean, label = "N_{0}", title = input_channel + "_calibration_pull_N0_mean", xmin = -0.05*xsec*lumi, xmax = 0.75*xsec*lumi, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	PlotCalPullDistribution(graph_NL_mean, label = "N_{L}", title = input_channel +  "_calibration_pull_NL_mean", xmin = 0.25*xsec*lumi, xmax = 0.7*xsec*lumi, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	PlotCalPullDistribution(graph_NR_mean, label = "N_{R}", title = input_channel +  "_calibration_pull_NR_mean", xmin = -0.05*xsec*lumi, xmax = 0.4*xsec*lumi, ymin=-0.1, ymax=0.1, input_channel, string_mode="Mean+MeanError");
	//RMS+RMSError
	PlotCalPullDistribution(graph_F0_RMS, label = "F_{0}", title = input_channel + "_calibration_pull_F0_RMS", xmin = -0.05, xmax = 0.75, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");
	PlotCalPullDistribution(graph_FL_RMS, label = "F_{L}", title = input_channel +  "_calibration_pull_FL_RMS", xmin = 0.25, xmax = 0.7, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");
	PlotCalPullDistribution(graph_FR_RMS, label = "F_{R}", title = input_channel +  "_calibration_pull_FR_RMS", xmin = -0.05, xmax = 0.4, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");
	PlotCalPullDistribution(graph_N0_RMS, label = "N_{0}", title = input_channel + "_calibration_pull_N0_RMS", xmin = -0.05*xsec*lumi, xmax = 0.75*xsec*lumi, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");
	PlotCalPullDistribution(graph_NL_RMS, label = "N_{L}", title = input_channel +  "_calibration_pull_NL_RMS", xmin = 0.25*xsec*lumi, xmax = 0.7*xsec*lumi, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");
	PlotCalPullDistribution(graph_NR_RMS, label = "N_{R}", title = input_channel +  "_calibration_pull_NR_RMS", xmin = -0.05*xsec*lumi, xmax = 0.4*xsec*lumi, ymin=0.7, ymax=1.1, input_channel, string_mode="RMS+RMSError");



	graph_F0->Delete();	
	graph_FL->Delete();
	graph_FR->Delete();
	graph_N0->Delete();	
	graph_NL->Delete();
	graph_NR->Delete();
	//graph_Wjets->Delete();	
	//graph_QCD->Delete();
	//graph_RemBkg->Delete();
	
	diff_F0->Delete();
	diff_FL->Delete();
	diff_FR->Delete();
	diff_N0->Delete();
	diff_NL->Delete();
	diff_NR->Delete();
	//diff_Wjets->Delete();
	//diff_QCD->Delete();
	//diff_RemBkg->Delete();

	graph_F0_comb->Delete();
	graph_FL_comb->Delete();
	graph_FR_comb->Delete();
	graph_N0_comb->Delete();
	graph_NL_comb->Delete();
	graph_NR_comb->Delete();
	graph_F0_mean->Delete();
	graph_FL_mean->Delete();
	graph_FR_mean->Delete();
	graph_N0_mean->Delete();
	graph_NL_mean->Delete();
	graph_NR_mean->Delete();
	graph_F0_RMS->Delete();
	graph_FL_RMS->Delete();
	graph_FR_RMS->Delete();
	graph_N0_RMS->Delete();
	graph_NL_RMS->Delete();
	graph_NR_RMS->Delete();

  	return 0 ;


}
