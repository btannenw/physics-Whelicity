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

#include "fill_histograms.h"

using namespace std;

int main(int argc, char *argv[])
{

  std::cout << "Number of input parameters   " << argc << std::endl;

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
    cout << "F0 value         "               << argv[3] << endl;
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

  std::string fInput            = cwd.str()+"/"+argv[1]+"/";
  std::string input_channel     = argv[2];
  std::string F0_value		= argv[3];
  std::string number_nui	= argv[4];
  std::string validation_mode   = argv[5];
  std::string fOutput           = cwd.str()+"/"+argv[6]+"/";

  std::string htmltitle         = "index_pull_"+input_channel+"_"+F0_value+"_"+number_nui+"_"+validation_mode+".html";

  //  int counter = atoi(scounter.c_str());
  int num_of_nui_parameters = atoi(number_nui.c_str());

  // find out if folder exists
  if ( !opendir(fInput.c_str()) ){
   
    std::cout << fInput << " does not exist!!!" << std::endl;
 
    return false;

  }

  histograms *fNewHisto = new histograms();

  fNewHisto -> SetInputFolder(fInput);

  if (!opendir(fOutput.c_str())){
    
    std::cout << "Folder does not exist, make new folder: " << fOutput << std::endl;
    
    mkdir(fOutput.c_str(), 0777);

  }

  fNewHisto -> SetOutputFolder(fOutput);

  fNewHisto -> SetValidationMode(validation_mode);

  cout << "make chain of files " << endl;

  fNewHisto -> MakeChainOfFiles(num_of_nui_parameters);

  histograms::Parameter nui[num_of_nui_parameters];

  // F0, FL, FR, N0, NL, NR, Wjets, QCD, RemBkg, nui0, nui1
  //as well as additional parameter for number of sugnal evenst Nges

  cout << "Define parameters " << endl;

  histograms::Parameter N0,NL,NR,Wjets,QCD,RemBkg,nui0;//, nui1;
  histograms::Parameter F0,FL,FR;
  histograms::Parameter Nges;

  F0.name     = "F0";      F0.label     = "F0";
  FL.name     = "FL";      FL.label     = "FL";
  FR.name     = "FR";      FR.label     = "FR";
  N0.name     = "N0";      N0.label     = "N0";
  NL.name     = "NL";      NL.label     = "NL";
  NR.name     = "NR";      NR.label     = "NR";
  Nges.name   = "Nges";    Nges.label   = "Nges";
  Wjets.name  = "Wjets";   Wjets.label  = "Wjets";
  QCD.name    = "QCD";     QCD.label    = "QCD";
  RemBkg.name = "RemBkg";  RemBkg.label = "RemBkg";

  for (int i=0; i != num_of_nui_parameters; i++)
  {

   // std::cout << "nui" << "\t" << i << std::endl;

	std::stringstream nui_label;
	nui_label << i;
	nui[i].name  = "nui_"+nui_label.str() ; 
	nui[i].label = "Nuisance parameter "+nui_label.str() ;

  }

  F0.Up     =  0.0;      F0.Low     = 1.0;
  FL.Up     =  0.0;      FL.Low     = 1.0;
  FR.Up     =  0.0;      FR.Low     = 1.0;
  N0.Up     =  0.0;      N0.Low     = 10e5;
  NL.Up     =  0.0;      NL.Low     = 10e5;
  NR.Up     =  0.0;      NR.Low     = 10e5;
  Nges.Up   =  0.0;      Nges.Low   = 10e5;
  Wjets.Up  =  0.0;      Wjets.Low  = 10e5;
  QCD.Up    =  0.0;      QCD.Low    = 10e5;
  RemBkg.Up =  0.0;      RemBkg.Low = 10e5;

  for (int i=0; i != num_of_nui_parameters; i++)
  {
	nui[i].Up = -0.001; nui[i].Low = 0.001;
  }

  F0.Up_unc     =  0.0;      F0.Low_unc     = 10.0;
  FL.Up_unc     =  0.0;      FL.Low_unc     = 10.0;
  FR.Up_unc     =  0.0;      FR.Low_unc     = 10.0;
  N0.Up_unc     =  0.0;      N0.Low_unc     = 10e5;
  NL.Up_unc     =  0.0;      NL.Low_unc     = 10e5;
  NR.Up_unc     =  0.0;      NR.Low_unc     = 10e5;
  Nges.Up_unc   =  0.0;      Nges.Low_unc   = 10e5;
  Wjets.Up_unc  =  0.0;      Wjets.Low_unc  = 10e5;
  QCD.Up_unc    =  0.0;      QCD.Low_unc    = 10e5;
  RemBkg.Up_unc =  0.0;      RemBkg.Low_unc = 10e5;

  for (int i=0; i != num_of_nui_parameters; i++)
  {
	nui[i].Up_unc = -0.001; nui[i].Low_unc = 0.001;
  }

  F0.Up_prof_unc     =  0.0;      F0.Low_prof_unc     = 10.0;
  FL.Up_prof_unc     =  0.0;      FL.Low_prof_unc     = 10.0;
  FR.Up_prof_unc     =  0.0;      FR.Low_prof_unc     = 10.0;
  N0.Up_prof_unc     =  0.0;      N0.Low_prof_unc     = 10e5;
  NL.Up_prof_unc     =  0.0;      NL.Low_prof_unc     = 10e5;
  NR.Up_prof_unc     =  0.0;      NR.Low_prof_unc     = 10e5;
  Nges.Up_prof_unc   =  0.0;      Nges.Low_prof_unc   = 10e5;
  Wjets.Up_prof_unc  =  0.0;      Wjets.Low_prof_unc  = 10e5;
  QCD.Up_prof_unc    =  0.0;      QCD.Low_prof_unc    = 10e5;
  RemBkg.Up_prof_unc =  0.0;      RemBkg.Low_prof_unc = 10e5;

  for (int i=0; i != num_of_nui_parameters; i++)
  {
	nui[i].Up_prof_unc = -0.001; nui[i].Low_prof_unc = 0.001;
  }


  std::cout << "Make parameter list" << std::endl;

  std::vector<histograms::Parameter> ParameterList;

  ParameterList.push_back(F0);
  ParameterList.push_back(FL);
  ParameterList.push_back(FR);
  ParameterList.push_back(N0);
  ParameterList.push_back(NL);
  ParameterList.push_back(NR);
  ParameterList.push_back(Wjets);
  ParameterList.push_back(QCD);
  ParameterList.push_back(RemBkg);
  ParameterList.push_back(Nges);
  
  for (int i=0; i != num_of_nui_parameters; i++)
  {
        //cout << "here " << i << endl;
	ParameterList.push_back(nui[i]);
  }

  cout << "Set parameter list " << endl;

  fNewHisto -> SetParameterList(ParameterList);

  cout << "Fill histograms " << endl;

  fNewHisto -> FillHistograms(input_channel, num_of_nui_parameters);
  fNewHisto -> MakeHTML(htmltitle);
  fNewHisto -> MakeOutputTrees(num_of_nui_parameters);

  delete fNewHisto;

  return 0 ;


}
