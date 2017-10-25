#include "CalibrationCurves.h"
#include "StatusLogbook.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

int main(int argc, char *argv[])
{


  if(argc != 7){

    WriteErrorStatus("ProduceCalibration", "Wrong number of input parameters");
    
  }

  std::string input_folder  = argv[1];
  std::string output_folder = argv[2];
  std::string channel       = argv[3];
  std::string type          = argv[4];
  std::string fNumNui       = argv[5];
  std::string fNuiVar       = argv[6];

  CalibrationCurves *fCalib = new CalibrationCurves();

  if(type == "CaliNui"){

    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_m1.5_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_m1_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_m0.5_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0.5_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p1_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p1.5_Plots/");

    fCalib -> SetXValueList(-1.5);
    fCalib -> SetXValueList(-1.0);
    fCalib -> SetXValueList(-0.5);
    fCalib -> SetXValueList(0.0);
    fCalib -> SetXValueList(0.5);
    fCalib -> SetXValueList(1.0);
    fCalib -> SetXValueList(1.5);
    
    fCalib -> SetOutputFolder(output_folder);
    
    // true = use XValueList as XValues
    fCalib -> FillGraphs("nui1_tree",   "Nui1",   channel, true, false);
    fCalib -> FillGraphs("F0_tree",     "F0",     channel, true, false);
    fCalib -> FillGraphs("FL_tree",     "FL",     channel, true, false);
    fCalib -> FillGraphs("FR_tree",     "FR",     channel, true, false);
    fCalib -> FillGraphs("N0_tree",     "N0",     channel, true, false);
    fCalib -> FillGraphs("NL_tree",     "NL",     channel, true, false);
    fCalib -> FillGraphs("NR_tree",     "NR",     channel, true, false);
    fCalib -> FillGraphs("Wjets_tree",  "Wjets",  channel, true, false);
    fCalib -> FillGraphs("QCD_tree",    "QCD",    channel, true, false);
    fCalib -> FillGraphs("RemBkg_tree", "RemBkg", channel, true, false);

    int num_nui = TString(argv[5]).Atoi();

    for(int i = 1; i < num_nui+1; ++i){

      std::stringstream oss;
      oss << i;

      fCalib -> FillGraphs("nui"+oss.str()+"_tree",   "Nui"+oss.str(),   channel, true, false);

    }
    
  }
  else if(type == "CaliFraction"){
    
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.4_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.5_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.6_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.7_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.8_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_0.9_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    fCalib -> SetInputFolderList(input_folder+"/Output_"+channel+"_F0_1_NrNui"+fNumNui+"_NuiVar"+fNuiVar+"_p0_Plots/");
    
    fCalib -> SetXValueList(0.4);
    fCalib -> SetXValueList(0.5);
    fCalib -> SetXValueList(0.6);
    fCalib -> SetXValueList(0.7);
    fCalib -> SetXValueList(0.8);
    fCalib -> SetXValueList(0.9);
    fCalib -> SetXValueList(1.0);
    
    fCalib -> SetOutputFolder(output_folder);
    
    // true = use XValueList as XValues
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for F0");
    fCalib -> FillGraphs("F0_tree",     "F0",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for FL");
    fCalib -> FillGraphs("FL_tree",     "FL",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for FR");
    fCalib -> FillGraphs("FR_tree",     "FR",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for N0");
    fCalib -> FillGraphs("N0_tree",     "N0",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for NL");
    fCalib -> FillGraphs("NL_tree",     "NL",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for NR");
    fCalib -> FillGraphs("NR_tree",     "NR",     channel, false, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for Wjets");
    fCalib -> FillGraphs("Wjets_tree",  "Wjets",  channel, true, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for QCD");
    fCalib -> FillGraphs("QCD_tree",    "QCD",    channel, true, true);
    WriteInfoStatus("ProduceCalibration", "Run Calicurve for RemBkg");
    fCalib -> FillGraphs("RemBkg_tree", "RemBkg", channel, true, true);
    
    int num_nui = TString(argv[5]).Atoi();
    
    for(int i = 1; i < num_nui+1; ++i){
      
      std::stringstream oss;
      oss << i;

      fCalib -> FillGraphs("nui"+oss.str()+"_tree",   "Nui"+oss.str(),   channel, true, true);
      
    }

    
    
  }
  else
    WriteErrorStatus("ProduceCalibration", "This type of calicurve does not exist!!!");
  
  return(0);

}
