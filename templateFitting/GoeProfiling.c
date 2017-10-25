#include "StatusLogbook.h"
#include "ProfilingClass.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
 
     if(argc != 8){

       WriteErrorStatus("GoeProfiling", "Number of input variables is wrong!!!");
       WriteErrorStatus("GoeProfiling", "EXIT   ");

       return 1;

     }
     else{
       
       WriteParameterStatus("GoeProfiling", "Input Channel: "               + string(argv[1]));
       WriteParameterStatus("GoeProfiling", "Fit method: "                  + string(argv[2]));
       WriteParameterStatus("GoeProfiling", "Validation mode: "             + string(argv[3]));
       WriteParameterStatus("GoeProfiling", "Number of PE: "                + string(argv[5]));
       WriteParameterStatus("GoeProfiling", "Output sample number: "        + string(argv[7]));
       WriteParameterStatus("GoeProfiling", "Nuisance parm. for variation: "+ string(argv[6]));
       WriteParameterStatus("GoeProfiling", "Background treatment: "        + string(argv[4]));
       //   WriteParameterStatus("GoeProfiling", "Number of histogram bins: "    + string(argv[7]));
     }
     
     //Obtain input arguments
     std::string InputChannel   = argv[1];
     std::string FitMethod      = argv[2];
     std::string ValidationMode = argv[3];
     std::string BackgroundMode = argv[4];
     int NumberOfPE             = TString(argv[5]).Atoi();
     int OutputSampleNumber     = TString(argv[7]).Atoi(); 
     int NuisanceParamVar       = TString(argv[6]).Atoi();

     //     int NumberOfBins  =   15;
     double lower_edge = -1.0;
     double upper_edge =  1.0;

     /*    //Globally defined histograms, needed for fit function
     hist      = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_sum  = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_help = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);*/


     double x[8]       = {-1.,-0.5,-0.25,0.0,0.25,0.5,0.75,1.};
     
     int    n    = 7;
     
     int NumberOfBins = 7;

     //Globally defined histograms, needed for fit function                                                                                                            
     hist      = TH1D("hist", "hist", n, x);
     hist_sum  = TH1D("hist", "hist", n, x);
     hist_help = TH1D("hist", "hist", n, x);


     //Define ProfilingClass object
     ProfilingClass fProf = ProfilingClass();

     //Add signal parameter with param name and hist name
     fProf.AddSignalParameter("N0",  "CosTheta_F0", "F0eff");
     fProf.AddSignalParameter("NL",  "CosTheta_FL", "FLeff");
     fProf.AddSignalParameter("NR",  "CosTheta_FR", "FReff");

     //Add bkg parameter with param name and hist name
     fProf.AddBackgroundParameter("Wjets",  "Wjets",  "WjetsNorm",  "WjetsUnc");
     fProf.AddBackgroundParameter("QCD",    "QCD",    "QCDNorm",    "QCDUnc");
     fProf.AddBackgroundParameter("RemBkg", "RemBkg", "RemBkgNorm", "RemBkgUnc");

     //Two different modes characterised by different treatment of bkg parameters
     if(BackgroundMode == "BkgFit"){ //bkg fitted via scale/normalisation parameters

       WriteParameterStatus("GoeProfiling", "Run profiling with fitted Bkg");

       fProf.SetBackgroundMode("BkgFit");

     }
     else if(BackgroundMode == "BkgNui"){//bkg fitted via nuisance parameters

       WriteParameterStatus("GoeProfiling", "Run profiling with Bkg as nuisance parameter");
       
       fProf.SetBackgroundMode("BkgNui");

       fProf.AddNuisanceParameter("WjetsNorm");
       fProf.AddNuisanceParameter("QCDNorm");
       fProf.AddNuisanceParameter("RemBkgNorm");

     }
     else{
       
       WriteErrorStatus("GoeProfiling", "No such Bkg Mode, possible values are BkgFit and BkgNui!!!");
       WriteErrorStatus("GoeProfiling", "EXIT   ");
       
       return 1;
       
     }

     //Needed to specify name of data histogram
     fProf.SetDataHist("Data");
     
     //     fProf.SetOutputFolder("Root_Files_Cora_neu");

     //Define different input templates
     // fProf.SetNominalInputFile("Root_Files_Cora_neu/TotalTemplatesFile_4inclJets_1inclTags_"+InputChannel+"_TS.root");

     fProf.SetNominalInputFile("Templates_Combined_"+InputChannel+".root");                                              



     //Add nuisance parameters to the fit

     //fProf.AddNuisanceParameter("BTAG");
     //fProf.AddNuisanceParameter("CTAG");
     //fProf.AddNuisanceParameter("MISTAG");
     //fProf.AddNuisanceParameter("JVF_SF");
     //fProf.AddNuisanceParameter("MUON_ID");
     //fProf.AddNuisanceParameter("MUON_RECO");
     //fProf.AddNuisanceParameter("MUON_TRIG");
     //fProf.AddNuisanceParameter("JES");
     //fProf.AddNuisanceParameter("JER");
     //fProf.AddNuisanceParameter("JEFF");
     //fProf.AddNuisanceParameter("MUSC");
     //fProf.AddNuisanceParameter("FLAVOR_COMB");

     //Specify the available up/down variations (2k, 3k, 4k)
     for(int iK = 1; iK < 2; ++iK) { //change number depending on input
       //fProf.AddSystematicFile("BTAG",      "root_files/Templates_" + InputChannel, "btagu",     "btagd",     iK, -iK);
       //fProf.AddSystematicFile("CTAG",      "root_files/Templates_" + InputChannel, "ctagu",     "ctagd",     iK, -iK);
       //fProf.AddSystematicFile("MISTAG",    "root_files/Templates_" + InputChannel, "mistagu",   "mistagd",   iK, -iK);
       //fProf.AddSystematicFile("JVF_SF", "root_files/Templates_" + InputChannel, "jvfsfu",   "jvfsfd",   iK, -iK);
       //fProf.AddSystematicFile("MUON_ID",   "root_files/Templates_" + InputChannel, "muonidu",   "muonidd",   iK, -iK);
       //fProf.AddSystematicFile("MUON_RECO", "root_files/Templates_" + InputChannel, "muonrecou", "muonrecod", iK, -iK);
       //fProf.AddSystematicFile("MUON_TRIG", "root_files/Templates_" + InputChannel, "muontrigu", "muontrigd", iK, -iK);
       //fProf.AddSystematicFile("FLAVOR_COMB", "root_files/Templates_" + InputChannel, "flavorcompu", "flavorcompd", iK, -iK);
	//fProf.AddSystematicFile("JES",   "root_files/Templates_" + InputChannel, "jesu",   "jesd",   iK, -iK);
       }

     for(int iK = 1; iK < 3; ++iK) { //change number depending on input
       //fProf.AddSystematicFile("JES",   "root_files/Templates_" + InputChannel, "jesu",   "jesd",   iK, -iK);

     }
     
     for(int iK = 1; iK < 4; ++iK) { //change number depending on input
       //fProf.AddSystematicFile("JES",   "root_files_05122012_4incl/Templates_" + InputChannel, "jesu",   "jesd",   iK, -iK);

     }

     //Parameters for which only up/down variation is available
     //last string indicates which kind of variation is in input!!!:
     //"up" or "down" 
     for(int iK = 1; iK < 2; ++iK) { //change number depending on input
       //fProf.SymmetriseSystematicFile("JER",  "root_files/Templates_" + InputChannel, "jeru", "jerd", iK, -iK, "up");
       //fProf.SymmetriseSystematicFile("JEFF", "root_files/Templates_" + InputChannel, "jeffu", "jeffd", iK, -iK, "up");
       //fProf.SymmetriseSystematicFile("MUSC", "root_files/Templates_" + InputChannel, "muscu", "muscd", iK, -iK, "up");
     }
     
     // lumi in pb^-1
     //fProf.SetLumi(4655.74); 
     fProf.SetLumi(20276.9); //mjk : //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopData12
 
     // xsec in pb
     //fProf.SetXsec(96.0);
     fProf.SetXsec(114.48); //mjk :: for signal 110404
     

     //Set other important variables, all in ProfilingClass.h
     fProf.SetInputChannel(InputChannel);
     fProf.SetFitMethod(FitMethod);
     fProf.SetValidationMode(ValidationMode);
     fProf.SetNumberOfPE(NumberOfPE);
     fProf.SetOutputSampleNumber(OutputSampleNumber);
     fProf.SetNuiVariation(NuisanceParamVar);
     fProf.SetNumberOfBins(NumberOfBins);
     fProf.SetLowerEdge(lower_edge);
     fProf.SetUpperEdge(upper_edge);

     cout << "fNBins " << fNBins << endl;

     //Set Interpolation method:
     // options are: QuadraticFit,QuadraticInterp,LinearInterp,PiecewiseLinear
     std::string InterpolationMethod = "QuadraticFit";
     //std::string InterpolationMethod = "LinearInterp";
     //std::string InterpolationMethod = "PiecewiseLinear";
     //std::string InterpolationMethod = "QuadraticInterp";

     fProf.SetInterpolationMethod(InterpolationMethod);

     fProf.ReadNominalInfo();

     WriteInfoStatus("GoeProfiling", "Read systematic info");
     
     fProf.ReadSystematicInfo();

     //Curves needed for interpolation step 
     WriteInfoStatus("GoeProfiling", "Make interpolation curves");
     WriteInfoStatus("GoeProfiling", "That can take a few seconds...");
       
     fProf.MakeInterpolationCurves();

     WriteInfoStatus("GoeProfiling", "Fill interpolation objects");
     
     fProf.FillInterpolationObjects();

     WriteInfoStatus("GoeProfiling", "Do Validation");

     ///Validation step; perform fit to data, pseudo-data or validate fitting method 
     fProf.DoValidation();

     WriteInfoStatus("GoeProfiling", "finalize");

     WriteInfoStatus("GoeProfiling", "finished");

    return 0;

}
