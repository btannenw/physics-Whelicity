#include "ProfilingClass.h"
#include "PlotInterpolationCurves.h"
#include "ExternalSystematics.h"
#include "ValidationClass.h"
#include "functions.h"
#include "StatusLogbook.h"

#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

std::string fInterpolationMethod;
std::vector<std::vector<std::vector<TF1> > > fFitFunc;
std::vector<std::vector<std::vector<TF1> > > fFitFunc_QuadFit;

int fNBins;
double fLowerEdge;
double fUpperEdge;

TH1D hist;   
TH1D hist_sum; 
TH1D hist_help;
TH1D hist_data;

ProfilingClass::ProfilingClass()
{

  fSystematics.clear();

}

ProfilingClass::~ProfilingClass()
{

}

void ProfilingClass::ProducePseudoData(double F0, double FL, double FR, std::vector<std::string> HistoNames, std::vector<std::string> VarNames, std::string FileName, std::string OutputFileName)
{
  
  TFile *fFile = new TFile(FileName.c_str(), "READ");

  TH1D histF0 = *(TH1D*) fFile -> Get(HistoNames[0].c_str());
  TH1D histFL = *(TH1D*) fFile -> Get(HistoNames[1].c_str());
  TH1D histFR = *(TH1D*) fFile -> Get(HistoNames[2].c_str());

  std::vector<TH1D> BkgHisto;
  
  int nFrac = fSignalParameter.size(); // = 3 for 3D  &  =2 for 2D

  for(int i = nFrac; i < HistoNames.size(); ++i){ // HistoNames contains 3 signal & 3 bkg histos 

    TH1D help = *(TH1D*) fFile -> Get(HistoNames[i].c_str());

    BkgHisto.push_back(help);

  }

  TTree *tree = (TTree*) fFile -> Get("InfoTree");

  double F0eff, FLeff, FReff; // float -> double (by MJK)

  std::vector<float> BackgroundNorm(HistoNames.size()-nFrac);

  tree->SetBranchAddress(VarNames[0].c_str(), &F0eff);
  tree->SetBranchAddress(VarNames[1].c_str(), &FLeff);
  if(fTemplateFittingMode == "3D")
    tree->SetBranchAddress(VarNames[2].c_str(), &FReff);

  for(int i = nFrac; i < HistoNames.size(); ++i){

    tree->SetBranchAddress(VarNames[i].c_str(), &(BackgroundNorm[i-nFrac]));
    
  }

  // Get the Event
  tree->GetEntry(0);
  //double fKfactor=1.1994; // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopMC12DiTopSamples#Powheg_Pythia_FS
  double fKfactor=1.1995; //for 110404
  //fKfactor=1.1994; // for 117050
  
  double NormF0 = fXsec*fLumi*F0*F0eff*fKfactor;
  double NormFL = fXsec*fLumi*FL*FLeff*fKfactor;
  double NormFR = fXsec*fLumi*FR*FReff*fKfactor;

  histF0.Scale(NormF0/histF0.Integral());
  histFL.Scale(NormFL/histFL.Integral());
  if(fTemplateFittingMode == "3D")
    histFR.Scale(NormFR/histFR.Integral());

  //  std::cout << histF0.Integral() << "\t" << histFL.Integral() << "\t" << histFR.Integral() << std::endl;
  
  for(int i = 0; i < BkgHisto.size(); ++i){
    
    BkgHisto[i].Scale(BackgroundNorm[i]/BkgHisto[i].Integral());
    
    //    std::cout << "Bkg Histo: " << BkgHisto[i].Integral() << std::endl;

  }

  // now add all Histos
  TH1D TotalHist = histF0;
 
  if(fTemplateFittingMode == "3D")
    TotalHist = TotalHist + histFL + histFR;
  else if(fTemplateFittingMode == "2D")
    TotalHist = TotalHist + histFL;

  std::cout<<"TotalHist.Integral (signal ONLY)= "<<TotalHist.Integral()<<std::endl;

  for(int i = 0; i < BkgHisto.size(); ++i){
    
    TotalHist = TotalHist + BkgHisto[i];

  }

  std::cout<<"TotalHist.Integral (ALL)= "<<TotalHist.Integral()<<std::endl;

  TFile *fFile2 = new TFile(OutputFileName.c_str(), "RECREATE");
  
  TotalHist.Write("PseudoData");

  fFile2 -> Close();

}

//Check whether a certain systematic effect already included 
bool ProfilingClass::CheckNewSystematic(std::string systname)
{

  bool isNew = true;

  for(int iSyst = 0; iSyst < fSystematics.size(); ++iSyst){

    if(systname == fSystematics[iSyst])
      isNew = false;

  }

  return isNew;

}

//Get index of a certain systematic effect
int ProfilingClass::GetSystIndex(std::string systname)
{

  int index = -1;

  for(int iSyst = 0; iSyst < fSystematics.size(); ++iSyst){

    if(systname == fSystematics[iSyst])
      index = iSyst;

  }

  return index;

}

//Function related to external evaluation
void ProfilingClass::AddSystematicPseudoData(std::string Syst, std::string input_file, std::string HistoName)
{
  
  TFile *fFile = new TFile(input_file.c_str(), "READ");

  SystematicInfo::MySystematicPD HelpType;

  HelpType.SystematicType = Syst;
  //std::cout << "HistoName.c_str(): "<< HistoName.c_str() << std::endl;
  HelpType.histoPD = *(TH1D*) fFile -> Get(HistoName.c_str());

  std::cout << "Add systematic:   " << Syst.c_str() << "  with file   " << input_file.c_str() << "   and Name  " << HistoName << "  and Integral:  " << HelpType.histoPD.Integral() << std::endl;

  fSystematicPD.push_back(HelpType);

}


//___Following functions ordered according to GoeProfiling.c___

//Add up- and down-variations of the systematic effects used in the fit and fill fSystematicFiles vector
void ProfilingClass::AddSystematicFile(std::string Syst, std::string Folder, std::string systu, std::string systd, int ku, int kd)
{

  // check if Systematic is new or already exists
  if(CheckNewSystematic(Syst)){
    // if new, add it to Systematic vector
    fSystematics.push_back(Syst);
    
    // and create new vector that should contain the input file path and the k-Value
    //Fill fSystematicsFiles vector
    std::vector<SystematicInfo::MySystematic> HelpVec;
    fSystematicFiles.push_back(HelpVec);
  }

  int SystIndex = GetSystIndex(Syst);

  SystematicInfo::MySystematic HelpType;
  
  std::stringstream kval_u, kval_d;
  kval_u << ku;
  kval_d << kd;

  WriteInfoStatus("ProfilingClass", "Add Systematic with type = "+Syst);

  HelpType.SystematicType = Syst;
  HelpType.Filename_up    = Folder+"_"+systu + "_" + kval_u.str() + "k.root";
  HelpType.Filename_down  = Folder+"_"+systd + "_" + kval_u.str() + "k.root"; // has to be kval_u !!!
  HelpType.kup            = ku;
  HelpType.kdown          = kd;

  fSystematicFiles[SystIndex].push_back(HelpType);

}


//Read signal and background parameters
void ProfilingClass::ReadNominalInfo()
{

  //  SetNumberOfSignalHist(fSignalParameter.size());

  //  SetNumberOfBkgHist(fBkgParameter.size());

  fNSignal = fSignalParameter.size(); // TIP: fNSignal declared as extern variable in functions.h and defined in functions.cxx, so it can be used anywhere by including the header file. 
  fNBkg    = fBkgParameter.size();

  std::stringstream oss1, oss2;
  oss1 << fNSignal;
  oss2 << fNBkg;

  WriteParameterStatus("ProfilingClass", "Number of Bkg histograms = "+oss2.str());
  WriteParameterStatus("ProfilingClass", "Number of Signal histograms = "+oss1.str());

  TFile* file = new TFile(fNomInputFile.c_str(), "open");

  // read out signal information...
  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){

    std::string histo     = fSignalParameter[iSig].HistName;
    std::string param_eff = fSignalParameter[iSig].EffName;
    std::string histoTS   = fSignalParameter[iSig].HistName+"_TS";

    fSignalParameter[iSig].hist   = *(TH1D*) file -> Get(histo.c_str());
    //    fSignalParameter[iSig].histTS = *(TH1D*) file -> Get(histoTS.c_str());

    TTree *tree = (TTree*) file -> Get("InfoTree");

    double eff; // float -> double (by MJK)

    // Link to Branch in Tree
    tree->SetBranchAddress(param_eff.c_str(), &eff);
    //std::cout<<"!!!! Branch: "<<param_eff.c_str()<<"\t Value: "<<eff<<std::endl;
    // Get the Event
    tree->GetEntry(0);

    fSignalParameter[iSig].SelEff = eff;

    std::cout << "Parameter " << iSig << " has efficiency " << eff << std::endl;

    delete tree;
  }

  // read out bkg information...
  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){

    std::string histo      = fBkgParameter[iBkg].HistName;
    std::string histoTS    = fBkgParameter[iBkg].HistName+"_TS";
    std::string param_exp  = fBkgParameter[iBkg].ExpName;
    std::string param_unc  = fBkgParameter[iBkg].UncName;

    std::cout << "Bkg param " << iBkg << "\t" << histo.c_str() << std::endl;

    fBkgParameter[iBkg].hist   = *(TH1D*) file -> Get(histo.c_str());
    //    fBkgParameter[iBkg].histTS = *(TH1D*) file -> Get(histoTS.c_str());

    std::cout << "Bkg param " << iBkg << "   Bins   " << fBkgParameter[iBkg].hist.GetNbinsX() << "\t" << fBkgParameter[iBkg].hist.GetBinLowEdge(1) << std::endl;

    TTree *tree = (TTree*) file -> Get("InfoTree");

    float exp, unc;

    // Link to Branch in Tree
    tree->SetBranchAddress(param_exp.c_str(), &exp);
    //tree->SetBranchAddress(param_unc.c_str(), &unc);

    // Get the Event
    tree->GetEntry(0);


    fBkgParameter[iBkg].Exp = (double)exp;

    /*    if(iBkg == 0)
      fBkgParameter[iBkg].Unc = (double)unc*3.0;
      else */
    
    //fBkgParameter[iBkg].Unc = (double)unc*1.0;

    // BT, Mar 23: Add hard-coded normalization uncertainties for backgrounds
    double normUnc = 0.0;
    if (param_unc.find("Fixed")!=string::npos) // Fixed Bkg
      normUnc = 0.005;
    else if (param_unc.find("WjetsUnc")!=string::npos) // W+jets (1W)
      normUnc = 0.48;
    else if (param_unc.find("QCDUnc")!=string::npos){ // QCD
      std::cout<<"\n\n###############\n###############\n bkg parameter uncertainty "<<param_unc<<" Set to 0.005 [FIXED] \n###############\n###############\n\n"<<std::endl;
      //normUnc = 0.30;
      normUnc = 0.005;
    }
    else if (param_unc.find("WLightUnc")!=string::npos){ // W+Light (3W)
      //normUnc = 0.05; // CF
      // normUnc = 0.10; // 2CF
      //normUnc = 0.48;
       std::cout<<"\n\n###############\n###############\n bkg parameter uncertainty "<<param_unc<<" Set to 0.005 [FIXED] \n###############\n###############\n\n"<<std::endl;
       normUnc = 0.005;
    }
    else if (param_unc.find("WcUnc")!=string::npos){ // W+c (3W)
      //normUnc = 0.25; // CF
      //normUnc = 0.50; // 2CF
      //normUnc = 0.48;
      std::cout<<"\n\n###############\n###############\n bkg parameter uncertainty "<<param_unc<<" Set to 0.005 [FIXED] \n###############\n###############\n\n"<<std::endl;
      normUnc = 0.005;
    }
    else if (param_unc.find("WbbccUnc")!=string::npos){ // W+bb/cc (3W)
      //normUnc = 0.07; // CF
      // normUnc = 0.14; // 2CF
      //normUnc = 0.48;              
      std::cout<<"\n\n###############\n###############\n bkg parameter uncertainty "<<param_unc<<" Set to 0.005 [FIXED] \n###############\n###############\n\n"<<std::endl;
      normUnc = 0.005;
    }
    else if (param_unc.find("RemBkgUnc")!=string::npos){ // RemBkg
      if (fInputChannel == "el_mu_lephad_bTag") // 8 - channel combination
	     normUnc = 0.125; //single top unc only cross-section
	     //normUnc = 0.41696*0.8061 + 0.48*0.176 + 0.48*0.0174; //single top unc x-sec + 24%/jet in quadrature
      if (fInputChannel == "el_mu_lephad" || fInputChannel == "el_mu_bTag") // 4 - channel combination
	     normUnc = 0.093; //single top unc only cross-section
	     //normUnc = 0.41696*0.8846 + 0.48*0.10 + 0.48*0.0564; //single top unc x-sec + 24%/jet in quadrature
    }
    else{ // Name not found.... BIG PROBLEM
      std::cout<<"\n\n###############\n###############\n bkg parameter uncertainty "<<param_unc<<" NOT FOUND! Setting Default to 0.48\n###############\n###############\n\n"<<std::endl;
      normUnc = 0.48; 
      exit(1);
    }
    std::cout<<"\n############### "<<param_unc<<" uncertainty "<<normUnc<<" ###############"<<std::endl;

    fBkgParameter[iBkg].Unc = fBkgParameter[iBkg].Exp * normUnc;
    
    
    //   std::cout << "Parameter " << iBkg << " has norm " << fBkgParameter[iBkg].Exp << " +- " << fBkgParameter[iBkg].Unc << std::endl;
    
    delete tree;
  }
  
  std::cout << fDataHist.c_str() << std::endl;
  
  //datahist
  hist_data = *(TH1D*) file -> Get(fDataHist.c_str());

  std::cout << hist_data.GetNbinsX() << std::endl;

  WriteInfoStatus("ProfilingClass", "Nominal Input has been read out");

}


//Read different systematic parameters
void ProfilingClass::ReadSystematicInfo() // mjk: used in GoeProfiling.c
{

  fNSyst = fSystematicFiles.size();//filled in AddSystematicFile

  if(fBkgMode == "BkgFit")
    fNNui  = fSystematicFiles.size();
  else
    fNNui  = fSystematicFiles.size()+fNBkg;

  for(int iSyst = 0; iSyst < fNSyst; ++iSyst){

    WriteInfoStatus("ProfilingClass", "Read Information for Systematic "+fSystematics[iSyst]);

    int NKfiles = fSystematicFiles[iSyst].size();

    for(int iK = 0; iK < NKfiles; ++iK){

      std::string file_up   = fSystematicFiles[iSyst][iK].Filename_up;
      std::string file_down = fSystematicFiles[iSyst][iK].Filename_down;

      TFile* file_u = new TFile(file_up.c_str(), "open");
      
      for(int iSig = 0; iSig < fNSignal; ++iSig){
	
	std::string histo     = fSignalParameter[iSig].HistName;
	
	fSystematicFiles[iSyst][iK].HistUp.push_back(*(TH1D*) file_u -> Get(histo.c_str()));
	fSystematicFiles[iSyst][iK].ParamName.push_back(fSignalParameter[iSig].ParamName);
	
      }
      for(int iBkg = 0; iBkg < fNBkg; ++iBkg){

	std::string histo     = fBkgParameter[iBkg].HistName;

        fSystematicFiles[iSyst][iK].HistUp.push_back(*(TH1D*) file_u -> Get(histo.c_str()));
        fSystematicFiles[iSyst][iK].ParamName.push_back(fBkgParameter[iBkg].ParamName);

      }

      TFile* file_d = new TFile(file_down.c_str(), "open");

      for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){

	std::string histo     = fSignalParameter[iSig].HistName;

        fSystematicFiles[iSyst][iK].HistDown.push_back(*(TH1D*) file_d -> Get(histo.c_str()));

      }
      for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){

	std::string histo     = fBkgParameter[iBkg].HistName;

        fSystematicFiles[iSyst][iK].HistDown.push_back(*(TH1D*) file_d -> Get(histo.c_str()));

      }
      //    delete file_u;
      //    delete file_d;
    }
    
  }

}


//Frame for creating the interpolation curves, see PlotInterpolationCurves.cxx
void ProfilingClass::MakeInterpolationCurves()
{
 
    //Only used if background fitted via nuisance parameters 
    if(fBkgMode == "BkgNui"){

      for(int iBkgNui = 0; iBkgNui < fBkgParameter.size(); ++iBkgNui){
	
	//Define PlotInterpolationCurves class object, see corresponding file
        PlotInterpolationCurves *fPlot = new PlotInterpolationCurves(fSignalParameter, fBkgParameter, iBkgNui, fOutputFolder+"/Templates_"+fInputChannel);
	
        fPlot -> FillVectorsBkgNorm(iBkgNui);
	
        fPlot -> MakePlots(fInterpolationMethod);
	
        fFitFunc.push_back(fPlot -> GetInterpolationFunctions());//globally defined
	
        delete fPlot;
	
      }
      
    }
    
    for(int iSyst = 0; iSyst < fSystematicFiles.size(); ++iSyst){

      //Define PlotInterpolationCurves class object, see corresponding file
      PlotInterpolationCurves *fPlot = new PlotInterpolationCurves(fSignalParameter, fBkgParameter, fSystematicFiles[iSyst]);

      fPlot -> FillVectorsNominal();

      fPlot -> FillVectorsSystematics();
     //cout << "fNBins " << fNBins << endl;
      fPlot -> MakePlots(fInterpolationMethod);

      //Returns fFitFunc vector as defined in plotInterpolationCurves.cxx
      //fFitFunc as vector of vectors globally defined here; 
      fFitFunc.push_back(fPlot -> GetInterpolationFunctions());
      fFitFunc_QuadFit.push_back(fPlot -> GetInterpolationFunctions_QuadFit());
      
      delete fPlot;
      
    }

  

}


//Defined in functions.cxx, according to TemplateInfo
void ProfilingClass::FillInterpolationObjects()
{
 
  for(int iSyst = 0; iSyst < fSystematicFiles.size(); ++iSyst){
    
    AddSignalSystematic(fSignalParameter, fSystematicFiles[iSyst]);
    AddBkgSystematic(fBkgParameter, fSystematicFiles[iSyst]);
    
  }
  
}


//External evaluation, described in EvaluateExternal.c
void ProfilingClass::DoExternalSystematicEvaluation(std::string SystType)
{

  ExternalSystematics *fExt = new ExternalSystematics(fTemplateFittingMode);

  WriteInfoStatus("ProfilingClass-DoExtSystEval", "Set Background parameters");

  fExt -> SetBkgParameters(fBkgParameter);

  WriteInfoStatus("ProfilingClass-DoExtSystEval", "Set Signal parameters");

  fExt -> SetSignalParameters(fSignalParameter);

  WriteInfoStatus("ProfilingClass-DoExtSystEval", "Set Systematic PseudoData");
  
  fExt -> SetSystematicPseudoData(fSystematicPD);

  fExt -> SetSystematicLabel(SystType);
  fExt -> SetOutputFolder(fOutputFolder);

  fExt -> SetLumi(fLumi);
  fExt -> SetXsec(fXsec);
  fExt -> SetKfactor(fKfactor); // added by MJK
  fExt -> SetNumberOfPE(fNPseudoExp);
  fExt -> SetChannel(fInputChannel);

  fExt -> SetFractions(fFraction);
 
  //std::cout << fLeptonLabel.size() << fJetBinLabel.size() << fBTagLabel.size() << std::endl;

  fExt -> SetSubPlotLabels(fLeptonLabel, fJetBinLabel, fBTagLabel);

  fExt -> SetFitParameters("");
  fExt -> FillTemplateVector(fSignalParameter, fBkgParameter);
  
  WriteInfoStatus("ProfilingClass", "Set fit parameters");

  fExt -> SetEvaluationMode(fEvaluationMode);

  std::cout << "SystematicsType: " << SystType.c_str() << std::endl;

  if(SystType == "Calibration")
    fExt -> CallEvaluation(fOutputTxtFile, "Cali");
  else if(SystType == "Datafit")
    fExt -> CallEvaluation(fOutputTxtFile, "Datafit");
  else if( SystType.find("CT10")!=string::npos || SystType.find("MSTW")!=string::npos || SystType.find("NNPDF")!=string::npos ) //
    fExt -> CallEvaluation(fOutputTxtFile, "PDF");
  else if( SystType.find("TopMass")!=string::npos) //
    fExt -> CallEvaluation(fOutputTxtFile, "mass");// for all other systs
  else
    fExt -> CallEvaluation(fOutputTxtFile, "Syst");// for all other systs

  WriteInfoStatus("ProfilingClass", "finalize");

  //delete fExt;

  WriteInfoStatus("ProfilingClass", "finished");


}

//Validation step; perform fit to data, pseudo-data or validate fitting method, described in GoeProfiling.c
void ProfilingClass::DoValidation()
{

  //Define validation object
  tvalidation *val_obj = new tvalidation();

  WriteInfoStatus("ProfilingClass-DoValidation", "Set Background parameters");

  val_obj -> SetBkgParameters(fBkgParameter);

  WriteInfoStatus("ProfilingClass-DoValidation", "Set Signal parameters");

  val_obj -> SetSignalParameters(fSignalParameter);

  WriteInfoStatus("ProfilingClass-DoValidation", "Set Nuisance parameters");

  val_obj -> SetNuisanceParameters(fNuisanceParameter);

  WriteInfoStatus("ProfilingClass-DoValidation", "Set output sample number");

  val_obj -> SetOutputSampleNumber(fNOutputSample);

  val_obj -> SetLumi(fLumi);
  val_obj -> SetXsec(fXsec);
  val_obj -> SetNumberOfPE(fNPseudoExp);
  
  WriteInfoStatus("ProfilingClass", "Set fit parameters");

  val_obj -> SetFitParameters(fBkgMode);

  WriteInfoStatus("ProfilingClass", "Fill histogram vector");

  if(fNuisanceParameter.size() > 0) {
    val_obj -> FillInterpolationVector();}
  else {
    val_obj -> FillTemplateVector(fSignalParameter, fBkgParameter); }


  WriteInfoStatus("ProfilingClass", "Validation mode "+fValidationMode);

  val_obj -> CallValidation(fValidationMode, fNuiVar, fInputChannel);
  
  WriteInfoStatus("ProfilingClass", "finalize");

  delete val_obj;

  WriteInfoStatus("ProfilingClass", "finished");


}



//Symmetrise one-sided systematic effects in case only up- OR down-variation is available, also in GoeProfiling.c
void ProfilingClass::SymmetriseSystematicFile(std::string Syst, std::string Folder, std::string systvar, std::string systvar_new, int kvar, int kvar_new, std::string sign)
{

  // check if Systematic is new or already exists
  if(CheckNewSystematic(Syst)){
    // if new, add it to Systematic vector
    fSystematics.push_back(Syst);
    
    // and create new vector that should contain the input file path and the k-Value
    std::vector<SystematicInfo::MySystematic> HelpVec;
    fSystematicFiles.push_back(HelpVec);
  }

  int SystIndex = GetSystIndex(Syst);

  std::stringstream kval_var, kval_var_new;
  kval_var << kvar;
  kval_var_new << kvar_new;
  


  //two separate histograms
  std::string filename_var = Folder+"_"+systvar+ "_" + kval_var.str() + "k.root";
  std::string outputname_var = Folder+"_"+systvar_new+"_"+kval_var.str() + "k.root";

  TFile* file_nom = new TFile(fNomInputFile.c_str(), "open");
  TFile* file_var = new TFile(filename_var.c_str(), "open");

  //Create new output file
  TFile output_file_var(outputname_var.c_str(), "recreate");

  // read out signal and bkg information...
  for(int iSig = 0; iSig < fSignalParameter.size(); ++iSig){

    std::string hist_name = fSignalParameter[iSig].HistName;

    TH1D hist_nom = *(TH1D*) file_nom -> Get(hist_name.c_str());
    TH1D hist_var = *(TH1D*) file_var -> Get(hist_name.c_str());

    TH1D hist_new = CalculateNewVariation (hist_nom, hist_var);

    hist_new.Write(hist_name.c_str() );
  }
  for(int iBkg = 0; iBkg < fBkgParameter.size(); ++iBkg){

    std::string hist_name  = fBkgParameter[iBkg].HistName;

    TH1D hist_nom = *(TH1D*) file_nom -> Get(hist_name.c_str());
    TH1D hist_var = *(TH1D*) file_var -> Get(hist_name.c_str());

    TH1D hist_new = CalculateNewVariation (hist_nom, hist_var);

    hist_new.Write(hist_name.c_str() );
  }

  output_file_var.Close();

  //Read new symmetrised input templates

  SystematicInfo::MySystematic HelpType;
  
  WriteInfoStatus("ProfilingClass", "Add Systematic with type = "+Syst);

  HelpType.SystematicType = Syst;

  if (sign == "up") {
    HelpType.Filename_up    = Folder+"_"+systvar + "_" + kval_var.str() + "k.root";
    HelpType.Filename_down  = Folder+"_"+systvar_new + "_" + kval_var.str() + "k.root"; // has to be kval_u !!!
  }
  else {
    HelpType.Filename_up    = Folder+"_"+systvar_new + "_" + kval_var.str() + "k.root";
    HelpType.Filename_down  = Folder+"_"+systvar + "_" + kval_var.str() + "k.root";
  }

  if (sign == "up") {
    HelpType.kup            = kvar;
    HelpType.kdown          = kvar_new;
  }
  else {
    HelpType.kup            = kvar_new;
    HelpType.kdown          = kvar;

  }
  fSystematicFiles[SystIndex].push_back(HelpType);

  file_nom->Close();
  file_var->Close();

}

//Used in SymmetriseSystematicFile to calculate the desired variation depending on the given one-sided systematic effect
TH1D ProfilingClass::CalculateNewVariation(TH1D hist_nom, TH1D hist_var)
{

  //  int NBins = hist_nom.GetNbinsX();

  TH1D hist_var_new = TH1D("", "", fNBins, fLowerEdge, fUpperEdge);

  for (int i = 1; i != fNBins+1; i++) {

    double value_nom = hist_nom.GetBinContent(i);
    double value_var = hist_var.GetBinContent(i);

    double value_new = 2*value_nom - value_var;

    hist_var_new.SetBinContent(i,value_new);
  }

  return hist_var_new;
}












