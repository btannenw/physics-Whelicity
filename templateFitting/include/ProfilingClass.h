#ifndef ProfilingClass_H_
#define ProfilingClass_H_

#include "TemplateInfo.h"
#include "SystematicInfo.h"
#include "StatusLogbook.h"
#include "ExternalSystematics.h"

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>

extern std::string fInterpolationMethod;
extern std::vector<std::vector<std::vector<TF1> > > fFitFunc;
extern std::vector<std::vector<std::vector<TF1> > > fFitFunc_QuadFit;
extern int    fNBins;
extern double fLowerEdge;
extern double fUpperEdge;

extern TH1D hist_help;
extern TH1D hist_sum;
extern TH1D hist;

class ProfilingClass{
  
 public:

  ProfilingClass();
  virtual ~ProfilingClass();
  
  //fProf.AddSignalParameter("N0",  "CosTheta_F0", "F0eff");
  void AddSignalParameter(std::string param, std::string hist, std::string eff)
  {

    TemplateInfo::MySignalTemplate HelpTemplate;

    HelpTemplate.ParamName = param;
    HelpTemplate.HistName  = hist;
    HelpTemplate.EffName   = eff;

    fSignalParameter.push_back(HelpTemplate);

  };

  void AddBackgroundParameter(std::string param, std::string hist, std::string exp, std::string unc)
  {

    TemplateInfo::MyBkgTemplate HelpTemplate;

    HelpTemplate.ParamName = param;
    HelpTemplate.HistName  = hist;
    HelpTemplate.ExpName   = exp;
    HelpTemplate.UncName   = unc;

    fBkgParameter.push_back(HelpTemplate);

  };

  void SetFractions(double F0, double FL, double FR){

    std::vector<double> help;
    help.push_back(F0);
    help.push_back(FL);
    help.push_back(FR);

    fFraction.push_back(help);

  };

  void ProducePseudoData(double , double , double , std::vector<std::string>, std::vector<std::string>, std::string,  std::string);

  void AddNuisanceParameter(std::string nui){fNuisanceParameter.push_back(nui);};

  void AddSystematicFile(std::string, std::string, std::string, std::string, int, int);
  
  void AddSystematicPseudoData(std::string, std::string, std::string);

  void SetTemplateFittingMode(std::string fittingMode){fTemplateFittingMode = fittingMode;};

  void SetOutputTxtFile(std::string file){fOutputTxtFile = file;};

  void SetEvaluationMode(std::string mode){fEvaluationMode = mode;};

  void SymmetriseSystematicFile(std::string, std::string, std::string, std::string, int, int, std::string);

  TH1D CalculateNewVariation(TH1D, TH1D);

  bool CheckNewSystematic(std::string);

  void DoValidation();

  void DoExternalSystematicEvaluation(std::string);

  void FillInterpolationObjects();

  int  GetSystIndex(std::string);

  void ReadNominalInfo();

  void ReadSystematicInfo();

  void MakeInterpolationCurves();

  void SetDataHist(std::string hist){fDataHist = hist;};

  void SetFitMethod(std::string method){fFitMethod = method;};
  void SetInputChannel(std::string channel){fInputChannel = channel;};

  void SetInterpolationMethod(std::string method){

    fInterpolationMethod = method;
    
    WriteParameterStatus("interpolation", "Using the interpolation method "+fInterpolationMethod);

  };

  void SetSubplotLabels(std::string Lepton, std::string Jet, std::string BTag){

    std::cout << Lepton.c_str() << std::endl;

    fLeptonLabel.push_back(Lepton);
    fJetBinLabel.push_back(Jet);
    fBTagLabel.push_back(BTag);

  };

  void SetOutputFolder(std::string folder){fOutputFolder = folder;};
  void SetBackgroundMode(std::string mode){fBkgMode = mode;}; //mjk: used in GoeProfiling.c
  void SetLumi(double lumi){fLumi = lumi;};
  void SetNominalInputFile(std::string file){fNomInputFile = file;};
  void SetNuiVariation(int nui){fNuiVar = nui;};
  void SetNumberOfPE(int pe){fNPseudoExp = pe;};
  void SetOutputSampleNumber(int nr){fNOutputSample = nr;};
  void SetValidationMode(std::string mode){fValidationMode = mode;};
  void SetXsec(double xsec){fXsec = xsec;};
  void SetNumberOfBins (int bins){fNBins = bins;};
  void SetLowerEdge(double low){fLowerEdge = low;};
  void SetUpperEdge(double up){fUpperEdge = up;};

 protected:

  std::vector<std::vector<double> > fFraction;
  
  int fNOutputSample;
  int fNPseudoExp;
  int fNSyst;
  int fNuiVar;
  //  int fNBins;

  //  double fLowerEdge;
  //  double fUpperEdge;

  //  int fNBkg;
  //  int fNSignal;

  double fLumi, fXsec, fKfactor;

  std::string fBkgMode;
  std::string fDataHist;
  std::string fFitMethod;
  std::string fOutputFolder;
  std::string fInputChannel;
  std::string fNomInputFile;
  std::string fValidationMode;
  std::string fEvaluationMode;
  std::string fOutputTxtFile;
  std::string fTemplateFittingMode;

  //  std::vector<std::vector<std::vector<TF1> > > fFitFunc;
  std::vector<std::string> fSystematics;
  std::vector<std::string> fNuisanceParameter;
  std::vector<TemplateInfo::MyBkgTemplate> fBkgParameter;
  std::vector<TemplateInfo::MySignalTemplate> fSignalParameter;
  std::vector<std::vector<SystematicInfo::MySystematic> > fSystematicFiles;
  std::vector<SystematicInfo::MySystematicPD> fSystematicPD;

  std::vector<std::string> fLeptonLabel,fJetBinLabel,fBTagLabel;

};

#endif
