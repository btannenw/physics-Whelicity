#ifndef FILL_HISTOGRAMS_H
#define FILL_HISTOGRAMS_H

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "OutputTreeReader.h"

// extern std::vector<Double_t> nui_value;//(num_nuisance);
// extern std::vector<Double_t> nui_nom;//(num_nuisance);
// extern std::vector<Double_t> nui_err;//(num_nuisance);
// extern std::vector<Double_t> nui_prof_err;//(num_nuisance);
// extern std::vector<TBranch*> b_nui_value;//(num_nuisance);
// extern std::vector<TBranch*> b_nui_nom;//(num_nuisance);
// extern std::vector<TBranch*> b_nui_err;//(num_nuisance);
// extern std::vector<TBranch*> b_nui_prof_err;//(num_nuisance);

class histograms
{
 public:
  
  histograms();      //constructor
  ~histograms();     //destructor
  
  struct Parameter{
    std::string name;
    std::string label;
    double      nom;
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
    double      pull_mean_prof;
    double      pull_mean_prof_err;
    double      pull_rms_prof;
    double      pull_rms_prof_err;
    double      gauss_pull_mean_prof;
    double      gauss_pull_mean_prof_err;
    double      gauss_pull_rms_prof;
    double      gauss_pull_rms_prof_err;
    double      Up;
    double      Low;
    double      Up_unc;
    double      Low_unc;
    double      Up_prof_unc;
    double      Low_prof_unc;
    TH1D*       hist_param;
    TH1D*       hist_pull;
    TH1D*       hist_pull_prof;
    TH1D*       hist_unc;
    TH1D*       hist_prof_unc;
    std::vector<TH2D> hist2D;
    std::vector<double> val_vec;
  };
  
  void SetInputFolder(std::string input){fInputFolder = input;};
  void SetOutputFolder(std::string output){fOutputFolder = output;};
  void SetParameterList(std::vector<Parameter> param){fParameter = param;};
  void SetValidationMode(std::string validation_mode){fValidationMode = validation_mode;};

  void FindRanges(int);
  void FillHistograms(std::string, int);
  void MakeChainOfFiles(int);
  void MakePlots(TH1D* , double, std::string, std::string, double, double, std::string, std::string, int, int);
  void MakeErrorDataPlots(TH1D* , double, double, std::string, std::string, double, double, std::string, std::string, int);

  std::vector<double> MakeFits(TH1D*, std::string, double, double);

  void MakeHTML(std::string);

  TFile *GetInputFile(std::string);
  std::vector<double> GetErrorVector(TFile *);

  TTree *InitializeOutputTree(std::string, int);

  void MakeOutputTrees(int);

 protected:

  OutputTreeReader *fInputTree;

  std::string fInputFolder;
  std::string fOutputFolder;
  std::string fValidationMode;

  std::vector<std::string> fFileNames;
  std::vector<Parameter>   fParameter;

  int fTotalEntries;

};

#endif
