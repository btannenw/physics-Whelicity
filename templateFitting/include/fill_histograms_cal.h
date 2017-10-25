#ifndef FILL_HISTOGRAMS_CAL_H
#define FILL_HISTOGRAMS_CAL_H

#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "OutputTreeReader_cal.h"

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
  
  typedef struct Parameter{
    std::string name;
    std::string label;
    double      nom;
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
  };

  void SetInputFolder(std::string input){fInputFolder = input;};
  void SetOutputFolder(std::string output){fOutputFolder = output;};
  void SetParameterList(std::vector<Parameter> param){fParameter = param;};

  void FindRanges(int);
  std::vector<double> FillHistograms(int, std::string, int, double, double);
  void MakeChainOfFiles(int);
  void MakePlots(TH1D* , double, std::string, std::string, double, double, std::string, std::string, int);

 protected:

  OutputTreeReader *fInputTree;

  std::string fInputFolder;
  std::string fOutputFolder;

  std::vector<Parameter> fParameter;

  int fTotalEntries;

};

#endif
