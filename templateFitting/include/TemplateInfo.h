#ifndef TemplateInfo_h
#define TemplateInfo_h

#include <TH1D.h>

#include <string>
#include <iostream>

#include <vector>

class TemplateInfo {

 public:

  TemplateInfo();
  ~TemplateInfo();

  struct MySignalTemplate{
    std::string ParamName;
    std::string HistName;
    std::string EffName;
    TH1D        hist;
    TH1D        histTS;
    double      SelEff;
  };

  struct MyBkgTemplate{
    std::string ParamName;
    std::string HistName;
    std::string ExpName;
    std::string UncName;
    TH1D        hist;
    TH1D        histTS;
    double      Exp;
    double      Unc;
  };

  struct SigInterpolObject{
    std::string SystematicType;
    std::vector<TH1D>  hist;
    std::vector<double> k;
    double SelEff;
  };
  
  struct BkgInterpolObject{
    std::string SystematicType;
    std::vector<TH1D>  hist;
    std::vector<double> k;
    double BkgExp;
    double BkgUnc;
  };

};

#endif
