#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "TF1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "TemplateInfo.h"
#include "SystematicInfo.h"

class tinterpolation
{
  public:
  
  tinterpolation(std::string, int); 	//constructor

  ~tinterpolation();	//destructor
	
  TH1D interpolate(std::string, std::vector<double> &);
  TH1D *interpolate(std::string, TH1D*, std::vector<double> &); 

 protected:
  
  int _hsize;   
  std::vector<TH1D> _hist;
  std::vector<TH1D> _hist_2k;
  std::vector<TH1D> _hist_3k;

  TH1D* result00_help;

  int fParamNumber;

  //  int fNBins;

};

extern double kfix;
extern int    fNBins;

TH1D* interpolate_single ( std::vector<TH1D*>&, TH1D*, double); 
double sign(double);
double max(double, double);

TH1D ReturnInterpolationHistogramScaled(int, std::string, std::vector<double> &, double, bool);

void CreateInterpFunction (const int);
std::vector<double> CreateSlopes (std::vector<double> &,   std::vector<TH1D*> &, std::vector<TH1D*> &, std::vector<TH1D*> &, const int);
std::vector<TF1>    CreateFunctions (std::vector<TH1D*> &, std::vector<TH1D*> &, std::vector<TH1D*> &, const int);

TGraphAsymmErrors * TemplateFitErrors (double*, const int, TH1D, std::vector<TH1D> &, std::vector<double> &, TGraphAsymmErrors *, int, int); 

TH1D TemplateFitErrors (double* , std::vector<TH1D>, std::vector<double>, int);


#endif
