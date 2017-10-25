#ifndef CalibrationCurves_H
#define CalibrationCurves_H

#include "TGraphErrors.h"

#include <vector>

class CalibrationCurves
{
 public:
  
  CalibrationCurves();      //constructor
  ~CalibrationCurves();     //destructor
  
  void SetOutputFolder(std::string output){fOutputFolder = output;};
  void SetInputFolderList(std::string input_folder){fInputFolderList.push_back(input_folder);};
  void SetXValueList(double i){fXValueList.push_back(i);};

  void FillGraphs(std::string, std::string, std::string, bool, bool);

  void ReadInputValues(std::string, int, bool, std::string);

 protected:

  std::vector<std::string> fInputFolderList;
  std::vector<double>      fXValueList;
  std::string fOutputFolder;

  TGraphErrors *fCalibCurve;
  TGraphErrors *fPullCurve;
  TGraphErrors *fRMSCurve;
  TGraphErrors *fPullCombCurve;
  TGraphErrors *fPullProfileCurve;
  TGraphErrors *fRMSProfileCurve;
  TGraphErrors *fPullCombProfileCurve;
  TGraphErrors *fPullCurveGauss;
  TGraphErrors *fRMSCurveGauss;
  TGraphErrors *fPullCombCurveGauss;
  TGraphErrors *fPullProfileCurveGauss;
  TGraphErrors *fRMSProfileCurveGauss;
  TGraphErrors *fPullCombProfileCurveGauss;


  TGraphErrors *fDiffCalibCurve;

  double  fCalibCurveVecX[7];
  double  fCalibCurveVecY[7];
  double  fPullCurveVecX[7];
  double  fPullCurveVecY[7];
  double  fRMSCurveVecX[7];
  double  fRMSCurveVecY[7];
  //double  fPullCombCurveVecX[7];
  //double  fPullCombCurveVecY[7];
  double  fPullProfileCurveVecX[7];
  double  fPullProfileCurveVecY[7];
  double  fRMSProfileCurveVecX[7];
  double  fRMSProfileCurveVecY[7];
  //double  fPullCombProfileCurveVecX[7];
  //double  fPullCombProfileCurveVecY[7];

};

#endif
