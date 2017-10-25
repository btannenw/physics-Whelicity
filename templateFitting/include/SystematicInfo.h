#ifndef SystematicInfo_h
#define SystematicInfo_h

#include <TH1D.h>

#include <string>
#include <iostream>

#include <vector>

class SystematicInfo {

 public:

  SystematicInfo();
  ~SystematicInfo();

  struct MySystematic{
    std::string SystematicType;
    std::string Filename_up;
    std::string Filename_down;
    double      kup;
    double      kdown;
    std::vector<TH1D>  HistUp;
    std::vector<TH1D>  HistDown;
    std::vector<std::string> ParamName;
  };

  // PD = Pseudo Data
  struct MySystematicPD{
    std::string SystematicType;
    TH1D histoPD;
  };

};

#endif
