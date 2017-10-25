// info: got colours from http://unitedsoft.ch/cpp/index.php/Console_colors

#include "StatusLogbook.h"

#include <iostream>
#include <cstring>

using namespace std;

int fMaxStrLen  = 110;

std::string fFinalStr        = " =======";
std::string fInitialStrInfo  = "======= INFO:::";
std::string fInitialStrError = "======= ERROR::";
std::string fInitialStrParam = "======= PARAM::";

void WriteParameterStatus(std::string classname, std::string info)
{

  int fNWhitespace = fMaxStrLen - 32 - strlen(info.c_str()) - strlen(classname.c_str());

  if(fNWhitespace <= 0) std::cout << "\033[31m" << "ERROR-PARAM::StatusLogbook: length of cout string is <= 0 !!!" << std::cout;

  std::string Whitespace   = "";

  for(int i = 0; i < fNWhitespace; ++i)
    Whitespace += " ";

  std::string outputstring = fInitialStrParam+classname+": "+info+Whitespace+fFinalStr;

  std::cout << "\033[34m" << "=======================================================================================================" << std::endl;
  std::cout << "\033[34m" << "=======                                                                                         =======" << std::endl; 
  std::cout << "\033[34m" << outputstring.c_str() << std::endl;
  std::cout << "\033[34m" << "=======                                                                                         =======" << std::endl;
  std::cout << "\033[34m" << "=======================================================================================================" << std::endl;

}


void WriteInfoStatus(std::string classname, std::string info)
{
  
  int fNWhitespace = fMaxStrLen - 32 - strlen(info.c_str()) - strlen(classname.c_str());

  if(fNWhitespace <= 0) std::cout << "\033[31m" << "ERROR-INFO::StatusLogbook: length of cout string is <= 0 !!!" << std::cout;

  std::string Whitespace = "";

  for(int i = 0; i < fNWhitespace; ++i)
    Whitespace += " ";
  
  std::string outputstring = fInitialStrInfo+classname+": "+info+Whitespace+fFinalStr;

  std::cout << "\033[32m"  << outputstring.c_str() << "\33[0m" << std::endl;

}


void WriteErrorStatus(std::string classname, std::string info)
{

  int fNWhitespace = fMaxStrLen - 32 - strlen(info.c_str()) - strlen(classname.c_str());

  if(fNWhitespace <= 0) std::cout << "\033[31m" << "ERROR-ERROR::StatusLogbook: length of cout string is <= 0 !!!" << std::cout;

  std::string Whitespace = "";

  for(int i = 0; i < fNWhitespace; ++i)
    Whitespace += " ";

  std::string outputstring = fInitialStrError+classname+": "+info+Whitespace+fFinalStr;

  std::cout << "\033[31m" << outputstring.c_str() << std::endl;

}

