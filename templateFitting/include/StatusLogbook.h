#ifndef StatusLogbook_H_
#define StatusLogbook_H_

#include <string>
#include <iostream>

extern int fMaxStrLen;

extern std::string fFinalStr;
extern std::string fInitialStrInfo;
extern std::string fInitialStrError;
extern std::string fInitialStrParam;

void WriteInfoStatus(std::string,  std::string);
void WriteErrorStatus(std::string, std::string);
void WriteParameterStatus(std::string, std::string);

#endif
