using namespace std;
#include "setInOut.h"

void InOut::contruct(std::string initFileName)
{
  ifstream initFile(initFileName);
  string fileName;
  initFile >> fileName;
  ifstream paramsFile(initFileName);
  this->paramsFile.open(fileName) ;
  initFile >> fileName;
  this->initialValueFile.open(fileName) ;
  initFile >> fileName;
  this->outPutFile.open(fileName);
  initFile >> fileName;
  this->sectionFile.open(fileName);
  initFile >> fileName;
  this->initialValuePoincareFile.open(fileName);
  initFile >> fileName;
  this->outPutSubSectionsFile.open(fileName,std::ios::in | std::ios::out);
  //this->outPutSubSectionsFile <<"eloo";
  initFile.close();
}
InOut::InOut(string initFileName)
{
  contruct(initFileName);
}
InOut::InOut()
{
  contruct("init.txt");
}
