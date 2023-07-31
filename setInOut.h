#ifndef SetInOut_H_
#define SetInOut_H_
#include <iostream>
#include <fstream>
#include <string>
struct InOut
{
     std::ifstream paramsFile;
     std::ifstream initialValueFile;
     std::ifstream initialValuePoincareFile;
     std::ifstream sectionFile;
     std::ofstream outPutFile;
     std::fstream outPutSubSectionsFile;
     InOut();
     InOut(std::string initFileName);
     void contruct(std::string initFileName);
};
#endif //_EXAMPLES_PROJECTSTARTER_OUTPUT_H_
