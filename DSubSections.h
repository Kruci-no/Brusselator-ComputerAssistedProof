#ifndef DSubSections_H_
#define DSubSection_H_
#include <iostream>
# include "capd/capdlib.h"
struct Section_Cordinates
{
	capd::DVector orgin;
	capd::DVector normal;
	capd::DMatrix cordinates;
	
};
typedef std::vector<Section_Cordinates> SectionCordinatesArray;

void saveSectionsToFile(std::fstream& inOut,SectionCordinatesArray sectionArray);
SectionCordinatesArray readSectionsFromFile(std::fstream& inOut);
//void saveSectionsToFile(std::fstream& inOut,SectionCordinatesArray sectionArray);
#endif //_EXAMPLES_PROJECTSTARTER_OUTPUT_H_
