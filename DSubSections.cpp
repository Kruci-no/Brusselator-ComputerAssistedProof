#include "DSubSections.h"
using namespace std;
using namespace capd ;
void saveSectionsToFile(std::fstream& file,SectionCordinatesArray sectionArray){
  int n = sectionArray.size();
	file.precision(18);
	file << sectionArray.size()<<endl;
	for(int i=0;i<n;i++){
		file << sectionArray[i].orgin << '\n';
		file << sectionArray[i].normal << '\n';
		file << sectionArray[i].cordinates << '\n';
	}


}
SectionCordinatesArray readSectionsFromFile(std::fstream& file){
  int n;
  file >> n ;
	SectionCordinatesArray sectionArray(n);
  for(int i=0;i<n;i++){
    file >> sectionArray[i].orgin;
		file >> sectionArray[i].normal;
		file >> sectionArray[i].cordinates;
  }
  return sectionArray;


}
