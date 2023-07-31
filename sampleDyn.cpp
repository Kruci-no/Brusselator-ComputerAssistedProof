/*start time duration time stepNum */
# include <iostream>
# include "capd/capdlib.h"
# include "BrusselatorTools.h"
# include "SetInOut.h"
using namespace capd ;
using namespace std ;
using capd::autodiff::Node;
#include <cmath>
#include <vector>
/*start time duration time stepNum */
int main (int argc, char* argv[]){
  InOut inOut;
	DVector  Params;
	inOut.paramsFile >> Params;
	DVector  initialValue;
	inOut.initialValueFile >> initialValue;
	int dim = initialValue.dimension();
	int dim_u=dim/2;
	int noParam=Params.dimension();
  double durationTime = 1;
  double startTime = 0;
  int stepNum = 1000;
  int n = 5;
  if(argc == 3){
     durationTime = stod(argv[1]);
     stepNum = stoi(argv[2]);
  }
  if(argc == 4){
     startTime = stod(argv[1]);
     durationTime = stod(argv[2]);
     stepNum = stoi(argv[3]);
  }
  if(argc == 5){
     startTime = stod(argv[1]);
     durationTime = stod(argv[2]);
     stepNum = stoi(argv[3]);
     n = stoi(argv[4]);
  }
  std::cout << "dimension"<< dim << '\n';
  cout << "d1: " << Params[0] <<'\n';
  cout << "d2: " << Params[1] <<'\n';
  cout << "B:  " << Params[2] <<'\n';
  cout << "A:  " << Params[3] <<'\n';
  std::cout << "dt="<< durationTime/stepNum << '\n';
  std::cout << initialValue << '\n';
  DMap Brusselator (BrusselatorVectorFieldSym(dim_u));
	Brusselator.setParameters(Params);
	DOdeSolver Dsolver(Brusselator,22);
  DTimeMap timeMap(Dsolver);
  string x="{";
  DVector v;
  vector<double> tabNormU(n, 0);
  vector<double> tabNormV(n, 0);
  vector<double> tabMaxNormU(n, 0);
  vector<double> tabMaxNormV(n, 0);
  timeMap(startTime,initialValue);
  for(int i = 0; i<=stepNum; i++){

     x=x+"{";
     v=timeMap((durationTime*1.0)/stepNum,initialValue);
     for(int j=0;j<dim-1;j++){
       x=x+to_string(v[j])+",";
     }
     for(int j=0;j<dim/2;j++){
       for(int k=0;k<n;k++)
       {
         tabNormU[k]+=pow ((2*j+1),2*k)*v[2*j]*v[2*j];
         tabNormV[k]+=pow ((2*j+1),2*k)*v[2*j+1]*v[2*j+1];

       }
     }

     x=x+to_string(v[dim-1]);
     x= x +"},\n";
     for(int k=0;k<n;k++){
       if(tabMaxNormU[k] < tabNormU[k]){
         tabMaxNormU[k]=tabNormU[k];
       }
       if(tabMaxNormV[k] < tabNormV[k]){
         tabMaxNormV[k]=tabNormV[k];
       }
       tabNormV[k] = 0;
       tabNormU[k] = 0;
     }
    //x = x +  timeMap((finalTime*1.0)/stepNum,initialValue) +",\n";
    //std::cout <<"At time:"<< finalTime * ((i * 1.0)/stepNum) <<" "
              //inOut.outPutFile << fixed <<timeMap((finalTime*1.0)/stepNum,initialValue)<< ","<< '\n';
  }
  x.pop_back();
  x.pop_back();
  x = x+"}";
  inOut.outPutFile << x;
  for(int k=0;k<n;k++){
    std::cout<<"Max L^2 norm u_x^"<<k<<": "<< sqrt(tabMaxNormU[k]) << '\n';
    std::cout<<"Max L^2 norm v_x^"<<k<<": "<< sqrt(tabMaxNormV[k]) << '\n';
  }

}
