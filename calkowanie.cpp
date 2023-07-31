#include <iostream>
#include "capd/capdlib.h"
#include "algebra.h"
#include <chrono>
#include <vector>
#include <cmath>
#include<unordered_map>
#include <limits>
#include <utility>
#include "solverDissPDE.h"
#include "BrusselatorTools.h"
#include "SetInOut.h"
#include "DSubSections.h"
using namespace std;
using namespace capd;
using capd::autodiff::Node;




int main (){
	InOut inOut;
	DVector  paramsDVector;
	inOut.paramsFile >> paramsDVector;
    ParamsMap params = {{"d1",paramsDVector[0]},{"d2",paramsDVector[1]},{"B",paramsDVector[2]},{"A",paramsDVector[3]}};
	
	int u_vMainSize = 18;
	int seriesesLen = 20;
	VectorField brusselatorVectorFieldPDESym = getBrusselatorVectorField(u_vMainSize,seriesesLen,params);
    
    Encloser encloser;
    Mover mover(brusselatorVectorFieldPDESym,encloser);
    
    IVector init(2*u_vMainSize) ;
    init[0]= interval(0.1);
    init[1] = interval(1);
    //init[2*u_vMainSize] = interval(1);
    InclRect2Set mainSet(init);
    
    auto seriesVector = SeriesVector(2);
    for(int i=0;i<2;i++){   
        	seriesVector[i] = Series(interval(0),interval(6),SeriesType::sin_odd);
        	seriesVector[i] = seriesVector[i].resize(seriesesLen);  
		   
    }
	//seriesVector[0].C = interval(-1,1) ;
	//seriesVector[1].C = interval(-1,1) ;

	
	//seriesVector[2].C= interval(0,1);
	//seriesVector[2].s = interval(-3);
	
	//cout << seriesVector[2].valueAt(seriesesLen+1);
	//cout << seriesVector[2].valueAt(seriesesLen+2);
    
    Set set(seriesVector,mainSet);
    
    set.makeCosistend(brusselatorVectorFieldPDESym.indexer);
    set.vector.print();
    for(int i=0;i<10;i++){
        mover.move(set,brusselatorVectorFieldPDESym,true);
        //set.vector.print();      
    }
    set.vector.print();
	cout<<"Czas: "<< set.getCurrentTime();
}
