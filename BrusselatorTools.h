#ifndef BrusselatorTools_H
#define BrusselatorTools_H
#include <iostream>
#include <vector>
#include "capd/capdlib.h"
#include "vectorFieldMaker.h"
#include "solverDissPDE.h"
#include "algebra.h"
using namespace capd;
using namespace std;

string BrusselatorVectorField(int size) ;

string BrusselatorVectorFieldSym(int size);

string BrusselatorVectorFieldSym2(int size);

string BrusselatorVectorFieldFRest(int size,int mainModesNum);

string BrusselatorVectorFieldFDiss(int size, int mainModesNum);

string BrusselatorVectorFieldFMain(int size, int mainModesNum);

SeriesVector linearBursselator(int size,ParamsMap params);

SeriesVector BursselatorNonlinearity(SeriesVector x,ParamsMap params);


Series u_squere_v(SeriesVector x,ParamsMap params);

IVector BrusselatorRest(SeriesVector x,Indexer indexer,ParamsMap params);

Indexer makeBrusselatorIndexer(int u_v_size);

VectorField getBrusselatorVectorField(int u_v_main_size,int serisesSizes,ParamsMap params);


#endif
