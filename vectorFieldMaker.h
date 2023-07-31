#ifndef vectorFieldMaker_H_
#define vectorFieldMaker_H_
#include <iostream>
#include <vector>

using namespace std;

int static getCoffSinTimesCos(int i1,int i2,int k);

int static getCoffSinTimesSin(int i1,int i2,int k);

void static printMat(std::vector< std::vector<int> > x);

void static printVect(std::vector<int>  x);

string static coffFromInt(int x);

string static analizeRowSinSquered(std::vector<string> u, std::vector<int>  x,int row_number);

int static ComputeRowPoints(std::vector<int> x,int num);

void static optimiseMatSinSquered (std::vector< std::vector<int> >&  x);

string static analizeMatSinSquered(std::vector<string> u, std::vector< std::vector<int> >  x);

std::vector<string> sinSquered(std::vector<string> u,int maxSize);

std::vector<string> sinSquered(std::vector<string> u);

string static analizeRow(std::vector<string> u,std::vector<string> v, std::vector<int>  x,int row_number);

string static analizeMatrix(std::vector<string> u,std::vector<string> v, std::vector< std::vector<int> > x);

std::vector<string> sinTimesCos(std::vector<string> sin,std::vector<string> cos,int maxSize);

std::vector<string> sinTimesCos(std::vector<string> sin,std::vector<string> cos);

std::vector<string> sinTimesSin(std::vector<string> sin1,std::vector<string> sin2,int maxSize);

std::vector<string> sinTimesSin(std::vector<string> sin1,std::vector<string> sin2);

std::vector<string> multiply(string mult ,std::vector<string> x);

std::vector<string> multiply(std::vector<string> x, std::vector<string> y);

std::vector<string> add(std::vector<string> x, std::vector<string> y);

std::vector<string> makeStringVector(int size,int begin,int end, string name);

std::vector<string> makeStringVectorWithGap(int size,int gap,string name);

std::vector<string> makeStringVector(int size,string name);

std::vector<string> merge(std::vector<string> x, std::vector<string> y);


std::vector<string> shuffle(std::vector<string> x, std::vector<string> y);

string toFormula(std::vector<string> parameters, std::vector<string> variables , std::vector<string> vectorField);

string toFormula(std::vector<string> variables,std::vector<string> vectorField);

void split(std::vector<string> u,std::vector<string>& uMain,std::vector<string>& uDiss,int mainModesNum);

#endif
