#ifndef __SOLVERDISSPDE_H__
#define __SOLVERDISSPDE_H__


#include <iostream>
#include "capd/capdlib.h"
#include "algebra.h"
#include <chrono>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <limits>
#include <utility>
typedef std::unordered_map<std::string,capd::interval> ParamsMap;

struct Indexer{
    std::vector<std::pair<int,int> > pairs;
    std::pair<SeriesVector, SeriesVector > splitVector (SeriesVector& x);
    int size();
    capd::IVector getIVector(SeriesVector& x);
    void intersectRepresetations(capd::InclRect2Set& mainModes,SeriesVector& x);
    void makeCosistend(capd::InclRect2Set& mainModes,SeriesVector& x);
};
struct VectorField{
    int nMain;
    int numOfNormalParams;
    std::function<SeriesVector(SeriesVector&, ParamsMap)> F_nonLinearity;
    std::function<capd::IVector (SeriesVector&, Indexer , ParamsMap)> F_rest;
	capd::IMap f;
    Indexer indexer;
    ParamsMap params;
    SeriesVector L;
    SeriesVector computeNonLinearity(SeriesVector x);
    capd::IVector computeInclusion(SeriesVector x);
};
struct Set{
    SeriesVector vector;
    capd::InclRect2Set mainModes;
    Set(SeriesVector vector,capd::InclRect2Set mainModes);
    void intersectRepresetations(Indexer indexer);
    void makeCosistend(Indexer indexer);
    capd::interval getCurrentTime();
};
struct Encloser{
    SeriesVector enclosureExtent; 
    SeriesVector enclosurePointWiseExtent; 
    capd::interval validatedTimeStep;
    void enclose(
                Set& x,
                VectorField& vectorField,
                capd::interval dt,
                int refineNum,
                bool constStep, 
                bool comPointWiseEnclose);
};

struct Mover{
    
    capd::interval step;
    Encloser encloser;
    capd::IMap perturbMap;
	capd::IMultiMap* multiMap;
	capd::CWDiffInclSolver* solver;
    capd::interval maxDecay;
    Mover(VectorField& vectorField,Encloser encloser);
    void setStep(capd::interval step);
    void move(Set& x,VectorField& vectorField,bool constStep);
    void static perturb(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, 
                        capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams);
};
struct Section{
	capd::IAffineSection sectionMainModes;
	Section (capd::IAffineSection sectionMainModes);
	capd::interval evalAt(Set& set);
};
struct PoincareMapDiss{
    capd::interval reachSectionTime;
	Encloser& encloser;
	Mover& mover;
	Section section;
	PoincareMapDiss(Encloser& encloser,Mover& mover ,Section section);
    std::pair<capd::IVector, SeriesVector> compute(Set& set,capd::IVector x0 ,capd::IMatrix A,VectorField& vectorField );
};
#endif // __SOLVERDISSPDE_H__