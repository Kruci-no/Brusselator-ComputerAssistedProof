#include <iostream>
#include "capd/capdlib.h"
#include "algebra.h"
#include <chrono>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <limits>
#include <utility>
#include "solverDissPDE.h"
using namespace capd;
using namespace std;


std::pair<SeriesVector, SeriesVector > Indexer::splitVector(SeriesVector& x) 
{
    SeriesVector x1 = x*interval(0);
    SeriesVector x2 = x;
    for(int i=0;i<this->pairs.size();i++){
            int l = pairs[i].first;
            int r = pairs[i].second; 
            x1[l].main[r] = x2[l].main[r];
            x2[l].main[r] = interval(0);
    }
    return std::pair<SeriesVector, SeriesVector>(x1,x2); 
}

int Indexer::size() 
{
    return pairs.size();    
}

capd::IVector Indexer::getIVector(SeriesVector& x) 
{
    IVector result(size());
        for(int i=0;i<size();i++){
            int l = pairs[i].first;
            int r = pairs[i].second; 
            result[i] = x[l].main[r];
        }
        return result;
}

void Indexer::intersectRepresetations(capd::InclRect2Set& mainModes,SeriesVector& x) 
{
    IVector main = (IVector) mainModes;
    for(int i=0;i<size();i++){
        int l = pairs[i].first;
        int r = pairs[i].second; 
        if(!intersection(main[i], x[l].main[r], x[l].main[r]) ){
            x.print();
            cout<< main<<"\n";
            throw std::runtime_error("Set intersectRepresetations - empty intersection\n");
        }
    }
}

void Indexer::makeCosistend(capd::InclRect2Set& mainModes,SeriesVector& x) 
{
    IVector main = (IVector) mainModes;
    for(int i=0;i<size();i++){
        int l = pairs[i].first;
        int r = pairs[i].second; 
        x[l].main[r] = main[i];
    }
}
SeriesVector VectorField::computeNonLinearity(SeriesVector x) 
{
    return F_nonLinearity(x,params);
}

capd::IVector VectorField::computeInclusion(SeriesVector x) 
{
    return F_rest(x,indexer,params);
}


Set::Set(SeriesVector vector,capd::InclRect2Set mainModes) 
:vector(vector),mainModes(mainModes)
{
}

void Set::intersectRepresetations(Indexer indexer) 
{
    indexer.intersectRepresetations(mainModes,vector);    
}

void Set::makeCosistend(Indexer indexer) 
{
    indexer.makeCosistend(mainModes,vector);
}

interval Set::getCurrentTime() 
{
    return mainModes.getCurrentTime();
}
void Encloser::enclose(Set& x,VectorField& vectorField,interval dt,int refineNum = 0,bool constStep = false,bool comPointWiseEnclose = false){
    SeriesVector L = vectorField.L;
    SeriesVector eps(x.vector.vec.size(),interval(-0.01,0.01)*1,x.vector.vec[0].s,SeriesType::sin_odd);
    enclosureExtent = eps;
    enum class Status {nonValidated,validated,end};
    Status status = Status::nonValidated;
    SeriesVector image;
    while(status != Status::end){
        if(dt < 1e-50){
            throw std::runtime_error("cannot find enlosure");
        }
        SeriesVector expLminusOne = expMinusOne(L*dt*interval(0,1));
        SeriesVector x_Extend = x.vector + enclosureExtent;
        SeriesVector nonLinearity = vectorField.computeNonLinearity(x_Extend);
        SeriesVector mainVectorField = dt * interval(0,1)*(elementWiseMult(L,x_Extend) + nonLinearity);
        SeriesVector rightBound = elementWiseMult(expLminusOne,x.vector.upperBound()+ elementWiseMult(nonLinearity.upperBound(),L.elementWiseInverse()));
        SeriesVector leftBound = elementWiseMult(expLminusOne,x.vector.lowerBound()+ elementWiseMult(nonLinearity.lowerBound(),L.elementWiseInverse()));            
        SeriesVector hull = leftBound.lowerBound()*interval(0,1)+ rightBound.upperBound() * interval(0,1);
        //hull = semiIntersection( elementWiseMult(expLminusOne,x.vector + elementWiseMult(nonLinearity,L.elementWiseInverse())) ,hull);
        image = semiIntersection(mainVectorField,hull);
        switch(status){
            case Status::nonValidated: 
                if(image.subsetInterior(enclosureExtent)){
                    enclosureExtent = image;
                    if(refineNum == 0 )
                
                        status = Status::end;
                    else
                        status = Status::validated;
                    
                }
                else{
                    if(!constStep){
                        dt = dt*interval(1./2);
                    }
                    enclosureExtent = (image + eps)*interval(-0.1,1.2);
                }
                break; 
            case Status::validated: 
                enclosureExtent = semiIntersection(enclosureExtent,image);
                refineNum = refineNum - 1;
                if(refineNum == 0 ){
                    status = Status::end;
                    if(comPointWiseEnclose){
                        SeriesVector dissPointWiseBounds = elementWiseMult(expLminusOne,x.vector + elementWiseMult(nonLinearity,L.elementWiseInverse()));
                        enclosurePointWiseExtent = semiIntersection(mainVectorField,dissPointWiseBounds) ; 
                        for(int i=0; i<enclosurePointWiseExtent.vec.size();i++){
                            enclosurePointWiseExtent[i] = enclosurePointWiseExtent[i].resize(x.vector[i].mainSize);
            
                        }
                    }
                }
                break; 
            case Status::end: 
                break;
        }  
        for(int i=0; i<enclosureExtent.vec.size();i++){
            enclosureExtent[i] = enclosureExtent[i].resize(x.vector[i].mainSize);
            
        }

    }
    validatedTimeStep = dt;
}
Mover::Mover(VectorField& vectorField,Encloser encloser) 
{
    this->step = interval(1./512);
    this->perturbMap = IMap(perturb, vectorField.nMain, vectorField.nMain, vectorField.nMain);
    this->multiMap = new IMultiMap(vectorField.f, this->perturbMap);
    int order = 7;
    solver = new  CWDiffInclSolver( *(this->multiMap) ,order , IMaxNorm() );
}

void Mover::setStep(interval step) 
{
    this->step = step;
}

void Mover::move(Set& x,VectorField& vectorField,bool constStep = false) 
{
    int order = 7;
     CWDiffInclSolver solver( *(this->multiMap) ,order , IMaxNorm() );
    
    //solver = new  CWDiffInclSolver( *(this->multiMap) ,order , IMaxNorm() );
    encloser.enclose(x,vectorField,step,4,constStep);
    SeriesVector x_Extend= x.vector + encloser.enclosureExtent; 
    //encloser.enclosureExtent.print();
    SeriesVector nonLinearity = vectorField.computeNonLinearity(x_Extend);
    interval dt = encloser.validatedTimeStep;
    SeriesVector dtL  = dt * vectorField.L;
    SeriesVector nonLinearDuhamelTerm = elementWiseMult(expMinusOne(dtL),
                                                        elementWiseMult( vectorField.L.elementWiseInverse(),nonLinearity));

    SeriesVector x_Eval1 = elementWiseMult(exp(dtL,interval(0)) ,x.vector) + nonLinearDuhamelTerm;
    SeriesVector x_Eval2 = elementWiseMult(exp(dtL, -vectorField.L[0].s*(1./4)) ,x.vector) + nonLinearDuhamelTerm;
    interval x_Eval1Sum = interval(0);
    interval x_Eval2Sum = interval(0);
    for(int i=0; i<x.vector.vec.size();i++){
        x_Eval1[i] = x_Eval1[i].resize(x.vector[i].mainSize);
        x_Eval2[i] = x_Eval2[i].resize(x.vector[i].mainSize);
        x_Eval1Sum += abs((x_Eval1[i].tailSum()).right()) + abs((x_Eval1[i].tailSum()).left());
        x_Eval2Sum += abs((x_Eval2[i].tailSum()).right()) + abs((x_Eval2[i].tailSum()).left());
    }
    if( (x_Eval1Sum < x_Eval2Sum) || ( (abs(x.vector[0].C)).right() >= 20)  ){
        x.vector = x_Eval1;     
    }
    else{
        x.vector = x_Eval2;
    }
    IVector inclusion = vectorField.computeInclusion(x_Extend); 
    //cout <<"incusion\n:" <<inclusion <<"\n";
    for(int i = 0;i<vectorField.nMain; i++){
        vectorField.f.setParameter(i+vectorField.numOfNormalParams,inclusion[i].mid());
        perturbMap.setParameter(i,inclusion[i] - inclusion[i].mid());
    }
    solver.setStep(dt);
    x.mainModes.move( solver );
    x.intersectRepresetations(vectorField.indexer);
}

void  Mover::perturb(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, 
                        capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams) 
{
    for(int k=0;k<dimIn;k++)
    {
        out[k] = params[k];
    }   
}

Section::Section(IAffineSection sectionMainModes)
:sectionMainModes(sectionMainModes) 
{
    
}

interval Section::evalAt(Set& set) 
{
    return set.mainModes.evalAt(sectionMainModes);
}
PoincareMapDiss::PoincareMapDiss(Encloser& encloser,Mover& mover ,Section section) 
:encloser(encloser), mover(mover), section(section)
{
    
}

std::pair<IVector, SeriesVector> PoincareMapDiss::compute(Set& set,IVector x0 ,IMatrix A,VectorField& vectorField) 
{
    interval moverStep = interval(1./1024);
    IVector initSet = (IVector) set.mainModes  ;
    Set setNext = set;
    mover.setStep(moverStep);
    interval prev= interval(0);
    bool IsMinPeriod = true;
    bool isFirst = true;
    encloser.enclose(set,vectorField,moverStep,4,true,true);
    SeriesVector x_Extend= set.vector + encloser.enclosureExtent; 
    IVector extentedMainModes = vectorField.indexer.getIVector(x_Extend); 
    IVector enclosureExtentMainModes = vectorField.indexer.getIVector(encloser.enclosurePointWiseExtent);
    IVector fEstimation  = vectorField.f(extentedMainModes);
    IVector inclusion = vectorField.computeInclusion(x_Extend); 
     if(subset(interval(0),section.sectionMainModes.getNormalVector()*(fEstimation+inclusion) )){
        IsMinPeriod = false;
    }
    while(1)
    {
        

        mover.move(setNext,vectorField,true);
        if(section.evalAt(set)<interval(0) && ( section.evalAt(setNext)>interval(0) || subset(interval(0),section.evalAt(setNext)) ) ){
                while(subset(interval(0),section.evalAt(setNext))){
                    mover.move(setNext,vectorField,true);
                }

                if(section.evalAt(setNext)>interval(0)){
                    break;
                }
        }

        else{
            encloser.enclose(set,vectorField,moverStep,4,true,true);
            enclosureExtentMainModes = vectorField.indexer.getIVector(encloser.enclosurePointWiseExtent);
            if(subset(interval(0),
            section.evalAt(set) + section.sectionMainModes.getNormalVector()*enclosureExtentMainModes) )
            {
                if(!(intersectionIsEmpty(initSet, (IVector) set.mainModes + enclosureExtentMainModes  )) )
                {
                    if(!isFirst){
                        IsMinPeriod = false;
                    }
                    
                }
            }
           // cout << IsMinPeriod<< " ";
            isFirst = false;
            set = setNext;
        }
    }


    interval returnTimeLeft = set.getCurrentTime();
    interval returnTimeRight = setNext.getCurrentTime();
    interval tempTimeRight = setNext.getCurrentTime();
    interval dt = setNext.getCurrentTime() - set.getCurrentTime();

        encloser.enclose(set,vectorField,dt,4,true,true);
     x_Extend= set.vector + encloser.enclosureExtent; 
    extentedMainModes = vectorField.indexer.getIVector(x_Extend); 
    enclosureExtentMainModes = vectorField.indexer.getIVector(encloser.enclosurePointWiseExtent);
    fEstimation  = vectorField.f(extentedMainModes);
    inclusion = vectorField.computeInclusion(x_Extend); 
     if(subset(interval(0),section.sectionMainModes.getNormalVector()*(fEstimation+inclusion) )){
        IsMinPeriod = false;
    }
    cout <<"Is the period minimal: " <<IsMinPeriod <<"\n";
    setNext = set;
    int counter = 0;
    while(1){
        mover.setStep(dt/2);
        mover.move(setNext,vectorField,true);
    if(section.evalAt(setNext)<interval(0)){
        set = setNext;
        returnTimeLeft = set.getCurrentTime();
        dt = dt * (1./2);
    }
    else if (section.evalAt(setNext)>interval(0)){
        tempTimeRight = setNext.getCurrentTime();
        setNext = set;
        dt = dt * (1./2);
    }
    else {
            counter = counter + 1;
            //tempTimeRight = setNext.getCurrentTime();
            if(counter >= 10)
                break;
            setNext = set;
            dt = dt * (1./2);
            
            
        }

    } 
    setNext = set;
    Set setNext2 = setNext;
    Set setTest = setNext;
    dt = 2*(tempTimeRight -returnTimeLeft).right() ;
    counter= 0 ;
    bool wasPositive = false;
    while(counter != 49){
        setNext = set;
        mover.setStep(dt);
        mover.move(setNext,vectorField,true);       
        if(section.evalAt(setNext) > interval(0) ){
            wasPositive = true;
            returnTimeRight = setNext.getCurrentTime();
            setTest = setNext;
            dt = dt * (1./2);
        }
        else if(interval(0).subset(section.evalAt(setNext)) && wasPositive ){
            counter = counter + 1;
           // dt = (setNext.getCurrentTime() - set.getCurrentTime()).right();
            dt = (1.25) * dt;
        } 
        else  {
            //dt = (setNext.getCurrentTime() - set.getCurrentTime()).right();
            dt = 16 * dt;
        } 
    }
    //cout << "Value left " <<section.evalAt(set) <<"\n"; 
    //cout << "Value right" <<section.evalAt(setTest) <<"\n"; 
    interval step = (returnTimeRight - returnTimeLeft).right();
    reachSectionTime = returnTimeLeft + step * interval(0,1);
    setNext = set;
    mover.setStep(step);
    mover.move(setNext,vectorField);
    encloser.enclose(set,vectorField,step,4,true,true);
     x_Extend= set.vector + encloser.enclosureExtent; 
     inclusion = vectorField.computeInclusion(x_Extend); 
    //cout<<"-----------------------------------\n";
    //encloser.enclosureExtent.print();
    //encloser.enclosurePointWiseExtent.print();
    //cout <<"enclosure subset PointWiseEnclosure"  <<encloser.enclosureExtent.subset(encloser.enclosurePointWiseExtent)<<"\n";
    //cout <<"PointWiseEnclosure subset enclosure"  <<encloser.enclosurePointWiseExtent.subset(encloser.enclosureExtent)<<"\n";
     extentedMainModes = vectorField.indexer.getIVector(x_Extend);
     enclosureExtentMainModes = vectorField.indexer.getIVector(encloser.enclosurePointWiseExtent);
    IMatrix DF = vectorField.f[extentedMainModes];
    IMatrix D = A*DF;        
    IVector v = set.mainModes.affineTransformation(A,x0);
    IVector rest = A*vectorField.f(set.mainModes.get_x()) + (D*set.mainModes.get_C())*set.mainModes.get_r0() +
                                (D*set.mainModes.get_B())*set.mainModes.get_r() + D*enclosureExtentMainModes + A*inclusion;
    IVector resultMain = v  +  (interval(0,1)* step) * rest;
    cout << "crossing time"<< step<<"\n";
     fEstimation = vectorField.f(set.mainModes.get_x()) + (DF*set.mainModes.get_C())*set.mainModes.get_r0() +
                                (DF*set.mainModes.get_B())*set.mainModes.get_r() + DF*enclosureExtentMainModes + inclusion;
    if(subset(interval(0),section.sectionMainModes.getNormalVector()*fEstimation )){
        cout << "Possible Non Transversal Crossing!!! "<<"\n";
    }
    else{
        cout << "Transversal Crossing"<<"\n";
    }
    cout <<"v" << v <<"\n";
    cout << "rest:" <<(interval(0,1)* step) * rest <<"\n";
    //cout <<"D*extentMainModes" << D*enclosureExtentMainModes  <<"\n";
    //cout <<"A*inclusion" << A*inclusion  <<"\n";
    return std::pair<IVector, SeriesVector> (resultMain, x_Extend);
    
}
