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



/*
Struktura, która zawiera wpolrzedne sekcji w postaci intervalowej
dodatkowo macierz invCordinates zawiera odwrocone wpolrzedne
tej sekcji potrzebne do przeliczania wartosci obrazu
*/
struct RigSectionCordinates
{
	IVector orgin;
	IVector normal;
	IMatrix cordinates;
	IMatrix invCordinates;

 	/*
	setPointer - wskaznik na zbior C0HOTripletonSet
	solver - scisly inegrator
	metoda zwraca wartosc obraz zbioru na sekcji w
	wpolrzednych sekcji
	*/
	std::pair<IVector, SeriesVector> transportToSection(Set set,VectorField& vectorField,Encloser& encloser,Mover& mover){
		Section section(IAffineSection(this->orgin, this->normal) );
		PoincareMapDiss pm(encloser,mover ,section);
		return pm.compute(set, this->orgin, this->invCordinates, vectorField);
	}
	std::pair<IVector, SeriesVector> transportToSection(Set set,VectorField& vectorField,Encloser& encloser,Mover& mover, interval& reachTime){
		Section section(IAffineSection(this->orgin, this->normal) );
		PoincareMapDiss pm(encloser,mover ,section);
		auto result = pm.compute(set, this->orgin, this->invCordinates, vectorField);
		reachTime = pm.reachSectionTime;
		return result;
	}
	/*
	r0 - wektor wpolrzednych na sekcji
	metoda zwraca wskaznik na set C0HOTripletonSet
	reprezentujacy ten wektor
	*/
	Set getSet(IVector r0,SeriesVector r1,Indexer indexer,bool intersect = true){
		r0[0] = interval(0);
		InclRect2Set mainModes = InclRect2Set(this->orgin, this->cordinates,r0);
		Set set(r1,mainModes);
		if(intersect)
        	set.intersectRepresetations(indexer);
		else 
			set.makeCosistend(indexer);
		return set;
	}

};

typedef std::vector<RigSectionCordinates> RigSectionCordinatesArray;
RigSectionCordinates setRigCordinates(Section_Cordinates section_Cordinates)
{
	IVector orgin = (IVector) section_Cordinates.orgin;
	IVector normal = (IVector) section_Cordinates.normal;
	IMatrix cordinates = (IMatrix) section_Cordinates.cordinates;
	int dim = orgin.dimension();
	for(int i=1;i<dim;i++)
	{
		interval w(1/ (normal* normal));
		cordinates.column(i) =  cordinates.column(i) -
														normal * (normal *cordinates.column(i)) * w ;
	}
	IMatrix invCordinates = matrixAlgorithms::krawczykInverse(cordinates);
	RigSectionCordinates rigSectionCordinates;
	rigSectionCordinates.orgin = orgin;
	rigSectionCordinates.normal = normal;
	rigSectionCordinates.cordinates = cordinates;
	rigSectionCordinates.invCordinates = invCordinates;
	return rigSectionCordinates;
}


RigSectionCordinatesArray setRigCordinatesArray(SectionCordinatesArray tabSections){
	RigSectionCordinatesArray tabRigSections(tabSections.size());
	for(int i=0;i<tabSections.size();i++){
		tabRigSections[i] = setRigCordinates(tabSections[i]);

	}
	return tabRigSections;


}

bool cheakInclusion(IVector x0,SeriesVector x1,IVector y0,SeriesVector y1,Indexer indexer){
	bool isInclusionMainModes = true;
	for(int i=1;i<x0.dimension();i++){

		std::cout <<"InitSet["<<i<< "] = "<< y0[i] << '\n';
		std::cout <<"Image[i] = "<< x0[i] << '\n';
		std::cout <<"inclusion on "<< i << " variable: "  << subset(x0[i],y0[i]) << '\n';
		std::cout <<"diam(Image[i])/diam(InitSet[i]) = "  << diam(x0[i])/ diam(y0[i]) << '\n';
		isInclusionMainModes = isInclusionMainModes && subset(x0[i],y0[i]) ; 
	}
	for(int i=0;i<indexer.size();i++){
		int l = indexer.pairs[i].first;
		int r = indexer.pairs[i].second;
		double inf = std::numeric_limits<double>::infinity();
		y1[l].main[r] =  interval((-1)*inf,inf);
	}
	cout << "Inclusion on rest of Series: " << x1.subset(y1);
	return x1.subset(y1) &&isInclusionMainModes;

}
std::pair<IVector, SeriesVector > TestSubsections (
RigSectionCordinatesArray tabRigSection,VectorField& vectorField,IVector r0, SeriesVector r1,interval& period,bool& isValidated){
	auto start_time = std::chrono::high_resolution_clock::now();
	int dimMain = tabRigSection[0].orgin.dimension();
	//int dimDiss = r1.dimension();
	period = interval(0);
	int numOfSections = tabRigSection.size();
	Encloser encloser;
    Mover mover(vectorField,encloser);
    std::pair<IVector, SeriesVector > result;
    Set set = tabRigSection[0].getSet(r0,r1,vectorField.indexer,false);
	SeriesVector initSeriesVector = set.vector;
	//cout << (IVector) set.mainModes<<"\n";
	//set.vector.print();
    int k;
    for(int i=1;i<=numOfSections;i++){
		interval reachTime = interval(0);
		std::cout <<"SectionsNumber: " <<i << '\n';
		k = i%numOfSections;
		result = tabRigSection[k].transportToSection(set,vectorField,encloser,mover,reachTime);
		period = reachTime+ period;
		set = tabRigSection[k].getSet(result.first,result.second,vectorField.indexer);

		for(int k=0;k<dimMain;k++){
			std::cout <<"Image[i] = "<< result.first[k] << '\n';
			}

	}
	//set.vector.print();
	cout<<endl <<"Period" <<period <<endl;
	isValidated = cheakInclusion(result.first,result.second,r0,initSeriesVector,vectorField.indexer);
	auto end_time = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "Czas dowodu: " << duration.count() << " mikrosekund." << std::endl;
	return result;
}


std::pair<IVector, SeriesVector> getInitialSetPoincare(InOut& inOut,int u_v_main_size,int seriesesLen){
	interval x;
	IVector r0(2*u_v_main_size);
	r0[0] = interval(0);
	IVector fileIntervals;
	inOut.initialValuePoincareFile >> fileIntervals;
	int numIntervals = fileIntervals.dimension();
	for(int i=0;i<2*u_v_main_size-1;i++){
		if(i<fileIntervals.dimension())
			r0[i+1] = fileIntervals[i];
		else
			r0[i+1] = fileIntervals[numIntervals - 1];
	}
	inOut.initialValuePoincareFile >> fileIntervals;
	interval C = fileIntervals[0];
	interval s = fileIntervals[1];
	SeriesVector r1(2, C ,s ,SeriesType::sin_odd );
	r1[0] = r1[0].resize(seriesesLen);
	r1[1] = r1[1].resize(seriesesLen);
	return std::make_pair(r0,r1);
}

void EstimateNorm(RigSectionCordinatesArray tabRigSection,VectorField& vectorField,IVector r0, SeriesVector r1,interval time){
	Encloser encloser;
    Mover mover(vectorField,encloser);
	Set set = tabRigSection[0].getSet(r0,r1,vectorField.indexer,false);
	interval maxL2uNorm = interval(0);
	interval maxL2uxNorm = interval(0);
	interval maxL2uxxNorm = interval(0);
	interval maxL2uxxxNorm = interval(0);
	interval maxL2uxxxxNorm = interval(0);
	interval maxL2vNorm = interval(0);
	interval maxL2vxNorm = interval(0);
	interval maxL2vxxNorm = interval(0);
	interval maxL2vxxxNorm = interval(0);
	interval maxL2vxxxxNorm = interval(0);
	interval maxDiamL2u = interval(0);
	interval maxDiamL2v = interval(0);
	interval maxDiamL2ux = interval(0);
	interval maxDiamL2vx = interval(0);
	interval uC = interval(0);
	interval vC = interval(0);

	cout <<"Size2" << mover.encloser.enclosureExtent.vec.size();

	
	while(time.right()>= (set.getCurrentTime()).right() ){
		mover.move(set,vectorField,true);
		
		auto stepEnclosure = set.vector + mover.encloser.enclosureExtent;
		maxL2uNorm = max(maxL2uNorm ,stepEnclosure[0].seminormHs(0));
		maxL2vNorm = max(maxL2vNorm ,stepEnclosure[1].seminormHs(0));
		maxL2uxNorm = max(maxL2uxNorm ,stepEnclosure[0].seminormHs(1));
		maxL2vxNorm = max(maxL2vxNorm ,stepEnclosure[1].seminormHs(1));
		maxL2uxxNorm = max(maxL2uxxNorm ,stepEnclosure[0].seminormHs(2));
		maxL2vxxNorm = max(maxL2vxxNorm ,stepEnclosure[1].seminormHs(2));

		
		maxDiamL2u =  max(maxDiamL2u ,diam(stepEnclosure[0].seminormHs(0)) );
		maxDiamL2v = max(maxDiamL2v ,diam(stepEnclosure[1].seminormHs(0)) );
		maxDiamL2ux = max(maxDiamL2ux ,diam(stepEnclosure[0].seminormHs(1)) );
		maxDiamL2vx = max(maxDiamL2vx ,diam(stepEnclosure[1].seminormHs(1)) );
		Series uS = stepEnclosure[0].refineTail(interval(4));
		Series vS = stepEnclosure[1].refineTail(interval(4));
		uC = max(uS.C,uC);
		vC = max(vS.C,vC);
		
	}
	cout << "maxL2uNorm" << maxL2uNorm <<"\n"; 
	cout << "maxL2uxNorm" << maxL2uxNorm <<"\n"; 
	cout << "maxL2uxxNorm" << maxL2uxxNorm <<"\n";
	//cout << "maxL2uxxxNorm" << maxL2uxxxNorm <<"\n";  
	//cout << "maxL2uxxxxNorm" << maxL2uxxxxNorm <<"\n";  


	cout << "maxL2vNorm" << maxL2vNorm <<"\n"; 
	cout << "maxL2vxNorm" << maxL2vxNorm <<"\n"; 
	cout << "maxL2vxxNorm" << maxL2vxxNorm <<"\n"; 
	//cout << "maxL2vxxxNorm" << maxL2vxxxNorm <<"\n";
	//cout << "maxL2vxxxxNorm" << maxL2vxxxxNorm <<"\n";


	cout << "maxDiamL2u" << maxDiamL2u <<"\n";
	cout << "maxDiamL2v" << maxDiamL2v <<"\n";

	cout << "maxDiamL2ux" << maxDiamL2ux <<"\n";
	cout << "maxDiamL2vx" << maxDiamL2vx <<"\n";  
	cout << "uC" << uC <<"\n"; 
	cout << "vC" << vC <<"\n"; 


}
int main (){
	InOut inOut;
	DVector  paramsDVector;
	inOut.paramsFile >> paramsDVector;
    ParamsMap params = {{"d1",paramsDVector[0]},{"d2",paramsDVector[1]},{"B",paramsDVector[2]},{"A",paramsDVector[3]}};
	RigSectionCordinatesArray tabRigSections = setRigCordinatesArray(readSectionsFromFile(inOut.outPutSubSectionsFile));
	
	int u_vMainSize =tabRigSections[0].orgin.dimension()/2;
	int seriesesLen = 30;
	cout << "serieses Finite Lenght: " << seriesesLen << "\n";
	cout << "u v main sizes: " << u_vMainSize <<"\n";
	auto inital = getInitialSetPoincare(inOut,u_vMainSize,seriesesLen);

	VectorField brusselatorVectorFieldPDESym = getBrusselatorVectorField(u_vMainSize,seriesesLen,params);
	interval period = interval(0.1);
	bool validatedPeriodicOrbit = false;
	int iter = 1;
	while(iter>0 && validatedPeriodicOrbit==false ){
		auto result = TestSubsections(tabRigSections,brusselatorVectorFieldPDESym,inital.first,inital.second,period,validatedPeriodicOrbit);
		inital.first = intersection(inital.first,interval(1.5)*result.first) ;
		inital.second = interval(4)*result.second ;
		for(int i=0;i<inital.second.vec.size();i++){
			inital.second[i] = inital.second[i].refineTail(inital.second[i].s-1);
		}
		iter = iter - 1;
	}
	

}