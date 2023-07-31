# include <iostream>
# include "capd/capdlib.h"
# include "BrusselatorTools.h"
# include "SetInOut.h"
# include "DSubSections.h"
using namespace capd ;
using namespace std ;
using capd::autodiff::Node;
#include <vector>
#include <cmath>

// Struktura do trzymania sekcji i jej wpółrzędnych w doublach
/*struct Section_Cordinates
{
	DVector orgin;
	DVector normal;
	DMatrix cordinates;
};
*/
/*funkcja szuka w macierzy monodromyMatrix
lewego wektora wlasnego odpowiadajacego 1
*/
DVector getLeftEigenVector(DMatrix monodromyMatrix){
	int dim = monodromyMatrix.column(0).dimension();
	DVector leftValue_2(dim), leftIValue_2(dim);
	DMatrix leftVec_2(dim,dim), leftIVec_2(dim,dim);
	capd::alglib::computeEigenvaluesAndEigenvectors(transpose(monodromyMatrix),leftValue_2,leftIValue_2,leftVec_2,leftIVec_2);
	std::cout << leftValue_2 << '\n';
	std::cout << leftIValue_2 << '\n';
	double eps = 1e-12;
	while(eps < 1e-4){
		for(int i = 0;i<dim;i++){
			if( abs(1- leftValue_2[i]) +abs(leftIValue_2[i]) <eps)
			{
				return leftVec_2.column(i)/leftVec_2.column(i).euclNorm();
			}
		}
		eps = eps*10;
	}
	throw std::runtime_error("Cannot find left eigenvector coresponding to 1\n");
}
/*
dmap - odwzorowanie,
initialValue - wartość początkowa do znalezienia punktu stałego odw. poincare,
DSection_1 - ustalona sekcja początkowa postaci x_i = a,
threshold - stała która determinuje jakie wektory własne wezmiemy do
						macierzy wpółrzędnych,
Funkcja zwraca wpółrzedne sekcji Section_Cordinates, która
- jest zcetrowana w punkcie stalym,
- wektor normalny odpowiada lewemu wektorowi wlasnymu
- wpolrzedne na sekcji sa zadane poprzez wektory wlasne
macierzy rozniczki
*/
Section_Cordinates setCordinates(DMap& dmap,DVector initialValue,DCoordinateSection DSection_1,double threshold){
	int dim = initialValue.dimension();
	DOdeSolver Dsolver(dmap,22);
	DPoincareMap pm_1(Dsolver, DSection_1, poincare::MinusPlus);
	double err = 1;
	DVector  u = initialValue;
	//Szukanie punktu stałego
	while(err>1e-12)
	{
		double returnTime = 0;
		DVector v = pm_1(u,returnTime);
		err = (v-u).euclNorm();
		u = v;
		cout << u << endl;
		std::cout << "returnTime: "<<returnTime << '\n';
	}

	DMatrix monodromyMatrix_1(dim,dim);
	double returnTime_1 = 0.;
	pm_1(u,monodromyMatrix_1,returnTime_1);
	cout <<"returnTime_1: " <<returnTime_1<<"\n";
	
	// Biore lewy wektor własny jako wektor normalny do tworzonej sekcji
	DVector normal = getLeftEigenVector(monodromyMatrix_1);
	/*int hej;
	std::cin >> hej;*/
	normal = normal/normal.euclNorm();
	std::cout <<"(normal * f(orgin))/(||normal|| * ||f(orgin)||) = "  <<
								(normal * dmap(u)) /(normal.euclNorm() * dmap(u).euclNorm() )<< '\n';
	std::cout << DSection_1.getNormalVector()	 << '\n';
	std::cout << normal*DSection_1.gradient(u)	 << '\n';
	// Ustawiam wektor normalny tak by kierunki przejscia były takie same
	if(normal*DSection_1.gradient(u) < 0)
		normal = - normal;
	/*Tworze nową sekcje zcentrowaną w punkie stałym u i z normalną normal
	 i dla niej wyliczam macierz różniczki w punkie stałym
	*/
	DAffineSection DSection_2(u,normal);
	DPoincareMap pm_2(Dsolver, DSection_2, poincare::MinusPlus);
	DMatrix monodromyMatrix_2(dim,dim);
	double returnTime_2 = 0.;
	DVector P = pm_2(u,monodromyMatrix_2,returnTime_2);
	DMatrix DP = pm_2.computeDP(P,monodromyMatrix_2,returnTime_2);
	DVector rightValue_2(dim), rightIValue_2(dim);
	DMatrix rightVec_2(dim,dim), rightIVec_2(dim,dim);
	DMaxNorm maxNorm;
	cout <<"maxNorm of DP: " <<maxNorm(DP)<<"\n";
	if(abs(returnTime_1 - returnTime_2)>0.001 ){
		cout << "returnTime_1"<< returnTime_1<<"\n ";
		cout << "returnTime_2" << returnTime_1<<"\n ";
		cout <<"WARNING TIMES FOR SECTIONS DONT MATCH" <<"\n";

	}
	cout <<"maxNorm(DP) = "<<maxNorm(DP)<<"\n";
	
	//Licze lewe i prawe wektory własne dla wyliczonej różniczki
	capd::alglib::computeEigenvaluesAndEigenvectors(DP,rightValue_2,rightIValue_2,rightVec_2,rightIVec_2);
	DVector leftValue_2(dim), leftIValue_2(dim);
	DMatrix leftVec_2(dim,dim), leftIVec_2(dim,dim);
	capd::alglib::computeEigenvaluesAndEigenvectors(transpose(DP),leftValue_2,leftIValue_2,leftVec_2,leftIVec_2);
	DMatrix CordinateMatrix = DMatrix::Identity(dim) ;
	cout <<"orgin:\n" << u<<"\n";
	cout <<"normal:\n" << normal  << "\n";
	std::cout <<"eigenvalues of DP\n" <<rightValue_2 << '\n';
	std::cout <<"imaginaryPart eigenvalues of DP\n" <<rightIValue_2 << '\n';
	//Pierwsza wpółrzędna to wartość pola w punkie stałym/przesunieca sekcji
	CordinateMatrix.column(0) = dmap(u)/ (dmap(u).euclNorm());
	int k = 0;
	double imaginaryEpsilon = 1e-11;
	DMatrix LeftVectorsMatrix = DMatrix::Identity(dim) ;
	/*Pierwszą kolumne dla macierzy LeftVectorsMatrix uzupełniam
	kierunkiem normalnym by dalsze wektory utworzone z ortogonalizacji
	tej macierzy były do niego ortogonalne
	*/
	LeftVectorsMatrix.column(0) = normal;
	/*
		Jeśli czesc rzeczywita wartość własna macierzy rozniczki
	  jest wieksza od threshold to daje ten wektor do
		macierzy wpółrzędnych,
		dodatkowo do macierzy LeftVectorsMatrix dodaje
		lewe wektory własne odpowiadające tym wartoscia wlasnym
	*/
	while(std::abs(rightValue_2[k]) +  std::abs(rightIValue_2[k]) > threshold && k < dim)
	{
		//dwie procedury w zaleznosci czy wartość w. jest rzeczywita czy zepolona
		if(std::abs(rightIValue_2[k])  < imaginaryEpsilon){
			CordinateMatrix.column(k+1) = rightVec_2.column(k) ;
			LeftVectorsMatrix.column(k+1) = leftVec_2.column(k);
			k=k+1;
		}
		else{
			CordinateMatrix.column(k+1) = rightVec_2.column(k) + rightIVec_2.column(k);
			CordinateMatrix.column(k+2) = rightVec_2.column(k) - rightIVec_2.column(k);
			LeftVectorsMatrix.column(k+1) = leftVec_2.column(k) + leftIVec_2.column(k);
			LeftVectorsMatrix.column(k+2) = leftVec_2.column(k) - leftIVec_2.column(k);
			k = k + 2;
		}
	}
	std::cout << "number of taken eigenvectors:"<< k << '\n';
	DMatrix Q,R;
	/*Ortogonalizuje macierz LeftVectorsMatrix
	Na dalszych kierunkach kierunkach ustawiam
	kierunki na sekcji, ktore sa ortogonalne
	do lewych wektorow wlasnych i do normalnej
	na sekcji. Mozna zobaczyc ze dla rozniczki stanowia
	one niezmienicza podprzestrzen
	*/
	matrixAlgorithms::QR_decompose(LeftVectorsMatrix,Q,R);

	k = k + 1;
	while(k < dim)
	{
		CordinateMatrix.column(k) = Q.column(k);
		k=k+1;
	}
	//normowanie
	for(int i = 0;i<dim;i++)
	{
		CordinateMatrix.column(i) = CordinateMatrix.column(i)/CordinateMatrix.column(i).euclNorm();
	}

	Section_Cordinates section_Cordinates;
	section_Cordinates.orgin = u;
	section_Cordinates.normal = normal;
	section_Cordinates.cordinates = CordinateMatrix;
	return section_Cordinates;
}





//Struktura do trzymania wielu sekcji
//typedef std::vector<Section_Cordinates> SectionCordinatesArray;
/*
CordinateMatrix - referencja do macierzy wpolrzednych
DP - macierz rozniczki w punkcie stalym na sekcji
normal - wektor normalny do sekcji
threshold - stała która determinuje jakie lewe wektory własne wezmiemy do
						liczenia podprzestrzeni niezmieniczej,

CordinateMatrix zostaje uzupelniona kierunkami ortogonalnymi
do lewych wektorow wlasnych macierzy DP i do normalnej normal
*/
void FillRestOfEigenVectors(DMatrix& CordinateMatrix,DMatrix DP, DVector normal,double threshold){
	int dim = DP.column(0).dimension();
	DVector leftValue_2(dim), leftIValue_2(dim);
	DMatrix leftVec_2(dim,dim), leftIVec_2(dim,dim);
	capd::alglib::computeEigenvaluesAndEigenvectors(transpose(DP),leftValue_2,leftIValue_2,leftVec_2,leftIVec_2);
	DMatrix LeftVectorsMatrix = DMatrix::Identity(dim) ;
	LeftVectorsMatrix.column(0) = normal/(normal.euclNorm());
	double imaginaryEpsilon = 1e-11;
	int k=0;

	while(std::abs(leftValue_2[k]) + std::abs(leftIValue_2[k])> threshold && k < dim)
	{
		if(std::abs(leftIValue_2[k])  < imaginaryEpsilon){
			LeftVectorsMatrix.column(k+1) = leftVec_2.column(k);
			k=k+1;
		}
		else{
			LeftVectorsMatrix.column(k+1) = leftVec_2.column(k) + leftIVec_2.column(k);
			LeftVectorsMatrix.column(k+2) = leftVec_2.column(k) - leftIVec_2.column(k);
			k = k + 2;
		}
	}
	DMatrix Q,R;
	matrixAlgorithms::QR_decompose(LeftVectorsMatrix,Q,R);
	k = k + 1;
	while(k < dim)
	{
		CordinateMatrix.column(k) = Q.column(k);
		k=k+1;
	}
	for(int i = 0;i<dim;i++)
	{
		CordinateMatrix.column(i) = CordinateMatrix.column(i)/CordinateMatrix.column(i).euclNorm();
	}
}

/*
dmap - odwzorowanie,
section - sekcja początkowa policzona wczesniej za pomoca setCordinates
n - liczba sekcji, które chcemy utworzyc
threshold - stała która determinuje jakie wektory własne wezmiemy do
						macierzy wpółrzędnych
Funkcja zwraca tablice (wektor) sekcji posredniej
*/
SectionCordinatesArray setSubSections(DMap& dmap,Section_Cordinates section,int n,double threshold){
	int dim = section.orgin.dimension();
	SectionCordinatesArray tabSections(n);
	//jako pierwsza sekcje ustawiamy ta policzona wczesniej
	tabSections[0] = section;
	DOdeSolver Dsolver(dmap,22);
	DAffineSection DSection0(section.orgin,section.normal);
	DPoincareMap pm_0(Dsolver, DSection0, poincare::MinusPlus);
	double period;
	pm_0(section.orgin,period);
	//Licze okres orbity
	std::cout << "/* message */"<<period << '\n';
	DTimeMap timeMap(Dsolver);
	for(int i=1;i<n;i++){
		/*
		Jako początki ukladu wpolrzednych biore wartość
		potoku w punktach posrednich z orbity
		*/

		std::cout << "Section:"<< i << '\n';
		std::cout << "Move step"<<(1.0)/n * period << '\n';
		tabSections[i].orgin = timeMap((1.0)/n * period,section.orgin);
		std::cout << "/* message ORGIN */"<<tabSections[i].orgin << '\n';
		std::cout << "VectorFieldValue:"<< dmap(tabSections[i].orgin) << '\n';
		std::cout << "VectorFieldValueNorm:"<< dmap(tabSections[i].orgin).euclNorm() << '\n';
		DPoincareMap pm(Dsolver, DSection0, poincare::MinusPlus);
		DMatrix monodromyMatrix_1(dim,dim);
		pm(tabSections[i].orgin,monodromyMatrix_1);
		/*
		Transportuje nowy srodek do sekcji wybranej wczesniej
		Transpozycja macierzy monodromi daje mi odwzorowanie dualne
		pomiedzy przeszczenia dualna starej sekcji i nowej sekcji
		i transportuje nią wektor normalny z pierwszej sekcji
		*/
		DVector v =  dmap(tabSections[i].orgin)/dmap(tabSections[i].orgin).euclNorm();
		tabSections[i].normal = transpose(monodromyMatrix_1)*tabSections[0].normal;
		tabSections[i].normal = tabSections[i].normal/tabSections[i].normal.euclNorm();
		/*if(std::abs(v*tabSections[i].normal)<0.3){
			tabSections[i].normal = v;
		}
		*/
		if(v*tabSections[i].normal<0){
			tabSections[i].normal = -tabSections[i].normal;
		}
		//tabSections[i].normal = -tabSections[i].normal;
		std::cout <<"normal * f(orgin)" <<v*tabSections[i].normal << '\n';

		/*
		Tworze odwzorowanie poincare do nowej sekcji
		*/
		DAffineSection DSection2(tabSections[i].orgin,tabSections[i].normal);


		DPoincareMap pm2(Dsolver, DSection2, poincare::MinusPlus);
		DMatrix monodromyMatrix_2(dim,dim);
		/*
		Licze macierz rozniczki odwzorowania poincare z
		pierwszej sekcji do teraz utworzonej sekcji i
		transportuje nia wpolrzedne ze starej sekcji
		do nowo utworzonej sekjci
		*/
		pm2(tabSections[0].orgin,monodromyMatrix_2);
		DMatrix DP = pm2.computeDP(tabSections[i].orgin,monodromyMatrix_2);
		tabSections[i].cordinates = DP*tabSections[0].cordinates;

		/*Pierwsza kolumna macierzy wpolrzednych w nowo utworzonej sekcji
		jest kierunkiem pola w jej srodku
		*/
		tabSections[i].cordinates.column(0) = v;

		/*
		Uzupelniam czesc wektorow w macierzy wpolrzednych korzystajac
		z lewych wektorow wlasnych rozniczki odwzorwania z nowej
		sekcji w samą w siebie
		*/
		DMatrix monodromyMatrix_3(dim,dim);
		pm2(tabSections[i].orgin,monodromyMatrix_3);
		DP = pm2.computeDP(tabSections[i].orgin,monodromyMatrix_3);
		FillRestOfEigenVectors(tabSections[i].cordinates,DP,tabSections[i].normal,threshold);

		for(int k=0;k<dim;k++){
			const DVector& v1 = tabSections[i].cordinates.column(k);
			tabSections[i].cordinates.column(k) = v1/v1.euclNorm();
		}

	}
	return tabSections;

}

/*
Struktura, która zawiera wpolrzedne sekcji w postaci intervalowej
dodatkowo macierz invCordinates zawiera odwrocone wpolrzedne
tej sekcji potrzebne do przeliczania wartosci obrazu
*/
/*void SaveSectionsToFile(InOut& inOut,SectionCordinatesArray sectionArray){
	int n = sectionArray.size();
	cout.precision(18);
	inOut.outPutSubSectionsFile << sectionArray.size()<<endl;
	for(int i=0;i<n;i++){
		inOut.outPutSubSectionsFile << sectionArray[i].orgin << '\n';
		inOut.outPutSubSectionsFile << sectionArray[i].normal << '\n';
		inOut.outPutSubSectionsFile << sectionArray[i].cordinates << '\n';
	}

}*/
int main ()
{
	InOut inOut;
	DVector  Params;
	inOut.paramsFile >> Params;
	DVector  initialValue;
	inOut.initialValueFile >> initialValue;
	int dim = initialValue.dimension();
	int dim_u=dim/2;
	int noParam=Params.dimension();
	std::cout << "dimension = "<< dim << '\n';
	cout << "d1: " << Params[0] <<'\n';
	cout << "d2: " << Params[1] <<'\n';
	cout << "B:  " << Params[2] <<'\n';
	cout << "A:  " << Params[3] <<'\n';
	std::cout << initialValue << '\n';
	DMap Brusselator(BrusselatorVectorFieldSym(dim_u));
	int cordinate ;
	double valueAtCordinate;
	int numOfSections;
	double threshold;
	inOut.sectionFile >> cordinate;
	inOut.sectionFile >> valueAtCordinate;
  inOut.sectionFile >> numOfSections;
	inOut.sectionFile >> threshold;
	DCoordinateSection section(dim, cordinate, valueAtCordinate);
	Brusselator.setParameters(Params);

	Section_Cordinates section_Cordinates;
	section_Cordinates = setCordinates(Brusselator,initialValue,section,threshold);
	saveSectionsToFile(inOut.outPutSubSectionsFile,setSubSections(Brusselator,section_Cordinates,numOfSections,threshold));
	return 0 ;
}
