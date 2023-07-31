#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <iostream>
#include "capd/capdlib.h"
#include <cmath>
#include <vector>

enum class SeriesType {sin, cos, sin_odd, cos_even};
struct Series
{
    SeriesType type;
    capd::IVector main;
    int mainSize;
    int n;
    capd::interval C;
    capd::interval s;
    int static get_n(SeriesType type,int mainSize);
    int static getSeriesIndex(SeriesType type,int mainIndex);
    Series();
    Series(capd::IVector main,capd::interval C,capd::interval s,SeriesType type);
    Series(capd::interval C,capd::interval s,SeriesType type);
    void print() const;
    capd::interval valueAt(int i) const;//działa bez zera
    Series resize(int newSize);//działa bez zera
    Series upperBound();//działa bez zera
    Series lowerBound();//działa bez zera
    capd::interval tailSum();
    capd::interval seminormHs(capd::interval q);
    Series refineTail(capd::interval sNew);//działa bez zera ale zwraca szereg z zerem
    bool subset(const Series& x);//nie działa bez zera
    bool subsetInterior(const Series& x);//nie działa bez zera
    Series xx();//działa bez zera
    Series elementWiseInverse();//działa bez zera
    Series operator-();// działa bez zera
    friend Series operator*(const capd::interval&  a,const Series& x);//działa bez zera
    friend Series operator*(const Series& x,const capd::interval&  a);//działa bez zera
    friend Series operator*(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series operator+(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series operator-(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series intersection(const Series& x,const Series& y);//działa dla szeregów z takim samym decayej
    friend Series semiIntersection(const Series& x,const Series& y);
    friend Series squere(const Series& x);//działa bez zera ale zwraca szereg z zerem
    friend Series mult(const Series& x,const Series& y);//działa bez zera ale zwraca szereg z zerem
    friend Series elementWiseMult(const Series &x,const Series &y);//działa bez zera
    friend Series expMinusOne(const Series &x);
    friend Series exp(const Series &x,const capd::interval& sNew);
    //friend Series exp(const Series &x,capd::interval sNew);
    
};
struct SeriesVector{
    std::vector<Series> vec;
    SeriesVector(int size,capd::interval C,capd::interval s,SeriesType type);
    SeriesVector(int size);
    SeriesVector();
    void print();
    bool subset(const SeriesVector& x);
    bool subsetInterior(const SeriesVector& x);
    capd::IVector tailSum();
    Series& operator[](int i);
    SeriesVector elementWiseInverse();
    SeriesVector upperBound();
    SeriesVector lowerBound();
    friend SeriesVector operator+(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector operator*(const capd::interval&  a,const SeriesVector& x);
    friend SeriesVector operator*(const SeriesVector& x,const capd::interval&  a);
    friend SeriesVector elementWiseMult(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector expMinusOne(const SeriesVector& x);
    friend SeriesVector semiIntersection(const SeriesVector& x,const SeriesVector& y);
    friend SeriesVector exp(const SeriesVector &x,const capd::interval& sNew);
     
     
};    




#endif // __ALGEBRA_H__