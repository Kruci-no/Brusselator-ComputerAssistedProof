#include "algebra.h"
#include <limits>
using namespace capd;
int  Series::get_n(SeriesType type,int mainSize) 
{
    switch(type){
      case SeriesType::sin  : 
        return mainSize; 
      case SeriesType::cos  : 
        return mainSize - 1; 
      case SeriesType::sin_odd  : 
        return 2*mainSize - 1;   
      case SeriesType::cos_even  : 
        return 2*(mainSize - 1); 
    }
}

int  Series::getSeriesIndex(SeriesType type,int mainIndex) 
{
    switch(type){
      case SeriesType::sin  : 
        return mainIndex + 1; 
      case SeriesType::cos  : 
        return mainIndex ; 
      case SeriesType::sin_odd  : 
        return 2*mainIndex + 1;   
      case SeriesType::cos_even  : 
        return 2*mainIndex; 
    }    
}

Series::Series() 
{   
}

Series::Series(IVector main,capd::interval C,interval s,SeriesType type)
:main(main),C(C),s(s),type(type)
  {
    this->mainSize = this->main.dimension();
    this->n = get_n(type,this->mainSize);
  }
  
Series::Series(capd::interval C,capd::interval s,SeriesType type) 
:C(C),s(s),type(type)
{
  if((type==SeriesType::sin )||(type==SeriesType::sin_odd)){
    this->mainSize = 0;
    this->n = get_n(type,this->mainSize);
  }
  else{
    this->main = IVector(1);
    this->mainSize = this->main.dimension();
    this->n = get_n(type,this->mainSize);
  }    
    
}

void Series::print() const
{
    switch(type){
      case SeriesType::sin  : std::cout << "Sin series\n";   break;
      case SeriesType::cos  : std::cout << "Cos series\n"; break;
      case SeriesType::sin_odd  : std::cout << "Sin series with only odd nonzero terms\n";   break;
      case SeriesType::cos_even  : std::cout << "Cos series with only even nonzero terms\n"; break;
    }
    std::cout <<"Finite Part:\n " <<main;
    std::cout <<"\n Tail:\n " <<C <<"/(k^"<<s<<")\n";
    std::cout <<" n: " <<n ;
    std::cout <<"\n mainSize: " <<mainSize<<std::endl;
}

interval Series::valueAt(int i) const
{
    switch(type){
        case SeriesType::sin : 
            if(i<= n){
            return main[i-1];
            }
            else return C /(interval(i)^s);
            
        case SeriesType::cos: 
            if(i<=n){
            return main[i];
            }
            else return C /(interval(i)^s);

        case SeriesType::sin_odd : 
            if(i%2 == 0) 
                return interval(0);
            if(i<= n){
            return main[ (i-1)/2];
            }
            else return C /(interval(i)^s);
        case SeriesType::cos_even: 
            if(i%2 == 1) 
                return interval(0);
            if(i<= n){
            return main[ i/2];
            }
            else return C /(interval(i)^s);
    }
}

Series Series::resize(int newSize) 
{
  Series z = *(this);
  IVector newMain(newSize);
  z.s = this->s;
  z.n = get_n(this->type,newSize);
  z.mainSize = newSize;
  for(int i=0; i<newSize; i++){
    newMain[i]= valueAt(getSeriesIndex(this->type,i));
  }
  z.main = newMain;
  if(newSize<mainSize){
    interval Cplus;
    interval Cminus;
    for(int i=newSize; i<this->mainSize; i++){
      int k = getSeriesIndex(this->type,i);
      Cplus = max(this->valueAt(k) * (interval(k)^s), z.C);
      Cminus = min(this->valueAt(k) * (interval(k)^s), z.C);
      z.C = intervalHull(Cplus,Cminus);
    }
  }
  return z;
}

Series Series::upperBound() 
{
  Series z = *(this);
  for(int i=0;i<z.mainSize;i++){
    z.main[i] = (this->main[i]).right();
  }
  z.C = (this->C).right();
  return z;
}

Series Series::lowerBound() 
{
  Series z = *(this);
  for(int i=0;i<z.mainSize;i++){
    z.main[i] = (this->main[i]).left()*interval(0,1);
  }
  z.C = (this->C).left();
  return z;
}

capd::interval Series::tailSum() 
{
  if(s<interval(-1)){
    double inf = std::numeric_limits<double>::infinity();
    return interval(0,inf)*this->C;
  }
  return interval(0,1)* this->C * (interval(n)^(-this->s + interval(1)))/(this->s - interval(1));
}
capd::interval Series::seminormHs(capd::interval q){
  interval sum = interval(0);
  for(int i=0;i<this->mainSize;i++){
    int k = Series::getSeriesIndex(this->type,i);
    sum = sum + this->valueAt(k)*(k^q)*this->valueAt(k)*(k^q) ; 
  }
  interval sNew = 2*(this->s-q);

  interval rest = interval(0,1)* this->C * (interval(n)^(-sNew + interval(1)))/(sNew - interval(1));
  return  ( (interval::pi()/(interval(2)) ) *(rest + sum))^(interval(1./2)) ;
}

Series Series::refineTail(interval sNew) 
{
  if(sNew < s ){
    Series z = *(this);
    int next_n = getSeriesIndex(this->type,this->mainSize);
    z.C = interval(0,1) * this->valueAt(next_n) * (interval(next_n)^sNew);
    z.s = sNew;
    return z;
  }
  else{
    throw std::runtime_error("Can't refine Tail\n");
  }
}

bool Series::subset(const Series& x) 
{
  if(! (interval(0).subset(this->C)) && (interval(0).subset(this->C)) ){
    throw std::runtime_error("One of the series does't contain zero\n");
  }
  if(this->type != x.type){
    throw std::runtime_error("Type error\n");
  }
  if(this->s< x.s &&  this->C!= 0)
    return false;
  for(int i=0;i<=max(x.mainSize,this->mainSize);i++){
    int k = getSeriesIndex(x.type,i);
    if(! this->valueAt(k).subset(x.valueAt(k)) ){
      //std::cout << k <<"\n";
      //std::cout << this->valueAt(k)<<"\n";
      //std::cout << x.valueAt(k)<<"\n";
      return false;
    }
  }
  return true;
}

bool Series::subsetInterior(const Series& x) 
{ 
  if(! (interval(0).subset(this->C)) && (interval(0).subset(this->C)) ){
    throw std::runtime_error("One of the series does't contain zero\n");
  }
  if(this->type != x.type){
    throw std::runtime_error("Type error\n");
  }
  if(this->s< x.s &&  this->C!= 0)
      return false;
  for(int i=0;i<=max(x.mainSize,this->mainSize);i++){
      int k = getSeriesIndex(x.type,i);
      if(! (this->valueAt(k)).subsetInterior(x.valueAt(k)) ){
        //std::cout << k <<"\n";
        //std::cout << this->valueAt(k)<<"\n";
        //std::cout << x.valueAt(k)<<"\n";
        return false;
      }
  }
  return true;
}



Series Series::xx() 
{
  Series z = *this;
  for(int i=0;i<this->mainSize;i++){
    int k = Series::getSeriesIndex(this->type,i);
    z.main[i] = (-1)*(interval(k)^interval(2))* z.valueAt(k);
  }
  z.C = -z.C;
  z.s = z.s - interval(2);
  return z;
}

Series Series::elementWiseInverse() 
{
  Series z = *(this);
  int k;
  for(int i=0;i<z.mainSize;i++){
    k = Series::getSeriesIndex(this->type,i);
    try{
      z.main[i] = interval(1)/this->valueAt(k); 
    }
    catch(capd::intervals::IntervalError<double>){
      double inf = std::numeric_limits<double>::infinity();
      z.main[i] = interval((-1)*inf,inf);
    }
  }
  //Dorobić jak C zawiera zero
  z.C = interval(1)/this->C;
  z.s = - this->s;
  return z;
}

Series Series::operator-() 
{
  return interval(-1)*(*this);
}

Series exp(const Series &x,const interval& sNew) 
{
  Series y = x;
  if(x.s >= 0 || x.C.right()>=0){
    std::cout<<x.C<<"\n";
    std::cout<<x.s<<"\n";
    x.print();
    throw std::runtime_error("exp - operation not supported for x Series\n");
  }
  for(int i=0;i<x.mainSize;i++){
    y.main[i] = exp(x.main[i]);
  }
  int next_n = Series::getSeriesIndex(x.type,x.mainSize);
  if(sNew<= 0){
    y.C = interval(0,1)*exp(x.valueAt(next_n))*(interval(next_n)^sNew);
    y.s = sNew;
  }
  else{
    
    interval argMax =((abs(x.C*x.s/sNew)) ^ (interval(1/x.s)));
    y.s = sNew;
    //std::cout<<argMax <<"\n";
    if(argMax < next_n){
      y.C = interval(0,1)*exp(x.valueAt(next_n))*(interval(next_n)^sNew);
    }
    else{
      y.C = interval(0,1)*exp(x.C*((argMax)^(-x.s)) )*(interval(argMax)^sNew);
      //std::cout << y.C<<"\n";
    }

  }
  return y;
}

Series expMinusOne(const Series &x) 
{
  Series y = x;
  for(int i=0;i<x.mainSize;i++){
    y.main[i] = exp(x.main[i])-interval(1);
  }
  int n_next = Series::getSeriesIndex(x.type,x.mainSize);
  if(x.s > 0){
    y.s = interval(0);
    y.C = exp(interval(0,1)*x.valueAt(n_next))-interval(1);
  } else if (x.s == interval(0) ){
    y.s = interval(0);
    y.C = exp(x.C) - interval(1);
  } else{
    y.s = interval(0);
    if(x.C > interval(0))
    {
      throw std::runtime_error("ExpMinusOne error series in not decaying and is positive");
    }
    else{
      y.s = interval(0);
      y.C = interval(-1,0);
    }

  }
  return y;
}



Series elementWiseMult(const Series &x,const Series &y) 
{
  if(x.type!=y.type){
    throw std::runtime_error("Type error\n");
  }
  int newSize = max(x.mainSize,y.mainSize);
  IVector newMain(newSize);
  int k;
  for(int i=0;i<newSize;i++){
    k = Series::getSeriesIndex(x.type,i);
    newMain[i] = x.valueAt(k) * y.valueAt(k); 
  }
  Series z;
  z.type = x.type;
  z.main = newMain;
  z.mainSize = newSize;
  z.n = Series::get_n(x.type,newSize);
  z.C = x.C * y.C;
  z.s = x.s + y.s;
  return z;
}
Series intersection(const Series&x ,const Series&y) 
{
  if(x.type!=y.type)
    throw std::runtime_error("Non maching types");
  if(x.s!=y.s)
    throw std::runtime_error("Series have Different decays");
  Series result;
  result.type = x.type;
  int newSize = std::max(x.mainSize,y.mainSize);
  IVector newMain(newSize);
  for(int i=0;i<newSize;i++){
    int k = Series::getSeriesIndex(result.type,i);
    if(!intersection(x.valueAt(k),y.valueAt(k),newMain[i]))
    {
      std::cout <<"at value "<<k<<x.valueAt(k) <<" "<<y.valueAt(k)<<std::endl;
      x.print();
      y.print();
      throw std::runtime_error("Series intersection is empty!");
      
    }
  }
  result.main = newMain;
  result.mainSize = newSize;
  result.s = x.s;
  if(!intersection(x.C,y.C,result.C))
  {
    throw std::runtime_error("Series  intersection is empty!");
  }
  result.n = Series::get_n(result.type, newSize); 
  return result;  
}
Series semiIntersection(const Series& x,const Series& y) 
{
  if(x.type!=y.type)
    throw std::runtime_error("Non maching types");
  Series result;
  result.type = x.type;
  int newSize = std::max(x.mainSize,y.mainSize);
  IVector newMain(newSize);
  for(int i=0;i<newSize;i++){
    int k = Series::getSeriesIndex(result.type,i);
    if(!intersection(x.valueAt(k),y.valueAt(k),newMain[i]))
    {
      std::cout <<"at value "<<k<<x.valueAt(k) <<" "<<y.valueAt(k)<<std::endl;
      x.print();
      y.print();
      throw std::runtime_error("Series intersection is empty!");
    }
  }
  if(x.s > y.s){
    result.s = x.s;
    result.C = x.C;
  }
  else if(x.s < y.s){
    result.s = y.s;
    result.C = y.C;
  }
  else{
    result.s = x.s;
    if(!intersection(x.C,y.C,result.C))
    {
      throw std::runtime_error("Series intersection is empty!");
    }
  }
  result.mainSize = newSize;
  result.n = Series::get_n(result.type, newSize);
  result.main = newMain;
  
  return result;
}
Series operator-(const Series& x,const Series& y) 
{
  return x + (interval(-1)*y) ;
}

Series operator*(const Series& x,const Series&  y) 
{
  return mult(x,y);
}







Series operator*(const Series&x,const interval&  a) 
{
    Series y = x;
    y.main = y.main * a;
    y.C = y.C * a;
    return y; 
    
}

Series operator*(const interval&  a,const Series&x) 
{
    Series y = x;
    y.main = y.main * a;
    y.C = y.C * a;
    return y; 
}
Series operator+(const Series& x,const Series& y) 
{
    if(x.type!=y.type)
        throw std::runtime_error("Non maching types");
    Series result;
    result.type = x.type;
    int newSize = std::max(x.mainSize,y.mainSize);
    IVector newMain(newSize);
    for(int i=0;i<newSize;i++){
        newMain[i] = x.valueAt(Series::getSeriesIndex(result.type,i))+
                     y.valueAt(Series::getSeriesIndex(result.type,i));
    }
    result.main = newMain;
    result.mainSize = newSize;
    if(x.s == y.s){
      result.s = x.s;
      result.C = x.C + y.C;
      
    }
    else if(x.s <y.s)
    {
      //
      result.s = x.s;
      int next_n = Series::getSeriesIndex(result.type,newSize);
      interval coef= next_n^(x.s - y.s);
      result.C = x.C +y.C*interval(0,1)*coef;
    }
    else{
      result.s = y.s;
      int next_n = Series::getSeriesIndex(result.type,newSize);
      interval coef= next_n^(y.s - x.s);
      result.C = y.C +x.C*interval(0,1)*coef;
    }
    result.n = Series::get_n(result.type, newSize);
    return result;  
}

Series squere(const Series& x) 
{
    if(x.type== SeriesType::cos || x.type == SeriesType::cos_even)
        throw std::runtime_error("squere operation not supported");
    int newSize;
    Series result;
    if(x.type == SeriesType::sin ){
        result.type = SeriesType::cos;
        newSize = 2*x.mainSize+1;
    }
    if(x.type == SeriesType::sin_odd){
        result.type = SeriesType::cos_even;
        newSize = 2*x.mainSize;
    }
    IVector newMain(newSize);
    int n = x.n;
    for(int i=1;i<=n;i++){
        newMain[0] += (1./2) *  (x.valueAt(i)^2);
    }
    
    for(int j=1;j<newSize;j++){
        int k = Series::getSeriesIndex(result.type,j);
        for(int i=1;i<=n;i++)
            newMain[j] += x.valueAt(i)*x.valueAt(i+k);

        for(int i=1;i<= k-1;i++)
            newMain[j] +=  (-1)*(1./2) * x.valueAt(i)*x.valueAt(k-i);
    }
    
    interval estimationsMain = ((abs(x.C)^2) / (interval(2)*x.s -interval(1) ) ) * (interval(n)^(interval(-2)*x.s + interval(1) ) );
    newMain[0] += interval(0,1./2) * estimationsMain;
    for(int j=1;j<newSize;j++){
        newMain[j] += interval(-1,1) * estimationsMain;
    }
    
    interval C_1(0);
    for(int i=1;i<=n;i++){
        C_1 += abs(x.valueAt(i)) ;
    }
    C_1 = x.C * C_1;
    interval C_2(0);
    for(int i=1;i<=n;i++){
        interval coef = (interval(1) + interval(i)/interval(2*n+1-i))^x.s;
        C_2 += abs(x.valueAt(i)) * coef ;
    }
    C_2 = x.C * C_2;
    interval estimationTail = (abs(x.C)^2) *(interval(2)^x.s + 1) * 
                              (interval(1) /(interval(1)*x.s - interval(1))) * (interval(n)^(interval(-1)*x.s + interval(1) ) );
    result.s = x.s; 
    result.C = C_1 * interval(-1,1)+ C_2*interval(-1,1) + estimationTail * interval(-1,1);  
    result.main = newMain;
    result.mainSize = newSize;
    result.n = Series::get_n(result.type,result.mainSize);    
    return result;
}
Series static sinTimesCos(const Series& xSin,const Series& yCos){
  int newSize;
  Series result;
  if(xSin.mainSize + 1 != yCos.mainSize){
      throw std::runtime_error("sinTimesCos Size error\n");
  }
  if(xSin.type == SeriesType::sin ){
    result.type = SeriesType::sin;
    newSize = 2*xSin.mainSize;
    
  }
  if(xSin.type == SeriesType::sin_odd){
    result.type = SeriesType::sin_odd;
    newSize = 2*xSin.mainSize;
  }
  result.n = Series::get_n(result.type,newSize);
  int n = yCos.n;
  interval s = yCos.s;
  IVector newMain(newSize);
  for(int j=0;j<newSize;j++){
    int k = Series::getSeriesIndex(result.type,j);
    for(int i=1;i<=n;i++){
       newMain[j] += interval(1./2)*xSin.valueAt(i+k)*yCos.valueAt(i) - interval(1./2)*xSin.valueAt(i)*yCos.valueAt(i+k);
    }
      newMain[j] += yCos.valueAt(0)*xSin.valueAt(k);
    for(int i=1;i<=k-1;i++){
      newMain[j]+= interval(1./2) * xSin.valueAt(i) * yCos.valueAt(k-i);
    }
  }
  interval estimationsMain = ((abs(xSin.C * yCos.C)) / (interval(2)*s -interval(1) ) ) * (interval(n)^(interval(-2)*s + interval(1) ) );
  for(int j=0;j<newSize;j++){
    newMain[j] += interval(-1,1) * estimationsMain;
  }
  interval C_1(0);
  for(int i=1;i<=n;i++){
    C_1 +=interval(1./2)*(abs(xSin.C)*(abs(yCos.valueAt(i)) )+abs(yCos.C)*abs(xSin.valueAt(i) ) );
  }
  interval C_2(0);
  for(int i=1;i<=n;i++){
    interval coef = (interval(1) + interval(i)/interval(2*n+1-i))^s;

    C_2 += interval(1./2)*(abs(yCos.C)*(abs(xSin.valueAt(i)) )+abs(xSin.C)*abs(yCos.valueAt(i) ) )*coef;
  }
  interval estimationTail = abs(xSin.C*yCos.C) *(interval(2)^s + 1) * (interval(1) /(interval(1)*s - interval(1))) * (interval(n)^(interval(-1)*s + interval(1) ) );
  result.main = newMain;
  result.mainSize = newSize;
  result.n = Series::get_n(result.type,result.mainSize); 
  interval C = interval(-1,1)*C_1+ interval(-1,1)*C_2 + interval(-1,1)*estimationTail + abs(yCos.valueAt(0)) * xSin.C * interval(-1,1);
  result.C = C;
  result.s = s;
  return result;
}

Series static sinTimesSin(const Series& x,const Series& y){
  if(x.mainSize != y.mainSize){
      throw std::runtime_error("Size error\n");
  }
  int newSize;
  Series result;
  if(x.type == SeriesType::sin ){
    result.type = SeriesType::cos;
    newSize = 2*x.mainSize+1;
  }
  if(x.type == SeriesType::sin_odd){
    result.type = SeriesType::cos_even;
    newSize = 2*x.mainSize;
  }
  result.n = Series::get_n(result.type,newSize);
  int n = x.n;
  interval s = x.s;
  IVector newMain(newSize);
  for(int i=1;i<=n;i++){
    newMain[0] += interval(1./2) *  (x.valueAt(i)*y.valueAt(i));
  }
   for(int j=1;j<newSize;j++){
    int k = Series::getSeriesIndex(result.type,j);
    for(int i=1; i<= n;i++)
      newMain[j] += interval(1./2)* (x.valueAt(i)*y.valueAt(i+k) +x.valueAt(i+k)*y.valueAt(i));
    for(int i=1;i<= k-1;i++)
      newMain[j] +=  (-1)*interval(1./2) * x.valueAt(i)*y.valueAt(k-i);
  }
  interval estimationsMain = ((abs(x.C * y.C)) / (interval(2)*s -interval(1) ) ) * (interval(n)^(interval(-2)*s + interval(1) ) );
  newMain[0] += interval(-1./2,1./2) *estimationsMain; 
  for(int j=1;j<newSize;j++){
    newMain[j] += interval(-1,1) * estimationsMain;
  }
  interval C_1(0);
  for(int i=1;i<=n;i++){
    C_1 +=interval(1./2)*(abs(x.C)*abs(y.valueAt(i)) +abs(y.C)*abs(x.valueAt(i) ) );
  }
  interval C_2(0);
  for(int i=1;i<=n;i++){
    interval coef = (interval(1) + interval(i)/interval(2*n+1-i))^s;
    C_2 += interval(1./2)*(abs(x.C)*abs(y.valueAt(i)) + abs(y.C)*abs(x.valueAt(i)) )*coef;
  }
  interval estimationTail = abs(x.C*y.C) * (interval(2)^s + 1) * 
                            (interval(1) /(interval(1)*s - interval(1))) * (interval(n)^(interval(-1)*s + interval(1) ) );
  result.main = newMain;
  interval C;
  C = C_1*interval(-1,1)+ C_2*interval(-1,1) + estimationTail *interval(-1,1);
  result.C = C;
  result.s = s;
  result.mainSize = newSize;
  return result;
}
Series mult(const Series& x,const Series& y) 
{
    if(x.s<interval(1) || y.s<interval(1)){
      throw std::runtime_error("Series usufficient decays\n");
    }
    if(x.s != y.s){
        throw std::runtime_error("Series have different decays\n");
    }
    if((x.type == SeriesType::cos && y.type == SeriesType::cos) || (x.type == SeriesType::cos_even && y.type == SeriesType::cos_even)){
        throw std::runtime_error("cos series times cos series is not supported\n");
    }
    if((x.type == SeriesType::sin && y.type == SeriesType::cos )||(x.type == SeriesType::sin_odd && y.type == SeriesType::cos_even)){
        
        return sinTimesCos(x,y);
    }
    if((y.type == SeriesType::sin && x.type == SeriesType::cos )||(y.type == SeriesType::sin_odd && x.type == SeriesType::cos_even)){
        return sinTimesCos(y,x);
    }
    
    if((x.type == SeriesType::sin && y.type == SeriesType::sin )||(x.type == SeriesType::sin_odd && y.type == SeriesType::sin_odd)){
        return sinTimesSin(x,y);
    }
    throw std::runtime_error("no maching types multiplication");
}


SeriesVector::SeriesVector(int size,capd::interval C,capd::interval s,SeriesType type) 
{
  this->vec = std::vector<Series> (size);
  for(int i=0;i<size;i++){
    this->vec[i] =Series(C,s,type);
  }
}

SeriesVector::SeriesVector(int size) 
{
  this->vec = std::vector<Series> (size);
}

SeriesVector::SeriesVector() 
{
  
}

void SeriesVector::print() 
{
  for(int i=0;i<this->vec.size();i++){
    std::cout <<"Vector Series index: " << i <<"\n";
    vec[i].print();
  }
}

bool SeriesVector::subset(const SeriesVector& x) 
{
  for(int i=0;i<x.vec.size();i++){
    if( !(this->vec[i]).subset(x.vec[i]) )
      return false;
  }
  return true; 
}

bool SeriesVector::subsetInterior(const SeriesVector& x) 
{
  for(int i=0;i<x.vec.size();i++){
    if( !(this->vec[i]).subsetInterior(x.vec[i]) )
      return false;
  }
  return true;
}

IVector SeriesVector::tailSum() 
{
  IVector result((this->vec).size());
  for(int i=0;i<(this->vec).size();i++){
    result[i] = (this->vec[i]).tailSum();
  }
  return result;
}

SeriesVector elementWiseMult(const SeriesVector& x,const SeriesVector& y) 
{
  if(x.vec.size() != y.vec.size())
    throw std::runtime_error("Size error");
  SeriesVector z(x.vec.size());
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = elementWiseMult(x.vec[i], y.vec[i]);
  } 
  return z;
}

SeriesVector SeriesVector::elementWiseInverse() 
{
  SeriesVector z = *(this);
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = z.vec[i].elementWiseInverse();
  } 
  return z;
}

SeriesVector SeriesVector::upperBound() 
{
  SeriesVector z = *(this);
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = z.vec[i].upperBound();
  } 
  return z;
}

SeriesVector SeriesVector::lowerBound() 
{
  SeriesVector z = *(this);
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = z.vec[i].lowerBound();
  } 
  return z;
}

SeriesVector exp(const SeriesVector &x,const capd::interval& sNew) 
{
  SeriesVector z = x;
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = exp(x.vec[i],sNew);
  } 
  return z;
}

SeriesVector expMinusOne(const SeriesVector& x) 
{
  SeriesVector z = x;
  for(int i=0;i<z.vec.size();i++){
    z.vec[i] = expMinusOne(x.vec[i]);
  } 
  return z;
}

SeriesVector semiIntersection(const SeriesVector& x,const SeriesVector& y) 
{
  if(x.vec.size() != y.vec.size())
    throw std::runtime_error("Size intersection error\n");
  SeriesVector z(x.vec.size());
  for(int i=0;i<x.vec.size();i++){
    z.vec[i] = semiIntersection(x.vec[i],y.vec[i]);
  }
  return z;
}

Series& SeriesVector::operator[](int i) 
{
  return this->vec[i];
}
SeriesVector operator*(const capd::interval&  a,const SeriesVector& x) 
{
  SeriesVector y(x.vec.size());
  for(int i=0;i<x.vec.size();i++){
    y.vec[i] = a * x.vec[i];
  }
  return y;
}
SeriesVector operator*(const SeriesVector& x,const capd::interval&  a) 
{
  SeriesVector y(x.vec.size());
  for(int i=0;i<x.vec.size();i++){
    y.vec[i] = a * x.vec[i];
  }
  return y;
}
SeriesVector operator+(const SeriesVector& x,const SeriesVector& y) 
{
  if(x.vec.size()!=y.vec.size())
    throw std::runtime_error("wrong sizes!\n");
  SeriesVector z(x.vec.size());
  for(int i=0;i<x.vec.size();i++){
    z.vec[i] = x.vec[i] + y.vec[i];
  }    
  return z;
}
