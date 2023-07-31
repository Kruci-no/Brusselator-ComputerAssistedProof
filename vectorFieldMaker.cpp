#include "vectorFieldMaker.h"

int getCoffSinTimesCos(int i1,int i2,int k)
{
		return (i1 - i2 == k) + (-1)*(i1 - i2 == -k) + (i1+i2 == k )  ;
}
int getCoffSinTimesSin(int i1,int i2,int k)
{
		return (i1 - i2 == k) + (i1 - i2 == -k) + (-1)*(i1+i2 == k )  ;
}

void printMat(std::vector< std::vector<int> > x){
	for(int i=0;i<x.size();i++){
		for (int j = 0; j < x[i].size(); j++) {
			std::cout << x[i][j] << " ";
		}
		std::cout  << '\n';
	}

}

void printVect(std::vector<int>  x){
	for(int i=0;i<x.size();i++){
			std::cout << x[i] << " ";
	}
	std::cout  << '\n';
}

string coffFromInt(int x){
	if(x>=0)
		return to_string(x);
	else
		return "(-"+ to_string(x*(-1))+")";
}

string analizeRowSinSquered(std::vector<string> u, std::vector<int>  x,int row_number){
	bool allZero = true;
	for(int i=0;i<x.size();i++){
		if(x[i] != 0)
			allZero = false;
	}
	if(allZero == true)
		return "";
	bool squere = true;
	for(int i=0;i<x.size();i++){
		if(x[i] != 0 && i!= row_number)
			squere = false;
	}
	if(squere == true)
		return coffFromInt(x[row_number])+"*"+u[row_number]+"^2";
	//std::cout  << '\n';
	string temp="";
	bool first = true;
	//string coef="";
	for(int i=0;i<x.size();i++){
		if(x[i] != 0){
			if(first){
				first = false;
				temp = coffFromInt(x[i])+"*"+u[i];
			}
			else{
				temp = temp + "+" +coffFromInt(x[i])+"*" + u[i];
			}
		}
	}
	return  u[row_number]+"*("+temp+")";
}

int ComputeRowPoints(std::vector<int> x,int num){

	bool onlySquere = true;
	if(x[num] == 0 )
		false;
		for(int i=0;i<x.size();i++){
			if(i!=num && x[i]!=0)
				onlySquere = false;
		}
	if(onlySquere)
		return 1;
	bool allZero = true;
	int points=-1;
	for(int i=0;i<x.size();i++){
			if(x[i]!=0){
				points = points + 1;
				allZero = false;
			}
	}
	if(allZero)
		return 1;
	else
		return points;

}

void optimiseMatSinSquered (std::vector< std::vector<int> >&  x){
	int pointsR1 = 0 ;
	int pointsR2 = 0;
	std::vector<int> r1(x.size());
	std::vector<int> r2(x.size());
	//std::cout << "before" << '\n';
	//printMat(x);
	for(int i=0;i<x.size();i++){
			for(int j=0;j<x.size();j++){
				if(i!=j){
					 int prevPoints = ComputeRowPoints(x[i],i) + ComputeRowPoints(x[j],j);
					 r1 = x[i];
					 r2 = x[j];
					 r1[j] = x[j][i];
					 r2[i] = x[i][j];
					 int Points = ComputeRowPoints(r1,i) + ComputeRowPoints(r2,j);
					 if(Points > prevPoints){
						 x[i] = r1;
						 x[j] = r2;
					 }
				}
			}
	}
	//std::cout << "after Lifting" << '\n';
	//printMat(x);
	//std::cout << "----------------------" << '\n';
}

string analizeMatSinSquered(std::vector<string> u, std::vector< std::vector<int> >  x){
	optimiseMatSinSquered(x);
	//printMat(x);
	string result = "";
	string temp = "";
	bool first =true;
	for(int i=0;i<x.size();i++){
		temp = analizeRowSinSquered(u,x[i],i);
		if(temp!="")
			if(first){
				first = false;
				result = temp;
			}
			else{
				result = result+ "+"+temp;
			}
	}
	if(result!="")
		result = "0.5*("+ result+ ")";
	else
		result = "";
	return result;
}

std::vector<string> sinSquered(std::vector<string> u,int maxSize){
  std::vector<string> result(maxSize);
  string temp;

  for(int k = 1;k<maxSize;k++){
		std::vector< std::vector<int> > CoefMatrix(u.size());
		for(int i=0;i<u.size();i++){
			CoefMatrix[i].resize(u.size());
		}
		//std::cout << "elo" << '\n';
		//auto coeffArray = new int[u.size()][u.size()];
    temp="";
		for(int i1 = 0; i1<u.size(); i1++){
      if(u[i1]=="")
			{continue;}
			for(int i2 = 0; i2<=i1; i2++){
				if(u[i2]=="")
				{continue;}
				int coef =getCoffSinTimesSin(i1+1,i2+1,k);
        if(i2!=i1)
          coef = 2*coef;
				CoefMatrix[i2][i1] += coef;
			}

  	}
		result[k] = analizeMatSinSquered(u,CoefMatrix);
	}
	std::vector< std::vector<int> > CoefMatrix(u.size());
	for(int i=0;i<u.size();i++){
		CoefMatrix[i].resize(u.size());
	}
  for(int i1 = 0;i1<u.size();i1++){
		if(u[i1]== "")
			continue;
		CoefMatrix[i1][i1] = 1;
  }
	result[0] = analizeMatSinSquered(u,CoefMatrix);
  //std::cout << analizeMat(u,CoefMatrix);
  return result;
}

std::vector<string> sinSquered(std::vector<string> u){
  return sinSquered(u,2*u.size()+1);
}

string analizeRow(std::vector<string> u,std::vector<string> v, std::vector<int>  x,int row_number){
	bool allZero = true;
	for(int i=0;i<x.size();i++){
		if(x[i] != 0)
			allZero = false;
	}
	if(allZero == true)
		return "";
	bool first = true;
	string temp="";
	for(int i=0;i<x.size();i++){
		if(x[i] != 0){
			if(first){
				first = false;
				temp = coffFromInt(x[i])+"*"+v[i];
			}
			else{
				temp = temp + "+" +coffFromInt(x[i])+"*" + v[i];
			}
		}
	}
	return  u[row_number]+"*("+temp+")";
}

string analizeMatrix(std::vector<string> u,std::vector<string> v, std::vector< std::vector<int> > x){
	string result = "";
	string temp = "";
	bool first =true;
	for(int i=0;i<x.size();i++){
		temp = analizeRow(u,v,x[i],i);
		if(temp!="")
			if(first){
				first = false;
				result = temp;
			}
			else{
				result = result+ "+"+temp;
			}
	}
	if(result!="")
		result = "0.5*("+ result+ ")";
	else
		result = "";
	return result;
}

std::vector<string> sinTimesCos(std::vector<string> sin,std::vector<string> cos,int maxSize){
  std::vector<string> result(maxSize);
  string temp;
  for(int k = 0;k<maxSize;k++){
		std::vector< std::vector<int> > CoefMatrix(sin.size());
		for(int i=0;i<sin.size();i++){
			CoefMatrix[i].resize(cos.size());
		}
		for(int i1 = 0; i1<sin.size(); i1++){
      if(sin[i1]=="")
			{continue;}
			for(int i2 = 0;i2<cos.size(); i2++){
				if(cos[i2]=="")
				{continue;}
				int coef =getCoffSinTimesCos(i1+1,i2,k+1);
				CoefMatrix[i1][i2] += coef;
			}
  	}
		//printMat(CoefMatrix);
		result[k] = analizeMatrix(sin,cos,CoefMatrix);
	}
  return result;
}

std::vector<string> sinTimesCos(std::vector<string> sin,std::vector<string> cos){
	return sinTimesCos(sin,cos, sin.size()+cos.size()-1);
}

std::vector<string> sinTimesSin(std::vector<string> sin1,std::vector<string> sin2,int maxSize){
  std::vector<string> result(maxSize);
  string temp;
  for(int k = 1;k<maxSize;k++){
		std::vector< std::vector<int> > CoefMatrix(sin1.size());
		for(int i=0;i<sin1.size();i++){
			CoefMatrix[i].resize(sin2.size());
		}
		for(int i1 = 0; i1<sin1.size(); i1++){
      if(sin1[i1]=="")
			{continue;}
			for(int i2 = 0;i2<sin2.size(); i2++){
				if(sin2[i2]=="")
				{continue;}
				int coef =getCoffSinTimesSin(i1+1,i2+1,k);
				CoefMatrix[i1][i2] += coef;
			}
  	}
		//printMat(CoefMatrix);
		result[k] = analizeMatrix(sin1,sin2,CoefMatrix);
	}
	std::vector< std::vector<int> > CoefMatrix(sin1.size());
	for(int i=0;i<sin1.size();i++){
		CoefMatrix[i].resize(sin2.size());
	}
  for(int i1 = 0;i1<std::min(sin1.size(),sin2.size());i1++){
		if(sin1[i1]== "" || sin2[i1]== "")
			continue;
		CoefMatrix[i1][i1] = 1;
  }
	result[0] = analizeMatrix(sin1,sin2,CoefMatrix);
  return result;
}

std::vector<string> sinTimesSin(std::vector<string> sin1,std::vector<string> sin2){
  return sinTimesSin(sin1,sin2,sin1.size()+sin2.size()+1);
}

std::vector<string> multiply(string mult ,std::vector<string> x){
	std::vector<string> y(x.size());
	for(int i=0;i<x.size();i++){
		if(x[i]!="0" && x[i]!="")
			y[i] = mult+ "*"+ x[i];
	}
	return y;
}

std::vector<string> multiply(std::vector<string> x, std::vector<string> y){
	if(x.size()!= y.size())
		throw std::runtime_error("Size Error");
	std::vector<string> z(x.size());
	for(int i=0;i<x.size();i++){
			if(x[i]!="" && y[i]!="")
				z[i] = "("+x[i] + "*"+ y[i]+")";
			else
				z[i] = "";

	}
	return z;
}

std::vector<string> add(std::vector<string> x, std::vector<string> y){
	if(x.size()!= y.size())
		throw std::runtime_error("Size Error");
	std::vector<string> z(x.size());
	for(int i=0;i<x.size();i++){
			if(x[i]!="" && y[i]!="")
				z[i] = "("+x[i] + "+"+ y[i]+")";
			else if(x[i]!= "")
				z[i] = x[i];
			else
				z[i] = y[i];
	}
	return z;
}

std::vector<string> makeStringVector(int size,int begin,int end, string name){
	std::vector<string> x(size);
	for(int i=begin;i<end;i++){
		x[i] = name+to_string(i);
	}
	return x;
}

std::vector<string> makeStringVectorWithGap(int size,int gap,string name){
	std::vector<string> x(gap*size+  1 - gap);
	for(int i=0;i<size;i++){
		x[gap*i] = name+to_string(i);
	}
	return x;
}

std::vector<string> makeStringVector(int size,string name){
	std::vector<string> x(size);
	for(int i=0;i<size;i++){
		x[i] = name+to_string(i);
	}
	return x;
}

std::vector<string> merge(std::vector<string> x, std::vector<string> y){
	std::vector<string> z = x;
	z.insert(z.end(), y.begin(), y.end());
	return z;
}

std::vector<string> shuffle(std::vector<string> x, std::vector<string> y){
	if(x.size()!=y.size())
	 	throw std::runtime_error("Size Error");
	std::vector<string> z(x.size() + y.size());
	for(int i=0;i<x.size();i++){
		z[2*i] = x[i];
		z[2*i+1] = y[i];
	}
	return z;
}

string toFormula(std::vector<string> parameters, std::vector<string> variables , std::vector<string> vectorField){
	string par ="par:"+ parameters[0];
	for(int i=1;i<parameters.size();i++){
			if(parameters[i]!="")
				par = par + "," + parameters[i];
	}
	par = par + ";\n";

	string var ="var:"+ variables[0];
	for(int i=1;i<variables.size();i++){
		 if(variables[i]!="")
			var = var + "," + variables[i];
	}
	var = var + ";\n";

	string fun ="fun:\n";
	bool wasFirst = false;
	for(int i=0;i<vectorField.size();i++){
		  if(vectorField[i]!="")
				if(wasFirst)
					fun = fun + ",\n" + vectorField[i];
				else{
					fun = fun + vectorField[i];
					wasFirst = true;
				}
	}
	fun = fun+ ";";
	return par + var + fun;
}

string toFormula(std::vector<string> variables,std::vector<string> vectorField){
	std::vector<string> parameters(1);
	parameters[0] = "";
	//std::cout << "/* message */" << '\n';
	return toFormula(parameters,variables,vectorField);
}

void split(std::vector<string> u,std::vector<string>& uMain,std::vector<string>& uDiss,int mainModesNum){
	uMain = u;
	uDiss = u;
	int counter = mainModesNum;
	for(int i=0;i<u.size();i++){
		if(u[i]!="")
			if(counter >0){
				counter -= 1;
				uDiss[i] = "";
			}
			else{
				uMain[i] = "";
			}
	}
}
