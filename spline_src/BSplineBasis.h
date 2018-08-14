#ifndef BSPLINEBASIS_H
#define BSPLINEBASIS_H

#include <vector>

using namespace std;

class BSplineBasis
{
public:
	BSplineBasis();
	BSplineBasis(int deg,const vector<double>& v);
	void Set(int deg,const vector<double>& v);
	bool Check(int deg,const vector<double>& v);
	void BasisFunction(int pos,double par,int deriv,vector<double>& val);
	
private:
	int p;//polynomial order
	int nbs;//# of basis functions
	vector<double> kv;//knot vector
};

#endif