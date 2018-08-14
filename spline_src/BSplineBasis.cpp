#include "BSplineBasis.h"
#include <iostream>
//#include <fstream>
//#include <string>

using namespace std;

typedef unsigned int uint;

BSplineBasis::BSplineBasis()
{
	p=3; nbs=1;
	double tmp[5]={0.,1.,2.,3.,4.};
	kv.assign(tmp,tmp+5);
}

BSplineBasis::BSplineBasis(int deg,const vector<double>& v)
{
	if(Check(deg,v))
	{
		p=deg; kv=v; nbs=v.size()-deg-1;
	}
	else
	{
		p=3; nbs=1;
		double tmp[5]={0.,1.,2.,3.,4.};
		kv.assign(tmp,tmp+5);
	}
}

void BSplineBasis::Set(int deg,const vector<double>& v)
{
	if(Check(deg,v))
	{
		p=deg; kv=v; nbs=v.size()-deg-1;
	}
	else
	{
		p=3; nbs=1;
		double tmp[5]={0.,1.,2.,3.,4.};
		kv.assign(tmp,tmp+5);
	}
}

bool BSplineBasis::Check(int deg,const vector<double>& v)
{
	if(deg<0)
	{
		cerr<<"Wrong degree!\n";
		return false;
	}
	if(v.size()<=p+1)
	{
		cerr<<"Length of knot vector is wrong!\n";
		return false;
	}
	else
	{
		int flag(0);
		for(uint i=0;i<v.size()-1;i++)
		{
			if(v[i]>v[i+1])
			{
				flag=1;
				break;
			}
		}
		if(flag==1)
		{
			cerr<<"Knot vector is not in increasing order!\n";
			return false;
		}
		else
		{
			int rep(0);
			//check first knot
			for(uint i=0;i<v.size();i++)
			{
				if(v[i]==v[0]) rep++;
			}
			if(rep>p+1)
			{
				cerr<<"First knot repeated more than "<<p<<"+1 times\n";
				return false;
			}
			//check last knot
			rep=0;
			for(int i=v.size()-1;i>=0;i--)
			{
				if(v[i]==v.back()) rep++;
			}
			if(rep>p+1)
			{
				cerr<<"Last knot repeated more than "<<p<<"+1 times\n";
				return false;
			}
			//check interior knots
			for(uint i=1;i<v.size()-1;i++)
			{
				rep=0;
				for(uint j=i;j<v.size()-1;j++)
				{
					if(v[j]==v[i]) rep++;
				}
				if(rep>p)
				{
					cerr<<"Interior knots repeated more than "<<p<<" times\n";
					return false;
				}
			}
		}
	}
	return true;
}

void BSplineBasis::BasisFunction(int pos,double par,int deriv,vector<double>& val)
{
	val.clear();
	val.resize(deriv+1,0.);
	if(pos<0 || pos>=nbs)
	{
		cerr<<"Wrong basis functions ID!\n";
		return;
	}
	if(par<kv.front() || par>kv.back())
	{
		//cerr<<"Wrong parametric value!\n";
		return;
	}
	if(deriv>1)
	{
		cerr<<"Higher order of derivatives not available now!\n";
		return;
	}

	vector<double> N0(p+1,0.);
	double left,right,d1=0.;
	int i,j,q,flag;
	for(i=pos;i<=pos+p;i++)
	{
		if(par>=kv[i] && par<kv[i+1] || (par>kv[i] && par==kv[i+1] && par==kv.back()))
		{
			N0[i-pos]=1.;
			break;
		}
	}
	//recursive
	for(q=1;q<=p;q++)
	{
		j=p+1-q;
		for(i=0;i<j;i++)
		{
			if(kv[pos+i+q]==kv[pos+i]) left=0.;
			else left=(par-kv[pos+i])/(kv[pos+i+q]-kv[pos+i]);
			if(kv[pos+i+q+1]==kv[pos+i+1]) right=0.;
			else right=(kv[pos+i+q+1]-par)/(kv[pos+i+q+1]-kv[pos+i+1]);
			if(deriv==1 && q==p)//first derivative
			{
				double left1,right1;
				if(kv[pos+i+q]==kv[pos+i]) left1=0.;
				else left1=q/(kv[pos+i+q]-kv[pos+i]);
				if(kv[pos+i+q+1]==kv[pos+i+1]) right1=0.;
				else right1=q/(kv[pos+i+q+1]-kv[pos+i+1]);
				d1=left1*N0[0]-right1*N0[1];
			}
			//cout<<left<<" "<<right<<"\n";
			//cout<<N0[i]<<" "<<N0[i+1]<<"\n";
			N0[i]=left*N0[i]+right*N0[i+1];
			//cout<<N0[i]<<"\n";
			//getchar();
		}
	}
	val[0]=N0[0];
	if(deriv==1)
		val[1]=d1;
}