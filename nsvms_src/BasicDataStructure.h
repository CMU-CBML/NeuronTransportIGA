#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <petsc.h>
#include <petscksp.h>

#include "petscsys.h"   
#include "petscmat.h"
using namespace std;
/////////////////////////////////
class Vertex3D
{
public:
	double coor[3];
	int label; //ID for inlet and outlet
	Vertex3D();	
};

class Element3D
{
public:
	int degree;
	int order;
	int nbf;
	int type;//0 for interior and 1 for boundary, for visualization purpose
	int bzflag;//0 for spline element, 1 for Bezier element
	int label;

	vector<int> IEN;
	vector<array<double, 3>> pts;//tmp
	vector<array<double,64>> cmat;
	vector<int> IDBC;
	double velocity[3];

	Element3D(int p = 3);
	void BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const;
	void Basis(double u, double v, double w, vector<double>& Nt, vector<array<double, 3>>& dNdt) const;
	void Para2Phys(double u, double v, double w, double pt[3]) const;

};

//mesh format conversion

void Raw2Vtk_hex(string fn);

void Rawn2Vtk_hex(string fn);

#endif