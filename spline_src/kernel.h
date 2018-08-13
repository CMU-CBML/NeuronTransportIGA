#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "BasicDataStructure.h"
#include <sstream>
#include <iomanip>
using namespace std;

class kernel
{
public:
	void run_coupling_comp(string fn_in);//comparison using B-splines-like	
	
private:
	vector<Vertex3D> cp;//control points
	vector<array<double, 3>> cpa;
	vector<Element3D> tmesh;//elements in T-mesh
	vector<Edge3D> tmedge;
	vector<Face3D> tmface;
	vector<array<double, 3>> bzcp;//Bezier control points

	void Run_Coupling_Comp(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void OutputMesh_Coupling(const vector<BezierElement3D>& bzmesh, string fn);
	void BezierExtract_Comp(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void BuildSplines_Unstruct();
	void BuildElementSplines_Interior(int eid);
	void BuildElementSplines_Boundary(int eid);
	void InitializeMeshLabel(string fn);
	void RescaleDomain();
	void InitialConnect();
	void BuildInitialEdges();
};

#endif