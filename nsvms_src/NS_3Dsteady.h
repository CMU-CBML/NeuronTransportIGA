#ifndef NS_3DSTEADY_H
#define NS_3DSTEADY_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "BasicDataStructure.h"
#include "time.h"

using namespace std;

const int degree = 3;
const int dim = 3;
const int bzpt_num = 64;

class NS_3Dsteady
{
private:
	vector<double> Gpt;
	vector<double> wght;
	const double PI = 4 * atan(1.0);

	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;

	int rstart, rend;
	int n_bzmesh;
	vector<int> ele_process;
	vector<Element3D> bzmesh_process;

	KSP ksp;
	PC pc;
	Mat GK;
	Vec GR;
	Vec temp_solution;
	
	double dt;
	double velocity_max;
	double nu;
	double rou;
	double alphaM;
	double alphaF;
	double Gama;
	double fx, fy, fz;//elemental body force

	vector<double> Vel;
	vector<double> Pre;
	vector<double> gh; //prescribed boundary velocities
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-
	

public:
	NS_3Dsteady();
private:
	void BodyForce(double x, double y, double z, double &Fx, double &Fy, double &Fz);
	void ReadBezierElementProcess(string fn);

	/*Analysis*/
	void GaussInfo(int ng);	
	void BasisFunction(double u, double v, double w, int nen, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, vector<array<array<double, 3>, 3>> &dN2dx2, double dudx[3][3], double& detJ);
	void PointFormValue(vector<double>& Nx, const vector<array<double, 4>>& U, double Value[4]);
	void PointFormGrad(vector<array<double, 3>>& dNdx, const vector<array<double, 4>>& U, double Value[4][3]);
	void PointFormHess(vector<array<array<double, 3>, 3>>& d2Ndx2, const vector<array<double, 4>> &U, double Value[4][3][3]);
	void Tau(double J[3][3], double u[4], double &tauM, double &tauC);
	void FineScale(double tauM, double tauC, double u[3], double u_x[3], double u_y[3], double u_z[3], double u_xx[3], double u_yy[3], double u_zz[3], double p, double p_x, double p_y,double p_z, double u_s[3], double &p_s);
	void Residual(vector<double>& Nx, vector<array<double, 3>>& dNdx, vector<array<array<double, 3>, 3>>& dN2dx2, double dudx[3][3], const double detJ, const vector<array<double, 4>> &U, vector<array<double, 4>> Re);
	void Tangent(vector<double> &Nx, vector<array<double, 3>>& dNdx, double dudx[3][3], const double detJ, const vector<array<double, 4>>& U, vector<array<vector<array<double, 4>>, 4>>& Ke);
	void BuildLinearSystemProcess(const vector<Vertex3D>& cpts, const vector<array<double, 3>>& velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, 4>>, 4>>& Ke, vector<array<double, 4>> &Re);
	void MatrixAssembly(vector<array<vector<array<double, 4>>, 4>>& Ke, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<array<double,4>> &Re, const vector<int>& IEN, Vec& GR);

	/*Postprocessing*/
	void ResultCal_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double result[4], double dudx[3], double& detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex3D>& pts, const vector<Element3D>& mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 8>> sele, int step, string fn);
public:
	/*Preprocessing*/
	void InitializeProblem(const int ndof, const int n_bz, const vector<double>& Vel0, const vector<double>& Pre0, const vector<double>& var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex3D>& cpts, const vector<Element3D>& tmesh, const vector<array<double, 3>>& velocity_bc, string fn);
};

#endif