#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <vector>
#include <array>
#include "BasicDataStructure.h"
#include "time.h"

using namespace std;

class Transport
{
private:

	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;
	int n_bzmesh;
	int rstart, rend;
	vector<int> ele_process;
	vector<Element3D> bzmesh_process;

	KSP ksp;
	PC pc;
	Mat GK;
	Vec GR;
	Vec temp_solution;

	int judge;// judge if the matrix have been solved
	vector<double> Gpt;
	vector<double> wght;
	vector<double> N_0;
	vector<double> N_plus;
	vector<double> N_minus;
	vector<array<double, 3>> vplus, vminus;
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-, dt, nstep, n0_bc, n+_bc, n-_bc
	double dt;
	int nstep;
public:
	Transport();
private:
	void ReadBezierElementProcess(string fn);

	void GaussInfo(int ng);
	void SUPGcoefficient(double s, double v[3], double dudx[3][3], vector<array<double, 3>>& dNdx, double &tau_supg);
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[bzpt_num], double dNdx[bzpt_num][dim], double dudx[3][3], double& detJ);
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, double dudx[3][3], double& detJ);
	void WeightingFunction(const double velocity[3], const double& s, const double& tau, const vector<double> &Nx, const vector<array<double, 3>> &dNdx, vector<double> &Wx);
	void ElementValue(const vector<double> &Nx, const vector<double> value_node, double &value);
	void ElementVelocity(const vector<double> &Nx, const vector<array<double, 3>>& v_node, double v_tmp[3]);
	void Tangent(const int nen, vector<double>& Nx, vector<double>& Npx, vector<double>& Nmx, vector<array<double, 3>>& dNdx, double vp[3], double vm[3], double detJ, vector<vector<double>>& EMatrixSolve);
	void Residual(const int nen, const double CA, const double Nplus, const double Nminus, const vector<double> &Nx, const vector<double> &Npx, const vector<double> &Nmx, const double detJ, vector<double> &EVectorSolve);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve);
	void MatrixAssembly(vector<vector<double>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<double>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void BuildLinearSystemProcess(const vector<Element3D>& tmesh, const vector<Vertex3D> &cpts, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus);
	
	void ConcentrationCal_Coupling_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex3D>& pts, const vector<Element3D>& mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 8>> sele, int step, string fn);
public:	
	void InitializeProblem(const int n_bz, vector<array<double, 3>> &velocity_node, const vector<double>& N0_ini, const vector<double>& Nplus_ini, const vector<double>& Nminus_ini, const vector<double>& var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex3D>& cpts, const vector<array<double, 3>> velocity_node, const vector<Element3D> &tmesh, string path_in, string path_out);

};

#endif