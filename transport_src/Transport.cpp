#include "Transport.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;
const double PI = 3.1415926;

double MatrixDet(double dxdt[3][3])
{
	double det = dxdt[0][0] * dxdt[1][1] * dxdt[2][2] + dxdt[0][1] * dxdt[1][2] * dxdt[2][0] + dxdt[0][2] * dxdt[2][1] * dxdt[1][0] -
		(dxdt[0][2] * dxdt[1][1] * dxdt[2][0] + dxdt[0][0] * dxdt[1][2] * dxdt[2][1] + dxdt[1][0] * dxdt[0][1] * dxdt[2][2]);

	return det;
}

void Matrix3DInverse(double dxdt[3][3], double dtdx[3][3])
{
	double det = MatrixDet(dxdt);
	dtdx[0][0] = 1 / det*(dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]);
	dtdx[0][1] = 1 / det*(dxdt[2][1] * dxdt[0][2] - dxdt[0][1] * dxdt[2][2]);
	dtdx[0][2] = 1 / det*(dxdt[0][1] * dxdt[1][2] - dxdt[1][1] * dxdt[0][2]);
	dtdx[1][0] = 1 / det*(dxdt[2][0] * dxdt[1][2] - dxdt[1][0] * dxdt[2][2]);
	dtdx[1][1] = 1 / det*(dxdt[0][0] * dxdt[2][2] - dxdt[0][2] * dxdt[2][0]);
	dtdx[1][2] = 1 / det*(dxdt[1][0] * dxdt[0][2] - dxdt[0][0] * dxdt[1][2]);
	dtdx[2][0] = 1 / det*(dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	dtdx[2][1] = 1 / det*(dxdt[0][1] * dxdt[2][0] - dxdt[0][0] * dxdt[2][1]);
	dtdx[2][2] = 1 / det*(dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0]);
}

Transport::Transport()
{
	comm = MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;
	judge = 0;
}

void Transport::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open()) {
		fin_cmat >> neles;
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i<neles; i++) {
			if (i == ele_process[add]) {
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (int j = 0; j < nfunctions; j++) {
					for (int k = 0; k < 64; k++) {
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else {
				fin_cmat >> stmp >> nfunctions >> itmp;
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (int j = 0; j < nfunctions; j++)
					for (int k = 0; k < 64; k++)
						fin_cmat >> stmp;
			}
		}
		fin_cmat.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	add = 0;
	if (fin_bzpt.is_open()) {
		fin_bzpt >> npts;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++) {
			if (e == ele_process[add]) {
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++) {
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
				}
				add++;
			}
			else {
				for (int i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}
		}
		fin_bzpt.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

void Transport::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch(ng)
	{
	case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;							wght[1]=1.;
			break;
		}
	case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.1127016653792583;			Gpt[1]=0.5;							Gpt[2]=0.8872983346207417;
			wght[0]=0.5555555555555556;			wght[1]=0.8888888888888889;			wght[2]=0.5555555555555556;
			break;
		}
	case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.06943184420297371;			Gpt[1]=0.33000947820757187;			Gpt[2]=0.6699905217924281;			Gpt[3]=0.9305681557970262;
			wght[0]=0.3478548451374539;			wght[1]=0.6521451548625461;			wght[2]=0.6521451548625461;			wght[3]=0.3478548451374539;
			break;
		}
	case 5:
	{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;		Gpt[2] = 0.5;						Gpt[3] = 0.7692346550528415;		Gpt[4] = 0.953089922969332;
			wght[0] = 0.2369268850561891;		wght[1] = 0.4786286704993665;		wght[2] = 0.5688888888888889;		wght[3] = 0.4786286704993665;		wght[4] = 0.2369268850561891;
			break;
	}
	default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;							wght[1]=1.;
			break;
		}
	}
}

void Transport::InitializeProblem(const int n_bz, vector<array<double, 3>> &velocity_node, const vector<double>& N0_ini, const vector<double>& Nplus_ini, const vector<double>& Nminus_ini, const vector<double>& var)
{
	int i, j, e;
	
	MPI_Barrier(comm);

	PetscPrintf(PETSC_COMM_WORLD, "Initializing...\n");

	/*Initialize parameters*/
	GaussInfo(4);
	n_bzmesh = n_bz;

	N_0.resize(N0_ini.size());
	N_plus.resize(Nplus_ini.size());
	N_minus.resize(Nminus_ini.size());
	vplus.resize(N0_ini.size());
	vminus.resize(N0_ini.size());
	if (var.size() != 0)
		par = var;
	else
		cerr << "0 variables!\n"; getchar();
	dt = var[7];
	nstep = var[8];
	N_0 = N0_ini;
	N_plus = Nplus_ini;
	N_minus = Nminus_ini;

	/*Initialize petsc vector, matrix*/
	PetscInt mat_dim = N0_ini.size() * 2;
	ierr = MatCreate(PETSC_COMM_WORLD, &GK); 
	ierr = MatSetSizes(GK, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK, MATMPIAIJ); 
	ierr = MatMPIAIJSetPreallocation(GK, 250, NULL, 250, NULL); 																
	ierr = MatSetOption(GK, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); 
	ierr = MatSetUp(GK); 
	MatGetOwnershipRange(GK, &rstart, &rend);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_solution); 
}

void Transport::SUPGcoefficient(double s, double v[3], double dudx[3][3], vector<array<double, 3>>& dNdx, double &tau_supg)
{
	/*Calculate SUPG stabilization parameter*/
	double tau1(0.), tau2(0.);
	
	for (int i = 0; i < dNdx.size(); i++)
		tau1 += abs(v[0] * dNdx[i][0] + v[1] * dNdx[i][1] + v[2] * dNdx[i][2]);
	tau1 = 1 / tau1;
	tau2 = dt / 2.;
	tau_supg = 1. / sqrt(1. / (tau1*tau1) + 1. / (tau2*tau2));

}

void Transport::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[bzpt_num], double dNdx[bzpt_num][dim], double dudx[3][3], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[64][3];
	int i, j, k, a, b, loc(0);
	for (i = 0; i<4; i++){
		for (j = 0; j<4; j++){
			for (k = 0; k < 4; k++)	{
				Nx[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	double dxdt[3][3] = { { 0 } };
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a<3; a++)
			for (b = 0; b<3; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];
	double dtdx[3][3] = { { 0 } };
	Matrix3DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i<bzpt_num; i++){
		dNdx[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0] + dNdt[i][2] * dtdx[2][0];
		dNdx[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1] + dNdt[i][2] * dtdx[2][1];
		dNdx[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.125*detJ;
}

void Transport::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, double dudx[3][3], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[bzpt_num][3];
	double Nx_bz[bzpt_num];
	double dNdx_bz[bzpt_num][3];

	Nx.clear();
	dNdx.clear();
	Nx.resize(cmat.size(), 0);
	dNdx.resize(cmat.size(), { 0, 0, 0 });

	int i, j, k, a, b, c, loc;
	loc = 0;
	for (i = 0; i<4; i++){
		for (j = 0; j<4; j++){
			for (k = 0; k < 4; k++)	{
				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}

	double dxdt[3][3] = { { 0 } };
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a<3; a++)
			for (b = 0; b<3; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

	double dtdx[3][3] = { { 0 } };
	Matrix3DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i<bzpt_num; i++){
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0] + dNdt[i][2] * dtdx[2][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1] + dNdt[i][2] * dtdx[2][1];
		dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.125*detJ;

	for (i = 0; i < cmat.size(); i++){
		for (j = 0; j < bzpt_num; j++){
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 3; m++)	{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
			}
		}
	}
}

void Transport::WeightingFunction(const double velocity[3], const double& s, const double& tau, const vector<double> &Nx, const vector<array<double, 3>> &dNdx, vector<double> &Wx)
{
	/*Calculate weighting function using SUPG parameter Tau*/
	double U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
	for (int loc = 0; loc < Nx.size(); loc++)
		Wx[loc] = Nx[loc] + tau*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
}

void Transport::ElementValue(const vector<double> &Nx, const vector<double> value_node, double &value)
{
	value = 0.;
	for (int i = 0; i < Nx.size(); i++)
		value += value_node[i] * Nx[i];
}

void Transport::ElementVelocity(const vector<double> &Nx, const vector<array<double, 3>> &v_node, double v_tmp[3])
{
	for (int j = 0; j < dim; j++)
		v_tmp[j] = 0;
	for (int i = 0; i < Nx.size(); i++)
		for (int j = 0; j < dim; j++)
			v_tmp[j] += v_node[i][j] * Nx[i];
}

void Transport::Tangent(const int nen, vector<double>& Nx, vector<double>& Npx, vector<double>& Nmx, vector<array<double, 3>>& dNdx, double vp[3], double vm[3], double detJ, vector<vector<double>>& EMatrixSolve)
{
	/*calculate tangent matrix*/
	int i, j, A, B;
	for (i = 0; i < nen; i++){
		for (j = 0; j < nen; j++){
			EMatrixSolve[i + 0][j + 0] += ((1+ dt*(par[3] + par[4]))*Nx[i]*Nx[j] + dt*par[0]*(dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2]))*detJ;
			EMatrixSolve[i + 0][j + nen] += -par[5] * dt *Nx[i] * Nx[j] * detJ;
			EMatrixSolve[i + nen][j + 0] += -par[3] * dt *Npx[i] * Nx[j] * detJ;
			EMatrixSolve[i + nen][j + nen] += (1+ dt*par[5])*Npx[i] * Nx[j] * detJ + dt*Npx[i] * (vp[0] * dNdx[j][0] + vp[1] * dNdx[j][1] + vp[2] * dNdx[j][2])*detJ;
			//EMatrixSolve[i + 0][j + 16] += -par[6] * dt * EM[i][j];
			//EMatrixSolve[i + 16][j + 0] += -par[4] * dt * EM_minus[i][j];
			//EMatrixSolve[i + 16][j + 16] += EM_minus[i][j] + dt*EC_minus[i][j] + dt*par[6] * EM_minus[i][j];
		}
	}
}

void Transport::Residual(const int nen, const double N0, const double Nplus, const double Nminus, const vector<double> &Nx, const vector<double> &Npx, const vector<double> &Nmx, const double detJ, vector<double> &EVectorSolve)
{
	/*calculate residual of the equation*/
	int i, j, A;
	for (i = 0; i < nen; i++){
		EVectorSolve[i] += N0 * Nx[i] * detJ;
		EVectorSolve[i + nen] += Nplus * Npx[i] * detJ;
		//EVectorSolve[i + nen * 2] += Nminus*Nmx[i] * detJ;
	}
}

void Transport::ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve)
{
	int j, k;
	int nen = EVectorSolve.size() / 2;
	for (j = 0; j < nen * 2; j++){
		EVectorSolve[j] -= bc_value*EMatrixSolve[j][pt_num + variable_num * nen];
	}
	for (j = 0; j < nen * 2; j++){
		EMatrixSolve[j][pt_num + variable_num *nen] = 0.0; EMatrixSolve[pt_num + variable_num * nen][j] = 0.0;
	}
	EMatrixSolve[pt_num + variable_num * nen][pt_num + variable_num * nen] = 1.0;
	EVectorSolve[pt_num + variable_num * nen] = bc_value;
}

void Transport::MatrixAssembly(vector<vector<double>>& EMatrixSolve, const vector<int>& IEN, Mat& GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();


	PetscInt *nodeList = new PetscInt[nen * 2];	
	PetscReal *tmpGK = new PetscReal[nen * 2 * nen * 2];
	
	for (m = 0; m<IEN.size(); m++){
		A = IEN[m];
		for (i = 0; i < 2; i++)	{
			nodeList[2 * m + i] = 2 * A + i;
			for (n = 0; n<IEN.size(); n++){
				for (j = 0; j < 2; j++)	{
					tmpGK[add] = EMatrixSolve[m + i*nen][n + j*nen];
					add++;
				}
			}
		}
	}	
	MatSetValues(GK, nen * 2, nodeList, nen * 2, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void Transport::ResidualAssembly(vector<double>& EVectorSolve, const vector<int>& IEN, Vec& GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 2];
	PetscReal *tmpGR = new PetscReal[nen * 2];
	
	for (i = 0; i<IEN.size(); i++){
		A = IEN[i];
		nodeList[i * 2 + 0] = A * 2 + 0;
		nodeList[i * 2 + 1] = A * 2 + 1;
		tmpGR[add] = EVectorSolve[i];
		tmpGR[add + 1] = EVectorSolve[i + nen];
		add+=2;
	}
	VecSetValues(GR, nen * 2, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void Transport::BuildLinearSystemProcess(const vector<Element3D>& tmesh, const vector<Vertex3D> &cpts, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus)
{
	/*Build linear system in each process*/
	int e;
	cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for " << bzmesh_process.size() << " elements.\n";
	for (e = 0; e<bzmesh_process.size(); e++)
	{
		double detJ;
		double tau[2];
		double dudx[3][3];
		int nen = bzmesh_process[e].IEN.size();

		vector<vector<double>> EMatrixSolve;
		vector<double>EVectorSolve;
		vector<double>Nx, Npx, Nmx;
		vector<array<double, 3>> dNdx;

		EMatrixSolve.clear();	EVectorSolve.clear();
		Nx.clear(); Npx.clear(); Nmx.clear();
		dNdx.clear();

		EMatrixSolve.resize(nen * 2);	EVectorSolve.resize(nen * 2);
		Nx.resize(nen); Npx.resize(nen); Nmx.resize(nen);
		dNdx.resize(nen);

		for (int i = 0; i < nen; i++){
			EMatrixSolve[i].resize(nen * 2, 0.0);
			EMatrixSolve[i + nen].resize(nen * 2, 0.0);
		}

		for (int i = 0; i < nen * 2; i++){
			for (int j = 0; j < nen * 2; j++)	{
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		double vplus_tmp[3], vminus_tmp[3];
		double N0_tmp, Nplus_tmp, Nminus_tmp;

		vector<double> N0_cpt, Nplus_cpt, Nminus_cpt;
		N0_cpt.resize(nen); Nplus_cpt.resize(nen); Nminus_cpt.resize(nen);
		vector<array<double, 3>> vplus_cpt, vminus_cpt;
		vplus_cpt.resize(nen); vminus_cpt.resize(nen);

		for (int i = 0; i < nen; i++){
			N0_cpt[i] = N_0[bzmesh_process[e].IEN[i]];
			Nplus_cpt[i] = N_plus[bzmesh_process[e].IEN[i]];
			Nminus_cpt[i] = N_minus[bzmesh_process[e].IEN[i]];
			for (int j = 0; j < dim; j++){
				vplus_cpt[i][j] = velocity_node[bzmesh_process[e].IEN[i]][j];
				vminus_cpt[i][j] = -velocity_node[bzmesh_process[e].IEN[i]][j];
			}
		}

		for (int i = 0; i<Gpt.size(); i++){
			for (int j = 0; j<Gpt.size(); j++){
				for (int k = 0; k < Gpt.size(); k++){
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
					ElementVelocity(Nx, vplus_cpt, vplus_tmp);
					ElementVelocity(Nx, vminus_cpt, vminus_tmp);
					ElementValue(Nx, N0_cpt, N0_tmp);
					ElementValue(Nx, Nplus_cpt, Nplus_tmp);
					ElementValue(Nx, Nminus_cpt, Nminus_tmp);
					SUPGcoefficient(par[5], vplus_tmp, dudx, dNdx, tau[0]);
					SUPGcoefficient(par[6], vminus_tmp, dudx, dNdx, tau[1]);
					WeightingFunction(vplus_tmp, par[5], tau[0], Nx, dNdx, Npx);
					WeightingFunction(vminus_tmp, par[6], tau[1], Nx, dNdx, Nmx);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					Tangent(nen, Nx, Npx, Nmx, dNdx, vplus_tmp, vminus_tmp, detJ, EMatrixSolve);
					Residual(nen, N0_tmp, Nplus_tmp, Nminus_tmp, Nx, Npx, Nmx, detJ, EVectorSolve);
				}
			}
		}

		/*Apply Boundary Condition*/
		double N0_bc, Nplus_bc, Nminus_bc;
		for (int i = 0; i < nen; i++){
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1)	{
				//inlet
				N0_bc = par[9];
				ApplyBoundaryCondition(N0_bc, i, 0, EMatrixSolve, EVectorSolve);
				Nplus_bc = par[10];
				ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			}

			//movie2 photoactivation BCs
			//if (A<201)
			//{
			//	N0_bc = 0.0;
			//	ApplyBoundaryCondition(N0_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0.0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//}
			//if (A>=2412 && A<2613)
			//{
			//	N0_bc = 0.0;
			//	ApplyBoundaryCondition(N0_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0.0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//}

			//////movie5 photoactivation BCs
			//if (A<201)
			//{
			//	N0_bc = 0.0;
			//	ApplyBoundaryCondition(N0_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0.0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//}
			//if ((A >= 25821 && A <= 26021) || (A >= 16374 && A <= 16574))
			//{
			//	N0_bc = 0.0;
			//	ApplyBoundaryCondition(N0_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0.0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//}

		}
		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR);
		if (judge == 0) //The matrix is the same, so only need to assembly once
			MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK);
	}
	cout << "Process " << comRank << " :complete build matrix and vector!\n";
	VecAssemblyBegin(GR);
	if (judge == 0)
		MatAssemblyBegin(GK, MAT_FINAL_ASSEMBLY);
}

void Transport::ConcentrationCal_Coupling_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ)
{
	double dUdx[3][3];
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 3>> dNdx(bzel.IEN.size());
	bzel.Para2Phys(u, v, w, pt);
	BasisFunction(u, v, w, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
	disp = 0.;
	dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		disp += Nx[i] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[0] += dNdx[i][0] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[1] += dNdx[i][1] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[2] += dNdx[i][2] * (N_0[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
	}
}

void Transport::VisualizeVTK_ControlMesh(const vector<Vertex3D>& spt, const vector<Element3D>& mesh, int step, string fn)
{
	string fname;
	stringstream ss;
	ss << step;
	ofstream fout;
	unsigned int i;
	fname = fn + "/controlmesh_allparticle_" + ss.str() + ".vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
				<< " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << N_plus.size() << "\nSCALARS AllParticles float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<N_plus.size(); i++)
		{
			fout << N_0[i]+N_plus[i]+N_minus[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Transport::VisualizeVTK_PhysicalDomain(int step, string fn)
{
	vector<array<double, 3>> spt_all;//sample points
	vector<double> sresult_all;
	vector<array<int, 8>> sele_all;
	double detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	double spt_proc[num_bzmesh_ele * 24];
	double sresult_proc[num_bzmesh_ele * 8];
	int sele_proc[num_bzmesh_ele * 8];
	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3], dudx[3];
					double result;
					ConcentrationCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh_process[e] , pt1, result, dudx, detJ);
					spt_proc[24 * e + loc * 3 + 0] = pt1[0];
					spt_proc[24 * e + loc * 3 + 1] = pt1[1];
					spt_proc[24 * e + loc * 3 + 2] = pt1[2];
					sresult_proc[8 * e + loc] = result;
					loc++;
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					sele_proc[8 * e + 0] = 8 * e + a*nns[1] + b*ns + c;
					sele_proc[8 * e + 1] = 8 * e + a*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 2] = 8 * e + a*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 3] = 8 * e + a*nns[1] + (b + 1)*ns + c;
					sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
					sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				}
			}
		}
	}


	double *spts = NULL;
	double *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (double*)malloc(sizeof(double) * 24 * n_bzmesh);
		sresults = (double*)malloc(sizeof(double) * 8* n_bzmesh);
		seles = (int*)malloc(sizeof(int) * 8 * n_bzmesh);

		displs_spts = (int*)malloc(nProcess * sizeof(int));
		displs_sresults = (int*)malloc(nProcess * sizeof(int));
		displs_seles = (int*)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i<nProcess; i++) {
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 24;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 8;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 8;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 24;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 8;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 8;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 8 * 3, MPI_DOUBLE, spts, recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 8 , MPI_DOUBLE, sresults, recvcounts_sresults, displs_sresults, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 8, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);


	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				array<double, 3> pt = { spts[i * 24 + j * 3 + 0], spts[i * 24 + j * 3 + 1], spts[i * 24 + j * 3 + 2] };
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 8 + j]);

			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 8> el;
				el[0] = pstart + seles[8 * sum_ele + 0];
				el[1] = pstart + seles[8 * sum_ele + 1];
				el[2] = pstart + seles[8 * sum_ele + 2];
				el[3] = pstart + seles[8 * sum_ele + 3];
				el[4] = pstart + seles[8 * sum_ele + 4];
				el[5] = pstart + seles[8 * sum_ele + 5];
				el[6] = pstart + seles[8 * sum_ele + 6];
				el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 8;
		}
		cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
	}
}

void Transport::WriteVTK(const vector<array<double, 3>> spt, const vector<double> sdisp, const vector<array<int, 8>> sele, int step, string fn)
{
	stringstream ss;
	ss << step;
	string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Transport::AssignProcessor(vector<vector<int>> &ele_proc)
{
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void Transport::Run(const vector<Vertex3D>& cpts, const vector<array<double, 3>> velocity_node, const vector<Element3D> &tmesh, string path_in, string path_out)
{
	PetscPrintf(PETSC_COMM_WORLD, "Diffusion...\n");
	ReadBezierElementProcess(path_in);
	for (int n = 0; n < nstep; ++n)
	{
		stringstream ss;
		ss << n;
		PetscPrintf(PETSC_COMM_WORLD, "Step: %d\n", n);
		
		/*Visualize the result from initial step*/
		if (comRank == 0)
			VisualizeVTK_ControlMesh(cpts, tmesh, n, path_out);	//solution on control points
		VisualizeVTK_PhysicalDomain(n, path_out);//solution on physical domain
		
		
	    /*Build Linear System*/
		time_t t0, t1;
		PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");
		if(judge==0)
			MatZeroEntries(GK);
		VecSet(GR, 0);
		t0 = time(NULL);
		BuildLinearSystemProcess(tmesh, cpts, velocity_node, par[1], par[2]);
		VecAssemblyEnd(GR);
		PetscPrintf(PETSC_COMM_WORLD, "Done Vector Assembly...\n");
		if (judge == 0)
			MatAssemblyEnd(GK, MAT_FINAL_ASSEMBLY);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 - t0);

		MPI_Barrier(comm);

		PetscPrintf(PETSC_COMM_WORLD, "Solving...\n");
		/*Petsc solver setting*/
		if (judge == 0)	{
			KSPCreate(PETSC_COMM_WORLD, &ksp);
			KSPSetOperators(ksp, GK, GK);
			KSPGetPC(ksp, &pc);
			PCSetType(pc, PCBJACOBI);
			KSPSetType(ksp, KSPGMRES);
			KSPGMRESSetRestart(ksp, 100);
			KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);			
			KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000);			
			KSPSetFromOptions(ksp);
			KSPSetUp(ksp);
			judge = 1;
		}

		MPI_Barrier(comm);

		t0 = time(NULL);
		KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
		KSPSolve(ksp, GR, temp_solution);
		PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "iterations %d\n", its);
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
		PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Solving with time: %d \n", t1 - t0);

		
		Vec temp_solution_seq;
		VecScatter ctx;
		PetscScalar    *_a;
		VecScatterCreateToAll(temp_solution, &ctx, &temp_solution_seq);
		VecScatterBegin(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecGetArray(temp_solution_seq, &_a);
		
		MPI_Barrier(comm);
		

		for (int i = 0; i < N_0.size(); i++){
			N_0[i] = PetscRealPart(_a[i * 2 + 0]);
			N_plus[i] = PetscRealPart(_a[i * 2 + 1]);
		}
		VecRestoreArray(temp_solution_seq, &_a);
		VecScatterDestroy(&ctx);
		VecDestroy(&temp_solution_seq);			
		PetscPrintf(PETSC_COMM_WORLD, "Step %d done!\n", n);
	}

	/*Visualize the final step*/
	if (comRank == 0)
		VisualizeVTK_ControlMesh(cpts, tmesh, nstep, path_out);	//solution on control points
	VisualizeVTK_PhysicalDomain(nstep, path_out);//solution on physical domain

	MatDestroy(&GK);
	VecDestroy(&GR);
	VecDestroy(&temp_solution);
	KSPDestroy(&ksp);
	MPI_Barrier(comm);
}

