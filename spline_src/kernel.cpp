#include "kernel.h"

typedef unsigned int uint;

void kernel::run_coupling_comp(string fn_in)
{
	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	Run_Coupling_Comp(fn_in, 0, bzmesh, IDBC, gh);
	OutputMesh_Coupling(bzmesh, fn_in);
}

void kernel::Run_Coupling_Comp(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	InitializeMeshLabel(fn);
	BuildSplines_Unstruct();
	BezierExtract_Comp(bzmesh, IDBC, gh);
}

void kernel::BezierExtract_Comp(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh)
{
	vector<int> aflag(cp.size(), 0);
	vector<int> aloc(cp.size(), -1);
	int count(0);
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].bzflag == 0)
		{
			for (uint j = 0; j < tmesh[i].IEN.size(); j++)
			{
				aflag[tmesh[i].IEN[j]] = 1;
			}
		}
		else
		{
			cerr << "Wrong: bzflag is 1 in the comparison case!\n"; getchar();
		}
	}
	for (uint i = 0; i < cp.size(); i++)
	{
		if (aflag[i] == 1)
		{
			aloc[i] = count++;
		}
	}
	cpa.clear();
	cpa.resize(count);
	count = 0;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (aflag[i] == 1)
		{
			cpa[count][0] = cp[i].coor[0];
			cpa[count][1] = cp[i].coor[1];
			cpa[count][2] = cp[i].coor[2];
			count++;
		}
	}

	IDBC.clear();
	IDBC.resize(count + bzcp.size());//bzcp.size()==0
	for (uint i = 0; i < IDBC.size(); i++)
	{
		IDBC[i] = i;//boundary weakly imposed
	}
	gh.clear();
	gh.resize(IDBC.size());//initial guess for iterative solvers
	count = 0;
	//double xcoor[3];
	array<double, 3> xcoor;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (aflag[i] == 1)
		{
			xcoor[0] = cp[i].coor[0];
			xcoor[1] = cp[i].coor[1];
			xcoor[2] = cp[i].coor[2];
			count++;
		}
	}

	cout << "Bezier extracting...\n";
	cout << "# Bezier: " << tmesh.size() << "\n";
	bzmesh.resize(tmesh.size());
#pragma omp parallel for
	for (int eid = 0; eid < tmesh.size(); eid++)
	{
		if (eid != 0 && eid % 500 == 0)
		{
			cout << eid << " ";
		}
		//find types of interfaces
		//if (tmesh[eid].bzflag == 1)
		if (tmesh[eid].bzflag == 0)
		{
			for (int k = 0; k < 6; k++)
			{
				if (tmface[tmesh[eid].face[k]].hex.size() == 1)//must be a boundary element
				{
					bzmesh[eid].bc[k] = 1;//boundary face for boundary condition
										  //bzmesh[eid].bcflag = 1;
				}
				//else if (tmface[tmesh[eid].face[k]].hex.size() == 2)
				//{
				//	int hxnb(tmface[tmesh[eid].face[k]].hex[0]);
				//	if (hxnb == eid) hxnb = tmface[tmesh[eid].face[k]].hex[1];
				//	if (tmesh[hxnb].bzflag == 0)
				//	{
				//		bzmesh[eid].bc[k] = 2;//coupling interface
				//		bzmesh[eid].bzcouple = 1;
				//	}
				//}
			}
		}

		double tmp;
		bzmesh[eid].type = tmesh[eid].type;//used for visualization
		bzmesh[eid].bzflag = tmesh[eid].bzflag;
		if (bzmesh[eid].bzcouple == 0)
		{
			if (tmesh[eid].bzflag == 0)
			{
				bzmesh[eid].IEN.resize(tmesh[eid].IEN.size());
				bzmesh[eid].cmat.resize(tmesh[eid].IEN.size(), vector<double>(64));
				for (uint k = 0; k < tmesh[eid].IEN.size(); k++)
				{
					bzmesh[eid].IEN[k] = aloc[tmesh[eid].IEN[k]];//some IEN may be -1
					for (int k1 = 0; k1 < 64; k1++)
					{
						tmp = tmesh[eid].bemat[k][k1];
						bzmesh[eid].cmat[k][k1] = tmp;
						bzmesh[eid].pts[k1][0] += tmp * cp[tmesh[eid].IEN[k]].coor[0];
						bzmesh[eid].pts[k1][1] += tmp * cp[tmesh[eid].IEN[k]].coor[1];
						bzmesh[eid].pts[k1][2] += tmp * cp[tmesh[eid].IEN[k]].coor[2];
					}
				}
			}
			else
			{
				bzmesh[eid].IENb.resize(tmesh[eid].IENb.size());//64
				for (uint k = 0; k < tmesh[eid].IENb.size(); k++)
				{
					bzmesh[eid].IENb[k] = tmesh[eid].IENb[k] + count;
					bzmesh[eid].pts[k][0] = bzcp[tmesh[eid].IENb[k]][0];
					bzmesh[eid].pts[k][1] = bzcp[tmesh[eid].IENb[k]][1];
					bzmesh[eid].pts[k][2] = bzcp[tmesh[eid].IENb[k]][2];
				}
			}
		}
		else//bzflag must be 1
		{
			bzmesh[eid].IEN.resize(tmesh[eid].IEN.size());
			bzmesh[eid].cmat.resize(tmesh[eid].IEN.size(), vector<double>(64));
			for (uint k = 0; k < tmesh[eid].IEN.size(); k++)
			{
				bzmesh[eid].IEN[k] = aloc[tmesh[eid].IEN[k]];
				for (int k1 = 0; k1 < 64; k1++)
				{
					bzmesh[eid].cmat[k][k1] = tmesh[eid].bemat[k][k1];
				}
			}
			bzmesh[eid].IENb.resize(tmesh[eid].IENb.size());
			for (uint k = 0; k < tmesh[eid].IENb.size(); k++)
			{
				bzmesh[eid].IENb[k] = tmesh[eid].IENb[k] + count;
				bzmesh[eid].pts[k][0] = bzcp[tmesh[eid].IENb[k]][0];
				bzmesh[eid].pts[k][1] = bzcp[tmesh[eid].IENb[k]][1];
				bzmesh[eid].pts[k][2] = bzcp[tmesh[eid].IENb[k]][2];
			}
		}
	}
}

void kernel::BuildSplines_Unstruct()
{
	cout << "# elements: " << tmesh.size() << "\n";
#pragma omp parallel for
	for (int i = 0; i < tmesh.size(); i++)
	{
		if (i != 0 && i % 500 == 0)
		{
			cout << i << " ";
		}
		if (tmesh[i].type != 1)
		{
			BuildElementSplines_Interior(i);
		}
		else
		{
			BuildElementSplines_Boundary(i);
		}
	}
}

void kernel::BuildElementSplines_Interior(int eid)
{
	//find IENtmp first
	tmesh[eid].IEN.clear();
	uint i, j, k, hxid;
	vector<int> loc(cp.size(), -1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[tmesh[eid].cnct[i]] = tmesh[eid].IEN.size();
		tmesh[eid].IEN.push_back(tmesh[eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
		{
			hxid = cp[tmesh[eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(), hx1r.end(), hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), tmesh[hxid].cnct[k]);
					if (it == tmesh[eid].IEN.end())
					{
						loc[tmesh[hxid].cnct[k]] = tmesh[eid].IEN.size();
						tmesh[eid].IEN.push_back(tmesh[hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();
	tmesh[eid].bemat.resize(tmesh[eid].IEN.size(), vector<double>(64, 0.));
	//8 body points, not consider boundary yet
	double w[2] = { 2. / 3., 1. / 3. };
	double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
		w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },{ 1, 2, 3, 0, 5, 6, 7, 4 },{ 2, 3, 0, 1, 6, 7, 4, 5 },{ 3, 0, 1, 2, 7, 4, 5, 6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 },{ 5, 6, 7, 4, 1, 2, 3, 0 },{ 6, 7, 4, 5, 2, 3, 0, 1 },{ 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		for (j = 0; j < 8; j++)//which body point, bezier
		{
			for (k = 0; k < 8; k++)//which local corner point, b-splines
			{
				bpm[i][j][k] = a[k];
				bpmap[i][j][k] = loc[tmesh[hx1r[i]].cnct[bpi[j][k]]];
			}
		}
	}
	int layer[4] = { 0, 16, 32, 48 };
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			tmesh[eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 },{ 1, 2 },{ 2, 3 },{ 3, 0 },{ 0, 4 },{ 1, 5 },{ 2, 6 },{ 3, 7 },{ 4, 5 },{ 5, 6 },{ 6, 7 },{ 7, 4 } };
	int edbz[12][2] = { { 1, 2 },{ 7, 11 },{ 14, 13 },{ 8, 4 },{ 0 + layer[1], 0 + layer[2] },{ 3 + layer[1], 3 + layer[2] },
	{ 15 + layer[1], 15 + layer[2] },{ 12 + layer[1], 12 + layer[2] },{ 1 + layer[3], 2 + layer[3] },{ 7 + layer[3], 11 + layer[3] },{ 14 + layer[3], 13 + layer[3] },{ 8 + layer[3], 4 + layer[3] } };
	int pos1, pos2;
	for (i = 0; i < 12; i++)
	{
		uint nhex = tmedge[tmesh[eid].edge[i]].hex.size();
		for (j = 0; j<tmedge[tmesh[eid].edge[i]].hex.size(); j++)
		{
			hxid = tmedge[tmesh[eid].edge[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[edi[i][0]]);
			pos2 = it1 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k] / nhex;
			}
			int* it2 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[edi[i][1]]);
			pos2 = it2 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3] };
	for (i = 0; i < 8; i++)
	{
		uint nhex = cp[tmesh[eid].cnct[i]].hex.size();
		for (j = 0; j<nhex; j++)
		{
			hxid = cp[tmesh[eid].cnct[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[i]);
			pos2 = it1 - tmesh[hxid].cnct;
			for (k = 0; k < 8; k++)
			{
				tmesh[eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
			}
		}
	}
	//4*6 face points
	int fci[6][4] = { { 0, 1, 3, 2 },{ 0, 1, 4, 5 },{ 1, 2, 5, 6 },{ 3, 2, 7, 6 },{ 0, 3, 4, 7 },{ 4, 5, 7, 6 } };
	int fcbz[6][4] = { { 5, 6, 9, 10 },{ 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] },{ 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] },{ 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] },{ 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	for (i = 0; i < 6; i++)
	{
		uint nhex = tmface[tmesh[eid].face[i]].hex.size();
		for (j = 0; j < nhex; j++)
		{
			hxid = tmface[tmesh[eid].face[i]].hex[j];
			vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
			pos1 = it - hx1r.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[fci[i][j1]]);
				pos2 = it1 - tmesh[hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					tmesh[eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
	}
}

void kernel::BuildElementSplines_Boundary(int eid)
{
	//find IEN first
	tmesh[eid].IEN.clear();
	uint i, j, k, hxid;
	vector<int> loc(cp.size(), -1);
	vector<int> hx1r(1, eid);
	for (i = 0; i < 8; i++)
	{
		loc[tmesh[eid].cnct[i]] = tmesh[eid].IEN.size();
		tmesh[eid].IEN.push_back(tmesh[eid].cnct[i]);
	}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
		{
			hxid = cp[tmesh[eid].cnct[i]].hex[j];
			vector<int>::iterator it1 = find(hx1r.begin(), hx1r.end(), hxid);
			if (it1 == hx1r.end())
			{
				hx1r.push_back(hxid);
				for (k = 0; k < 8; k++)
				{
					vector<int>::iterator it = find(tmesh[eid].IEN.begin(), tmesh[eid].IEN.end(), tmesh[hxid].cnct[k]);
					if (it == tmesh[eid].IEN.end())
					{
						loc[tmesh[hxid].cnct[k]] = tmesh[eid].IEN.size();
						tmesh[eid].IEN.push_back(tmesh[hxid].cnct[k]);
					}
				}
			}
		}
	}
	for (i = 0; i < tmesh[eid].bemat.size(); i++)
	{
		tmesh[eid].bemat[i].clear();
	}
	tmesh[eid].bemat.clear();
	tmesh[eid].bemat.resize(tmesh[eid].IEN.size(), vector<double>(64, 0.));
	//8 body points
	double w[2] = { 2. / 3., 1. / 3. };
	double a[8] = { w[0] * w[0] * w[0], w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[0] * w[0],
		w[1] * w[0] * w[0], w[1] * w[1] * w[0], w[1] * w[1] * w[1], w[1] * w[1] * w[0] };
	int bpi[8][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },{ 1, 2, 3, 0, 5, 6, 7, 4 },{ 2, 3, 0, 1, 6, 7, 4, 5 },{ 3, 0, 1, 2, 7, 4, 5, 6 },
	{ 4, 5, 6, 7, 0, 1, 2, 3 },{ 5, 6, 7, 4, 1, 2, 3, 0 },{ 6, 7, 4, 5, 2, 3, 0, 1 },{ 7, 4, 5, 6, 3, 0, 1, 2 } };
	vector<array<array<double, 8>, 8>> bpm(hx1r.size());
	vector<array<array<int, 8>, 8>> bpmap(hx1r.size());
	for (i = 0; i < hx1r.size(); i++)//which element
	{
		for (j = 0; j < 8; j++)//which body point, bezier
		{
			for (k = 0; k < 8; k++)//which local corner point, b-splines
			{
				bpm[i][j][k] = a[k];
				bpmap[i][j][k] = loc[tmesh[hx1r[i]].cnct[bpi[j][k]]];
			}
		}
	}
	int layer[4] = { 0, 16, 32, 48 };
	int bpbz[8] = { 5 + layer[1], 6 + layer[1], 10 + layer[1], 9 + layer[1], 5 + layer[2], 6 + layer[2], 10 + layer[2], 9 + layer[2] };
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			tmesh[eid].bemat[bpmap[0][i][j]][bpbz[i]] = bpm[0][i][j];
		}
	}
	int pos1, pos2;
	double af[4] = { w[0] * w[0], w[1] * w[0], w[1] * w[1], w[1] * w[0] };
	int f_cnct[4][4] = { { 0, 1, 2, 3 },{ 1, 2, 3, 0 },{ 2, 3, 0, 1 },{ 3, 0, 1, 2 } };
	int fci[6][4] = { { 0, 1, 3, 2 },{ 0, 1, 4, 5 },{ 1, 2, 5, 6 },{ 3, 2, 7, 6 },{ 0, 3, 4, 7 },{ 4, 5, 7, 6 } };
	int fcbz[6][4] = { { 5, 6, 9, 10 },{ 1 + layer[1], 2 + layer[1], 1 + layer[2], 2 + layer[2] },{ 7 + layer[1], 11 + layer[1], 7 + layer[2], 11 + layer[2] },
	{ 13 + layer[1], 14 + layer[1], 13 + layer[2], 14 + layer[2] },{ 4 + layer[1], 8 + layer[1], 4 + layer[2], 8 + layer[2] },{ 5 + layer[3], 6 + layer[3], 9 + layer[3], 10 + layer[3] } };
	//find all boundary faces
	vector<int> fc_b;
	vector<array<array<double, 4>, 4>> fpm;
	vector<array<array<int, 4>, 4>> fpmap;//reference is face, note that reference is not solid
	for (i = 0; i < hx1r.size(); i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (tmface[tmesh[hx1r[i]].face[j]].hex.size() == 1)
			{
				int fcid(tmesh[hx1r[i]].face[j]);
				fc_b.push_back(fcid);
				array<array<double, 4>, 4> fpm_tmp;
				array<array<int, 4>, 4> fpmap_tmp;
				for (int j1 = 0; j1 < 4; j1++)//Bezier
				{
					for (k = 0; k < 4; k++)//B-splines
					{
						fpm_tmp[j1][k] = af[k];
						fpmap_tmp[j1][k] = loc[tmface[fcid].cnct[f_cnct[j1][k]]];
					}
				}
				fpm.push_back(fpm_tmp);
				fpmap.push_back(fpmap_tmp);
			}
		}
	}
	//determine coefs
	for (i = 0; i < 6; i++)
	{
		uint nhex = tmface[tmesh[eid].face[i]].hex.size();
		if (nhex == 2)
		{
			for (j = 0; j < nhex; j++)
			{
				hxid = tmface[tmesh[eid].face[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				for (int j1 = 0; j1 < 4; j1++)
				{
					int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[fci[i][j1]]);
					pos2 = it1 - tmesh[hxid].cnct;
					for (k = 0; k < 8; k++)
					{
						tmesh[eid].bemat[bpmap[pos1][pos2][k]][fcbz[i][j1]] += bpm[pos1][pos2][k] / nhex;
					}
				}
			}
		}
		else if (nhex == 1)//new, for boundary face
		{
			int fcid = tmesh[eid].face[i];
			vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
			pos1 = it - fc_b.begin();
			for (int j1 = 0; j1 < 4; j1++)
			{
				int* it1 = find(tmface[fcid].cnct, tmface[fcid].cnct + 4, tmesh[eid].cnct[fci[i][j1]]);
				pos2 = it1 - tmface[fcid].cnct;
				for (k = 0; k < 4; k++)
				{
					tmesh[eid].bemat[fpmap[pos1][pos2][k]][fcbz[i][j1]] = fpm[pos1][pos2][k];
				}
			}
		}
	}
	//2*12 edge points
	int edi[12][2] = { { 0, 1 },{ 1, 2 },{ 2, 3 },{ 3, 0 },{ 0, 4 },{ 1, 5 },{ 2, 6 },{ 3, 7 },{ 4, 5 },{ 5, 6 },{ 6, 7 },{ 7, 4 } };
	int edbz[12][2] = { { 1, 2 },{ 7, 11 },{ 14, 13 },{ 8, 4 },{ 0 + layer[1], 0 + layer[2] },{ 3 + layer[1], 3 + layer[2] },
	{ 15 + layer[1], 15 + layer[2] },{ 12 + layer[1], 12 + layer[2] },{ 1 + layer[3], 2 + layer[3] },{ 7 + layer[3], 11 + layer[3] },{ 14 + layer[3], 13 + layer[3] },{ 8 + layer[3], 4 + layer[3] } };
	for (i = 0; i < 12; i++)
	{
		if (tmedge[tmesh[eid].edge[i]].type != 1)//non-boundary
		{
			uint nhex = tmedge[tmesh[eid].edge[i]].hex.size();
			for (j = 0; j<tmedge[tmesh[eid].edge[i]].hex.size(); j++)
			{
				hxid = tmedge[tmesh[eid].edge[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[edi[i][0]]);
				pos2 = it1 - tmesh[hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][0]] += bpm[pos1][pos2][k] / nhex;
				}
				int* it2 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[edi[i][1]]);
				pos2 = it2 - tmesh[hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					tmesh[eid].bemat[bpmap[pos1][pos2][k]][edbz[i][1]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
		else if (tmedge[tmesh[eid].edge[i]].type == 1 && tmedge[tmesh[eid].edge[i]].sharp == 0)//boundary, non-sharp
		{
			int nfc_b(0);
			for (j = 0; j < tmedge[tmesh[eid].edge[i]].face.size(); j++)
			{
				if (tmface[tmedge[tmesh[eid].edge[i]].face[j]].type == 1) nfc_b++;
			}
			for (j = 0; j<tmedge[tmesh[eid].edge[i]].face.size(); j++)
			{
				int fcid(tmedge[tmesh[eid].edge[i]].face[j]);
				if (tmface[fcid].type == 1)//boundary face
				{
					vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
					pos1 = it - fc_b.begin();
					for (int j1 = 0; j1 < 2; j1++)
					{
						int* it1 = find(tmface[fcid].cnct, tmface[fcid].cnct + 4, tmesh[eid].cnct[edi[i][j1]]);
						pos2 = it1 - tmface[fcid].cnct;
						for (k = 0; k < 4; k++)
						{
							tmesh[eid].bemat[fpmap[pos1][pos2][k]][edbz[i][j1]] += fpm[pos1][pos2][k] / nfc_b;
						}
					}
				}
			}
		}
		else if (tmedge[tmesh[eid].edge[i]].type == 1 && tmedge[tmesh[eid].edge[i]].sharp == 1)//boundary, sharp edge
		{
			int bid[2] = { loc[tmesh[eid].cnct[edi[i][0]]], loc[tmesh[eid].cnct[edi[i][1]]] };
			tmesh[eid].bemat[bid[0]][edbz[i][0]] = w[0];
			tmesh[eid].bemat[bid[0]][edbz[i][1]] = w[1];
			tmesh[eid].bemat[bid[1]][edbz[i][0]] = w[1];
			tmesh[eid].bemat[bid[1]][edbz[i][1]] = w[0];
		}
	}
	//8 corner points
	int cnbz[8] = { 0, 3, 15, 12, 0 + layer[3], 3 + layer[3], 15 + layer[3], 12 + layer[3] };
	for (i = 0; i < 8; i++)
	{
		if (cp[tmesh[eid].cnct[i]].type != 1)//non-boundary vertex points
		{
			uint nhex = cp[tmesh[eid].cnct[i]].hex.size();
			for (j = 0; j<nhex; j++)
			{
				hxid = cp[tmesh[eid].cnct[i]].hex[j];
				vector<int>::iterator it = find(hx1r.begin(), hx1r.end(), hxid);
				pos1 = it - hx1r.begin();
				int* it1 = find(tmesh[hxid].cnct, tmesh[hxid].cnct + 8, tmesh[eid].cnct[i]);
				pos2 = it1 - tmesh[hxid].cnct;
				for (k = 0; k < 8; k++)
				{
					tmesh[eid].bemat[bpmap[pos1][pos2][k]][cnbz[i]] += bpm[pos1][pos2][k] / nhex;
				}
			}
		}
		else if (cp[tmesh[eid].cnct[i]].type == 1 && cp[tmesh[eid].cnct[i]].sharp == 0)//boundary vertex points, non-sharp
		{
			int nfc_b(0);
			for (j = 0; j < cp[tmesh[eid].cnct[i]].face.size(); j++)
			{
				if (tmface[cp[tmesh[eid].cnct[i]].face[j]].type == 1) nfc_b++;
			}
			for (j = 0; j<cp[tmesh[eid].cnct[i]].face.size(); j++)
			{
				int fcid(cp[tmesh[eid].cnct[i]].face[j]);
				if (tmface[fcid].type == 1)//boundary face
				{
					vector<int>::iterator it = find(fc_b.begin(), fc_b.end(), fcid);
					pos1 = it - fc_b.begin();
					int* it1 = find(tmface[fcid].cnct, tmface[fcid].cnct + 4, tmesh[eid].cnct[i]);
					pos2 = it1 - tmface[fcid].cnct;
					for (k = 0; k < 4; k++)
					{
						tmesh[eid].bemat[fpmap[pos1][pos2][k]][cnbz[i]] += fpm[pos1][pos2][k] / nfc_b;
					}
				}
			}
		}
		else if (cp[tmesh[eid].cnct[i]].type == 1 && cp[tmesh[eid].cnct[i]].sharp == 1)//boundary vertex points, sharp edge
		{
			//tmesh[eid].bemat[loc[tmesh[eid].cnct[i]]][cnbz[i]] = 1.;
			vector<int> shp_ed;
			for (j = 0; j < cp[tmesh[eid].cnct[i]].edge.size(); j++)
			{
				if (tmedge[cp[tmesh[eid].cnct[i]].edge[j]].sharp == 1)
				{
					shp_ed.push_back(cp[tmesh[eid].cnct[i]].edge[j]);
				}
			}
			if (shp_ed.size() == 2)
			{
				for (j = 0; j < shp_ed.size(); j++)
				{
					if (tmedge[shp_ed[j]].pt[0] == tmesh[eid].cnct[i])
					{
						tmesh[eid].bemat[loc[tmedge[shp_ed[j]].pt[0]]][cnbz[i]] = 2. / 3.;
						tmesh[eid].bemat[loc[tmedge[shp_ed[j]].pt[1]]][cnbz[i]] = 1. / 6.;
					}
					else
					{
						tmesh[eid].bemat[loc[tmedge[shp_ed[j]].pt[1]]][cnbz[i]] = 2. / 3.;
						tmesh[eid].bemat[loc[tmedge[shp_ed[j]].pt[0]]][cnbz[i]] = 1. / 6.;
					}
				}
			}
			else
			{
				cerr << "# of sharp edges of the point is wrong!\n";
				cout << shp_ed.size() << "\n";
				getchar();
			}
		}
		else if (cp[tmesh[eid].cnct[i]].type == 1 && cp[tmesh[eid].cnct[i]].sharp == 2)//boundary vertex points, sharp corner
		{
			tmesh[eid].bemat[loc[tmesh[eid].cnct[i]]][cnbz[i]] = 1.;
		}
	}
}

void kernel::OutputMesh_Coupling(const vector<BezierElement3D>& bzmesh, string fn)
{
	int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	string fname = fn + "bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nBezier mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8 * bzmesh.size() << " float\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << bzmesh.size() << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < bzmesh.size(); i++)
		{
			fout << bzmesh[i].type << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	string fname3(fn + "bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (uint l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l] + 1;
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname3 << '\n';
	}

	string fname1(fn + "cmat.txt");
	//ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{

			//fout << i << " " << bzmesh[i].IENb.size() << " " << bzmesh[i].IEN.size() << "\n";
			//for (uint l = 0; l < bzmesh[i].IENb.size(); l++)
			//{
			//	fout << bzmesh[i].IENb[l];
			//	if (l == bzmesh[i].IENb.size() - 1)
			//	{
			//		fout << "\n";
			//	}
			//	else
			//	{
			//		fout << " ";
			//	}
			//}
			fout << i << " " << bzmesh[i].IEN.size() << " " << bzmesh[i].type << "\n";
			for (uint l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l];
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
			for (uint j = 0; j < bzmesh[i].cmat.size(); j++)
			{
				for (uint k = 0; k < bzmesh[i].cmat[j].size(); k++)
				{
					fout << bzmesh[i].cmat[j][k];
					if (k == bzmesh[i].cmat[j].size() - 1)
					{
						fout << "\n";
					}
					else
					{
						fout << " ";
					}
				}
			}

		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << '\n';
	}

	string fname2(fn + "bzpt.txt");
	//ofstream fout;
	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() * 64 << "\n";
		for (uint i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname2 << '\n';
	}
	//string fname1(fn + "_mp.vtk");
	////ofstream fout;
	//fout.open(fname1.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nMeshfree nodes\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << mp.size() << " float\n";
	//	for (uint i = 0; i<mp.size(); i++)
	//	{
	//		fout << mp[i].coor[0] << " " << mp[i].coor[1] << " " << mp[i].coor[2] << "\n";
	//	}
	//	fout << "\nCELLS " << mp.size() << " " << 2 * mp.size() << '\n';
	//	for (uint i = 0; i<mp.size(); i++)
	//	{
	//		fout << "1 " << i << "\n";
	//	}
	//	fout << "\nCELL_TYPES " << mp.size() << '\n';
	//	for (uint i = 0; i < mp.size(); i++)
	//	{
	//		fout << "1\n";
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cerr << "Can't open " << fname1 << '\n';
	//}
}

void kernel::InitializeMeshLabel(string fn)
{
	//read hex vtk
	string fname(fn + "controlmesh.vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			cp[i].act = 1;
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2] ;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		for (int i = 0; i<neles+4; i++) getline(fin, stmp);//skip lines
		for (int i = 0; i<npts; i++)	fin >> cp[i].label;
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}


	RescaleDomain();

	InitialConnect();
	//SetSharpFeature_1();

	//string fn1("../io/global3/rod1");
	//OutputEdge(fn1);
	//OutputCM(fn1);
	//cout << "done setting sharp feature\n";
	//getchar();
}

void kernel::RescaleDomain()
{
	double tmp(1.e6);
	double x_range[3][2] = { { tmp, -tmp },{ tmp, -tmp },{ tmp, -tmp } };
	uint i, j;
	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (cp[i].coor[j] < x_range[j][0]) x_range[j][0] = cp[i].coor[j];
			if (cp[i].coor[j] > x_range[j][1]) x_range[j][1] = cp[i].coor[j];
		}
	}
	double xh[3] = { x_range[0][1] - x_range[0][0], x_range[1][1] - x_range[1][0], x_range[2][1] - x_range[2][0] };
	double dim_min(1.e6);
	for (i = 0; i < 3; i++)
	{
		if (xh[i] < dim_min) dim_min = xh[i];
	}
	//xh[0] /= dim_min; xh[1] /= dim_min; xh[2] /= dim_min;

	for (i = 0; i <cp.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			//cp[i].coor[j] = (cp[i].coor[j] - x_range[j][0]) / xh[j];
			cp[i].coor[j] = (cp[i].coor[j] - x_range[j][0]) / dim_min;
		}
	}
}

void kernel::InitialConnect()
{
	uint i, j;
	tmedge.clear();
	tmface.clear();
	BuildInitialEdges();

	//find BC face, edge, vertex

	//int ed0[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };//order could be wrong, but doesn't matter
	for (i = 0; i<tmface.size(); i++)
	{
		if (tmface[i].hex.size() == 1)
		{
			tmface[i].type = 1;
			tmesh[tmface[i].hex[0]].type = 1;
			for (j = 0; j<4; j++)
			{
				cp[tmface[i].cnct[j]].type = 1;
				tmedge[tmface[i].edge[j]].type = 1;
			}
		}
	}
	//additional boundary elements
	for (i = 0; i<tmesh.size(); i++)
	{
		if (tmesh[i].type != 1)
		{
			for (j = 0; j<8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 1)
				{
					tmesh[i].type = 1;
					break;
				}
			}
		}
	}
	//find extraordinary edges and vertices
	for (i = 0; i<tmedge.size(); i++)
	{
		if (tmedge[i].type != 1 && tmedge[i].hex.size() != 4)
		{
			tmedge[i].type = 2;
			if (cp[tmedge[i].pt[0]].type != 1)
				cp[tmedge[i].pt[0]].type = 3;
			if (cp[tmedge[i].pt[1]].type != 1)
				cp[tmedge[i].pt[1]].type = 3;
		}
	}
	
	//find irregular elements
	for (i = 0; i<tmesh.size(); i++)
	{
		if (tmesh[i].type != 1)
		{
			for (j = 0; j<12; j++)
			{
				if (tmedge[tmesh[i].edge[j]].type == 2)
				{
					tmesh[i].type = 2;
					break;
				}
			}
			//additional
			for (j = 0; j<8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 3)
				{
					tmesh[i].type = 2;
					break;
				}
			}
		}
	}

	//boundry extraordinary points
	for (i = 0; i<cp.size(); i++)
	{
		if (cp[i].type == 1)
		{
			int count(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].type == 2) count++;
			}
			if (count == 1) cp[i].bcxp = 1;
			else if (count>1) cp[i].bcxp = 2;
		}
	}

	//check 3D EP
	int n3d(0);
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 3)
		{
			int ned(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].type == 2)
				{
					ned++;
				}
			}
			if (ned>2)
			{
				n3d++;
			}
		}
	}
	cout << "# 3D EP: " << n3d << "\n";
	//getchar();
}

void kernel::BuildInitialEdges()
{
	int edloc[12][2] = { { 0,1 },{ 1,2 },{ 2,3 },{ 3,0 },{ 0,4 },{ 1,5 },{ 2,6 },{ 3,7 },{ 4,5 },{ 5,6 },{ 6,7 },{ 7,4 } };
	int fcloc[6][4] = { { 0,3,2,1 },{ 0,1,5,4 },{ 1,2,6,5 },{ 2,3,7,6 },{ 0,4,7,3 },{ 4,5,6,7 } };
	int fced[6][4] = { { 3,2,1,0 },{ 0,5,8,4 },{ 1,6,9,5 },{ 2,7,10,6 },{ 4,11,7,3 },{ 8,9,10,11 } };

	uint i, j, k;
	tmedge.clear();
	tmface.clear();
	//point-hex relation
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 8; j++)
		{
			cp[tmesh[i].cnct[j]].hex.push_back(i);
		}
	}
	//construct edges
	for (i = 0; i<tmesh.size(); i++)
	{
		vector<int> nb;
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < cp[tmesh[i].cnct[j]].hex.size(); k++)
			{
				int eid(cp[tmesh[i].cnct[j]].hex[k]);
				if (eid < i)
				{
					vector<int>::iterator it = find(nb.begin(), nb.end(), eid);
					if (it == nb.end())
					{
						nb.push_back(eid);
					}
				}
			}
		}
		for (j = 0; j < 12; j++)//edge
		{
			Edge3D edtmp;
			edtmp.pt[0] = tmesh[i].cnct[edloc[j][0]];
			edtmp.pt[1] = tmesh[i].cnct[edloc[j][1]];
			int flag(-1);
			for (k = 0; k < nb.size(); k++)
			{
				for (int k0 = 0; k0 < 12; k0++)
				{
					if (edtmp == tmedge[tmesh[nb[k]].edge[k0]])
					{
						flag = tmesh[nb[k]].edge[k0]; break;
					}
				}
				if (flag != -1) break;
			}
			if (flag != -1)
			{
				tmesh[i].edge[j] = flag;
			}
			else
			{
				tmedge.push_back(edtmp);
				tmesh[i].edge[j] = tmedge.size() - 1;
			}
		}
		for (j = 0; j < 6; j++)//face
		{
			Face3D fctmp;
			for (k = 0; k < 4; k++)
			{
				fctmp.cnct[k] = tmesh[i].cnct[fcloc[j][k]];
				fctmp.edge[k] = tmesh[i].edge[fced[j][k]];
			}
			int flag(-1);
			for (k = 0; k < nb.size(); k++)
			{
				for (int k0 = 0; k0 < 6; k0++)
				{
					if (fctmp == tmface[tmesh[nb[k]].face[k0]])
					{
						flag = tmesh[nb[k]].face[k0]; break;
					}
				}
				if (flag != -1) break;
			}
			if (flag != -1)
			{
				tmesh[i].face[j] = flag;
			}
			else
			{
				tmface.push_back(fctmp);
				tmesh[i].face[j] = tmface.size() - 1;
			}
		}
	}

	for (i = 0; i<tmesh.size(); i++)
	{
		for (j = 0; j<12; j++)
		{
			tmedge[tmesh[i].edge[j]].hex.push_back(i);
		}
		for (j = 0; j<6; j++)
		{
			tmface[tmesh[i].face[j]].hex.push_back(i);
		}
	}
	for (i = 0; i<tmface.size(); i++)
	{
		for (j = 0; j<4; j++)
		{
			cp[tmface[i].cnct[j]].face.push_back(i);
			tmedge[tmface[i].edge[j]].face.push_back(i);
		}
	}
	for (i = 0; i<tmedge.size(); i++)
	{
		for (j = 0; j<2; j++)
		{
			cp[tmedge[i].pt[j]].edge.push_back(i);
		}
	}
}