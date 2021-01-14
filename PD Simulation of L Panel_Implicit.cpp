// This program is used to simulate the L-shaped panel test in the implicit framework of peridynamic 
// Author：Chen XU (1536343106@qq.com)
// Version: Early Beta (too terrible，to much parameter tuning)
// Environment configurations：
//		Please firstly install Visual studio2017，then install Inter Fortran Parallel Studio2018 for integration
//		Project→Properties→Inter Performance Library→Use Inter MKL：Parallel（注意配置为所有配置、平台为所有平台）
//		Project→Properties→VC++ Directories→Include Directories：add the file path of Eigen
//		Project→Properties→C / C++→Code generation→Runtime library：multithreading
//		Project→Properties→C / C++→Language→OpenMP support：Yes
//		Configuration Manager:    Release    x64

#include "pch.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <Eigen/PardisoSupport>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "mkl.h"
#include <omp.h>
#include <windows.h>
//not used, i was intended to use these libraries to comfirm whether the tangential stiffness matrix is negative definite
//#include <armadillo>
//#include <Spectra/SymEigsSolver.h>
//#include <Spectra/MatOp/SparseSymMatProd.h>


//inter MKL
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
//π
#define pi 3.1415926536
//"NUM_BOND - 1"（Do not ignore -1）: the number of bonds terminated in every iterative step
#define NUM_BOND 5
#define Infi 100000.0
//tolerence
#define etol 0.00001
//to describe the bilinear bond force-stretch relationship
#define Constitutive 3.0


//////////////////////////////////////////////////////////////////////////////////////////////
//****************************************Stiffness*************************************
//Output: element stiffness matrix k(6×6)
//Input: coordinate of point i and point j，PD parameter--bc、bd
Eigen::Matrix <double, 6, 6> Stiffness(double xi, double yi, double xj, double yj, double bc, double bd)
{
	Eigen::Matrix<double, 6, 6> k1;
	double idist = sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi));
	k1 << bc, 0, 0, -bc, 0, 0,
		0, 12.0 * bd / idist / idist, 6.0 * bd / idist, 0, -12.0 * bd / idist / idist, 6.0 * bd / idist,
		0, 6.0 * bd / idist, 4.0 * bd, 0, -6.0 * bd / idist, 2.0 * bd,
		-bc, 0, 0, bc, 0, 0,
		0, -12.0 * bd / idist / idist, -6.0 * bd / idist, 0, 12.0 * bd / idist / idist, -6.0 * bd / idist,
		0, 6.0 * bd / idist, 2.0 * bd, 0, -6.0 * bd / idist, 4.0 * bd;
	//rotation matrix R
	Eigen::Matrix<double, 6, 6> R;
	double cosalpha = (xj - xi) / idist;
	double sinalpha = (yj - yi) / idist;
	R << cosalpha, sinalpha, 0, 0, 0, 0,
		-sinalpha, cosalpha, 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 0,
		0, 0, 0, cosalpha, sinalpha, 0,
		0, 0, 0, -sinalpha, cosalpha, 0,
		0, 0, 0, 0, 0, 1.0;
	Eigen::Matrix<double, 6, 6> k;
	k = 1.0 / idist * (R.transpose()) * k1 * R;
	return k;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//**************************************PartialStiffness*********************************
//Output: PD element stiffness matrix k(6×6) in the soften stage
//Input: coordinate of point i and point j, PD parameter--bc、bd, bond stretch, damage of bond--μ, critical stretch--set
Eigen::Matrix <double, 6, 6> PartialStiffness(double xi, double yi, double xj, double yj, double bc, double bd, double stretch, double miu, double set)
{
	Eigen::Matrix<double, 6, 6> k1;
	double idist = sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi));
	k1 << bc, 0, 0, -bc, 0, 0,
		0, 12.0 * bd / idist / idist, 6.0 * bd / idist, 0, -12.0 * bd / idist / idist, 6.0 * bd / idist,
		0, 6.0 * bd / idist, 4.0 * bd, 0, -6.0 * bd / idist, 2.0 * bd,
		-bc, 0, 0, bc, 0, 0,
		0, -12.0 * bd / idist / idist, -6.0 * bd / idist, 0, 12.0 * bd / idist / idist, -6.0 * bd / idist,
		0, 6.0 * bd / idist, 2.0 * bd, 0, -6.0 * bd / idist, 4.0 * bd;
	k1 = k1 * miu / idist;
	double pmiu_pu1 = Constitutive * set / (Constitutive - 1) / stretch / stretch / idist;
	k1(0, 0) = k1(0, 0) - bc / idist * pmiu_pu1*idist*stretch;
	k1(0, 3) = k1(0, 3) + bc / idist * pmiu_pu1*idist*stretch;
	k1(3, 0) = k1(3, 0) + bc / idist * pmiu_pu1*idist*stretch;
	k1(3, 3) = k1(3, 3) - bc / idist * pmiu_pu1*idist*stretch;
	Eigen::Matrix<double, 6, 6> R;
	double cosalpha = (xj - xi) / idist;
	double sinalpha = (yj - yi) / idist;
	R << cosalpha, sinalpha, 0, 0, 0, 0,
		-sinalpha, cosalpha, 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 0,
		0, 0, 0, cosalpha, sinalpha, 0,
		0, 0, 0, -sinalpha, cosalpha, 0,
		0, 0, 0, 0, 0, 1.0;
	Eigen::Matrix<double, 6, 6> t;
	t = (R.transpose()) * k1 * R;
	return t;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//****************************************EQsolver*************************************
//Solve quadratic equation:  a1*x*x + a2*x + a3 = 0
double EQsolver(double a1, double a2, double a3)
{
	double delta_lambda;
	double lambda1;
	double lambda2;
	double fac;
	if (a1 == 0 && a2 == 0)
	{
		std::cout << "No roots found" << std::endl;
		exit(0);
	}
	else if (a1 == 0 && a2 != 0)
	{
		delta_lambda = -a3 / a2;
	}
	else
	{
		fac = a2 * a2 - 4 * a1*a3;
		if (fac < 0)
		{
			std::cout << "No roots found" << std::endl;
			exit(0);
		}
		lambda1 = (-a2 + sqrt(fac)) / (2 * a1);
		lambda2 = (-a2 - sqrt(fac)) / (2 * a1);
		if (lambda2 >= lambda1)
		{
			delta_lambda = lambda2;
		}
		else
		{
			delta_lambda = lambda1;
		}
	}
	return delta_lambda;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//****************************************QuickSort************************************
//Matrix s: first column--element number, first column--corresponding bond stretch
void QuickSort(double s[NUM_BOND][2], int l, int r)
{
	if (l < r)
	{
		int i = l, j = r;
		double x_val = s[l][1];
		int x_row = s[l][0];
		while (i < j)
		{
			while (i < j && s[j][1] < x_val)
			{
				j--;
			}
			if (i < j)
			{
				s[i][1] = s[j][1];
				s[i][0] = s[j][0];
				i = i + 1;
			}
			while (i < j && s[i][1] >= x_val) 
			{
				i++;
			}
			if (i < j)
			{
				s[j][1] = s[i][1];
				s[j][0] = s[i][0];
				j = j - 1;
			}
		}
		s[i][1] = x_val;
		s[i][0] = x_row;
		QuickSort(s, l, i - 1); // recursive
		QuickSort(s, i + 1, r);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//************************************Main Function***********************************
int main()
{
	double start = GetTickCount();


	//pre-processing
	//number of points in the X-direction
	const int ndivx = 100;
	//number of points in the Y-direction
	const int ndivy = 100;
	//thickness of boundary
	const int nbnd = 3;
	//total number of material points
	const int totnode = ndivx * ndivy / 4.0 *3.0 + nbnd * ndivx / 2 + 1;
	//total load step
	const int num_loadstep = 70;
	//max number of points in a horizon
	const int maxfam = 40;
	//length of structure
	double length = 0.5;
	//width of structure
	double width = 0.5;
	//thickness of structure
	double thick = 0.1;
	//distance between points
	double dx = length / ndivx;
	//horizon radius
	double delta = 3.015 * dx;
	double area = dx * dx;
	double vol = area * thick;


	//Mesh generation
	//global coordinate matrix
	Eigen::MatrixXd coord = Eigen::MatrixXd::Zero(totnode, 2);
	int nnum = 0;
	for (int i = 1; i <= ndivx; i = i + 1)
	{
		for (int j = 1; j <= ndivy; j = j + 1)
		{
			double coordx = (-1.0 * length / 2.0) + (dx / 2.0) + (i - 1) * dx;
			double coordy = (-1.0 * width / 2.0) + (dx / 2.0) + (j - 1) * dx;
			if (coordx <= 0 || coordy <= 0)
			{
				coord(nnum, 0) = coordx;
				coord(nnum, 1) = coordy;
				nnum = nnum + 1;
			}
		}
	}
	int totint = nnum - 1;
	for (int i = 1; i <= ndivx / 2; i = i + 1)
	{
		for (int j = 1; j <= nbnd; j = j + 1)
		{
			coord(nnum, 0) = -1.0 / 2.0 * length + dx / 2.0 + (i - 1) * dx + 0 * dx;
			coord(nnum, 1) = 1.0 / 2.0 * width + dx / 2.0 + (j - 1) * dx;
			nnum = nnum + 1;
		}
	}
	int totboundary = nnum - 1;
	for (int i = 1; i <= 1; i = i + 1)
	{
		for (int j = 1; j <= 1; j = j + 1)
		{
			coord(nnum, 0) = 1.0 / 2.0 * length - dx / 2.0 - (i - 1) * dx - 1 * dx;
			coord(nnum, 1) = 0.0 / 2.0 * width + dx / 2.0 + (j - 1) * dx;
			nnum = nnum + 1;
		}
	}
	int totload = nnum - 1;


	//determine the horizon of each point
	Eigen::VectorXi numfam = Eigen::VectorXi::Zero(totnode);
	Eigen::VectorXi pointfam = Eigen::VectorXi::Zero(totnode);
	Eigen::VectorXi nodefam = Eigen::VectorXi::Zero(maxfam*totnode);
	Eigen::VectorXi uninumfam = Eigen::VectorXi::Zero(totnode);
	for (int i = 0; i < totnode; i = i + 1)
	{
		if (i == 0)
		{
			pointfam(i) = 0;
		}
		else
		{
			pointfam(i) = pointfam(i - 1) + numfam(i - 1);
		}
		for (int j = 0; j < totnode; j = j + 1)
		{
			double idist = sqrt(pow(coord(j, 0) - coord(i, 0), 2) + pow(coord(j, 1) - coord(i, 1), 2));
			if (i != j)
			{
				if (idist <= delta)
				{
					numfam(i) = numfam(i) + 1;
					nodefam(pointfam(i) + numfam(i) - 1) = j;
				}
			}
		}
	}
	int nelem = 0; //total number of elements
	for (int i = 0; i < totnode; i = i + 1)
	{
		nelem = nelem + numfam(i);
	}
	nelem = nelem / 2;
	Eigen::MatrixXd element = Eigen::MatrixXd::Zero(nelem, 2);
	int jj = 0;
	double tempnode;
	for (int i = 0; i < totnode; i = i + 1)
	{
		for (int j = 0; j < numfam(i); j = j + 1)
		{
			tempnode = nodefam(pointfam(i) + j);
			if (tempnode > i)
			{
				element(jj, 0) = i;
				element(jj, 1) = tempnode;
				jj = jj + 1;
			}
		}
	}
	if (jj != nelem)
	{
		std::cout << "Element matrix initialization error！" << std::endl;
		exit(0);
	}
	std::cout << "Mesh generation completed！" << std::endl;


	//concrete : material parameter
	//c_emod : Elastic modulus
	double c_emod = 25.85e9;
	//c_pratio : Poisson's ratio
	double c_pratio = 0.18;
	//c_bc : Bond constant 1
	double c_bc = 6.0 * c_emod / (pi * thick * (1.0 - c_pratio) * pow(delta, 3));
	//c_bd : Bond constant 2
	double c_bd = (1.0 - 3 * c_pratio) * c_emod / (6.0 * pi * thick * delta * (1.0 - c_pratio * c_pratio));
	//c_set : Critical stretch
	double c_set = 2.7 / 25850.0;
	//c_sut : Ultimate stretch
	double c_sut = Constitutive * c_set;


	//initial bond length
	Eigen::VectorXd oridist = Eigen::VectorXd::Zero(nelem);
	for (int i = 0; i < nelem; i = i + 1)
	{
		oridist(i) = sqrt(pow(coord((element(i, 1)), 0) - coord((element(i, 0)), 0), 2) + pow(coord((element(i, 1)), 1) - coord((element(i, 0)), 1), 2));
	}
	

	Eigen::VectorXd hori_correct = Eigen::VectorXd::Zero(nelem);
	for (int i = 0; i < nelem; i = i + 1)
	{
		if (oridist(i) <= (delta - dx / 2.0))
		{
			hori_correct(i) = 1.0;
		}
		else if (oridist(i) <= delta)
		{
			hori_correct(i) = (delta + dx / 2.0 - oridist(i)) / dx;
		}
		else
		{
			hori_correct(i) = 0.0;
		}
	}


	//arc-length
	double deltaL = 200;
	//global force vector FFG
	Eigen::VectorXd FFG = Eigen::VectorXd::Zero(3 * totnode);
	Eigen::VectorXd FFGREF = Eigen::VectorXd::Zero(3 * totnode);
	for (int i = totboundary + 1; i <= totload; i = i + 1)
	{
		FFGREF(3 * i + 1) = -deltaL;
	}
	//unbalanced force vector RRG
	Eigen::VectorXd RRG = Eigen::VectorXd::Zero(3 * totnode);
	//RRGOLD: to record RRG
	Eigen::VectorXd RRGOLD = Eigen::VectorXd::Zero(3 * totnode);
	//global displacement vector UUG，iterative increment ΔU、ΔΔU1、ΔΔU2
	Eigen::VectorXd UUG = Eigen::VectorXd::Zero(3 * totnode);
	Eigen::VectorXd DUUG = Eigen::VectorXd::Zero(3 * totnode);
	Eigen::VectorXd DDUUG1 = Eigen::VectorXd::Zero(3 * totnode);
	Eigen::VectorXd DDUUG2 = Eigen::VectorXd::Zero(3 * totnode);
	double Lambda = 0.0;
	double DLambda;
	double DDLambda;
	//bond stretch
	Eigen::VectorXd stretch = Eigen::VectorXd::Zero(nelem);
	//global stiffness matrix KKG
	Eigen::SparseMatrix<double> KKG(3 * totnode, 3 * totnode);
	Eigen::SparseMatrix<double> KKGOLD(3 * totnode, 3 * totnode);
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(6, 6);
	//damage of material points
	Eigen::VectorXd dmg = Eigen::VectorXd::Zero(totnode);
	//damage of bonds
	Eigen::VectorXd fail = Eigen::VectorXd::Zero(nelem);
	for (int i = 0; i < nelem; i = i + 1)
	{
		fail(i) = 1.0;
	}


	//volume correction
	Eigen::VectorXd volume = Eigen::VectorXd::Zero(totnode);
	double idist;
	double cnode;
	double surf_correct;
	for (int i = 0; i < totnode; i = i + 1)
	{
		for (int j = 0; j < numfam(i); j = j + 1)
		{
			cnode = nodefam(pointfam(i) + j);
			idist = sqrt(pow(coord(cnode, 0) - coord(i, 0), 2) + pow(coord(cnode, 1) - coord(i, 1), 2));
			if (idist <= (delta - dx / 2.0))
			{
				surf_correct = 1.0;
			}
			else if (idist <= delta)
			{
				surf_correct = (delta + dx / 2.0 - idist) / dx;
			}
			else
			{
				surf_correct = 0.0;
			}
			volume(i) = volume(i) + vol * surf_correct;
		}
		volume(i) = volume(i) + vol;
	}
	//surface correction
	Eigen::VectorXd vol_correct = Eigen::VectorXd::Zero(nelem);
	for (int i = 0; i < nelem; i = i + 1)
	{
		vol_correct(i) = 2 * pi*delta*delta*thick / (volume(element(i, 0)) + volume(element(i, 1)));
	}


	Eigen::VectorXd stiffpara = Eigen::VectorXd::Zero(nelem);
	for (int i = 0; i < nelem; i = i + 1)
	{
		stiffpara(i) = hori_correct(i) * vol_correct(i) * vol * vol;
	}
	std::cout << "Pre-processing completed！" << std::endl;


	//************************************Computation************************************
	std::cout << "Computation started！" << std::endl;
	//mark of unloading, 0--loading, 1--unloading
	int flag = 0;
	//apply load step
	for (int loadstep = 1; loadstep <= num_loadstep; loadstep = loadstep + 1)
	{
		//reset ΔU and Δλ
		DUUG = Eigen::VectorXd::Zero(3 * totnode);
		DLambda = 1.0;
		if (loadstep >= 30 && flag == 0)
		{
			//Spectra and Armadillo can also be helpful to judge whether a matrix is negative definite
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> lltOfA(KKGOLD);
			if (lltOfA.info() == Eigen::NumericalIssue)
			{
				std::cout << "Possibly non semi-positive definitie matrix!" << std::endl;
				flag = 1;
			}
		}
		if (flag == 1)
			DLambda = -1.0;
		Lambda += DLambda;


		int dof[6];
		double tol = 1.0;
		//apply iterative step
		for (int iterstep = 1; tol > etol && iterstep <= 7; iterstep = iterstep + 1)
		{
			//reset ΔΔU and ΔΔλ
			DDLambda = 0.0;
			DDUUG1 = Eigen::VectorXd::Zero(3 * totnode);
			DDUUG2 = Eigen::VectorXd::Zero(3 * totnode);
			//assembly global stiffness matrix
			std::vector < Eigen::Triplet < double > > tripletList;
			for (int i = 0; i < nelem; i = i + 1)
			{
				dof[0] = 3 * element(i, 0) + 0;
				dof[1] = 3 * element(i, 0) + 1;
				dof[2] = 3 * element(i, 0) + 2;
				dof[3] = 3 * element(i, 1) + 0;
				dof[4] = 3 * element(i, 1) + 1;
				dof[5] = 3 * element(i, 1) + 2;
				if (fail(i) != 0.0)
				{
					if (fail(i) == 1.0)
					{
						k = stiffpara(i) * Stiffness(coord((element(i, 0)), 0), coord((element(i, 0)), 1), coord((element(i, 1)), 0), coord((element(i, 1)), 1), c_bc, c_bd);
					}
					else
					{
						k = stiffpara(i) * PartialStiffness(coord((element(i, 0)), 0), coord((element(i, 0)), 1), coord((element(i, 1)), 0), coord((element(i, 1)), 1), c_bc, c_bd, stretch(i), fail(i), c_set);
					}

					for (int p = 0; p < 6; p = p + 1)
					{
						for (int q = 0; q < 6; q = q + 1)
						{
							if (k(p, q) != 0.0)
							{

								tripletList.push_back(Eigen::Triplet<double>(dof[p], dof[q], k(p, q)));
							}
						}
					}
				}
			}
			KKG.setFromTriplets(tripletList.begin(), tripletList.end());
			KKGOLD = KKG;


			//calculate the unbalanced force
			RRG = Eigen::VectorXd::Zero(3 * totnode);
			FFG = Eigen::VectorXd::Zero(3 * totnode);
			for (int i = 0; i < nelem; i = i + 1)
			{
				dof[0] = 3 * element(i, 0) + 0;
				dof[1] = 3 * element(i, 0) + 1;
				dof[2] = 3 * element(i, 0) + 2;
				dof[3] = 3 * element(i, 1) + 0;
				dof[4] = 3 * element(i, 1) + 1;
				dof[5] = 3 * element(i, 1) + 2;
				if (fail(i) != 0.0)
				{
					k = fail(i)*stiffpara(i) * Stiffness(coord((element(i, 0)), 0), coord((element(i, 0)), 1), coord((element(i, 1)), 0), coord((element(i, 1)), 1), c_bc, c_bd);
				}
				else
				{
					k = Eigen::MatrixXd::Zero(6, 6);
				}
				Eigen::VectorXd u = Eigen::VectorXd::Zero(6);
				u(0) = UUG(dof[0]);
				u(1) = UUG(dof[1]);
				u(2) = UUG(dof[2]);
				u(3) = UUG(dof[3]);
				u(4) = UUG(dof[4]);
				u(5) = UUG(dof[5]);
				Eigen::VectorXd fi = Eigen::VectorXd::Zero(6);
				fi = k * u;
				RRG(dof[0]) = RRG(dof[0]) + fi(0);
				RRG(dof[1]) = RRG(dof[1]) + fi(1);
				RRG(dof[2]) = RRG(dof[2]) + fi(2);
				RRG(dof[3]) = RRG(dof[3]) + fi(3);
				RRG(dof[4]) = RRG(dof[4]) + fi(4);
				RRG(dof[5]) = RRG(dof[5]) + fi(5);
			}
			RRGOLD = -(RRG - FFG);
			FFG = FFGREF * Lambda;
			RRG = -(RRG - FFG);

			//introduce the boundary condition
#pragma omp parallel for
			for (int j = 0; j < KKG.outerSize(); ++j)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(KKG, j); it; ++it)
				{
					for (int i = totint + 1; i <= totboundary; i = i + 1)
					{
						if (it.col() == (3 * i + 0))
						{
							RRG(it.row()) = RRG(it.row()) - (0.0 - DUUG(3 * i + 0)) * it.value();
							FFG(it.row()) = FFG(it.row()) - (0.0 - DUUG(3 * i + 0)) * it.value();
						}
						if (it.col() == (3 * i + 1))
						{
							RRG(it.row()) = RRG(it.row()) - (0.0 - DUUG(3 * i + 1)) * it.value();
							FFG(it.row()) = FFG(it.row()) - (0.0 - DUUG(3 * i + 1)) * it.value();
						}
					}
				}
			}


#pragma omp parallel for
			for (int j = 0; j < KKG.outerSize(); ++j)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(KKG, j); it; ++it)
				{
					for (int i = totint + 1; i <= totboundary; i = i + 1)
					{
						if (it.row() == (3 * i + 0) && it.col() != (3 * i + 0))
						{
							KKG.coeffRef(it.row(), it.col()) = 0.0;
						}
						if (it.row() != (3 * i + 0) && it.col() == (3 * i + 0))
						{
							KKG.coeffRef(it.row(), it.col()) = 0.0;
						}
						if (it.row() == (3 * i + 0) && it.col() == (3 * i + 0))
						{
							KKG.coeffRef(it.row(), it.col()) = 1.0;
						}

						if (it.row() == (3 * i + 1) && it.col() != (3 * i + 1))
						{
							KKG.coeffRef(it.row(), it.col()) = 0.0;
						}
						if (it.row() != (3 * i + 1) && it.col() == (3 * i + 1))
						{
							KKG.coeffRef(it.row(), it.col()) = 0.0;
						}
						if (it.row() == (3 * i + 1) && it.col() == (3 * i + 1))
						{
							KKG.coeffRef(it.row(), it.col()) = 1.0;
						}
					}
				}
			}

			
			for (int i = totint + 1; i <= totboundary; i = i + 1)
			{
				RRG(3 * i + 0) = 0.0 - DUUG(3 * i + 0);
				RRG(3 * i + 1) = 0.0 - DUUG(3 * i + 1);
				FFG(3 * i + 0) = 0.0 - DUUG(3 * i + 0);
				FFG(3 * i + 1) = 0.0 - DUUG(3 * i + 1);
			}


			//solved by Pardiso
			Eigen::PardisoLU <Eigen::SparseMatrix<double> > solver;
			solver.compute(KKG);
			DDUUG1 = solver.solve(FFG);
			DDUUG2 = solver.solve(RRG);
			//spherical arc-length method----Not used
			//double a1 = DDUUG1.dot(DDUUG1) + FFGREF.dot(FFGREF);
			//double a2 = 2.0 * DDUUG1.dot(DUUG + DDUUG2) + 2.0 * DLambda*FFGREF.dot(FFGREF);
			//double a3 = (DUUG + DDUUG2).dot(DUUG + DDUUG2) - deltaL * deltaL + DLambda * DLambda*FFGREF.dot(FFGREF);
			//DDLambda = EQsolver(a1, a2, a3);
			DDLambda = -DUUG.dot(DDUUG2) / (DUUG.dot(DDUUG1) + DLambda * FFGREF.dot(FFGREF));
			DUUG = DUUG + DDUUG1 * DDLambda + DDUUG2;
			UUG = UUG + DDUUG1 * DDLambda + DDUUG2;
			DLambda += DDLambda;
			Lambda += DDLambda;


			//number of broken bonds in each iterative step
			int broken = 0;
			//record the limited number of bonds which have the maximum stretch
			double bondcut[NUM_BOND][2] = { 0 };
			for (int i = 0; i < nelem; i = i + 1)
			{
				double nlength = sqrt(pow(coord((element(i, 1)), 0) + UUG(3 * element(i, 1) + 0) - coord((element(i, 0)), 0) - UUG(3 * element(i, 0) + 0), 2) + pow(coord((element(i, 1)), 1) + UUG(3 * element(i, 1) + 1) - coord((element(i, 0)), 1) - UUG(3 * element(i, 0) + 1), 2));
				//set broken bonds' stretch to ∞
				if (stretch(i) < (Infi - 1.0))
				{
					stretch(i) = (nlength - oridist(i)) / oridist(i);
					bondcut[NUM_BOND - 1][1] = stretch(i);
					bondcut[NUM_BOND - 1][0] = i;
					//divide and conquer
					QuickSort(bondcut, 0, NUM_BOND - 1);
				}
			}


			//terminate the recorded bonds
			for (int i = 0; i < NUM_BOND - 1; i = i + 1)
			{
				if (bondcut[i][1] >= c_sut)
				{
					fail(bondcut[i][0]) = 0.0;
					stretch(bondcut[i][0]) = Infi;
					broken = broken + 1;
				}
				if (bondcut[i][1] >= c_set && bondcut[i][1] < c_sut)
				{
					fail(bondcut[i][0]) = c_set / stretch(bondcut[i][0])*(stretch(bondcut[i][0]) - c_sut) / (c_set - c_sut);
				}
			}


			//calculate the tolerence
			if (iterstep >= 2)
			{
				if (broken == 0)
					tol = etol;
				if (loadstep >= 40 && iterstep >= 5)
					tol = etol;
				//in the soften stage, the number of iterative steps should be decreased, or the stiffness will be overestimated.
				//it may be weird but works, if you have better strategy, please email me!
			}
			std::cout << "loadstep" << loadstep << " iterstep" << iterstep << " tolerence is " << tol << " broken " << broken << std::endl;
		}


		//caculate the damage of material points
		dmg = Eigen::VectorXd::Zero(totnode);
		for (int i = 0; i < nelem; i = i + 1)
		{
			if (fail(i) == 0.0)
			{
				dmg(element(i, 0)) = dmg(element(i, 0)) + 1.0 / numfam(element(i, 0));
				dmg(element(i, 1)) = dmg(element(i, 1)) + 1.0 / numfam(element(i, 1));
			}
		}
		

		//output VTK files which can be post-processing by paraview (please change the output file path by yourself)
		std::string filename = "C:\\Users\\mis\\Desktop\\PD Simulation of L Panel_Implicit\\PD_";
		filename = filename + std::to_string(loadstep) + ".vtk";
		std::ofstream outfile;
		outfile.open(filename);
		outfile << "# vtk DataFile Version 3.0" << std::endl;
		outfile << "PD_Implicit" << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		outfile << std::endl;
		outfile << "POINTS  " << (totload + 1) << "  double" << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << std::setiosflags(std::ios::left) << std::setw(20) << coord(i, 0) << std::setw(20) << coord(i, 1) << std::setw(20) << 0.0 << std::endl;
		}
		outfile << std::endl;
		outfile << "CELLS  " << (totload + 1) << "  " << 2 * (totload + 1) << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << std::setiosflags(std::ios::left) << std::setw(20) << 1 << std::setw(20) << i << std::endl;
		}
		outfile << std::endl;
		outfile << "CELL_TYPES  " << (totload + 1) << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << 1 << std::endl;
		}
		outfile << std::endl;
		outfile << "POINT_DATA  " << (totload + 1) << std::endl;
		outfile << "SCALARS dmg double 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << dmg(i, 0) << std::endl;
		}
		outfile << "VECTORS disp double" << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << std::setiosflags(std::ios::left) << std::setw(20) << UUG(3 * i + 0) << std::setw(20) << UUG(3 * i + 1) << std::setw(20) << 0.0 << std::endl;
		}
		outfile << "VECTORS force double" << std::endl;
		for (int i = 0; i <= totload; i = i + 1)
		{
			outfile << std::setiosflags(std::ios::left) << std::setw(20) << RRGOLD(3 * i + 0) << std::setw(20) << RRGOLD(3 * i + 1) << std::setw(20) << 0.0 << std::endl;
		}


		//load-displacement curve
		double displacement = 0.0;
		double force = 0.0;
		for (int i = totboundary + 1; i <= totload; i = i + 1)
		{
			displacement -= UUG(3 * i + 1);
			force += RRGOLD(3 * i + 1);
		}
		std::ofstream in;
		in.open("C:\\Users\\mis\\Desktop\\PD Simulation of L Panel_Implicit\\f-d-curve.txt", std::ios_base::app);
		in << std::setiosflags(std::ios::left) << std::setw(20) << displacement * 1000.0 << std::setw(20) << force / 1000.0 << std::endl;
		in.close();


		std::cout << "loadstep" << loadstep << " calculation completed！" << std::endl;
	}


	double  end = GetTickCount();
	std::cout << "Computational time (min):    " << (end - start) / 60000 << std::endl;
	return 0;
}