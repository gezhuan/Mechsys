/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/* J. N. Reddy's Finite Element Method:   Example 8.3.1   */

// STL
#include <iostream>
#include <fstream>
#include <cmath>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/tri3diffusion.h"
#include "fem/elems/quad4diffusion.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using boost::make_tuple;
using Util::_4;
using Util::_6;
using Util::_8s;

double u_correct(double s, double k, double x, double y)
{
	// s: source, k:conductivity
	// x,y: position
	size_t nterms = 200;
	double sum    = 0.0;
	for (size_t i=1; i<=nterms; ++i)
	{
		double alp = 0.5*(2.0*i-1.0)*Util::PI;
		sum += pow(-1.0,i)*cos(alp*y)*cosh(alp*x)/(pow(alp,3.0)*cosh(alp));
	}
	return 0.5*s*(1.0-y*y+4.0*sum)/k;
}

int main(int argc, char **argv) try
{
	// Constants
	int    ndiv = 2;     // number of divisions along x and y
	String linsol("LA"); // LAPACK
	Array<double> source(1); source[0] = 1.0;
	
	// Input
	cout << "Input: " << argv[0] << "  ndiv  linsol(LA,UM,SLU)\n";
	if (argc==2) ndiv        = atoi(argv[1]);
	if (argc==3) linsol.Printf("%s",argv[2]);

	double err_Ke_1 = 0.0;
	{
		cout << "\n[0;1m================================================================== Coarse triangular mesh[0m\n";

		////////////////////////////////////////////////////////////////////////////////////////// FEM /////

		// Geometry
		FEM::Geom g(2); // 2D

		// Nodes
		g.SetNNodes (6);
		g.SetNode   (0, 0.0, 0.0);
		g.SetNode   (1, 0.5, 0.0);
		g.SetNode   (2, 0.5, 0.5);
		g.SetNode   (3, 1.0, 0.0);
		g.SetNode   (4, 1.0, 0.5);
		g.SetNode   (5, 1.0, 1.0);

		// Elements
		g.SetNElems (4);
		g.SetElem   (0, "Tri3Diffusion");
		g.SetElem   (1, "Tri3Diffusion");
		g.SetElem   (2, "Tri3Diffusion");
		g.SetElem   (3, "Tri3Diffusion");

		// Set connectivity
		g.Ele(0)->Connect(0, g.Nod(0))->Connect(1, g.Nod(1))->Connect(2, g.Nod(2));
		g.Ele(1)->Connect(0, g.Nod(4))->Connect(1, g.Nod(2))->Connect(2, g.Nod(1));
		g.Ele(2)->Connect(0, g.Nod(1))->Connect(1, g.Nod(3))->Connect(2, g.Nod(4));
		g.Ele(3)->Connect(0, g.Nod(2))->Connect(1, g.Nod(4))->Connect(2, g.Nod(5));

		// Parameters and initial values
		g.Ele(0)->SetModel("LinDiffusion", "k=1.0", "");
		g.Ele(1)->SetModel("LinDiffusion", "k=1.0", "");
		g.Ele(2)->SetModel("LinDiffusion", "k=1.0", "");
		g.Ele(3)->SetModel("LinDiffusion", "k=1.0", "");
		
		// Properties (heat source)
		g.Ele(0)->SetProps(source);
		g.Ele(1)->SetProps(source);
		g.Ele(2)->SetProps(source);
		g.Ele(3)->SetProps(source);

		// Boundary conditions (must be after connectivity)
		g.Ele(0)->EdgeBry("q", 0.0, 0)->EdgeBry("q", 0.0, 2);
		g.Ele(2)->EdgeBry("q", 0.0, 0)->EdgeBry("u", 0.0, 1);
		g.Ele(3)->EdgeBry("u", 0.0, 1)->EdgeBry("q", 0.0, 2);
		g.Nod(3)->Bry("u",0.0);
		g.Nod(5)->Bry("u",0.0);

		// Stiffness
		LinAlg::Matrix<double> Ke0, Ke1, Ke2, Ke3;
		LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(3,3);
		g.Ele(0)->Order1Matrix(0,Ke0);
		g.Ele(1)->Order1Matrix(0,Ke1);
		g.Ele(2)->Order1Matrix(0,Ke2);
		g.Ele(3)->Order1Matrix(0,Ke3);
		Ke_correct =  0.5, -0.5,  0.0,
					 -0.5,  1.0, -0.5,
					  0.0, -0.5,  0.5;

		// Check conductivity matrices
		for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			err_Ke_1 += fabs(Ke0(i,j)-Ke_correct(i,j));
			err_Ke_1 += fabs(Ke1(i,j)-Ke_correct(i,j));
			err_Ke_1 += fabs(Ke2(i,j)-Ke_correct(i,j));
			err_Ke_1 += fabs(Ke3(i,j)-Ke_correct(i,j));
		}

		// Solve
		FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
		sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
		sol -> Solve();
		double norm_resid = LinAlg::Norm(sol->Resid());
		cout << "\n[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
		cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n\n";

		// Output: Nodes
		cout << _6<<"Node #" << _8s<<"u" << _8s<<"q" << endl;
		for (size_t i=0; i<g.NNodes(); ++i)
			cout << _6<<i << _8s<<g.Nod(i)->Val("u") << _8s<<g.Nod(i)->Val("q") << endl;

		//////////////////////////////////////////////////////////////////////////////////////// Check /////

		std::ofstream of("tex831_1.cal", std::ios::out);
		of << _8s<<"x" << _8s<<"u" << _8s<<"ucorr" << endl;
		Array<double> err_u;
		for (size_t i=0; i<g.NNodes(); ++i)	
		{
			double x     = g.Nod(i)->X();
			double y     = g.Nod(i)->Y();
			double u     = g.Nod(i)->Val("u");
			double ucorr = u_correct(1.0,1.0,x,y);
			err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
			if (fabs(y)<=1e-5) of << _8s<<x << _8s<<u << _8s<<ucorr << endl;
		}
		of.close();

		// Error summary
		double tol_u     = 1.0e-2;
		double min_err_u = err_u[err_u.Min()];
		double max_err_u = err_u[err_u.Max()];
		cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
		cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
		cout << endl;
	}

	double err_Ke_2 = 0.0;
	{
		cout << "\n[0;1m================================================================== Fine quadrangular mesh[0m\n";

		///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

		// Blocks
		Mesh::Block b;
		b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
		b.SetCoords (false, 4,           // Is3D, NNodes
					 0., 1., 1., 0.,     // x coordinates
					 0., 0., 1., 1.);    // y coordinates
		b.SetNx     (ndiv);                // x weights and num of divisions along x
		b.SetNy     (ndiv);                // y weights and num of divisions along y
		b.SetETags  (4,-10,-20,-30,-40); // edge tags
		Array<Mesh::Block*> blocks;
		blocks.Push (&b);

		// Generate
		Mesh::Structured mesh;
		clock_t start = std::clock();           // Initial time
		size_t  ne    = mesh.Generate (blocks); // Discretize domain
		clock_t total = std::clock() - start;   // Time elapsed
		cout << "\nNumber of quadrangles   = " << ne << endl;
		cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

		////////////////////////////////////////////////////////////////////////////////////////// FEM /////

		// Geometry
		FEM::Geom g(2);

		// Edges brys (the order matters!)
		FEM::EBrys_T ebrys;
		ebrys.Push (make_tuple(-10, "q", 0.0));
		ebrys.Push (make_tuple(-30, "q", 0.0));
		ebrys.Push (make_tuple(-20, "u", 0.0));
		ebrys.Push (make_tuple(-40, "u", 0.0));

		// Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4Diffusion", "LinDiffusion", "k=1.0", ""));

		// Set geometry: nodes, elements, attributes, and boundaries
		FEM::SetNodesElems (&mesh, &eatts, &g);
		FEM::SetBrys       (&mesh, NULL, &ebrys, NULL, &g);

		// Set heat source
		for (size_t i=0; i<g.NElems(); ++i)
			g.Ele(i)->SetProps(source);

		// Check conductivity matrices
		LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(4,4);
		Ke_correct =  4.0, -1.0, -2.0, -1.0,
		             -1.0,  4.0, -1.0, -2.0,
		             -2.0, -1.0,  4.0, -1.0,
		             -1.0, -2.0, -1.0,  4.0;
		Ke_correct = (1.0/6.0)*Ke_correct;
		for (size_t i=0; i<g.NElems(); ++i)
		{
			LinAlg::Matrix<double> Ke;
			g.Ele(i)->Order1Matrix(0,Ke);
			for (int i=0; i<4; ++i)
			for (int j=0; j<4; ++j)
				err_Ke_2 += fabs(Ke(i,j)-Ke_correct(i,j));
		}

		// Solve
		FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
		sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
		start = std::clock();
		sol -> Solve();
		total = std::clock() - start;
		double norm_resid = LinAlg::Norm(sol->Resid());
		cout << "Time elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
		cout << "[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
		cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";

		// Output: VTU
		Output o; o.VTU (&g, "tex831.vtu");
		cout << "[1;34mFile <tex831.vtu> saved.[0m\n\n";

		// Output: Nodes
		if (ndiv<3)
		{
			cout << _6<<"Node #" << _8s<<"u" << _8s<<"q" << endl;
			for (size_t i=0; i<g.NNodes(); ++i)
				cout << _6<<i << _8s<<g.Nod(i)->Val("u") << _8s<<g.Nod(i)->Val("q") << endl;
			cout << endl;
		}

		//////////////////////////////////////////////////////////////////////////////////////// Check /////

		// Check
		std::ofstream of("tex831_2.cal", std::ios::out);
		of << _8s<<"x" << _8s<<"u" << _8s<<"ucorr" << endl;
		Array<double> err_u;
		for (size_t i=0; i<g.NNodes(); ++i)	
		{
			double x     = g.Nod(i)->X();
			double y     = g.Nod(i)->Y();
			double u     = g.Nod(i)->Val("u");
			double ucorr = u_correct(1.0,1.0,x,y);
			err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
			if (fabs(y)<=1e-5) of << _8s<<x << _8s<<u << _8s<<ucorr << endl;
		}
		of.close();

		// Error summary
		double tol_u     = 1.0e-2;
		double min_err_u = err_u[err_u.Min()];
		double max_err_u = err_u[err_u.Max()];
		cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
		cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
		cout << endl;


	}

}
catch (Exception * e)
{
	e->Cout();
	if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}
catch (char const * m)
{
	std::cout << "Fatal: " << m << std::endl;
	exit (1);
}
catch (...)
{
	std::cout << "Some exception (...) ocurred\n";
}
