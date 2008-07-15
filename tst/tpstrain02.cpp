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

// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/quad4pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	double H     = 2.0;   // height
	double L     = 2.0;   // length
	double E     = 207.0; // Young
	double nu    = 0.3;   // Poisson
	double q     = 1.0;   // Load
	int    ndivx = 2;     // number of divisions along x
	int    ndivy = 2;     // number of divisions along y

	/*        | | | | | | | | | | | | | | | | |  q
	          V V V V V V V V V V V V V V V V V 
	  -+-     @-------@-------@-------@-------@
	   |      |20     |21     |22     |23     |24
	   |      |   e12 |   e13 |   e14 |   e15 |
	   |      |15     |16     |17     |18     |19
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e8  |   e9  |   e10 |   e11 |
	          |10     |11     |12     |13     |14
	   H      @-------@-------@-------@-------@
	          |       |       |       |       |
	   |      |   e4  |   e5  |   e6  |   e7  |
	   |      |5      |6      |7      |8      |9
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e0  |   e1  |   e2  |   e3  |
	   |      |0      |1      |2      |3      |4
	  -+-     @-------@-------@-------@-------@
	         /_\     /_\     /_\     /_\     /_\
	         o o     o o     ///     o o     o o
	
	          |-------------- L --------------|
	 */

	// Input
	String linsol("LA");
	if (argc==2)   linsol.Printf("%s",argv[1]);
	if (argc==3) { linsol.Printf("%s",argv[1]); ndivx = ndivy = atoi(argv[2]); }
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol nDiv\n  where LinSol:\n \tLA   => LAPACK_T  : DENSE\n \tUM   => UMFPACK_T : SPARSE\n \tSLU  => SuperLU_T : SPARSE\n \tnDiv => Number of division along x and y (must be even)\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	String wx; for (int i=0; i<ndivx; ++i) wx.Printf("%s %f",wx.CStr(),1.0);
	String wy; for (int i=0; i<ndivy; ++i) wy.Printf("%s %f",wy.CStr(),1.0);
	Mesh::Block b;
	b.SetTag (-1); // tag to be replicated to all generated elements inside this block
	b.Set2D  ();   // 2D
	b.C      () = 0.,  L, L, 0.,    L/2.,    L, L/2.,   0., // x coordinates
	              0., 0., H,  H,      0., H/2.,    H, H/2.; // y coordinates
	b.SetWx  (wx.CStr());                                   // x weights and num of divisions along x
	b.SetWy  (wx.CStr());                                   // y weights and num of divisions along y
	b.ETags  () = 0, 0, -10, -20;                           // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	cout << "\nMesh Generation: --------------------------------------------------------------" << endl;
	Mesh::Structured ms;
	clock_t start = std::clock(); // Initial time
	size_t  ne    = ms.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "[1;33m"<<ne<<" elements[0m. Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

	cout << ms << endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// 0) Geometry
	FEM::Geom g(2); // 2D

	// 1) Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (make_tuple(L/2., 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val

	// 2) Faces brys
	FEM::FBrys_T fbrys;
	fbrys.Push (make_tuple(-10, "uy", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-20, "fy",  -q)); // tag, key, val

	// 3) Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0")); // tag, type, model, prms, inis

	// 4) Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetGeom (&ms, nbrys, fbrys, eatts, &g);

	// 5) Solve
	cout << "\nSolution: ---------------------------------------------------------------------" << endl;
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	start = std::clock(); // Initial time
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> SetCte("DTOL", 1.0e-10);
	sol -> Solve();
	total = std::clock() - start; // Time elapsed
	//cout << "NormResid = "<<sol->GetVar("NormResid")<<". Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;
	cout << "RelError  = "<<sol->GetVar("RelError" )<<". Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    double errors = 0.0;

	double Sy = q;
	double Ex = -nu*(1.0+nu)*Sy/E;
	double Ey =  (1.0-nu*nu)*Sy/E;
	double Sz = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey);

	// Stress and strains
	for (size_t i=0; i<g.nElems(); ++i)
	{
		errors += fabs(g.Ele(i)->Val("Ex" ) - (Ex));
		errors += fabs(g.Ele(i)->Val("Ey" ) - (Ey));
		errors += fabs(g.Ele(i)->Val("Exy") - (0.0));
		errors += fabs(g.Ele(i)->Val("Sx" ) - (0.0));
		errors += fabs(g.Ele(i)->Val("Sy" ) - (Sy ));
		errors += fabs(g.Ele(i)->Val("Sz" ) - (Sz ));
		errors += fabs(g.Ele(i)->Val("Sxy") - (0.0));
	}

	// Displacements
	for (size_t i=0; i<g.nNodes(); ++i)
	{
		if (fabs(g.Nod(i)->Y())<1.0e-5) // bottom nodes
		{
			errors += fabs(g.Nod(i)->Val("uy") - (0.0));
			if (fabs(g.Nod(i)->X()-(L/2.0))<1.0e-5) // central node
				errors += fabs(g.Nod(i)->Val("ux") - (0.0));
		}
		else if (fabs(g.Nod(i)->Y()-(H/2.0))<1.0e-5) // mid nodes
			errors += fabs(g.Nod(i)->Val("uy") - (-0.5*H*Ey));
		else if (fabs(g.Nod(i)->Y()-(H))<1.0e-5) // top nodes
			errors += fabs(g.Nod(i)->Val("uy") - (-H*Ey));
		if (fabs(g.Nod(i)->X())<1.0e-5) // left nodes
			errors += fabs(g.Nod(i)->Val("ux") - (0.5*L*Ex));
		if (fabs(g.Nod(i)->X()-L)<1.0e-5) // right nodes
			errors += fabs(g.Nod(i)->Val("ux") - (-0.5*L*Ex));
	}

	if (fabs(errors)>1.0e-10) cout << "[1;31m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Write file
	cout << "Write output file: ------------------------------------------------------------" << endl;
	start = std::clock(); // Initial time
	FEM::WriteVTUEquilib (g, "tpstrain02.vtu");
	total = std::clock() - start; // Time elapsed
	cout << "[1;34mFile <tpstrain02.vtu> saved.[0m" << endl;
	cout << "Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;
	cout << endl;

	return 0;
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
