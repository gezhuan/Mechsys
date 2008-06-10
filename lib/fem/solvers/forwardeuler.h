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

#ifndef MECHSYS_FEM_FORWARDEULER_H
#define MECHSYS_FEM_FORWARDEULER_H

// MechSys
#include "fem/solver.h"

namespace FEM
{

class ForwardEuler: public Solver
{
private:
	// Data
	LinAlg::Vector<double> _resid;

	// Private methods
	void _do_solve_for_an_increment(double dTime);

}; // class ForwardEuler


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* private */

inline void ForwardEuler::_do_solve_for_an_increment(double dTime)
{
	// Allocate and fill dF_ext and dU_ext
	int ndofs = _dF_ext.Size();
	LinAlg::Vector<double> dF_ext(ndofs);
	LinAlg::Vector<double> dU_ext(ndofs);
	LinAlg::CopyScal(1.0/GFE_nSI,_dF_ext, dF_ext); // dF_ext <- _dF_ext/GFE_nSI
	LinAlg::CopyScal(1.0/GFE_nSI,_dU_ext, dU_ext); // dU_ext <- _dU_ext/GFE_nSI
	double h = dTime/GFE_nSI;

	// Allocate auxiliar vector
	if (_resid.Size()==0) _resid.Resize(ndofs);

	// Start
	for (int i=0; i<GFE_nSI; ++i)
	{
		// Assemble G matrix and calculate dU_ext
		_inv_G_times_dF_minus_hKU(h, dF_ext, dU_ext); // dU_ext <- inv(G)*(dF_ext - hKU)

		// Update nodes and elements state
		_update_nodes_and_elements(h, dU_ext);

		// Calculate residual (internal)
		if (_has_hKU)
		LinAlg::Axpy      (+1.0,_hKU, dF_ext);                // dF_ext <- dF_ext + hKU
		LinAlg::AddScaled (1.0,dF_ext, -1.0,_dF_int, _resid); // _resid <- dF_ext - dF_int
		double denom = 0.0;                                   // Normalizer
		for (int i=0; i<ndofs; i++) denom += pow((dF_ext(i)+_dF_int(i))/2.0, 2.0);
		GFE_Resid = LinAlg::Norm(_resid)/(sqrt(denom)+1.0);
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new ForwardEuler solver
Solver * ForwardEulerMaker()
{
	return new ForwardEuler();
}

// Register an ForwardEuler solver into SolverFactory array map
int ForwardEulerRegister()
{
	SolverFactory["ForwardEuler"] = ForwardEulerMaker;
	return 0;
}

// Execute the autoregistration
int __ForwardEuler_dummy_int = ForwardEulerRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_FORWARDEULER_H
