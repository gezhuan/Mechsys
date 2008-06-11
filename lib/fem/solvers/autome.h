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

#ifndef MECHSYS_FEM_AUTOME_H
#define MECHSYS_FEM_AUTOME_H

// MechSys
#include "fem/solver.h"
#include "util/string.h"
#include "util/lineparser.h"

namespace FEM
{

class AutoME: public Solver
{
private:
	// Data
	LinAlg::Vector<double> _dF_1;
	LinAlg::Vector<double> _dU_1;
	LinAlg::Vector<double> _dF_2;
	LinAlg::Vector<double> _dU_2;
	LinAlg::Vector<double> _dU_ME;
	LinAlg::Vector<double> _Err_U;
	LinAlg::Vector<double> _Err_F;

	// Private methods
	void _do_solve_for_an_increment(double dTime);

}; // class AutoME


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* private */

inline void AutoME::_do_solve_for_an_increment(double dTime)
{
	// Allocate and fill local dF_ext and dU_ext
	LinAlg::Vector<double> dF_ext(_dF_ext);
	LinAlg::Vector<double> dU_ext(_dU_ext);

	// Allocate auxiliar vectors
	if (_dF_1.Size()==0)
	{
		int ndofs = dF_ext.Size();
		_dF_1 .Resize (ndofs);
		_dU_1 .Resize (ndofs);
		_dF_2 .Resize (ndofs);
		_dU_2 .Resize (ndofs);
		_dU_ME.Resize (ndofs);
		_Err_U.Resize (ndofs);
		_Err_F.Resize (ndofs);
	}

	// Start substeps
	double  T = 0.0;
	double dT = GME_dTini;
	for (int k=0; k<GME_maxSI; ++k)
	{
		if (T>=1.0) return;

		// Sub-divide timestep
		double h = dT*dTime;

		// Calculate scaled increment vectors for this sub-step
		LinAlg::CopyScal(dT,dF_ext, _dF_1); // _dF_1 <- dT*dF_ext
		LinAlg::CopyScal(dT,dU_ext, _dU_1); // _dU_1 <- dT*dU_ext;
		LinAlg::CopyScal(dT,dF_ext, _dF_2); // _dF_2 <- dT*dF_ext
		LinAlg::CopyScal(dT,dU_ext, _dU_2); // _dU_2 <- dT*dU_ext;

		// Backup element state (needed when updating disp. state for a ME increment)
		_backup_nodes_and_elements();

		// Forward-Euler: Assemble G matrix and calculate _dU_1
		_inv_G_times_dF_minus_hKU(h, _dF_1, _dU_1); // _dU_1 <- inv(G)*(dF_ext - hKU)
	
		// Forward-Euler: update nodes and elements state
		_update_nodes_and_elements(h, _dU_1);

		// Modified-Euler: Assemble G matrix and calculate dU_2
		_inv_G_times_dF_minus_hKU(h, _dF_2, _dU_2); // _dU_2 <- inv(G)*(dF_ext - hKun)
	
		// Save the norm of essential and natural vectors
		double normU = _norm_essential_vector();
		double normF = _norm_natural_vector();

		// Calculate local error
		if (normF<=GME_ZTOL) throw new Message(_("AutoME::_do_solve_for_an_increment: k=%d: normF=%e cannot be equal to ZTOL (%e)"),k,normF,GME_ZTOL);
		LinAlg::AddScaled (0.5,_dU_2, -0.5,_dU_1, _Err_U); // Error on U (%)
		LinAlg::AddScaled (0.5,_dF_2, -0.5,_dF_1, _Err_F); // Error on F
		double Rerr_U = LinAlg::Norm(_Err_U)/(1.0+normU);  // (R)elative error on U
		double Rerr_F = LinAlg::Norm(_Err_F)/normF;        // (R)elative error on F
		GME_Rerr = (Rerr_U>Rerr_F ? Rerr_U : Rerr_F);      // (R)elative error
		double m = GME_mCoef*sqrt(GME_DTOL/GME_Rerr);

		// Restore nodes and element for initial disp. state given at the start of the increment
		_restore_nodes_and_elements();

		if (GME_Rerr<=GME_DTOL)
		{
			// Calculate Modified-Euler displacement increment vector
			LinAlg::AddScaled(0.5,_dU_1, 0.5,_dU_2, _dU_ME);

			// Update nodes and elements state for a Modified-Euler evaluation of displacements
			_update_nodes_and_elements(h, _dU_ME);

			// Next pseudo time
			T = T + dT;
			if (m>GME_mMax) m=GME_mMax;
		}
		else
			if (m<GME_mMin) m=GME_mMin;

		// Next substep size
		dT = m*dT;
		if (dT>1.0-T) dT=1.0-T;
	}
	if (GME_Cconv) throw new Fatal(_("AutoME::_do_solve_for_an_increment: did not converge for %d substeps"), GME_maxSI);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new AutoME solver
Solver * AutoMEMaker()
{
	return new AutoME();
}

// Register an AutoME solver into SolverFactory array map
int AutoMERegister()
{
	SolverFactory["AutoME"] = AutoMEMaker;
	return 0;
}

// Execute the autoregistration
int __AutoME_dummy_int = AutoMERegister();

}; // namespace FEM

#endif // MECHSYS_FEM_AUTOME_H
