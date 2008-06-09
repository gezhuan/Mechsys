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

#ifndef MECHSYS_LINELASTIC2D_H
#define MECHSYS_LINELASTIC2D_H

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "models/equilibmodel2d.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"
#include "linalg/laexpr.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class LinElastic2D : public EquilibModel2D
{
public:
	// Typedefs
	typedef blitz::TinyVector<double,3>    Tensor2_2d;
	typedef blitz::TinyMatrix<double,3,3>  Tensor4_2d;

	// Constructor
	LinElastic2D ();
	
	// Destructor
	virtual ~LinElastic2D () {}

	// Parameters and initial values
	void SetPrms (String const & Prms);
	void SetInis (String const & Inis);

	// Derived Methods
	void TgStiffness  (Matrix<double> & D) const;
	int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig);
	void BackupState  ();
	void RestoreState ();

	// Access Methods
	void Sig (Vector<double> & Stress ) const { Stress .Resize(3); Stress = _sig(0), _sig(1), _sig(2); }
	void Eps (Vector<double> & Strain ) const { Strain .Resize(3); Strain = _eps(0), _eps(1), _eps(2); }
	void Ivs (Array<double>  & IntVals) const { IntVals.Resize(0); }

private:
	// Data
	Tensor2_2d _sig;
	Tensor2_2d _eps;
	Tensor2_2d _sig_bkp;
	Tensor2_2d _eps_bkp;
	Tensor4_2d _De;

}; // class LinElastic2D


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline LinElastic2D::LinElastic2D()
{
	// Parameters
	SetPrms ("E=10000.0 nu=0.25");

	// Data
	SetInis ("Sx=0.0 Sy=0.0 Sxy=0.0");
}

inline void LinElastic2D::SetPrms(String const & Prms)
{
	/* "E=20000.0 nu=0.2" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	if (names.Size()==2 && values.Size()==2)
	{
		if (names[0]=="E" && names[1]=="nu")
		{
			double E  = values[0];
			double nu = values[1];
			double c  = E/((1.0+nu)*(1.0-2.0*nu));
			_De = c*(1.0-nu), c*(    nu), c*(0.0),
			      c*(    nu), c*(1.0-nu), c*(0.0),
			      c*(   0.0), c*(   0.0), c*((1.0-2.0*nu)/2.0);
			return;
		}
	}
	else throw Fatal("LinElastic2D::SetPrms: Parameters definition is incorrect. The syntax must be as in:\n\t E=10000.0 nu=0.25\n");
}

inline void LinElastic2D::SetInis(String const & Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	if (names.Size()==3 && values.Size()==3)
	{
		if (names[0]=="Sx" && names[1]=="Sy" && names[2]=="Sxy")
		{
			_sig  = values[0], values[1], values[2];
			_eps  =       0.0,       0.0,       0.0;
			return;
		}
	}
	else throw Fatal("LinElastic2D::SetInis: Initial values definition is incorrect. The syntax must be as in:\n\t Sx=0.0 Sy=0.0 Sxy=0.0\n");
}

inline void LinElastic2D::TgStiffness(Matrix<double> & De) const
{
	De.Resize(3,3);
	De = _De(0,0), _De(0,1), _De(0,2),
	     _De(1,0), _De(1,1), _De(1,2),
	     _De(2,0), _De(2,1), _De(2,2);
}

inline int LinElastic2D::StressUpdate(Vector<double> const & DEps, Vector<double> & DSig)
{
	Tensor2_2d deps;  deps = DEps(0), DEps(1), DEps(2);
	Tensor2_2d dsig;
	blitz::product(_De, deps);
	_eps += deps;
	_sig += dsig;
	DSig.Resize(3);
	DSig = dsig(0), dsig(1), dsig(2);
	return 1;
}

inline void LinElastic2D::BackupState()
{
	_sig_bkp = _sig;
	_eps_bkp = _eps;
}

inline void LinElastic2D::RestoreState()
{
	_sig = _sig_bkp;
	_eps = _eps_bkp;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new LinElastic2D model
EquilibModel2D * LinElastic2DMaker()
{
	return new LinElastic2D();
}

// Register an LinElastic2D model into ModelFactory array map
int LinearElastic2DRegister()
{
	EquilibModel2DFactory["LinElastic2D"] = LinElastic2DMaker;
	return 0;
}

// Execute the autoregistration
int __LinearElastic2D_dummy_int = LinearElastic2DRegister();


#endif // MECHSYS_LINELASTIC2D_H
