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

#ifndef MECHSYS_FEM_HEX8EQUILIB_H
#define MECHSYS_FEM_HEX8EQUILIB_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/hex8.h"
#include "util/exception.h"

namespace FEM
{

class Hex8Equilib : public Hex8, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int  _geom     () const { return 3; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize()
	{
		if (_ndim!=3) throw new Fatal("Hex8Equilib::_initialize: For this element, _ndim must be equal to 3 (%d is invalid)",_ndim);
		_d  = _ndim-1;
		_nd = EquilibElem::ND[_d];
		_nl = EquilibElem::NL[_geom()-1];
	}

}; // class Hex8Equilib

// Hex8Equilib constants
char const * Hex8Equilib::NAME = "Hex8Equilib";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Hex8Equilib element
Element * Hex8EquilibMaker()
{
	return new Hex8Equilib();
}

// Register a Hex8Equilib element into ElementFactory array map
int Hex8EquilibRegister()
{
	ElementFactory[Hex8Equilib::NAME] = Hex8EquilibMaker;
	return 0;
}

// Execute the autoregistration
int __Hex8Equilib_dummy_int  = Hex8EquilibRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_HEX8EQUILIB_H
