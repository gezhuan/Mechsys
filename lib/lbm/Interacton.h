/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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


#ifndef MECHSYS_LBM_INTERACTON_H
#define MECHSYS_LBM_INTERACTON_H

// Mechsys
#include <mechsys/lbm/Dem.h>
#include <mechsys/dem/basic_functions.h>


class Interacton
{
public:
    //Constructor
    Interacton () {};
    Interacton(Disk * D1, Disk * D2);

    //Methods
    void CalcForce      (double dt);
    bool UpdateContacts (double Alpha);

    //Data
    Disk * D1;       //Pointer to first particle
    Disk * D2;       //Pointer to second particle
    double Kn;       //Spring constant 
};

Interacton::Interacton(Disk * Dp1, Disk * Dp2)
{
    D1 = Dp1;
    D2 = Dp2;
    Kn = 2.0*ReducedValue(D1->Kn,D2->Kn);
}

void Interacton::CalcForce(double dt)
{
    double dist  = norm(D2->X - D1->X);
    double delta = D1->R + D2->R - dist;
    if (delta>0)
    {
        Vec3_t n = (D2->X - D1->X)/dist;
        D1->F   -= Kn*delta*n;
        D2->F   += Kn*delta*n;
    }
}

bool Interacton::UpdateContacts(double Alpha)
{
    if (norm(D1->X-D2->X) <= D1->R + D2->R + 2*Alpha) return true;
    else                                              return false;
}

#endif
