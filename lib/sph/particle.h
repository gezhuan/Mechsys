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

#ifndef MECHSYS_SPH_PARTICLE_H
#define MECHSYS_SPH_PARTICLE_H

// Std lib
#include <iostream>

// MechSys
#include <mechsys/linalg/matvec.h>

class SPHParticle
{
public:

    // Constructor
    SPHParticle(Vec3_t const & x0, Vec3_t const & v0, double density0, double h0, bool Fixed);


    // Data
    bool   IsFree;                  ///< Check is the particle is free to move or not
    Vec3_t x;                       ///< Position of the particle
    Vec3_t v;                       ///< Velocity of the particle
    Vec3_t vb;                      ///< Previous velocity of the particle for the leap frog algorithm
    Vec3_t a;                       ///< Acceleration of the particle
    double Pressure;                ///< Pressure at the position of the particle
    double Density;                 ///< Density at the position of the particle
    double h;                       ///< Smoothing length of the particle


    // Methods
    void Move (double dt);                                                  ///< Update the important quantities of te simulation
    void StartAccel (Vec3_t acc = Vec3_t(0.0,0.0,0.0)) {a = acc;};          ///< Start the acceleration of the particle with one predefined value

};

inline SPHParticle::SPHParticle(Vec3_t const & x0, Vec3_t const & v0, double density0, double h0,bool Fixed)
{
    x = x0;
    v = v0;
    Density = density0;
    IsFree = !Fixed;
    h = h0;
}

inline void SPHParticle::Move (double dt)
{
}


#endif // MECHSYS_SPH_PARTICLE_H
