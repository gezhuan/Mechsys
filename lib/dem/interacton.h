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

#ifndef MECHSYS_DEM_INTERACTON_H
#define MECHSYS_DEM_INTERACTON_H

// Std Lib
#include <math.h>
#include <map>
#include <vector>
#include <utility>

// MechSys
#include <mechsys/dem/particle.h>

class Interacton
{
public:
    // typedefs
    typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;
    typedef Array<pair<int,int> > ListContacts_t;

    // Constructor and destructor
    Interacton () {};                            ///< Default constructor
    Interacton (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles

    // Methods
    virtual bool UpdateContacts (double alpha);    ///< Update contacts by verlet algorithm
    virtual void CalcForce      (double dt = 0.0); ///< Calculates the contact force between particles

    // Data
    Particle     * P1;        ///< First particle
    Particle     * P2;        ///< Second particle
    size_t         I1;        ///< Index of the first particle
    size_t         I2;        ///< Index of the second particle
    double         Kn;        ///< Normal stiffness
    double         Kt;        ///< Tengential stiffness
    double         Gn;        ///< Normal viscous coefficient
    double         Gt;        ///< Tangential viscous coefficient
    double         Mu;        ///< Microscpic coefficient of friction
    double         Epot;      ///< Potential elastic energy
    double         dEvis;     ///< Energy dissipated in viscosity at time step
    double         dEfric;    ///< Energy dissipated by friction at time step
    size_t         Nc;        ///< Number of contacts
    size_t         Nsc;       ///< Number of sliding contacts
    Vec3_t         Fn;        ///< Normal force between elements
    Vec3_t         Fnet;      ///< Net normal force
    Vec3_t         Ftnet;     ///< Net tangential force
    ListContacts_t Lee;       ///< List of edge-edge contacts 
    ListContacts_t Lvf;       ///< List of vertex-face contacts 
    ListContacts_t Lfv;       ///< List of face-vertex contacts
    FrictionMap_t  Fdee;      ///< Static friction displacement for pair of edges
    FrictionMap_t  Fdvf;      ///< Static friction displacement for pair of vertex-face
    FrictionMap_t  Fdfv;      ///< Static friction displacement for pair of face-vertex
protected:
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt);
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_contacts        (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha);
};


class InteractonSphere: public Interacton
{
public:
    // Methods 
    InteractonSphere (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    bool UpdateContacts (double alpha);                ///< Update contacts by verlet algorithm
    void CalcForce (double dt = 0.0);                  ///< Calculates the contact force between particles
    
    // Data
    FrictionMap_t  Fdvv;                                ///< Static Friction displacement for the vertex vertex pair
    ListContacts_t Lvv;                                 ///< List of vertices
    Vec3_t         Fdr;                                 ///< Rolling displacement 
    double         beta;                                ///< Rolling stiffness coefficient
    double         eta;                                 ///< Plastic moment coefficient

protected:
    void _update_rolling_resistance(double dt);               ///< Calculates the rolling resistance torque
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Interacton::Interacton (Particle * Pt1, Particle * Pt2)
    : P1(Pt1), P2(Pt2), Kn(2*ReducedValue(P1->Kn,P2->Kn)), Kt(2*ReducedValue(P1->Kt,P2->Kt)), Gn(2*ReducedValue(P1->Gn,P2->Gn)), 
      Gt(2*ReducedValue(P1->Gt,P2->Gt)), Mu(2*ReducedValue(P1->Mu,P2->Mu)), Epot(0.0)
{
    I1 = P1->Index;
    I2 = P2->Index;
    CalcForce(0.0);
}


inline bool Interacton::UpdateContacts (double alpha)
{
    Lee.Resize(0);
    Lvf.Resize(0);
    Lfv.Resize(0);
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha)
    {
        _update_contacts (P1->Edges,P2->Edges,Lee,alpha);
        _update_contacts (P1->Verts,P2->Faces,Lvf,alpha);
        _update_contacts (P1->Faces,P2->Verts,Lfv,alpha);
        if (Lee.Size()>0||Lvf.Size()>0||Lfv.Size()>0) return true;
        else return false;
    }
    else return false;
}

inline void Interacton::CalcForce (double dt)
{
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Fnet   = 0.0;
    Ftnet   = 0.0;
    //if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax)
    //{
        //_update_disp_calc_force (P1->Edges,P2->Edges,Fdee,Lee,dt);
        //_update_disp_calc_force (P1->Verts,P2->Faces,Fdvf,Lvf,dt);
        //_update_disp_calc_force (P1->Faces,P2->Verts,Fdfv,Lfv,dt);
    //}
    _update_disp_calc_force (P1->Edges,P2->Edges,Fdee,Lee,dt);
    _update_disp_calc_force (P1->Verts,P2->Faces,Fdvf,Lvf,dt);
    _update_disp_calc_force (P1->Faces,P2->Verts,Fdfv,Lfv,dt);

    //If there is at least a contact, increase the coordination number of the particles
    if (Nc>0) 
    {
        P1->Cn++;
        P2->Cn++;
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline void Interacton::_update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt)
{
    // update
    for (size_t k=0; k<L.Size(); ++k)
    {
        size_t i = L[k].first;
        size_t j = L[k].second;
        Vec3_t xi, xf;
        Distance ((*A[i]), (*B[j]), xi, xf);
        double dist  = norm(xf-xi);
        double delta = P1->R + P2->R - dist;
        if (delta>0)
        {
            // Count a contact
            Nc++;

            // update force
            Vec3_t n = (xf-xi)/dist;
            Vec3_t x = xi+n*((P1->R*P1->R-P2->R*P2->R+dist*dist)/(2*dist));
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;
            pair<int,int> p;
            p = make_pair(i,j);
            Fn = Kn*delta*n;
            Fnet += Fn;

            FMap[p] += vt*dt;
            FMap[p] -= dot(FMap[p],n)*n;
            Vec3_t tan = FMap[p];
            if (norm(tan)>1.0e-22) tan/=norm(tan);
            if (norm(FMap[p])>Mu*norm(Fn)/Kt)
            {
                // Count a sliding contact
                Nsc++;

                FMap[p] = Mu*norm(Fn)/Kt*tan;
                dEfric += Kt*dot(FMap[p],vt)*dt;
            }
            Ftnet += Kt*FMap[p];
            Vec3_t F = Fn + Kt*FMap[p] + Gn*dot(n,vrel)*n + Gt*vt;
            P1->F += -F;
            P2->F +=  F;
            dEvis += (Gn*dot(vrel-vt,vrel-vt)+Gt*dot(vt,vt))*dt;

            // torque
            Vec3_t T, Tt, temp;
            temp = x - P1->x;
            Tt = cross (temp,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            P1->T -= T;
            temp = x - P2->x;
            Tt = cross (temp,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            P2->T += T;

            //Transfering the branch vector information
            for (size_t m=0;m<3;m++)
            {
                for (size_t n=0;n<3;n++)
                {
                    P1->M(m,n)-=F(m)*x1(n);
                    P2->M(m,n)+=F(m)*x2(n);
                    P1->B(m,n)+=x1(m)*x1(n);
                    P2->B(m,n)+=x2(m)*x2(n);
                }
            }

            // potential energy
            Epot += 0.5*Kn*delta*delta+0.5*Kt*dot(FMap[p],FMap[p]);
        }
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline void Interacton::_update_contacts (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha)
{
    for (size_t i=0; i<A.Size(); ++i)
    for (size_t j=0; j<B.Size(); ++j)
    {
        if (Distance ((*A[i]), (*B[j]))<=P1->R+P2->R+2*alpha)
        {
            pair<int,int> p;
            p = make_pair(i,j);
            L.Push(p);
        }
    }

}

inline InteractonSphere::InteractonSphere (Particle * Pt1, Particle * Pt2)
{
    P1   = Pt1;
    P2   = Pt2;
    I1 = P1->Index;
    I2 = P2->Index;
    Kn   = 2*ReducedValue(Pt1->Kn,Pt2->Kn);
    Kt   = 2*ReducedValue(Pt1->Kt,Pt2->Kt);
    Gn   = 2*ReducedValue(Pt1->Gn,Pt2->Gn);
    Gt   = 2*ReducedValue(Pt1->Gt,Pt2->Gt);
    Mu   = 2*ReducedValue(Pt1->Mu,Pt2->Mu);
    beta = 2*ReducedValue(Pt1->Beta,Pt2->Beta);
    eta  = 2*ReducedValue(Pt1->Eta,Pt2->Eta);
    Nc = 0;
    Nsc = 0;

    Epot = 0.0;
    Fdr  = 0.0, 0.0, 0.0;
    Lvv.Push(make_pair(0,0));

    CalcForce(0.0);
}

inline void InteractonSphere::_update_rolling_resistance(double dt)
{
    Vec3_t t1,t2;
    Rotation(P1->w,P1->Q,t1);
    Rotation(P2->w,P2->Q,t2);
    Vec3_t Normal = Fn/norm(Fn);
    Vec3_t Vr = P1->R*P2->R*cross(Vec3_t(t1 - t2),Normal)/(P1->R+P2->R);
    Fdr += Vr*dt;
    Fdr -= dot(Fdr,Normal)*Normal;
    Vec3_t tan = Fdr;
    if (norm(tan)>1.0e-22) tan/=norm(tan);
    double Kr = beta*Kt;
    if (norm(Fdr)>eta*Mu*norm(Fn)/Kr)
    {
        Fdr = eta*Mu*norm(Fn)/Kr*tan;
    }
    Vec3_t Ft = -Kr*Fdr;

    Vec3_t Tt = P1->R*cross(Normal,Ft);
    Vec3_t T;
    Quaternion_t q;
    Conjugate (P1->Q,q);
    Rotation  (Tt,q,T);
    P1->T += T;

    Tt = P2->R*cross(Normal,Ft);
    Conjugate (P2->Q,q);
    Rotation  (Tt,q,T);
    P2->T -= T;
}

inline void InteractonSphere::CalcForce(double dt)
{
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Fnet   = 0.0;
    Ftnet  = 0.0;
    _update_disp_calc_force (P1->Verts,P2->Verts,Fdvv,Lvv,dt);
    if (Epot>0.0) _update_rolling_resistance(dt);


    //If there is at least a contact, increase the coordination number of the particles
    if (Nc>0) 
    {
        P1->Cn++;
        P2->Cn++;
    }
}

inline bool InteractonSphere::UpdateContacts (double alpha)
{
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha) return true;
    else return false;
}
#endif //  MECHSYS_DEM_INTERACTON_H
