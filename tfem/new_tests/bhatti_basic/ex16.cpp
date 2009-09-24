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

/*  Bhatti (2005): Example 1.6, p32  *
 *  ===============================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/elems/tri3.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (6/*verts*/, 4/*cells*/);
    mesh.SetVert   (0, -100, 0.0, 0.0, 0);
    mesh.SetVert   (1, -100, 0.0, 2.0, 0);
    mesh.SetVert   (2,    0, 2.0, 0.0, 0);
    mesh.SetVert   (3,    0, 2.0, 1.5, 0);
    mesh.SetVert   (4,    0, 4.0, 0.0, 0);
    mesh.SetVert   (5,    0, 4.0, 1.0, 0);
    mesh.SetCell   (0,   -1, /*NVerts*/3, 0,2,3);
    mesh.SetCell   (1,   -1, /*NVerts*/3, 3,1,0);
    mesh.SetCell   (2,   -1, /*NVerts*/3, 2,4,5);
    mesh.SetCell   (3,   -1, /*NVerts*/3, 5,3,2);
    mesh.SetBryTag (1, 0, -10);
    mesh.SetBryTag (3, 0, -10);
    //mesh.WriteVTU  ("ex16");
    //cout << mesh << endl;

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active  h pse", PROB("Equilib"), GEOM("Tri3"), TRUE, 0.25, TRUE);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu pse", MODEL("LinElastic"), 1.0e+4, 0.2, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(/*NDim*/2, prps, mdls, inis);
    dom.SetMesh    (mesh);
    dom.SetOutNods ("ex16",/*NNod*/1,/*ID*/5);
    dom.SetOutEles ("ex16",/*NEle*/1,/*ID*/3);

    // check matrices
    {
        double tol   = 1.0e-10;
        double error = 0.0;
        Mat_t K0c(6,6),K1c(6,6),K2c(6,6),K3c(6,6);
        K0c =
          9.7656250000000011e+02,  0.0000000000000000e+00, -9.7656250000000011e+02,  2.6041666666666669e+02,  0.0000000000000000e+00, -2.6041666666666669e+02,
          0.0000000000000000e+00,  3.9062500000000000e+02,  5.2083333333333337e+02, -3.9062500000000000e+02, -5.2083333333333337e+02,  0.0000000000000000e+00,
         -9.7656250000000011e+02,  5.2083333333333337e+02,  1.6710069444444448e+03, -7.8125000000000000e+02, -6.9444444444444434e+02,  2.6041666666666669e+02,
          2.6041666666666669e+02, -3.9062500000000000e+02, -7.8125000000000000e+02,  2.1267361111111113e+03,  5.2083333333333337e+02, -1.7361111111111111e+03,
          0.0000000000000000e+00, -5.2083333333333337e+02, -6.9444444444444434e+02,  5.2083333333333337e+02,  6.9444444444444434e+02,  0.0000000000000000e+00,
         -2.6041666666666669e+02,  0.0000000000000000e+00,  2.6041666666666669e+02, -1.7361111111111111e+03,  0.0000000000000000e+00,  1.7361111111111111e+03;
        K1c =
          1.3020833333333335e+03,  0.0000000000000000e+00, -9.7656250000000011e+02,  2.6041666666666669e+02, -3.2552083333333337e+02, -2.6041666666666669e+02,
          0.0000000000000000e+00,  5.2083333333333337e+02,  5.2083333333333337e+02, -3.9062500000000000e+02, -5.2083333333333337e+02, -1.3020833333333334e+02,
         -9.7656250000000011e+02,  5.2083333333333337e+02,  1.2532552083333335e+03, -5.8593750000000000e+02, -2.7669270833333337e+02,  6.5104166666666657e+01,
          2.6041666666666669e+02, -3.9062500000000000e+02, -5.8593750000000000e+02,  1.5950520833333335e+03,  3.2552083333333331e+02, -1.2044270833333335e+03,
         -3.2552083333333337e+02, -5.2083333333333337e+02, -2.7669270833333337e+02,  3.2552083333333337e+02,  6.0221354166666674e+02,  1.9531250000000000e+02,
         -2.6041666666666669e+02, -1.3020833333333334e+02,  6.5104166666666657e+01, -1.2044270833333335e+03,  1.9531250000000000e+02,  1.3346354166666667e+03;
        K2c =
          6.5104166666666674e+02,  0.0000000000000000e+00, -6.5104166666666674e+02,  2.6041666666666669e+02,  0.0000000000000000e+00, -2.6041666666666669e+02,
          0.0000000000000000e+00,  2.6041666666666669e+02,  5.2083333333333337e+02, -2.6041666666666669e+02, -5.2083333333333337e+02,  0.0000000000000000e+00,
         -6.5104166666666674e+02,  5.2083333333333337e+02,  1.6927083333333335e+03, -7.8125000000000000e+02, -1.0416666666666667e+03,  2.6041666666666669e+02,
          2.6041666666666669e+02, -2.6041666666666669e+02, -7.8125000000000000e+02,  2.8645833333333335e+03,  5.2083333333333337e+02, -2.6041666666666670e+03,
          0.0000000000000000e+00, -5.2083333333333337e+02, -1.0416666666666667e+03,  5.2083333333333337e+02,  1.0416666666666667e+03,  0.0000000000000000e+00,
         -2.6041666666666669e+02,  0.0000000000000000e+00,  2.6041666666666669e+02, -2.6041666666666670e+03,  0.0000000000000000e+00,  2.6041666666666670e+03;
        K3c =
          9.7656250000000011e+02,  0.0000000000000000e+00, -6.5104166666666674e+02,  2.6041666666666669e+02, -3.2552083333333337e+02, -2.6041666666666669e+02,
          0.0000000000000000e+00,  3.9062500000000000e+02,  5.2083333333333337e+02, -2.6041666666666669e+02, -5.2083333333333337e+02, -1.3020833333333334e+02,
         -6.5104166666666674e+02,  5.2083333333333337e+02,  1.1284722222222222e+03, -5.2083333333333337e+02, -4.7743055555555554e+02,  0.0000000000000000e+00,
          2.6041666666666669e+02, -2.6041666666666669e+02, -5.2083333333333337e+02,  1.9097222222222224e+03,  2.6041666666666669e+02, -1.6493055555555557e+03,
         -3.2552083333333337e+02, -5.2083333333333337e+02, -4.7743055555555554e+02,  2.6041666666666669e+02,  8.0295138888888880e+02,  2.6041666666666669e+02,
         -2.6041666666666669e+02, -1.3020833333333334e+02,  0.0000000000000000e+00, -1.6493055555555557e+03,  2.6041666666666669e+02,  1.7795138888888889e+03;
        Mat_t K0,K1,K2,K3;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        dom.Eles[3]->CalcK(K3);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        error += CompareMatrices (K2,K2c);
        error += CompareMatrices (K3,K3c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    // solver
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "qn",   -20.0)
       .Set(-100, "ux uy", 0.0,0.0);
    dom.SetBCs (bcs);
    sol.Solve  (10,1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_12_6);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("                   ux                       uy", /*NRows*/ 6,
                 0.000000000000000e+00,   0.000000000000000e+00,
                 0.000000000000000e+00,   0.000000000000000e+00,
                -1.035527877607004e-02,  -2.552969847657423e-02,
                 4.727650463081949e-03,  -2.473565538172127e-02,
                -1.313941349422282e-02,  -5.549310752960183e-02,
                 8.389015766816341e-05,  -5.556637423271112e-02);

    Table ele_sol;
    ele_sol.Set("sx  sy  sz  sxy   ex  ey  ez  exy", /*NRows*/4,
                -5.283090599362460e+01, -5.272560566371797e+00, 0.000000000000000e+00, -1.128984616188524e+01, -5.177639388035024e-03,  5.293620632353122e-04,  1.162069331199928e-03, -2.709563078852457e-03/2.0,
                 2.462317949521848e+01,  4.924635899043697e+00, 0.000000000000000e+00, -5.153261537858599e+01,  2.363825231540974e-03,  0.000000000000000e+00, -5.909563078852436e-04, -1.236782769086064e-02/2.0,
                -1.465334062185674e+01, -3.663335155464233e+00, 0.000000000000000e+00, -7.326670310928396e+00, -1.392067359076390e-03, -7.326670310928846e-05,  3.663335155464196e-04, -1.758400874622815e-03/2.0,
                 3.102227081237862e+00,  5.914066048600676e+00, 0.000000000000000e+00, -2.178221979271434e+01,  1.919413871517726e-04,  5.293620632353103e-04, -1.803258625967707e-04, -5.227732750251441e-03/2.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-15, 1.0e-15);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-12,1.0e-12,1.0e-15,1.0e-12, 1.0e-15,1.0e-15,1.0e-15,1.0e-15);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
