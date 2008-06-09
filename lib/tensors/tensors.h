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

#ifndef MECHSYS_TENSORS_H
#define MECHSYS_TENSORS_H

// STL
#include <cfloat> // for DBL_EPSILON

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

/** \namespace Tensors %Tensors with components according to the Mandel's basis.
 
    2nd and 4th order tensors are reprsented according to the Mandel basis
    see Brannon, 2003, http://www.mech.utah.edu/~brannon/gobag.html
    \verbatim
        @CONFERENCE{brannon:00,
         author    = {Rebecca M. Brannon},
         title     = {Geometrical interpretation of radial and oblique return methods},
         booktitle = {Plasticity 2000: The Eighth International Symposium on Plasticity and its Current Applications},
         year      = {2000},
         pages     = {30}}
      \endverbatim

    Examples
	 - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/tensors/tst/tfuns.cpp?view=markup">tfuns.cpp</a>
	 - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/tensors/tst/tops.cpp?view=markup">tops.cpp</a>
	 - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/tensors/tst/ttij.cpp?view=markup">ttij.cpp</a>
 */
namespace Tensors
{

/** First order tensor (3D vector).
 * \f[ \{\TeFi{v}\} = \VeExpV{v} \f]
 */
typedef blitz::TinyVector<double,3> Tensor1;

/** Second-order symmetric tensor represented as a 6x1 vector (6D vector).
 * Introducing a base of nine second order tensors, known as \em Mandel's \em basis
 * (<a href="http://www.mech.utah.edu/~brannon/gobag.html">Brannon, 2000</a>),
 * second order tensors may be viewed as 9D vectors.
 * If a tensor \f$ \TeSe{a} \f$ is symmetric, i.e. \f$\TeSe{a}^T=\TeSe{a}\f$,
 * then it may be represented by a \em 6D vector:
 * 
 * \f[
 *   [\TeSe{a}] = \MaExpZ{a}
 * 			 \longrightarrow \begin{Bmatrix}
 * 								  a_{0}=a_{00}         \\
 * 								  a_{1}=a_{11}         \\
 * 								  a_{2}=a_{22}         \\
 * 								  a_{3}=\sqrt{2}a_{01} \\
 * 								  a_{4}=\sqrt{2}a_{12} \\
 * 								  a_{5}=\sqrt{2}a_{02}
 * 							 \end{Bmatrix}
 * \f]
 * 
 * Therefore, Tensor2 is a vector with six (6) components,
 * in which the last three components are equal to the \em orthonormal \em Cartesian \em components
 * multiplied by \f$ \sqrt{2} \f$.
 */
typedef blitz::TinyVector<double,6> Tensor2;

/** Fourth-order symmetric tensor represented as a 6x6 matrix (6D second order tensor).
 * Introducing a base of nine second order tensors, known as \em Mandel's \em basis
 * (<a href="http://www.me.unm.edu/~rmbrann/gobag.html">Brannon, 2000</a>),
 * completely symmetric (\f$T_{ijkl}=T_{jikl}=T_{jilk}=T_{klij}\f$) fourth order tensors
 * can be represented as 6D \em second \em order \em tensors:
 *  
 * \f[
 *   [\TeFo{T}] =
 *   \begin{bmatrix}
 * 			  T_{0000}         & T_{0011}         & T_{0022} & \sqrt{2}T_{0001} & \sqrt{2}T_{0012} & \sqrt{2}T_{0002} \\
 * 			  T_{1100}         & T_{1111}         & T_{1122} & \sqrt{2}T_{1101} & \sqrt{2}T_{1112} & \sqrt{2}T_{1102} \\
 * 			  T_{2200}         & T_{2211}         & T_{2222} & \sqrt{2}T_{2201} & \sqrt{2}T_{2212} & \sqrt{2}T_{2202} \\
 * 	  \sqrt{2}T_{0100} & \sqrt{2}T_{0111} & \sqrt{2}T_{0122} &        2T_{0101}        & 2T_{0112}       &  2T_{0102} \\
 * 	  \sqrt{2}T_{1200} & \sqrt{2}T_{1211} & \sqrt{2}T_{1222} &        2T_{1201}        & 2T_{1212}       &  2T_{1202} \\
 * 	  \sqrt{2}T_{0200} & \sqrt{2}T_{0211} & \sqrt{2}T_{0222} &        2T_{0201}        & 2T_{0212}       &  2T_{0202} \\
 *   \end{bmatrix}
 * \f]
 */
typedef blitz::TinyMatrix<double,6,6> Tensor4;

/** Second order identity tensor (symmetric/Mandel's basis) \f$ \TeSe{I} \f$.
 * \f[ \{\TeSe{I}\} =
 *     \begin{Bmatrix} 1 \\ 1 \\ 1 \\ 0 \\ 0 \\ 0 \end{Bmatrix}
 * \f]
 */
Tensor2 I;

/** Forth order identity tensor (symmetric/Mandel's basis) \f$ \Dc{\TeFo{I}^{sym}}{\TeSe{a}}=\TeSe{a}^{sym} \f$.
 * Definition:
 * \f[ \IIsym = \IIsymX \f]
 *
 * Cartesian components:
 * \f[ \left[\IIsym\right] =
 * \begin{bmatrix}
 *      1 & 0 & 0 & 0 & 0 & 0 \\
 *      0 & 1 & 0 & 0 & 0 & 0 \\
 *      0 & 0 & 1 & 0 & 0 & 0 \\
 *      0 & 0 & 0 & 1 & 0 & 0 \\
 *      0 & 0 & 0 & 0 & 1 & 0 \\
 *      0 & 0 & 0 & 0 & 0 & 1
 * \end{bmatrix}
 * \f]
 *
 * Properties:
 * \f[ \IIsym:\TeSe{a} = \TeSe{a} \f]
 */
Tensor4 IIsym;

/** Forth order tensor given by \f$ \Dy{\TeSe{I}}{\TeSe{I}} \f$ (symmetric/Mandel's basis).
 * Cartesian components:
 * \f[ \left[\Dy{\TeSe{I}}{\TeSe{I}}\right] =
 * \begin{bmatrix}
 *     1 &  1 &  1 &  0 &  0 &  0 \\ 
 *     1 &  1 &  1 &  0 &  0 &  0 \\ 
 *     1 &  1 &  1 &  0 &  0 &  0 \\ 
 *     0 &  0 &  0 &  0 &  0 &  0 \\ 
 *     0 &  0 &  0 &  0 &  0 &  0 \\ 
 *     0 &  0 &  0 &  0 &  0 &  0 
 * \end{bmatrix}
 * \f]
 */
Tensor4 IdyI;

/** Forth order \em symmetric-deviatoric tensor (symmetric/Mandel's basis) \f$ \TeFo{P}^{symdev}=\TeFo{I}^{sym}-\frac{I\otimes I}{3} \f$.
 * Definition:
 * \f[ \Psd = \TeFo{I}^{sym}-\frac{I\otimes I}{3} \f]
 * or:
 * \f[ \Psd = \nicefrac{1}{2}(\Leaf{\I}{\I}+\Palm{\I}{\I}) - \nicefrac{1}{3}\Dy{\I}{\I} \f]
 *
 * Cartesian components:
 * \f[ \left[\Psd\right] =
 * \begin{bmatrix}
 *      \nicefrac{2}{3} &  -\nicefrac{1}{3} &  -\nicefrac{1}{3} &  0 &  0 &  0 \\
 *     -\nicefrac{1}{3} &   \nicefrac{2}{3} &  -\nicefrac{1}{3} &  0 &  0 &  0 \\
 *     -\nicefrac{1}{3} &  -\nicefrac{1}{3} &   \nicefrac{2}{3} &  0 &  0 &  0 \\
 *                    0 &                 0 &                 0 &  1 &  0 &  0 \\
 *                    0 &                 0 &                 0 &  0 &  1 &  0 \\
 *                    0 &                 0 &                 0 &  0 &  0 &  1
 * \end{bmatrix}
 * \f]
 *
 * Properties:
 * \f[ \deviator{\TeSe{a}} =\Dc{\Psd}{\TeSe{a}} =\TeSe{a}^{sym} - \TeSe{a}^{iso} =\frac{1}{2}(\TeSe{a}+\TeSe{a}^T) - \frac{\tr{\TeSe{a} }}{3}\I \f]
 */
Tensor4 Psd;

/** Forth order \em isotropic tensor (symmetric/Mandel's basis) \f$ \TeFo{P}^{iso}=\frac{I\otimes I}{3} \f$.
 * Definition:
 * \f[ \Piso = \nicefrac{1}{3}\Dy{\I}{\I} \f]
 *
 * Cartesian components:
 * \f[ \Piso =
 * \begin{bmatrix}
 *     \nicefrac{1}{3} & \nicefrac{1}{3} & \nicefrac{1}{3} & 0 & 0 & 0 \\
 *     \nicefrac{1}{3} & \nicefrac{1}{3} & \nicefrac{1}{3} & 0 & 0 & 0 \\
 *     \nicefrac{1}{3} & \nicefrac{1}{3} & \nicefrac{1}{3} & 0 & 0 & 0 \\
 *                   0 &               0 &               0 & 0 & 0 & 0 \\
 *                   0 &               0 &               0 & 0 & 0 & 0 \\
 *                   0 &               0 &               0 & 0 & 0 & 0
 * \end{bmatrix}
 * \f]
 *
 * Properties:
 * \f[ \isotropic{\TeSe{a}} =\Dc{\Piso}{\TeSe{a}} =\TeSe{a}^{iso} =\frac{\tr{\TeSe{a} }}{3}\I \f]
 */
Tensor4 Piso;

// Constants
const double ZERO   = sqrt(DBL_EPSILON); ///< Machine epsilon (smaller positive)
const double SQ2    = sqrt(2.0);         ///< \f$ \sqrt{2} \f$
const double SQ3    = sqrt(3.0);         ///< \f$ \sqrt{3} \f$
const double SQ6    = sqrt(6.0);         ///< \f$ \sqrt{6} \f$
const double SQ2BY3 = sqrt(2.0/3.0);     ///< \f$ \sqrt{2/3} \f$


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Ignored by Doxygen \cond

/** Initialize (global) second order identity tensor. */
int __initialize_the_I(Tensor2 & theI)
{
	theI = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
	return 0;
}

/** Initialize (global) fourth order identity tensor. */
int __initialize_the_IIsym(Tensor4 & theIIsym)
{
	theIIsym = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
	return 0;
}

/** Initialize (global) IdyI fourth order tensor. */
int __initialize_the_IdyI(Tensor4 & theIdyI)
{
	theIdyI = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	          1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	          1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	          0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	return 0;
}

/** Initialize (global) fourth order symmetric-deviatoric tensor. */
int __initialize_the_Psd(Tensor4 & thePsd)
{
	thePsd =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
	         -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
	         -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0,
	              0.0,      0.0,      0.0, 1.0, 0.0, 0.0,
	              0.0,      0.0,      0.0, 0.0, 1.0, 0.0,
	              0.0,      0.0,      0.0, 0.0, 0.0, 1.0;
	return 0;
}

/** Initialize (global) fourth order isotropic tensor. */
int __initialize_the_Piso(Tensor4 & thePiso)
{
	thePiso = 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	          1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	          1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	              0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
	              0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
	              0.0,     0.0,     0.0, 0.0, 0.0, 0.0;
	return 0;
}

int __dummy1=__initialize_the_I(I);         ///< Dummy variable used only to initialize I
int __dummy2=__initialize_the_IIsym(IIsym); ///< Dummy variable used only to initialize IIsym
int __dummy3=__initialize_the_IdyI(IdyI);   ///< Dummy variable used only to initialize IdyI
int __dummy4=__initialize_the_Psd(Psd);     ///< Dummy variable used only to initialize Psd
int __dummy5=__initialize_the_Piso(Piso);   ///< Dummy variable used only to initialize Piso

// \endcond

}; // namespace Tensors

#endif // MECHSYS_TENSORS_H
