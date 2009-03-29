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

#define USE_BOOST_PYTHON

// Std Lib
#include <iostream>
#include <cfloat>    // for DBL_EPSILON

// Boost-Python
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module_init.hpp>
#include <boost/python/def.hpp>
#include <boost/python/call_method.hpp>
#include <boost/ref.hpp>
#include <boost/utility.hpp>
#include <boost/python.hpp> // this includes everything

// MechSys
#include "mechsys.h"

using namespace boost::python;

// Overloadings                                                        minargs  maxargs
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MB_SetNx,          SetNx,            1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MB_SetNy,          SetNy,            1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MB_SetNz,          SetNz,            1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_SetVert,        SetVert,          4, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_Generate,       Generate,         0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_GenBox,         GenBox,           3, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_SetPolySize,    SetPolySize,      2, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_SetPolyPoint,   SetPolyPoint,     3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_SetPolySegment, SetPolySegment,   3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_SetPolyRegion,  SetPolyRegion,    5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_SetPolyHole,    SetPolyHole,      3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_Generate,       Generate,         0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MA_AddCloudPoint,  AddCloudPoint,    2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MA_Generate,       Generate,         0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (GE_SetOnlyFrame,   SetOnlyFrame,     0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_Solve,          Solve,            0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_SolveWithInfo,  PySolveWithInfo,  4, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DA_OutElems,       OutElems,         0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DA_OutNodes,       OutNodes,         0, 1)

BOOST_PYTHON_MODULE (mechsys)
{
	// --------------------------------------------------------------------------- Mesh

	class_<Mesh::Generic>("mesh_generic","Class to hold mesh data",init<bool>())
	    .def("is_3d",         &Mesh::Generic::Is3D)
	    .def("set_o2",        &Mesh::Generic::SetO2)
	    .def("set_nverts",    &Mesh::Generic::SetNVerts)
	    .def("set_nelems",    &Mesh::Generic::SetNElems)
	    .def("set_vert",      &Mesh::Generic::SetVert, MG_SetVert())
	    .def("set_elem",      &Mesh::Generic::SetElem)
	    .def("set_elem_con",  &Mesh::Generic::SetElemCon)
	    .def("set_elem_etag", &Mesh::Generic::SetElemETag)
	    .def("set_elem_ftag", &Mesh::Generic::SetElemFTag)
	    .def("write_vtu",     &Mesh::Generic::PyWriteVTU)
	    .def("get_verts",     &Mesh::Generic::PyGetVerts)
	    .def("get_verts_bry", &Mesh::Generic::PyGetVertsBry)
	    .def("get_edges",     &Mesh::Generic::PyGetEdges)
	    .def("get_etags",     &Mesh::Generic::PyGetETags)
	    .def("get_ftags",     &Mesh::Generic::PyGetFTags)
	    .def("get_elems",     &Mesh::Generic::PyGetElems)
	    .def(self_ns::str(self))
	    ;

	class_<Mesh::Block>("mesh_block")
	    .def("set",    &Mesh::Block::PySet)
	    .def("set_nx", &Mesh::Block::SetNx, MB_SetNx())
	    .def("set_ny", &Mesh::Block::SetNy, MB_SetNy())
	    .def("set_nz", &Mesh::Block::SetNz, MB_SetNz())
	    ;

	class_<Mesh::Structured, bases<Mesh::Generic> >("mesh_structured","Class to generate structured meshes",init<bool>())
	    .def("set_blocks", &Mesh::Structured::PySetBlocks)
	    .def("generate",   &Mesh::Structured::Generate, MS_Generate())
	    .def("gen_box",    &Mesh::Structured::GenBox,   MS_GenBox())
	    .def(self_ns::str(self))
	    ;

	class_<Mesh::Unstructured, bases<Mesh::Generic> >("mesh_unstructured","Class to generate unstructured meshes",init<bool>())
	    .def("set_poly_size",        &Mesh::Unstructured::SetPolySize,    MU_SetPolySize())
	    .def("set_poly_point",       &Mesh::Unstructured::SetPolyPoint,   MU_SetPolyPoint())
	    .def("set_poly_segment",     &Mesh::Unstructured::SetPolySegment, MU_SetPolySegment())
	    .def("set_poly_region",      &Mesh::Unstructured::SetPolyRegion,  MU_SetPolyRegion())
	    .def("set_poly_hole",        &Mesh::Unstructured::SetPolyHole,    MU_SetPolyHole())
	    .def("set_max_area_global",  &Mesh::Unstructured::SetMaxAreaGlobal)
	    .def("set_min_angle_global", &Mesh::Unstructured::SetMinAngleGlobal)
	    .def("generate",             &Mesh::Unstructured::Generate,       MU_Generate())
	    .def(self_ns::str(self))
	    ;

	class_<Mesh::AlphaShape, bases<Mesh::Generic> >("mesh_alpha_shape","Class to generate meshes based on alpha-shape",init<bool>())
	    .def("reset_cloud",     &Mesh::AlphaShape::ResetCloud)
	    .def("add_cloud_point", &Mesh::AlphaShape::AddCloudPoint, MA_AddCloudPoint())
	    .def("generate",        &Mesh::AlphaShape::Generate,      MA_Generate())
	    .def(self_ns::str(self))
	    ;

	// ---------------------------------------------------------------------------- FEM
	
	class_<FEM::Node>("node")
	    .def("id",      &FEM::Node::ID)
	    .def("tag",     &FEM::Node::Tag)
	    .def("x",       &FEM::Node::X)
	    .def("y",       &FEM::Node::Y)
	    .def("z",       &FEM::Node::Z)
	    .def("bry",     &FEM::Node::PyBry, return_internal_reference<>())
	    .def("val",     &FEM::Node::PyVal)
	    .def(self_ns::str(self))
	    ;

	class_<FEM::Element>("elem")
	    .def("id",                &FEM::Element::ID)
	    .def("tag",               &FEM::Element::Tag)
	    .def("nnodes",            &FEM::Element::NNodes)
	    .def("nod",               &FEM::Element::PyNod1, return_internal_reference<>())
	    .def("nod",               &FEM::Element::PyNod2, return_internal_reference<>())
	    .def("is_active",         &FEM::Element::IsActive)
	    .def("calc_deps",         &FEM::Element::CalcDeps)
	    .def("val",               &FEM::Element::PyVal1)
	    .def("val",               &FEM::Element::PyVal2)
	    .def("edge_bry",          &FEM::Element::PyEdgeBry1)
	    .def("edge_bry",          &FEM::Element::PyEdgeBry2)
	    .def("face_bry",          &FEM::Element::PyFaceBry)
	    .def("calc_deps",         &FEM::Element::CalcDeps)
	    .def("has_extra",         &FEM::Element::HasExtra)
	    .def("out_extra",         &FEM::Element::PyOutExtra)
	    .def(self_ns::str(self))
	    ;

	class_<FEM::Data>("data", init<int>())
	    .def("set_tol",             &FEM::Data::SetTol)
	    .def("set_only_frame",      &FEM::Data::SetOnlyFrame, GE_SetOnlyFrame())
		.def("set_nodes_elems",     &FEM::Data::PySetNodesElems)
		.def("set_brys",            &FEM::Data::PySetBrys)
	    .def("add_lin_elems",       &FEM::Data::PyAddLinElems)
	    .def("add_vol_forces",      &FEM::Data::AddVolForces)
	    .def("clear_disp",          &FEM::Data::ClearDisp)
	    .def("activate",            &FEM::Data::Activate)
	    .def("deactivate",          &FEM::Data::Deactivate)
	    .def("nnodes",              &FEM::Data::NNodes)
	    .def("nelems",              &FEM::Data::NElems)
	    .def("nod",                 &FEM::Data::PyNod, return_internal_reference<>())
	    .def("ele",                 &FEM::Data::PyEle, return_internal_reference<>())
	    .def("elems_with_tag",      &FEM::Data::PyElemsWithTag)
	    .def("bounds_2d",           &FEM::Data::PyBounds2D)
	    .def("bounds_3d",           &FEM::Data::PyBounds3D)
	    .def("set_out_elems",       &FEM::Data::PySetOutElems)
	    .def("set_out_nodes",       &FEM::Data::PySetOutNodes)
	    .def("out_elems",           &FEM::Data::OutElems, DA_OutElems())
	    .def("out_nodes",           &FEM::Data::OutNodes, DA_OutNodes())
	    .def(self_ns::str(self))
	    ;

	class_<FEM::Solver>("solver", init<FEM::Data &>())
		.def(init<FEM::Data &, BPy::str const &>())
	    .def("set_type",        &FEM::Solver::PySetType,        return_internal_reference<>())
	    .def("set_cte",         &FEM::Solver::PySetCte,         return_internal_reference<>())
	    .def("set_lin_sol",     &FEM::Solver::PySetLinSol,      return_internal_reference<>())
	    .def("solve",           &FEM::Solver::Solve,            SO_Solve())
	    .def("solve_with_info", &FEM::Solver::PySolveWithInfo,  SO_SolveWithInfo())
	    .def("time",            &FEM::Solver::Time)
	    .def("ndof",            &FEM::Solver::nDOF)
	    .def("out",             &FEM::Solver::Out,              return_internal_reference<>())
	    ;

	class_<Output>("output", init<FEM::Data const &, BPy::str const &>())
	    .def("write",      &Output::Write)
	    .def("val",        &Output::PyVal)
	    .def("get_labels", &Output::PyGetLabels)
	    ;

	// ---------------------------------------------------------------------- Exceptions
	
	register_exception_translator<Exception *>(&PyExceptTranslator);
}
