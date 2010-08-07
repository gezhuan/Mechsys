########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

SET (FLAGS  "${FLAGS}")
SET (LIBS    ${LIBS})
SET (LFLAGS "${LFLAGS}")
SET (DEPS_OK TRUE)

SET     (Boost_USE_STATIC_LIBS ON)
INCLUDE (FindBoost)
if(Boost_FOUND)
	INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    SET (DEPS_OK FALSE)
endif(Boost_FOUND)

INCLUDE (FindThreads)
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindBLITZ.cmake    )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindGSL.cmake      )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindHDF5.cmake     )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindMTL.cmake      )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTETGEN.cmake   )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTRIANGLE.cmake )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindUMFPACK.cmake  )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindMUMPS.cmake    )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSCALAPACK.cmake)
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindVORO.cmake     )

if(BLITZ_FOUND)
	INCLUDE_DIRECTORIES (${BLITZ_INCLUDE_DIRS})
else(BLITZ_FOUND)
    SET (DEPS_OK FALSE)
endif(BLITZ_FOUND)

ENABLE_LANGUAGE (Fortran)
INCLUDE         (FindLAPACK)
if(LAPACK_FOUND)
    SET (LIBS ${LIBS} ${LAPACK_LIBRARIES})
else(LAPACK_FOUND)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindLocLAPACK.cmake)
    if(LocLAPACK_FOUND)
        SET (LIBS ${LIBS} ${LocLAPACK_LIBRARIES} "gfortran")
    else(LocLAPACK_FOUND)
        SET (DEPS_OK FALSE)
    endif(LocLAPACK_FOUND)
endif(LAPACK_FOUND)

if(GSL_FOUND)
	INCLUDE_DIRECTORIES (${GSL_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${GSL_LIBRARIES})
else(GSL_FOUND)
    SET (DEPS_OK FALSE)
endif(GSL_FOUND)

if(HDF5_FOUND)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
else(HDF5_FOUND)
    SET (DEPS_OK FALSE)
endif(HDF5_FOUND)

if(MTL_FOUND)
	INCLUDE_DIRECTORIES (${MTL_INCLUDE_DIRS})
else(MTL_FOUND)
    SET (DEPS_OK FALSE)
endif(MTL_FOUND)

if(TETGEN_FOUND)
	INCLUDE_DIRECTORIES (${TETGEN_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TETGEN_LIBRARIES})
else(TETGEN_FOUND)
    SET (DEPS_OK FALSE)
endif(TETGEN_FOUND)

if(TRIANGLE_FOUND)
	INCLUDE_DIRECTORIES (${TRIANGLE_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TRIANGLE_LIBRARIES})
else(TRIANGLE_FOUND)
    SET (DEPS_OK FALSE)
endif(TRIANGLE_FOUND)

if(UMFPACK_FOUND)
	INCLUDE_DIRECTORIES (${UMFPACK_INCLUDE_DIRS})
	SET(LIBS ${LIBS} ${UMFPACK_LIBRARIES})
else(UMFPACK_FOUND)
    SET (DEPS_OK FALSE)
endif(UMFPACK_FOUND)

if(VORO_FOUND)
	INCLUDE_DIRECTORIES (${VORO_INCLUDE_DIRS})
else(VORO_FOUND)
    SET (DEPS_OK FALSE)
endif(VORO_FOUND)

# Options ##############################################################

if(MUMPS_FOUND)
	INCLUDE_DIRECTORIES (${MUMPS_INCLUDE_DIRS})
	SET(LIBS ${LIBS} ${MUMPS_LIBRARIES} -lmpi_f77)
else(MUMPS_FOUND)
    MESSAGE ("Could not find MUMPS. OK, compiling without it")
endif(MUMPS_FOUND)

if(SCALAPACK_FOUND)
	SET(LIBS ${LIBS} ${SCALAPACK_LIBRARIES})
else(SCALAPACK_FOUND)
    MESSAGE ("Could not find SCALAPACK. OK, compiling without it")
endif(SCALAPACK_FOUND)

if(A_USE_MPI OR A_USE_PMETIS)
    INCLUDE(FindMPI)
    if(MPI_FOUND)
        INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})
        SET (FLAGS  "${FLAGS}  ${MPI_COMPILE_FLAGS}")
        SET (LIBS    ${LIBS}   ${MPI_LIBRARIES}     )
        SET (LFLAGS "${LFLAGS} ${MPI_LINK_FLAGS}"   )
        ADD_DEFINITIONS (-DUSE_MPI)
    else(MPI_FOUND)
        MESSAGE ("Could not find OpenMPI. OK, compiling without it")
    endif(MPI_FOUND)
endif(A_USE_MPI OR A_USE_PMETIS)

if(A_USE_VTK)
    INCLUDE (FindVTK)
    if(VTK_FOUND)
        INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIRS})
        SET (LIBS  ${LIBS} vtkRendering)
        SET (FLAGS "${FLAGS} -DVTK_EXCLUDE_STRSTREAM_HEADERS")
        ADD_DEFINITIONS (-DUSE_VTK)
    else(VTK_FOUND)
        MESSAGE ("Could not find VTK. OK, compiling without it")
    endif(VTK_FOUND)
endif(A_USE_VTK)

if (A_USE_WXW)
    if(A_MAKE_WX_MONO)
        FIND_PACKAGE (wxWidgets COMPONENTS base core)
    else(A_MAKE_WX_MONO)
        FIND_PACKAGE (wxWidgets COMPONENTS base core aui)
    endif(A_MAKE_WX_MONO)
    if(wxWidgets_FOUND)
        INCLUDE (${wxWidgets_USE_FILE})
        SET (LIBS ${LIBS} ${wxWidgets_LIBRARIES})
        ADD_DEFINITIONS (-DUSE_WXW)
    else(wxWidgets_FOUND)
        MESSAGE ("Could not find wxWidgets. OK, compiling without it")
    endif(wxWidgets_FOUND)
endif (A_USE_WXW)

if(A_USE_PMETIS)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindPARMETIS.cmake)
    if(PARMETIS_FOUND)
        INCLUDE_DIRECTORIES (${PARMETIS_INCLUDE_DIRS})
        SET (LIBS ${LIBS} ${PARMETIS_LIBRARIES})
        ADD_DEFINITIONS (-DUSE_PMETIS)
    else(PARMETIS_FOUND)
        MESSAGE ("Could not find ParMETIS. OK, compiling without it")
    endif(PARMETIS_FOUND)
endif(A_USE_PMETIS)
