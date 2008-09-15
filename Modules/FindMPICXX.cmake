#####################################################################################
# MechSys - A C++ library to simulate (Continuum) Mechanical Systems                #
# Copyright (C) 2005 Dorival de Moraes Pedroso <dorival.pedroso at gmail.com>       #
# Copyright (C) 2005 Raul Dario Durand Farfan  <raul.durand at gmail.com>           #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

SET(MPICXX_PROGRAM_SEARCH_PATH
  $ENV{HOME}/opt/bin
  /usr/bin
  /usr/local/bin)

IF(WIN32)
  SET(MPICXX_PROGRAM_SEARCH_PATH
    ${MPICXX_PROGRAM_SEARCH_PATH}
    C:/MinGW/bin
    D:/MinGW/bin)
ENDIF(WIN32)

FIND_PROGRAM(MPICXX_MPICXX NAMES mpic++ PATHS ${MPICXX_PROGRAM_SEARCH_PATH})

SET(MPICXX_FOUND 1)
FOREACH(var MPICXX_MPICXX)
  IF(NOT ${var})
	SET(MPICXX_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MPICXX_FOUND)
  SET(MPICXX ${MPICXX_MPICXX})
ENDIF(MPICXX_FOUND)
