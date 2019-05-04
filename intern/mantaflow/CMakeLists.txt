# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# The Original Code is Copyright (C) 2016, Blender Foundation
# All rights reserved.
#
# The Original Code is: all of this file.
#
# Contributor(s): Sebastian Barschkis (sebbas)
#
# ***** END GPL LICENSE BLOCK *****

set(MANTAVERSION "0.12")

if(WITH_OPENMP)
  set(MANTA_PP
    intern/manta_pp/omp
  )
  add_definitions(-DOPENMP=1)
else()
  set(MANTA_PP
    intern/manta_pp/tbb
  )
  add_definitions(-DTBB=1)
endif()

if(WITH_OPENVDB)
  add_definitions(-DOPENVDB=1)
endif()

set(INC
  extern
  intern/strings
  ${MANTA_PP}
  ${MANTA_PP}/fileio
  ${MANTA_PP}/python
  ${MANTA_PP}/plugin
  ${MANTA_PP}/pwrapper
  ${MANTA_PP}/util
  ../../source/blender/makesdna
  ../../source/blender/blenlib
)

if(WITH_PYTHON)
  add_definitions(-DWITH_PYTHON)
endif()

set(INC_SYS
  ${PYTHON_INCLUDE_DIRS}
  ${ZLIB_INCLUDE_DIRS}
)

if(NOT WITH_OPENMP)
  list(APPEND INC_SYS
    ${TBB_INCLUDE_DIRS}
  )
endif()

if(WITH_OPENVDB)
  list(APPEND INC_SYS
    ${BOOST_INCLUDE_DIR}
    ${OPENEXR_INCLUDE_DIRS}
    ${OPENVDB_INCLUDE_DIRS}
  )
endif()

set(SRC
  intern/manta_python_API.cpp
  intern/manta_fluid_API.cpp
  intern/FLUID.cpp

  extern/manta_python_API.h
  extern/manta_fluid_API.h
  intern/FLUID.h
  intern/strings/fluid_script.h
  intern/strings/smoke_script.h
  intern/strings/liquid_script.h

  ${MANTA_PP}/commonkernels.h
  ${MANTA_PP}/commonkernels.h.reg
  ${MANTA_PP}/commonkernels.h.reg.cpp
  ${MANTA_PP}/conjugategrad.cpp
  ${MANTA_PP}/conjugategrad.h
  ${MANTA_PP}/conjugategrad.h.reg
  ${MANTA_PP}/conjugategrad.h.reg.cpp
  ${MANTA_PP}/edgecollapse.cpp
  ${MANTA_PP}/edgecollapse.h
  ${MANTA_PP}/edgecollapse.h.reg
  ${MANTA_PP}/edgecollapse.h.reg.cpp
  ${MANTA_PP}/fastmarch.cpp
  ${MANTA_PP}/fastmarch.h
  ${MANTA_PP}/fastmarch.h.reg
  ${MANTA_PP}/fastmarch.h.reg.cpp
  ${MANTA_PP}/fileio/iogrids.cpp
  ${MANTA_PP}/fileio/iomeshes.cpp
  ${MANTA_PP}/fileio/ioparticles.cpp
  ${MANTA_PP}/fileio/mantaio.h
  ${MANTA_PP}/fileio/mantaio.h.reg
  ${MANTA_PP}/fileio/mantaio.h.reg.cpp
  ${MANTA_PP}/fluidsolver.cpp
  ${MANTA_PP}/fluidsolver.h
  ${MANTA_PP}/fluidsolver.h.reg
  ${MANTA_PP}/fluidsolver.h.reg.cpp
  ${MANTA_PP}/general.cpp
  ${MANTA_PP}/general.h
  ${MANTA_PP}/general.h.reg
  ${MANTA_PP}/general.h.reg.cpp
  ${MANTA_PP}/gitinfo.h
  ${MANTA_PP}/grid.cpp
  ${MANTA_PP}/grid.h
  ${MANTA_PP}/grid.h.reg
  ${MANTA_PP}/grid.h.reg.cpp
  ${MANTA_PP}/grid4d.cpp
  ${MANTA_PP}/grid4d.h
  ${MANTA_PP}/grid4d.h.reg
  ${MANTA_PP}/grid4d.h.reg.cpp
  ${MANTA_PP}/kernel.cpp
  ${MANTA_PP}/kernel.h
  ${MANTA_PP}/kernel.h.reg
  ${MANTA_PP}/kernel.h.reg.cpp
  ${MANTA_PP}/levelset.cpp
  ${MANTA_PP}/levelset.h
  ${MANTA_PP}/levelset.h.reg
  ${MANTA_PP}/levelset.h.reg.cpp
  ${MANTA_PP}/mesh.cpp
  ${MANTA_PP}/mesh.h
  ${MANTA_PP}/mesh.h.reg
  ${MANTA_PP}/mesh.h.reg.cpp
  ${MANTA_PP}/movingobs.cpp
  ${MANTA_PP}/movingobs.h
  ${MANTA_PP}/movingobs.h.reg
  ${MANTA_PP}/movingobs.h.reg.cpp
  ${MANTA_PP}/multigrid.cpp
  ${MANTA_PP}/multigrid.h
  ${MANTA_PP}/multigrid.h.reg
  ${MANTA_PP}/multigrid.h.reg.cpp
  ${MANTA_PP}/noisefield.cpp
  ${MANTA_PP}/noisefield.h
  ${MANTA_PP}/noisefield.h.reg
  ${MANTA_PP}/noisefield.h.reg.cpp
  ${MANTA_PP}/particle.cpp
  ${MANTA_PP}/particle.h
  ${MANTA_PP}/particle.h.reg
  ${MANTA_PP}/particle.h.reg.cpp
  ${MANTA_PP}/plugin/advection.cpp
  ${MANTA_PP}/plugin/apic.cpp
  ${MANTA_PP}/plugin/extforces.cpp
  ${MANTA_PP}/plugin/fire.cpp
  ${MANTA_PP}/plugin/flip.cpp
  ${MANTA_PP}/plugin/fluidguiding.cpp
  ${MANTA_PP}/plugin/initplugins.cpp
  ${MANTA_PP}/plugin/kepsilon.cpp
  ${MANTA_PP}/plugin/meshplugins.cpp
#   ${MANTA_PP}/plugin/numpyconvert.cpp
  ${MANTA_PP}/plugin/pressure.cpp
  ${MANTA_PP}/plugin/ptsplugins.cpp
  ${MANTA_PP}/plugin/secondaryparticles
  ${MANTA_PP}/plugin/sndparticles.cpp
  ${MANTA_PP}/plugin/surfaceturbulence.cpp
#   ${MANTA_PP}/plugin/tfplugins.cpp
  ${MANTA_PP}/plugin/vortexplugins.cpp
  ${MANTA_PP}/plugin/waveletturbulence.cpp
  ${MANTA_PP}/plugin/waves.cpp
  ${MANTA_PP}/pwrapper/manta.h
#   ${MANTA_PP}/pwrapper/numpyWrap.cpp
#   ${MANTA_PP}/pwrapper/numpyWrap.h
  ${MANTA_PP}/pwrapper/pclass.cpp
  ${MANTA_PP}/pwrapper/pclass.h
  ${MANTA_PP}/pwrapper/pconvert.cpp
  ${MANTA_PP}/pwrapper/pconvert.h
  ${MANTA_PP}/pwrapper/pvec3.cpp
  ${MANTA_PP}/pwrapper/pythonInclude.h
  ${MANTA_PP}/pwrapper/registry.cpp
  ${MANTA_PP}/pwrapper/registry.h
  ${MANTA_PP}/python/defines.py
  ${MANTA_PP}/python/defines.py.reg
  ${MANTA_PP}/python/defines.py.reg.cpp
  ${MANTA_PP}/registration.cpp
  ${MANTA_PP}/shapes.cpp
  ${MANTA_PP}/shapes.h
  ${MANTA_PP}/shapes.h.reg
  ${MANTA_PP}/shapes.h.reg.cpp
  ${MANTA_PP}/test.cpp
  ${MANTA_PP}/timing.cpp
  ${MANTA_PP}/timing.h
  ${MANTA_PP}/timing.h.reg
  ${MANTA_PP}/timing.h.reg.cpp
  ${MANTA_PP}/turbulencepart.cpp
  ${MANTA_PP}/turbulencepart.h
  ${MANTA_PP}/turbulencepart.h.reg
  ${MANTA_PP}/turbulencepart.h.reg.cpp
  ${MANTA_PP}/util/integrator.h
  ${MANTA_PP}/util/interpol.h
  ${MANTA_PP}/util/interpolHigh.h
  ${MANTA_PP}/util/matrixbase.h
  ${MANTA_PP}/util/mcubes.h
  ${MANTA_PP}/util/quaternion.h
  ${MANTA_PP}/util/randomstream.h
  ${MANTA_PP}/util/rcmatrix.h
  ${MANTA_PP}/util/simpleimage.cpp
  ${MANTA_PP}/util/simpleimage.h
  ${MANTA_PP}/util/solvana.h
  ${MANTA_PP}/util/vector4d.cpp
  ${MANTA_PP}/util/vector4d.h
  ${MANTA_PP}/util/vectorbase.cpp
  ${MANTA_PP}/util/vectorbase.h
  ${MANTA_PP}/vortexpart.cpp
  ${MANTA_PP}/vortexpart.h
  ${MANTA_PP}/vortexpart.h.reg
  ${MANTA_PP}/vortexpart.h.reg.cpp
  ${MANTA_PP}/vortexsheet.cpp
  ${MANTA_PP}/vortexsheet.h
  ${MANTA_PP}/vortexsheet.h.reg
  ${MANTA_PP}/vortexsheet.h.reg.cpp
)

set(LIB
)

blender_add_lib(bf_intern_mantaflow "${SRC}" "${INC}" "${INC_SYS}" "${LIB}")
