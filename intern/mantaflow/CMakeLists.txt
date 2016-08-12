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

set(INC
	extern
	intern/strings
	intern/manta_pp
	intern/manta_pp/python
	intern/manta_pp/plugin
	intern/manta_pp/pwrapper
	intern/manta_pp/util
	../../source/blender/makesdna
	../../source/blender/blenlib
)

set(INC_SYS
	${PYTHON_INCLUDE_DIRS}
	${ZLIB_INCLUDE_DIRS}
)

set(SRC
	intern/manta_python_API.cpp
	intern/manta_smoke_API.cpp
	intern/SMOKE.cpp
	intern/spectrum.cpp

	extern/manta_python_API.h
	extern/manta_smoke_API.h
	intern/SMOKE.h
	intern/spectrum.h
	intern/strings/shared_script.h
	intern/strings/smoke_script.h
	intern/strings/liquid_script.h

	intern/manta_pp/commonkernels.h
	intern/manta_pp/commonkernels.h.reg
	intern/manta_pp/commonkernels.h.reg.cpp
	intern/manta_pp/conjugategrad.cpp
	intern/manta_pp/conjugategrad.h
	intern/manta_pp/conjugategrad.h.reg
	intern/manta_pp/conjugategrad.h.reg.cpp
	intern/manta_pp/edgecollapse.cpp
	intern/manta_pp/edgecollapse.h
	intern/manta_pp/edgecollapse.h.reg
	intern/manta_pp/edgecollapse.h.reg.cpp
	intern/manta_pp/fastmarch.cpp
	intern/manta_pp/fastmarch.h
	intern/manta_pp/fastmarch.h.reg
	intern/manta_pp/fastmarch.h.reg.cpp
	intern/manta_pp/fileio.cpp
	intern/manta_pp/fileio.h
	intern/manta_pp/fileio.h.reg
	intern/manta_pp/fileio.h.reg.cpp
	intern/manta_pp/fluidsolver.cpp
	intern/manta_pp/fluidsolver.h
	intern/manta_pp/fluidsolver.h.reg
	intern/manta_pp/fluidsolver.h.reg.cpp
	intern/manta_pp/general.cpp
	intern/manta_pp/general.h
	intern/manta_pp/general.h.reg
	intern/manta_pp/general.h.reg.cpp
	intern/manta_pp/gitinfo.h
	intern/manta_pp/grid.cpp
	intern/manta_pp/grid.h
	intern/manta_pp/grid.h.reg
	intern/manta_pp/grid.h.reg.cpp
	intern/manta_pp/grid4d.cpp
	intern/manta_pp/grid4d.h
	intern/manta_pp/grid4d.h.reg
	intern/manta_pp/grid4d.h.reg.cpp
	intern/manta_pp/kernel.cpp
	intern/manta_pp/kernel.h
	intern/manta_pp/kernel.h.reg
	intern/manta_pp/kernel.h.reg.cpp
	intern/manta_pp/levelset.cpp
	intern/manta_pp/levelset.h
	intern/manta_pp/levelset.h.reg
	intern/manta_pp/levelset.h.reg.cpp
	intern/manta_pp/mesh.cpp
	intern/manta_pp/mesh.h
	intern/manta_pp/mesh.h.reg
	intern/manta_pp/mesh.h.reg.cpp
	intern/manta_pp/movingobs.cpp
	intern/manta_pp/movingobs.h
	intern/manta_pp/movingobs.h.reg
	intern/manta_pp/movingobs.h.reg.cpp
	intern/manta_pp/noisefield.cpp
	intern/manta_pp/noisefield.h
	intern/manta_pp/noisefield.h.reg
	intern/manta_pp/noisefield.h.reg.cpp
	intern/manta_pp/particle.cpp
	intern/manta_pp/particle.h
	intern/manta_pp/particle.h.reg
	intern/manta_pp/particle.h.reg.cpp
	intern/manta_pp/plugin/advection.cpp
	intern/manta_pp/plugin/extforces.cpp
	intern/manta_pp/plugin/fire.cpp
	intern/manta_pp/plugin/flip.cpp
	intern/manta_pp/plugin/initplugins.cpp
	intern/manta_pp/plugin/kepsilon.cpp
	intern/manta_pp/plugin/meshplugins.cpp
	intern/manta_pp/plugin/pressure.cpp
	intern/manta_pp/plugin/surfaceturbulence.cpp
	intern/manta_pp/plugin/vortexplugins.cpp
	intern/manta_pp/plugin/waveletturbulence.cpp
	intern/manta_pp/plugin/waves.cpp
	intern/manta_pp/pwrapper/manta.h
	intern/manta_pp/pwrapper/pclass.cpp
	intern/manta_pp/pwrapper/pclass.h
	intern/manta_pp/pwrapper/pconvert.cpp
	intern/manta_pp/pwrapper/pconvert.h
	intern/manta_pp/pwrapper/pvec3.cpp
	intern/manta_pp/pwrapper/pymain.cpp
	intern/manta_pp/pwrapper/pythonInclude.h
	intern/manta_pp/pwrapper/registry.cpp
	intern/manta_pp/pwrapper/registry.h
	intern/manta_pp/python/defines.py
	intern/manta_pp/python/defines.py.reg
	intern/manta_pp/python/defines.py.reg.cpp
	intern/manta_pp/shapes.cpp
	intern/manta_pp/shapes.h
	intern/manta_pp/shapes.h.reg
	intern/manta_pp/shapes.h.reg.cpp
	intern/manta_pp/test.cpp
	intern/manta_pp/timing.cpp
	intern/manta_pp/timing.h
	intern/manta_pp/timing.h.reg
	intern/manta_pp/timing.h.reg.cpp
	intern/manta_pp/turbulencepart.cpp
	intern/manta_pp/turbulencepart.h
	intern/manta_pp/turbulencepart.h.reg
	intern/manta_pp/turbulencepart.h.reg.cpp
	intern/manta_pp/util/integrator.h
	intern/manta_pp/util/interpol.h
	intern/manta_pp/util/interpolHigh.h
	intern/manta_pp/util/mcubes.h
	intern/manta_pp/util/quaternion.h
	intern/manta_pp/util/randomstream.h
	intern/manta_pp/util/simpleimage.cpp
	intern/manta_pp/util/simpleimage.h
	intern/manta_pp/util/solvana.h
	intern/manta_pp/util/vector4d.cpp
	intern/manta_pp/util/vector4d.h
	intern/manta_pp/util/vectorbase.cpp
	intern/manta_pp/util/vectorbase.h
	intern/manta_pp/vortexpart.cpp
	intern/manta_pp/vortexpart.h
	intern/manta_pp/vortexpart.h.reg
	intern/manta_pp/vortexpart.h.reg.cpp
	intern/manta_pp/vortexsheet.cpp
	intern/manta_pp/vortexsheet.h
	intern/manta_pp/vortexsheet.h.reg
	intern/manta_pp/vortexsheet.h.reg.cpp
)

if(WITH_OPENMP)
	add_definitions(-DOPENMP)
else()
	message("Error. Mantaflow pp files require an OpenMP enabled compiler")
endif()

blender_add_lib(bf_intern_mantaflow "${SRC}" "${INC}" "${INC_SYS}")

