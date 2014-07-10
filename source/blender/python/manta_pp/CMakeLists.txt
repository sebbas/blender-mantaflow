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
# The Original Code is Copyright (C) 2012, Blender Foundation
# All rights reserved.
#
# The Original Code is: all of this file.
#
# Contributor(s): .
#
# ***** END GPL LICENSE BLOCK *****

set(INC
	cuda
	preprocessor
	pwrapper
	python
	plugin
	util
)

set(INC_SYS
	${PYTHON_INCLUDE_DIRS}
)
#blender_include_dirs_sys("${PYTHON_INCLUDE_DIRS}")
add_definitions(-DWITH_PYTHON)

if(WITH_PYTHON_MODULE)
	add_definitions(-DWITH_PYTHON_MODULE)
endif()

if(WITH_PYTHON_SAFETY)
	add_definitions(-DWITH_PYTHON_SAFETY)
endif()

set(SRC
commonkernels.h
commonkernels.h.reg
commonkernels.h.reg.cpp
conjugategrad.cpp
conjugategrad.h
conjugategrad.h.reg
conjugategrad.h.reg.cpp
edgecollapse.cpp
edgecollapse.h
edgecollapse.h.reg
edgecollapse.h.reg.cpp
fastmarch.cpp
fastmarch.h
fastmarch.h.reg
fastmarch.h.reg.cpp
fileio.cpp
fileio.h
fileio.h.reg
fileio.h.reg.cpp
flip.cpp
flip.h
flip.h.reg
flip.h.reg.cpp
fluidsolver.cpp
fluidsolver.h
fluidsolver.h.reg
fluidsolver.h.reg.cpp
general.cpp
general.h
general.h.reg
general.h.reg.cpp
grid.cpp
grid.h
grid.h.reg
grid.h.reg.cpp
hginfo.h
kernel.cpp
kernel.h
kernel.h.reg
kernel.h.reg.cpp
levelset.cpp
levelset.h
levelset.h.reg
levelset.h.reg.cpp
manta.h
mesh.cpp
mesh.h
mesh.h.reg
mesh.h.reg.cpp
movingobs.cpp
movingobs.h
movingobs.h.reg
movingobs.h.reg.cpp
noisefield.cpp
noisefield.h
noisefield.h.reg
noisefield.h.reg.cpp
particle.cpp
particle.h
particle.h.reg
particle.h.reg.cpp
shapes.cpp
shapes.h
shapes.h.reg
shapes.h.reg.cpp
test.cpp
turbulencepart.cpp
turbulencepart.h
turbulencepart.h.reg
turbulencepart.h.reg.cpp
##vortexfilament.cpp
##vortexfilament.h
vortexpart.cpp
vortexpart.h
vortexpart.h.reg
vortexpart.h.reg.cpp
vortexsheet.cpp
vortexsheet.h
vortexsheet.h.reg
vortexsheet.h.reg.cpp

cuda/buoyancy.cu
cuda/cudatools.h
cuda/curlnoise.cu
cuda/curlnoise.h
cuda/meshtools.cu
cuda/particle.cu
cuda/turbulence.cu

#gui/customctrl.cpp
#gui/customctrl.h
#gui/customctrl.h.reg.cpp
#gui/glwidget.cpp
#gui/glwidget.h
#gui/glwidget.h.reg
#gui/glwidget.h.reg.cpp
#gui/mainwindow.cpp
#gui/mainwindow.h
#gui/mainwindow.h.reg
#gui/mainwindow.h.reg.cpp
#gui/meshpainter.cpp
#gui/meshpainter.h
#gui/meshpainter.h.reg
#gui/meshpainter.h.reg.cpp
#gui/moc_customctrl.cxx
#gui/moc_customctrl.cxx_parameters
#gui/moc_glwidget.cxx
#gui/moc_glwidget.cxx_parameters
#gui/moc_mainwindow.cxx
#gui/moc_mainwindow.cxx_parameters
#gui/moc_meshpainter.cxx
#gui/moc_meshpainter.cxx_parameters
#gui/moc_painter.cxx
#gui/moc_painter.cxx_parameters
#gui/moc_particlepainter.cxx
#gui/moc_particlepainter.cxx_parameters
#gui/moc_qtmain.cxx
#gui/moc_qtmain.cxx_parameters
#gui/painter.cpp
#gui/painter.h
#gui/painter.h.reg
#gui/painter.h.reg.cpp
#gui/particlepainter.cpp
#gui/particlepainter.h
#gui/particlepainter.h.reg
#gui/particlepainter.h.reg.cpp
#gui/qtmain.cpp
#gui/qtmain.h
#gui/qtmain.h.reg
#gui/qtmain.h.reg.cpp

plugin/advection.cpp
plugin/extforces.cpp
plugin/initplugins.cpp
plugin/kepsilon.cpp
plugin/meshplugins.cpp
plugin/pressure.cpp
plugin/vortexplugins.cpp
plugin/waveletturbulence.cpp

preprocessor/code.cpp
preprocessor/code.h
preprocessor/codegen_kernel.cpp
preprocessor/codegen_python.cpp
preprocessor/main.cpp
preprocessor/main.h
preprocessor/merge.cpp
preprocessor/parse.cpp
preprocessor/prep.h
preprocessor/tokenize.cpp
preprocessor/tokenize.h
preprocessor/util.cpp
preprocessor/util.h

pwrapper/manta.h
pwrapper/pclass.cpp
pwrapper/pclass.h
pwrapper/pconvert.cpp
pwrapper/pconvert.h
pwrapper/pvec3.cpp
pwrapper/pymain.cpp
pwrapper/pythonInclude.h
pwrapper/registry.cpp
pwrapper/registry.h
pwrapper/manta_api.h
pwrapper/manta_api.cpp

python/defines.py.reg.cpp

util/integrator.h
util/interpol.h
util/interpolHigh.h
util/mcubes.h
util/quaternion.h
util/randomstream.h
util/solvana.h
util/vectorbase.cpp
util/vectorbase.h
)

blender_add_lib(bf_python_manta "${SRC}" "${INC}" "${INC_SYS}")