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

SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)


if(WITH_PYTHON_MODULE)
	add_definitions(-DWITH_PYTHON_MODULE)
endif()

if(WITH_PYTHON_SAFETY)
	add_definitions(-DWITH_PYTHON_SAFETY)
endif()

#=========================================
#original standalone manta setup 
#=========================================
#SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
#SET(VERBOSE 1) 
#
#SET(SOURCES 
#	preprocessor/main.cpp
#	preprocessor/code.cpp
#	preprocessor/tokenize.cpp
#	preprocessor/parse.cpp
#	preprocessor/util.cpp
#	preprocessor/merge.cpp
#	preprocessor/codegen_python.cpp
#	preprocessor/codegen_kernel.cpp    
#)
#add_executable(prep ${SOURCES})
#if (NOT WIN32)
#	set_target_properties(prep PROPERTIES COMPILE_FLAGS "-Wall -O2")
#endif()
#
##******************************************************************************
## Setup main project
#
#SET(F_LIBS "" )
#SET(F_LIB_PATHS)
#SET(F_LINKADD "") # additional linker flags, not a list
#set(PP_PATH "pp")
#SET(SILENT_SOURCES)
#SET(EXECCMD manta)
#
#
## need pre-processing
#SET(PP_SOURCES
#	fluidsolver.cpp
#	conjugategrad.cpp
#	grid.cpp
#	levelset.cpp
#	fastmarch.cpp
#	shapes.cpp
#	mesh.cpp
#	general.cpp
#	particle.cpp
#	flip.cpp
#	movingobs.cpp
#	fileio.cpp   
#	noisefield.cpp
#	kernel.cpp
#	vortexsheet.cpp
##	vortexfilament.cpp
#	vortexpart.cpp
#	turbulencepart.cpp
#	edgecollapse.cpp
#	plugin/pressure.cpp
#	plugin/initplugins.cpp
#	plugin/advection.cpp
#	plugin/extforces.cpp
#	plugin/kepsilon.cpp
#	plugin/meshplugins.cpp
#	plugin/vortexplugins.cpp
#	plugin/waveletturbulence.cpp
#	#python/defines.py
#)
#
#SET(PP_HEADERS
#	commonkernels.h
#	conjugategrad.h
#	fastmarch.h
#	fluidsolver.h
#	grid.h
#	mesh.h
#	particle.h
#	levelset.h
#	shapes.h
#	noisefield.h
#	vortexsheet.h
##	vortexfilament.h
#	general.h
#	kernel.h
#	flip.h
#	movingobs.h
#	fileio.h
#	edgecollapse.h
#	vortexpart.h
#	turbulencepart.h
#)
#
## no pre-processing needed
#set(NOPP_SOURCES
#	pwrapper/pymain.cpp
#	pwrapper/pclass.cpp
#	pwrapper/pvec3.cpp
#	pwrapper/pconvert.cpp
#	pwrapper/registry.cpp
#	pwrapper/manta_api.cpp
#	util/vectorbase.cpp
#)
# 
#SET(NOPP_HEADERS
#	pwrapper/pythonInclude.h
#	pwrapper/pclass.h
#	pwrapper/registry.h
#	pwrapper/pconvert.h
#	pwrapper/manta_api.h
#	util/integrator.h
#	util/vectorbase.h
#	util/quaternion.h
#	util/interpol.h
#	util/mcubes.h
#	util/randomstream.h
#	util/solvana.h    
#)
#
## include dirs
#SET(INCLUDE_PATHS
#   pwrapper
#   cuda
#   util
#   ${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source   
#)
#
#list(APPEND INCLUDE_PATHS ${PYTHON_INCLUDE_DIR})
#list(APPEND F_LIBS ${PYTHON_LIBRARIES})
#
#
##******************************************************************************
## apply preprocessor
#
#set(SOURCES ${NOPP_SOURCES} ${SILENT_SOURCES})
#list(APPEND SOURCES ${NOPP_HEADERS})
#set(HEADERS ${NOPP_HEADERS})
#set(PP_REGCPP)
#set(PP_REGS)
#set(PP_PREPD "0")
#MESSAGE("PR_prepdebug")
#
#SET(PP_REGCPP
#commonkernels.h.reg.cpp
#conjugategrad.h.reg.cpp
#edgecollapse.h.reg.cpp
#fastmarch.h.reg.cpp
#fileio.h.reg.cpp
#flip.h.reg.cpp
#fluidsolver.h.reg.cpp
#general.h.reg.cpp
#grid.h.reg.cpp
#kernel.h.reg.cpp
#levelset.h.reg.cpp
#mesh.h.reg.cpp
#movingobs.h.reg.cpp
#noisefield.h.reg.cpp
#particle.h.reg.cpp
#shapes.h.reg.cpp
#turbulencepart.h.reg.cpp
#vortexpart.h.reg.cpp
#vortexsheet.h.reg.cpp
#defines.py.reg.cpp
#)
#
#SET(PP_REGS
#commonkernels.h.reg
#conjugategrad.h.reg
#edgecollapse.h.reg
#fastmarch.h.reg
#fileio.h.reg
#flip.h.reg
#fluidsolver.h.reg
#general.h.reg
#grid.h.reg
#kernel.h.reg
#levelset.h.reg
#mesh.h.reg
#movingobs.h.reg
#noisefield.h.reg
#particle.h.reg
#shapes.h.reg
#turbulencepart.h.reg
#vortexpart.h.reg
#vortexsheet.h.reg
#defines.py.reg
#)
#MESSAGE("PR_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}}")
## link reg files
#add_custom_command(OUTPUT ${PP_REGCPP}
#				  COMMAND prep link ${PP_REGS}
#				  DEPENDS prep ${PP_REGS}
#				  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#				  COMMENT "Linking reg files")
##set_source_files_properties(${PP_REGS} ${PP_REGCPP} PROPERTIES GENERATED 1)
#MESSAGE("PR_prepdebug2")
#
#list(APPEND SOURCES ${PP_SOURCES})
#list(APPEND SOURCES ${PP_HEADERS})
#list(APPEND SOURCES ${PP_REGCPP})
#list(APPEND SOURCES ${PP_REGS})
#
#MESSAGE("PR_prepdebug3")
#
#
##******************************************************************************
## setup executable
#
#SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${F_LINKADD} ")
#SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${F_LINKADD} ")
#MESSAGE("PR_13")
#
#include_directories( ${INCLUDE_PATHS}) 
#link_directories( ${F_LIB_PATHS} ) 
#MESSAGE("PR_14")
#
#
##add_executable(${EXECCMD} ${SOURCES} ${PP_REGISTRY})
##target_link_libraries( ${EXECCMD} ${F_LIBS} )
#
#MESSAGE("PR_15")
#
## debug
#MESSAGE(STATUS "PR_DEBUG Flag-Summary - Includes: '${INCLUDE_PATHS}' | Libs: '${F_LIBS}' | LibPaths: '${F_LIB_PATHS}' ")
#
#=========================================
#~~~original standalone manta setup 
#=========================================

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
python/defines.py.reg

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

MESSAGE ("MANTA_SOURCES ${SOURCES} \n PATHS__:${INCLUDE_PATHS} \n SYSTEM__:${INC_SYS}")
#blender_add_lib(bf_python_manta "${SOURCES}" "${INCLUDE_PATHS}" "${INC_SYS}")
blender_add_lib(bf_python_manta "${SRC}" "${INC}" "${INC_SYS}")