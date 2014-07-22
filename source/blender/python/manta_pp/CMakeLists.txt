#add_definitions(-DGUI=1)

set(INC
	.
	python
	plugin
	pwrapper
	util
	#gui
)

set(INC_SYS
	${PYTHON_INCLUDE_DIRS}
)
set(SRC
#PP SOURCES and HEADERS
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
plugin/advection.cpp
plugin/extforces.cpp
plugin/initplugins.cpp
plugin/kepsilon.cpp
plugin/meshplugins.cpp
plugin/pressure.cpp
plugin/vortexplugins.cpp
plugin/waveletturbulence.cpp
python/defines.py
python/defines.py.reg
python/defines.py.reg.cpp
shapes.cpp
shapes.h
shapes.h.reg
shapes.h.reg.cpp
test.cpp
turbulencepart.cpp
turbulencepart.h
turbulencepart.h.reg
turbulencepart.h.reg.cpp
vortexpart.cpp
vortexpart.h
vortexpart.h.reg
vortexpart.h.reg.cpp
vortexsheet.cpp
vortexsheet.h
vortexsheet.h.reg
vortexsheet.h.reg.cpp
#NOPP SOURCES
pwrapper/pymain.cpp
pwrapper/pclass.cpp
pwrapper/pvec3.cpp
pwrapper/pconvert.cpp
pwrapper/registry.cpp
pwrapper/pythonInclude.h
pwrapper/pclass.h
pwrapper/registry.h
pwrapper/pconvert.h
pwrapper/manta_api.cpp
pwrapper/manta_api.h
util/vectorbase.cpp
util/integrator.h
util/vectorbase.h
util/quaternion.h
util/interpol.h
util/mcubes.h
util/randomstream.h
util/solvana.h

#gui/customctrl.cpp
#gui/customctrl.h
#gui/customctrl.h.reg
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
)

#linking QT lbraries

cmake_policy(SET CMP0020 NEW)
find_package(Qt5Core QUIET)
if (Qt5Core_FOUND)
	message("Using Qt5")
	find_package(Qt5Widgets REQUIRED)
	find_package(Qt5OpenGL REQUIRED)

	qt5_wrap_cpp(MOC_OUTFILES ${QT_REMAP} )
	qt5_add_resources(QT_RES resources/res.qrc )
	add_definitions(-DGUI=1)
	add_definitions(${Qt5Widgets_DEFINITIONS})
	list(APPEND INCLUDE_PATHS ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS})
	list(APPEND INC_SYS ${Qt5Widgets_LIBRARIES} ${Qt5OpenGL_LIBRARIES})
	list(APPEND SOURCES ${MOC_OUTFILES} ${QT_RES})		
else()
	message("Using Qt4")
	find_package(Qt4 REQUIRED)
	set(QT_USE_QTOPENGL TRUE)

#	qt4_wrap_cpp(MOC_OUTFILES ${QT_REMAP} )
#	qt4_add_resources(QT_RES resources/res.qrc )
	include(${QT_USE_FILE})
	add_definitions(${QT_DEFINITIONS})
	list(APPEND INC_SYS ${QT_LIBRARIES})
	list(APPEND SOURCES ${MOC_OUTFILES} ${QT_RES})
	endif()

	if (APPLE)
	# mac opengl framework
#	SET(INC_SYS "${INC_SYS} -framework OpenGL ")
	else()
	find_package(OpenGL REQUIRED)
	list(APPEND INC_SYS ${OPENGL_LIBRARIES})    
	endif()



blender_add_lib(bf_python_manta "${SRC}" "${INC}" "${INC_SYS}")