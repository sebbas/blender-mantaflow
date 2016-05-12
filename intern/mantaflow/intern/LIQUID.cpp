/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2016 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sebastian Barschkis (sebbas)
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file mantaflow/intern/LIQUID.cpp
 *  \ingroup mantaflow
 */

#include <sstream>
#include <fstream>
#include <iostream>

#include "LIQUID.h"
#include "registry.h"
#include "liquid_script.h"
#include "shared_script.h"

bool LIQUID::mantaInitialized = false;

LIQUID::LIQUID()
{
	std::cout << "LIQUID" << std::endl;
	
	// Only start Mantaflow once. No need to start whenever new MANTA objected is allocated
	if (!mantaInitialized)
		startMantaflow();

	initSetup();
}

void LIQUID::initSetup()
{
	std::string tmpString = manta_import
		+ solver_low
		+ adaptive_time_stepping
		+ alloc_liquid
		+ liquid_variables
		+ prep_domain
		+ mesh_loading
		+ manta_step
		+ liquid_step;
//	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(tmpString);
	
	runPythonString(mCommands);
}

void LIQUID::step(int currentFrame)
{
	// Run manta step and handover current frame number
	mCommands.clear();
	std::ostringstream manta_step;
	manta_step <<  "manta_step(" << currentFrame << ")";
	mCommands.push_back(manta_step.str());
	
	runPythonString(mCommands);
}

LIQUID::~LIQUID()
{
	std::cout << "~LIQUID()" << std::endl;
}

void LIQUID::runPythonString(std::vector<std::string> commands)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	for (std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it) {
		std::string command = *it;
		PyRun_SimpleString(command.c_str());
	}
	PyGILState_Release(gilstate);
}

void LIQUID::startMantaflow()
{
	std::cout << "Starting mantaflow" << std::endl;
	std::string filename = "manta_scene.py";
	std::vector<std::string> fill = std::vector<std::string>();
	
	// Initialize extension classes and wrappers
	srand(0);
	PyGILState_STATE gilstate = PyGILState_Ensure();
	Pb::setup(filename, fill);  // Namespace from Mantaflow (registry)
	PyGILState_Release(gilstate);
	mantaInitialized = true;
}

