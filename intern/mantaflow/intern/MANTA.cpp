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

/** \file mantaflow/intern/MANTA.cpp
 *  \ingroup mantaflow
 */

#include <sstream>
#include <fstream>
#include <iostream>

#include "MANTA.h"
#include "registry.h"
#include "smoke.h"

#include "BLI_path_util.h"
#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_modifier_types.h"
#include "DNA_smoke_types.h"

bool MANTA::mantaInitialized = false;

MANTA::MANTA(int *res, SmokeModifierData *smd)
{
	std::cout << "MANTA" << std::endl;
	smd->domain->fluid = this;
	smd->domain->manta_solver_res = 3; // Why do we need to set this explicitly? When not set, fluidsolver throws exception (occurs when loading a new .blend file)
	
	// General variables used for low and high res
	std::string tmpScript = "";
	std::string finalScript = "";

	mUsingHeat = smd->domain->active_fields & SM_ACTIVE_HEAT;
	mUsingFire = smd->domain->active_fields & SM_ACTIVE_FIRE;
	mUsingColors = smd->domain->active_fields & SM_ACTIVE_COLORS;
	mUsingHighRes = smd->domain->flags & MOD_SMOKE_HIGHRES;
	
	// Make sure that string vector does not contain any previous commands
	mCommands.clear();

	// Simulation constants
	mTempAmb            = 0; // TODO: Maybe use this later for buoyancy calculation
	mResX               = res[0];
	mResY               = res[1];
	mResZ               = res[2];
	int maxRes          = MAX3(mResX, mResY, mResZ);
	mConstantScaling    = 64.0f / maxRes;
	mConstantScaling    = (mConstantScaling < 1.0f) ? 1.0f : mConstantScaling;
	mTotalCells         = mResX * mResY * mResZ;
	
	// Low res grids
	mDensity        = NULL;
	mHeat           = NULL;
	mVelocityX      = NULL;
	mVelocityY      = NULL;
	mVelocityZ      = NULL;
	mForceX         = NULL;
	mForceY         = NULL;
	mForceZ         = NULL;
	mFlame          = NULL;
	mFuel           = NULL;
	mReact          = NULL;
	mColorR         = NULL;
	mColorG         = NULL;
	mColorB         = NULL;
	mDensityInflow  = NULL;
	mFuelInflow     = NULL;
	mMantaFlags     = NULL;
	mObVelocityX    = new float[mTotalCells];               // TODO in Mantaflow
	mObVelocityY    = new float[mTotalCells];               // TODO in Mantaflow
	mObVelocityZ    = new float[mTotalCells];               // TODO in Mantaflow
	mObstacles      = new unsigned char[mTotalCells];       // TODO in Mantaflow
	mObstaclesAnim  = new unsigned char[mTotalCells];       // TODO in Mantaflow
	
	// High res grids
	mDensityHigh    = NULL;
	mFlameHigh      = NULL;
	mFuelHigh       = NULL;
	mReactHigh      = NULL;
	mColorRHigh     = NULL;
	mColorGHigh     = NULL;
	mColorBHigh     = NULL;
	
	// TODO: Obstacle grid not mantaflow optimized
	for (int x = 0; x < mTotalCells; x++)
		mObstacles[x] = false;

	if (!mantaInitialized)
		startMantaflow();
	
	// Initialize Mantaflow variables in Python
	initSetup(smd);
	if (mUsingHeat)   initHeat(smd);
	if (mUsingFire)   initFire(smd);
	if (mUsingColors) initColors(smd);
	
	updatePointers(smd); // Needs to be after heat, fire, color init
	
	if (mUsingHighRes)
	{		
		// Make sure that string vector does not contain any previous commands
		mCommands.clear();

		// simulation constants
		int amplify     = smd->domain->amplify + 1;
		mResXHigh       = amplify * mResX;
		mResYHigh       = amplify * mResY;
		mResZHigh       = amplify * mResZ;
		mTotalCellsHigh	= mResXHigh * mResYHigh * mResZHigh;
		
		mTextureU = new float[mTotalCells];                 // TODO in Mantaflow
		mTextureV = new float[mTotalCells];                 // TODO in Mantaflow
		mTextureW = new float[mTotalCells];                 // TODO in Mantaflow
		
		const float dx = 1.0f/(float)(mResX);               // TODO in Mantaflow
		const float dy = 1.0f/(float)(mResY);               // TODO in Mantaflow
		const float dz = 1.0f/(float)(mResZ);               // TODO in Mantaflow
		int index = 0;
		for (int z = 0; z < mResZ; z++)
			for (int y = 0; y < mResY; y++)
				for (int x = 0; x < mResX; x++, index++)
				{
					mTextureU[index] = x*dx;                // TODO in Mantaflow
					mTextureV[index] = y*dy;                // TODO in Mantaflow
					mTextureW[index] = z*dz;                // TODO in Mantaflow
				}
		
		// Initialize Mantaflow variables in Python
		initSetupHigh(smd);
		if (mUsingFire)   initFireHigh(smd);
		if (mUsingColors) initColorsHigh(smd);

		updatePointersHigh(smd); // Needs to be after fire, color init
	}
}

void MANTA::initSetup(SmokeModifierData *smd)
{
	std::string tmpString =
		manta_import +
		solver_setup_low +
		alloc_base_grids_low +
		prep_domain_low +
		flags +
		manta_step +
		smoke_step_low;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
	
	runPythonString(mCommands);
}

void MANTA::initSetupHigh(SmokeModifierData *smd)
{
	std::string tmpString =
		solver_setup_high +
		uv_setup +
		alloc_base_grids_high +
		prep_domain_high +
		wavelet_turbulence_noise +
		smoke_step_high;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
		
	runPythonString(mCommands);
	mUsingHighRes = true;
}

void MANTA::initHeat(SmokeModifierData *smd)
{
	if (!mHeat) {
		mCommands.clear();
		mCommands.push_back(alloc_heat_low);
		mCommands.push_back(with_heat);
		
		runPythonString(mCommands);
		mUsingHeat = true;
	}
}

void MANTA::initFire(SmokeModifierData *smd)
{
	if (!mFuel) {
		mCommands.clear();
		mCommands.push_back(alloc_fire_low);
		mCommands.push_back(with_fire);

		runPythonString(mCommands);
		mUsingFire = true;
	}
}

void MANTA::initFireHigh(SmokeModifierData *smd)
{
	if (!mFuelHigh) {
		mCommands.clear();
		mCommands.push_back(alloc_fire_high);
		mCommands.push_back(with_fire);

		runPythonString(mCommands);
		mUsingFire = true;
	}
}

void MANTA::initColors(SmokeModifierData *smd)
{
	if (!mColorR) {
		mCommands.clear();
		std::string colorCodes = parseScript(set_color_codes, smd);
		mCommands.push_back(colorCodes);
		mCommands.push_back(alloc_colors_low);
		mCommands.push_back(init_colors_low);
		mCommands.push_back(with_colors);

		runPythonString(mCommands);
		mUsingColors = true;
	}
}

void MANTA::initColorsHigh(SmokeModifierData *smd)
{
	if (!mColorRHigh) {
		mCommands.clear();
		std::string colorCodes = parseScript(set_color_codes, smd);
		mCommands.push_back(colorCodes);
		mCommands.push_back(alloc_colors_high);
		mCommands.push_back(init_colors_high);
		mCommands.push_back(with_colors);

		runPythonString(mCommands);
		mUsingColors = true;
	}
}

void MANTA::step(SmokeModifierData *smd)
{
	// manta_write_effectors(this);                         // TODO in Mantaflow

	mCommands.clear();
	mCommands.push_back("manta_step()");
	
	runPythonString(mCommands);
}

MANTA::~MANTA()
{
	std::cout << "~MANTA()" << std::endl;

	// Destruction in Python
	mCommands.clear();
	mCommands.push_back(del_base_grids_low);
	mCommands.push_back(del_vars_low);
	if (mUsingHeat)          mCommands.push_back(del_heat_low);
	if (mUsingFire)          mCommands.push_back(del_fire_low);
	if (mUsingColors)        mCommands.push_back(del_colors_low);
	mCommands.push_back(del_solver_low);
	
	if (mUsingHighRes)                 mCommands.push_back(del_base_grids_high);
	if (mUsingHighRes)                 mCommands.push_back(del_vars_high);
	if (mUsingColors && mUsingHighRes) mCommands.push_back(del_colors_high);
	if (mUsingFire && mUsingHighRes)   mCommands.push_back(del_fire_high);
	if (mUsingHighRes)                 mCommands.push_back(del_solver_high);
	
	mCommands.push_back(gc_collect);
	runPythonString(mCommands);
	
	// Reset pointers to avoid dangling pointers
	mDensity        = NULL;
	mHeat           = NULL;
	mVelocityX      = NULL;
	mVelocityY      = NULL;
	mVelocityZ      = NULL;
	mForceX         = NULL;
	mForceY         = NULL;
	mForceZ         = NULL;
	mFlame          = NULL;
	mFuel           = NULL;
	mReact          = NULL;
	mColorR         = NULL;
	mColorG         = NULL;
	mColorB         = NULL;
	mDensityInflow  = NULL;
	mFuelInflow     = NULL;
	mMantaFlags     = NULL;

	if (mObVelocityX)   delete[] mObVelocityX;              // TODO in Mantaflow
	if (mObVelocityY)   delete[] mObVelocityY;              // TODO in Mantaflow
	if (mObVelocityZ)   delete[] mObVelocityZ;              // TODO in Mantaflow
	if (mObstacles)     delete[] mObstacles;                // TODO in Mantaflow
	if (mObstaclesAnim) delete[] mObstaclesAnim;            // TODO in Mantaflow
	
	if (mUsingHighRes)
	{
		mDensityHigh    = NULL;
		mFlameHigh      = NULL;
		mFuelHigh       = NULL;
		mReactHigh      = NULL;
		mColorRHigh     = NULL;
		mColorGHigh     = NULL;
		mColorBHigh     = NULL;
	
		if (mTextureU) delete[] mTextureU;                  // TODO in Mantaflow
		if (mTextureV) delete[] mTextureV;                  // TODO in Mantaflow
		if (mTextureW) delete[] mTextureW;                  // TODO in Mantaflow
	}
	
	// Reset flags
	mUsingHeat = false;
	mUsingFire = false;
	mUsingColors = false;
	mUsingHighRes = false;	
}

void MANTA::runPythonString(std::vector<std::string> commands)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	for (std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it) {
		std::string command = *it;
		PyRun_SimpleString(command.c_str());
	}
	PyGILState_Release(gilstate);
}

void MANTA::startMantaflow()
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

std::string MANTA::getRealValue(const std::string& varName,  SmokeModifierData *smd)
{
	std::ostringstream ss;
	bool is2D = (smd->domain->manta_solver_res == 2);
	ModifierData *md = ((ModifierData*) smd);
	
	if (varName == "USING_COLORS")
		ss << (smd->domain->active_fields & SM_ACTIVE_COLORS ? "True" : "False");
	else if (varName == "USING_HEAT")
		ss << (smd->domain->active_fields & SM_ACTIVE_HEAT ? "True" : "False");
	else if (varName == "USING_FIRE")
		ss << (smd->domain->active_fields & SM_ACTIVE_FIRE ? "True" : "False");
	else if (varName == "USE_WAVELETS")
		ss << (smd->domain->flags & MOD_SMOKE_HIGHRES ? "True" : "False");
	else if (varName == "SOLVER_DIM")
		ss << smd->domain->manta_solver_res;
	else if (varName == "DO_OPEN")
		ss << (smd->domain->border_collisions == SM_BORDER_CLOSED ? "False" : "True");
	else if (varName == "BOUNDCONDITIONS") {
		if (smd->domain->manta_solver_res == 2) {
			if (smd->domain->border_collisions == SM_BORDER_OPEN) ss << "xXyY";
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "yY";
			else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "";
		}
		if (smd->domain->manta_solver_res == 3) {
			if (smd->domain->border_collisions == SM_BORDER_OPEN) ss << "xXyYzZ";
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "zZ";
			else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "";
		}
	} else if (varName == "RES")
		ss << smd->domain->maxres;
	else if (varName == "RESX")
		ss << mResX;
	else if (varName == "RESY")
		if (is2D) {	ss << mResZ;}
		else { 		ss << mResY;}
	else if (varName == "RESZ") {
		if (is2D){	ss << 1;}
		else { 		ss << mResZ;}
	} else if (varName == "DT_FACTOR")
		ss << smd->domain->time_scale;
	else if (varName == "FPS")
		ss << md->scene->r.frs_sec / md->scene->r.frs_sec_base;
	else if (varName == "VORTICITY")
		ss << smd->domain->vorticity / mConstantScaling;
	else if (varName == "UPRES")
		ss << smd->domain->amplify;
	else if (varName == "HRESX")
		ss << mResXHigh;
	else if (varName == "HRESY") {
		if (is2D) {	ss << mResZHigh;}
		else { 		ss << mResYHigh;}
	} else if (varName == "HRESZ") {
		if (is2D) {	ss << 1;}
		else { 		ss << mResZHigh;}
	} else if (varName == "WLT_STR")
		ss << smd->domain->strength;
	else if (varName == "NOISE_POSSCALE")
		ss << smd->domain->noise_pos_scale;
	else if (varName == "NOISE_TIMEANIM")
		ss << smd->domain->noise_time_anim;
	else if (varName == "COLOR_R")
		ss << smd->domain->active_color[0];
	else if (varName == "COLOR_G")
		ss << smd->domain->active_color[1];
	else if (varName == "COLOR_B")
		ss << smd->domain->active_color[2];
	else if (varName == "ADVECT_ORDER")
		ss << 2;
	else if (varName == "ALPHA")
		ss << smd->domain->alpha;
	else if (varName == "BETA")
		ss << smd->domain->beta;
	else if (varName == "BURNING_RATE")
		ss << smd->domain->burning_rate;
	else if (varName == "FLAME_SMOKE")
		ss << smd->domain->flame_smoke;
	else if (varName == "IGNITION_TEMP")
		ss << smd->domain->flame_ignition;
	else if (varName == "MAX_TEMP")
		ss << smd->domain->flame_max_temp;
	else if (varName == "FLAME_SMOKE_COLOR_X")
		ss << smd->domain->flame_smoke_color[0];
	else if (varName == "FLAME_SMOKE_COLOR_Y")
		ss << smd->domain->flame_smoke_color[1];
	else if (varName == "FLAME_SMOKE_COLOR_Z")
		ss << smd->domain->flame_smoke_color[2];
	else if (varName == "MANTA_EXPORT_PATH") {
		char parent_dir[1024];
		BLI_split_dir_part(smd->domain->manta_filepath, parent_dir, sizeof(parent_dir));
		ss << parent_dir;
	} else
		std::cout << "ERROR: Unknown option:" << varName << std::endl;
	return ss.str();
}

std::string MANTA::parseLine(const std::string& line, SmokeModifierData *smd)
{
	if (line.size() == 0) return "";
	std::string res = "";
	int currPos = 0, start_del = 0, end_del = -1;
	bool readingVar = false;
	const char delimiter = '$';
	while (currPos < line.size()){
		if (line[currPos] == delimiter && ! readingVar) {
			readingVar  = true;
			start_del   = currPos + 1;
			res        += line.substr(end_del + 1, currPos - end_del -1);
		}
		else if (line[currPos] == delimiter && readingVar){
			readingVar  = false;
			end_del     = currPos;
			res        += getRealValue(line.substr(start_del, currPos - start_del), smd);
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

std::string MANTA::parseScript(const std::string& setup_string, SmokeModifierData *smd)
{
	std::istringstream f(setup_string);
	std::ostringstream res;
	std::string line = "";
	while(getline(f, line)) {
		res << parseLine(line, smd) << "\n";
	}
	return res.str();
}

void MANTA::exportScript(SmokeModifierData *smd)
{
	// Setup low
	std::string manta_script =
		manta_import +
		solver_setup_low +
		uv_setup +
		alloc_base_grids_low;
	
	// Add heat grid low if needed
	if (smd->domain->active_fields & SM_ACTIVE_HEAT) {
		manta_script += alloc_heat_low;
	}
	
	// Add color grids low if needed
	if (smd->domain->active_fields & SM_ACTIVE_COLORS) {
		manta_script += alloc_colors_low;
	}
	
	// Add fire grids low if needed
	if (smd->domain->active_fields & SM_ACTIVE_FIRE) {
		manta_script += alloc_fire_low;
	}
	
	// Rest of low res setup
	manta_script += prep_domain_low + flags;
	
	// Setup high
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script +=
			solver_setup_high +
			alloc_base_grids_high;
	}
	
	// Add color grids high if needed
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->active_fields & SM_ACTIVE_COLORS) {
		manta_script += alloc_colors_high;
	}
	
	// Add fire grids high if needed
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->active_fields & SM_ACTIVE_FIRE) {
		manta_script += alloc_fire_high;
	}

	// Rest of high res setup
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += prep_domain_high + wavelet_turbulence_noise;
	}
	
	// Import low
	manta_script += smoke_import_low;
	
	// Import high
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += smoke_import_high;
	}
	
	// Inflow low
	manta_script += smoke_inflow_low;
	
	// Inflow High
	// TODO
	
	// Step low functions
	manta_script += smoke_step_low;
	
	// Step high functions
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += smoke_step_high;
	}
	
	// Step wrapper function
	manta_script += manta_step;
	
	// Fill in missing variables in script
	std::string final_script = MANTA::parseScript(manta_script, smd);
	
	// Add standalone mode (loop, gui, ...)
	final_script += standalone;
	
	// Write script
	std::ofstream myfile;
	myfile.open(smd->domain->manta_filepath);
	myfile << final_script;
	myfile.close();
}

void MANTA::exportGrids(SmokeModifierData *smd)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	
	// Export low res grids
	PyRun_SimpleString(MANTA::parseScript(smoke_export_low, smd).c_str());

	// Export high res grids
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		PyRun_SimpleString(MANTA::parseScript(smoke_export_high, smd).c_str());
	}
	PyGILState_Release(gilstate);
}

PyObject* MANTA::getPythonObject(std::string pyVariableName)
{
	if (pyVariableName == "") return NULL;
	
	PyGILState_STATE gilstate = PyGILState_Ensure();

	PyObject* main = PyImport_AddModule("__main__");
	PyObject* pyObject = PyObject_GetAttrString(main, pyVariableName.c_str());

	Py_DECREF(pyObject);

	PyGILState_Release(gilstate);
	return pyObject;
}

std::string MANTA::getGridPointer(std::string gridName, std::string solverName)
{
	if ((gridName == "") && (solverName == "")) return "";

	PyGILState_STATE gilstate = PyGILState_Ensure();

	PyObject* main = PyImport_AddModule("__main__");
	PyObject* gridObject = PyObject_GetAttrString(main, gridName.c_str());

	PyObject* func = PyObject_GetAttrString(gridObject, (char*) "getDataPointer");
	PyObject* returnedValue = PyObject_CallObject(func, NULL);
	PyObject* encoded = PyUnicode_AsUTF8String(returnedValue);

	std::string res = PyBytes_AsString(encoded);
	std::cout << "Pointer on "<< gridName << " " << res << std::endl;

	Py_DECREF(gridObject);
	Py_DECREF(func);
	Py_DECREF(returnedValue);
	Py_DECREF(encoded);

	PyGILState_Release(gilstate);
	return res;
}

void* MANTA::pointerFromString(const std::string& s)
{
	std::stringstream ss(s);
	void *gridPointer = NULL;
	ss >> gridPointer;
	return gridPointer;
}

void MANTA::updatePointers(SmokeModifierData *smd)
{
	std::cout << "Updating pointers low res" << std::endl;
	mDensity        = (float*) pointerFromString( getGridPointer("density",    "s") );
	mVelocityX      = (float*) pointerFromString( getGridPointer("x_vel",      "s") );
	mVelocityY      = (float*) pointerFromString( getGridPointer("y_vel",      "s") );
	mVelocityZ      = (float*) pointerFromString( getGridPointer("z_vel",      "s") );
	mForceX         = (float*) pointerFromString( getGridPointer("x_force",    "s") );
	mForceY         = (float*) pointerFromString( getGridPointer("y_force",    "s") );
	mForceZ         = (float*) pointerFromString( getGridPointer("z_force",    "s") );
	mDensityInflow  = (float*) pointerFromString( getGridPointer("inflow_grid","s") );
	mFuelInflow     = (float*) pointerFromString( getGridPointer("fuel_inflow","s") );
	
	if (mUsingHeat) {
		mHeat       = (float*) pointerFromString(getGridPointer("heat",        "s") );
	}
	if (mUsingFire) {
		mFlame      = (float*) pointerFromString( getGridPointer("flame",      "s") );
		mFuel       = (float*) pointerFromString( getGridPointer("fuel",       "s") );
		mReact      = (float*) pointerFromString( getGridPointer("react",      "s") );
	}
	if (mUsingColors) {
		mColorR     = (float*) pointerFromString( getGridPointer("color_r",    "s") );
		mColorG     = (float*) pointerFromString( getGridPointer("color_g",    "s") );
		mColorB     = (float*) pointerFromString( getGridPointer("color_b",    "s") );
	}
}

void MANTA::updatePointersHigh(SmokeModifierData *smd)
{
	std::cout << "Updating pointers high res" << std::endl;
	mDensityHigh    = (float*) pointerFromString( getGridPointer("xl_density", "xl") );

	if (mUsingFire) {
		mFlameHigh  = (float*) pointerFromString( getGridPointer("xl_flame",   "xl") );
		mFuelHigh   = (float*) pointerFromString( getGridPointer("xl_fuel",    "xl") );
		mReactHigh  = (float*) pointerFromString( getGridPointer("xl_react",   "xl") );
	}
	if (mUsingColors) {
		mColorRHigh = (float*) pointerFromString( getGridPointer("xl_color_r", "xl") );
		mColorGHigh = (float*) pointerFromString( getGridPointer("xl_color_g", "xl") );
		mColorBHigh = (float*) pointerFromString( getGridPointer("xl_color_b", "xl") );
	}
}


