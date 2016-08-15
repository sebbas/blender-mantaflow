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

/** \file mantaflow/intern/SMOKE.cpp
 *  \ingroup mantaflow
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <zlib.h>

#include "SMOKE.h"
#include "registry.h"
#include "shared_script.h"
#include "smoke_script.h"
#include "liquid_script.h"

#include "BLI_path_util.h"
#include "BLI_utildefines.h"
#include "BLI_fileops.h"

#include "DNA_scene_types.h"
#include "DNA_modifier_types.h"
#include "DNA_smoke_types.h"

bool SMOKE::mantaInitialized = false;

SMOKE::SMOKE(int *res, SmokeModifierData *smd)
{
	std::cout << "SMOKE" << std::endl;
	smd->domain->fluid = this;
	smd->domain->manta_solver_res = 3; // Why do we need to set this explicitly? When not set, fluidsolver throws exception (occurs when loading a new .blend file)
	
	mUsingHeat    = smd->domain->active_fields & SM_ACTIVE_HEAT;
	mUsingFire    = smd->domain->active_fields & SM_ACTIVE_FIRE;
	mUsingColors  = smd->domain->active_fields & SM_ACTIVE_COLORS;
	mUsingHighRes = smd->domain->flags & MOD_SMOKE_HIGHRES;
	mUsingLiquid  = smd->domain->type == MOD_SMOKE_DOMAIN_TYPE_LIQUID;
	mUsingSmoke   = smd->domain->type == MOD_SMOKE_DOMAIN_TYPE_GAS;
	
	// Make sure that string vector does not contain any previous commands
	mCommands.clear();

	// Simulation constants
	mTempAmb            = 0; // TODO: Maybe use this later for buoyancy calculation
	mResX               = res[0];
	mResY               = res[1];
	mResZ               = res[2];
	mMaxRes             = MAX3(mResX, mResY, mResZ);
	mConstantScaling    = 64.0f / mMaxRes;
	mConstantScaling    = (mConstantScaling < 1.0f) ? 1.0f : mConstantScaling;
	mTotalCells         = mResX * mResY * mResZ;
	
	// Smoke low res grids
	mDensity        = NULL;
	mFlags          = NULL;
	mHeat           = NULL;
	mVelocityX      = NULL;
	mVelocityY      = NULL;
	mVelocityZ      = NULL;
	mObVelocityX    = NULL;
	mObVelocityY    = NULL;
	mObVelocityZ    = NULL;
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
	mObstacles      = NULL;
	
	// Smoke high res grids
	mDensityHigh    = NULL;
	mFlameHigh      = NULL;
	mFuelHigh       = NULL;
	mReactHigh      = NULL;
	mColorRHigh     = NULL;
	mColorGHigh     = NULL;
	mColorBHigh     = NULL;
	mTextureU       = NULL;
	mTextureV       = NULL;
	mTextureW       = NULL;
	mTextureU2      = NULL;
	mTextureV2      = NULL;
	mTextureW2      = NULL;
	
	// Liquid low res grids
	mPhi            = NULL;
	mPhiInit        = NULL;
	
	// Liquid high res grids
	mPhiHigh        = NULL;
	
	mNumVertices  = 0;
	mNumNormals   = 0;
	mNumTriangles = 0;

	// Only start Mantaflow once. No need to start whenever new SMOKE objected is allocated
	if (!mantaInitialized)
		startMantaflow();
	
	// Initialize Mantaflow variables in Python
	// Liquid
	if (mUsingLiquid) {
		initDomain(smd);
		initLiquid(smd);

		updatePointers(smd);
		
		if (mUsingHighRes) {
			// Make sure that string vector does not contain any previous commands
			mCommands.clear();

			// simulation constants
			int amplify     = smd->domain->amplify + 1;
			mResXHigh       = amplify * mResX;
			mResYHigh       = amplify * mResY;
			mResZHigh       = amplify * mResZ;
			mTotalCellsHigh	= mResXHigh * mResYHigh * mResZHigh;
			
			// Initialize Mantaflow variables in Python
			initDomainHigh(smd);
			initLiquidHigh(smd);

			updatePointersHigh(smd);
		}

		return;
	}
	
	// Smoke
	if (mUsingSmoke) {
		initDomain(smd);
		initSmoke(smd);
		if (mUsingHeat)   initHeat(smd);
		if (mUsingFire)   initFire(smd);
		if (mUsingColors) initColors(smd);

		updatePointers(smd); // Needs to be after heat, fire, color init

		if (mUsingHighRes) {
			// Make sure that string vector does not contain any previous commands
			mCommands.clear();

			// simulation constants
			int amplify     = smd->domain->amplify + 1;
			mResXHigh       = amplify * mResX;
			mResYHigh       = amplify * mResY;
			mResZHigh       = amplify * mResZ;
			mTotalCellsHigh	= mResXHigh * mResYHigh * mResZHigh;
			
			// Initialize Mantaflow variables in Python
			initDomainHigh(smd);
			initSmokeHigh(smd);
			if (mUsingFire)   initFireHigh(smd);
			if (mUsingColors) initColorsHigh(smd);

			updatePointersHigh(smd); // Needs to be after fire, color init
		}
	}
}

void SMOKE::initDomain(SmokeModifierData *smd)
{
	std::string tmpString = manta_import
		+ fluid_solver_low
		+ fluid_adaptive_time_stepping_low;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
	
	runPythonString(mCommands);
}

void SMOKE::initDomainHigh(SmokeModifierData *smd)
{
	std::string tmpString = fluid_solver_high
		+ fluid_adaptive_time_stepping_high;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
	
	runPythonString(mCommands);
}

void SMOKE::initSmoke(SmokeModifierData *smd)
{
	std::string tmpString = smoke_alloc_low
		+ smoke_variables_low
		+ smoke_bounds_low
		+ smoke_adaptive_step
		+ smoke_step_low;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
	
	runPythonString(mCommands);
}

void SMOKE::initSmokeHigh(SmokeModifierData *smd)
{
	std::string tmpString = smoke_alloc_high
		+ smoke_variables_high
		+ smoke_uv_setup
		+ smoke_bounds_high
		+ smoke_wavelet_turbulence_noise
		+ smoke_step_high;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
		
	runPythonString(mCommands);
	mUsingHighRes = true;
}

void SMOKE::initHeat(SmokeModifierData *smd)
{
	if (!mHeat) {
		mCommands.clear();
		mCommands.push_back(smoke_alloc_heat_low);
		mCommands.push_back(smoke_with_heat);
		
		runPythonString(mCommands);
		mUsingHeat = true;
	}
}

void SMOKE::initFire(SmokeModifierData *smd)
{
	if (!mFuel) {
		mCommands.clear();
		mCommands.push_back(smoke_alloc_fire_low);
		mCommands.push_back(smoke_with_fire);

		runPythonString(mCommands);
		mUsingFire = true;
	}
}

void SMOKE::initFireHigh(SmokeModifierData *smd)
{
	if (!mFuelHigh) {
		mCommands.clear();
		mCommands.push_back(smoke_alloc_fire_high);
		mCommands.push_back(smoke_with_fire);

		runPythonString(mCommands);
		mUsingFire = true;
	}
}

void SMOKE::initColors(SmokeModifierData *smd)
{
	if (!mColorR) {
		mCommands.clear();
		std::string colorCodes = parseScript(smoke_set_color_codes, smd);
		mCommands.push_back(colorCodes);
		mCommands.push_back(smoke_alloc_colors_low);
		mCommands.push_back(smoke_init_colors_low);
		mCommands.push_back(smoke_with_colors);

		runPythonString(mCommands);
		mUsingColors = true;
	}
}

void SMOKE::initColorsHigh(SmokeModifierData *smd)
{
	if (!mColorRHigh) {
		mCommands.clear();
		std::string colorCodes = parseScript(smoke_set_color_codes, smd);
		mCommands.push_back(colorCodes);
		mCommands.push_back(smoke_alloc_colors_high);
		mCommands.push_back(smoke_init_colors_high);
		mCommands.push_back(smoke_with_colors);

		runPythonString(mCommands);
		mUsingColors = true;
	}
}

void SMOKE::initLiquid(SmokeModifierData *smd)
{
	if (!mPhi) {
		std::string tmpString = liquid_alloc_low
			+ liquid_variables_low
			+ liquid_bounds_low
			+ liquid_init_phi
			+ liquid_save_mesh
			+ liquid_export_low
			+ liquid_import_low
			+ liquid_adaptive_step
			+ liquid_step_low;
		std::string finalString = parseScript(tmpString, smd);
		mCommands.clear();
		mCommands.push_back(finalString);

		runPythonString(mCommands);
		mUsingLiquid = true;
	}
}

void SMOKE::initLiquidHigh(SmokeModifierData *smd)
{
	std::string tmpString = liquid_alloc_high
		+ liquid_variables_high
		+ liquid_bounds_high
		+ liquid_step_high;
	std::string finalString = parseScript(tmpString, smd);
	mCommands.clear();
	mCommands.push_back(finalString);
		
	runPythonString(mCommands);
	mUsingHighRes = true;
}

void SMOKE::step(SmokeModifierData *smd)
{
	// manta_write_effectors(this);                         // TODO in Mantaflow

	// Get the frame number for this step
	ModifierData *md = ((ModifierData*) smd);
	int startFrame = md->scene->r.cfra - 1; // Current frame is always one ahead
	
	// Run manta step and handover current frame number
	mCommands.clear();
	std::ostringstream manta_step;
	manta_step <<  "manta_step(" << startFrame << ")";
	mCommands.push_back(manta_step.str());
	
	runPythonString(mCommands);
}

SMOKE::~SMOKE()
{
	std::cout << "~SMOKE()" << std::endl;

	// Destruction in Python
	mCommands.clear();
	
	// Liquid
	if (mUsingLiquid) {
		mCommands.push_back(liquid_delete_grids_low);
		mCommands.push_back(liquid_delete_variables_low);
		
		if (mUsingHighRes) mCommands.push_back(liquid_delete_grids_high);
		if (mUsingHighRes) mCommands.push_back(liquid_delete_variables_high);
	}
	
	// Smoke
	if (mUsingSmoke) {
		mCommands.push_back(smoke_delete_grids_low);
		mCommands.push_back(smoke_delete_variables_low);
		if (mUsingHeat)          mCommands.push_back(smoke_delete_heat_low);
		if (mUsingFire)          mCommands.push_back(smoke_delete_fire_low);
		if (mUsingColors)        mCommands.push_back(smoke_delete_colors_low);
		
		if (mUsingHighRes)                 mCommands.push_back(smoke_delete_grids_high);
		if (mUsingHighRes)                 mCommands.push_back(smoke_delete_variables_high);
		if (mUsingColors && mUsingHighRes) mCommands.push_back(smoke_delete_colors_high);
		if (mUsingFire && mUsingHighRes)   mCommands.push_back(smoke_delete_fire_high);
	}
	
	// Make sure that everything is garbage collected
	mCommands.push_back(gc_collect);

	// Solvers always have to be the last objects to be deleted
	mCommands.push_back(fluid_delete_solver_low);
	if (mUsingHighRes) mCommands.push_back(fluid_delete_solver_high);
	
	// Just in case: gc again
	mCommands.push_back(gc_collect);
	runPythonString(mCommands);
	
	// Reset pointers to avoid dangling pointers
	mDensity        = NULL;
	mFlags          = NULL;
	mHeat           = NULL;
	mVelocityX      = NULL;
	mVelocityY      = NULL;
	mVelocityZ      = NULL;
	mObVelocityX    = NULL;
	mObVelocityY    = NULL;
	mObVelocityZ    = NULL;
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
	mObstacles      = NULL;
	
	if (mUsingHighRes)
	{
		mDensityHigh    = NULL;
		mFlameHigh      = NULL;
		mFuelHigh       = NULL;
		mReactHigh      = NULL;
		mColorRHigh     = NULL;
		mColorGHigh     = NULL;
		mColorBHigh     = NULL;
		mTextureU       = NULL;
		mTextureV       = NULL;
		mTextureW       = NULL;
		mTextureU2      = NULL;
		mTextureV2      = NULL;
		mTextureW2      = NULL;
	}
	
	// Liquid
	mPhi     = NULL;
	mPhiInit = NULL;
	mPhiHigh = NULL;
	
	// Reset flags
	mUsingHeat    = false;
	mUsingFire    = false;
	mUsingColors  = false;
	mUsingHighRes = false;	
}

void SMOKE::runPythonString(std::vector<std::string> commands)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	for (std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it) {
		std::string command = *it;

#ifdef WIN32
		// special treatment for windows when running python code
		size_t cmdLength = command.length();
		char* buffer = new char[cmdLength+1];
		memcpy(buffer, command.data(), cmdLength);

		buffer[cmdLength] = '\0';
		PyRun_SimpleString(buffer);
		delete[] buffer;
#else
		PyRun_SimpleString(command.c_str());
#endif
	}
	PyGILState_Release(gilstate);
}

void SMOKE::startMantaflow()
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

std::string SMOKE::getRealValue(const std::string& varName,  SmokeModifierData *smd)
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
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "xX";
			else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "";
		}
		if (smd->domain->manta_solver_res == 3) {
			if (smd->domain->border_collisions == SM_BORDER_OPEN) ss << "xXyYzZ";
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "zZ";
			else if (smd->domain->border_collisions == SM_BORDER_HORIZONTAL) ss << "xXyY";
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
		ss << smd->domain->amplify + 1;
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
	else if (varName == "CURRENT_FRAME")
		ss << md->scene->r.cfra;
	else if (varName == "RANDOMNESS")
		ss << smd->domain->particle_randomness;
	else if (varName == "GRAVITY_X")
		ss << smd->domain->gravity[0];
	else if (varName == "GRAVITY_Y")
		ss << smd->domain->gravity[1];
	else if (varName == "GRAVITY_Z")
		ss << smd->domain->gravity[2];
	else if (varName == "MANTA_EXPORT_PATH") {
		char parent_dir[1024];
		BLI_split_dir_part(smd->domain->manta_filepath, parent_dir, sizeof(parent_dir));
		ss << parent_dir;
	} else
		std::cout << "ERROR: Unknown option:" << varName << std::endl;
	return ss.str();
}

std::string SMOKE::parseLine(const std::string& line, SmokeModifierData *smd)
{
	if (line.size() == 0) return "";
	std::string res = "";
	int currPos = 0, start_del = 0, end_del = -1;
	bool readingVar = false;
	const char delimiter = '$';
	while (currPos < line.size()) {
		if (line[currPos] == delimiter && !readingVar) {
			readingVar  = true;
			start_del   = currPos + 1;
			res        += line.substr(end_del + 1, currPos - end_del -1);
		}
		else if (line[currPos] == delimiter && readingVar) {
			readingVar  = false;
			end_del     = currPos;
			res        += getRealValue(line.substr(start_del, currPos - start_del), smd);
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

std::string SMOKE::parseScript(const std::string& setup_string, SmokeModifierData *smd)
{
	std::istringstream f(setup_string);
	std::ostringstream res;
	std::string line = "";
	while(getline(f, line)) {
		res << parseLine(line, smd) << "\n";
	}
	return res.str();
}

void SMOKE::exportScript(SmokeModifierData *smd)
{
	// Setup low
	std::string manta_script =
		manta_import +
		fluid_solver_low +
		smoke_alloc_low;
	
	// Add heat grid low if needed
	if (smd->domain->active_fields & SM_ACTIVE_HEAT) {
		manta_script += smoke_alloc_heat_low;
	}
	
	// Add color grids low if needed
	if (smd->domain->active_fields & SM_ACTIVE_COLORS) {
		manta_script += smoke_alloc_colors_low;
	}
	
	// Add fire grids low if needed
	if (smd->domain->active_fields & SM_ACTIVE_FIRE) {
		manta_script += smoke_alloc_fire_low;
	}
	
	// Rest of low res setup
	manta_script += smoke_bounds_low + smoke_variables_low;
	
	// Setup high
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += fluid_solver_high
			+ smoke_variables_high
			+ smoke_uv_setup
			+ smoke_alloc_high;
	}
	
	// Add color grids high if needed
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->active_fields & SM_ACTIVE_COLORS) {
		manta_script += smoke_alloc_colors_high;
	}
	
	// Add fire grids high if needed
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->active_fields & SM_ACTIVE_FIRE) {
		manta_script += smoke_alloc_fire_high;
	}

	// Rest of high res setup
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += smoke_bounds_high + smoke_wavelet_turbulence_noise;
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
	manta_script += smoke_adaptive_step;
	
	// Fill in missing variables in script
	std::string final_script = SMOKE::parseScript(manta_script, smd);
	
	// Add standalone mode (loop, gui, ...)
	final_script += smoke_standalone;
	
	// Write script
	std::ofstream myfile;
	myfile.open(smd->domain->manta_filepath);
	myfile << final_script;
	myfile.close();
}

void SMOKE::exportGrids(SmokeModifierData *smd)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	
	// Export low res grids
	PyRun_SimpleString(SMOKE::parseScript(smoke_export_low, smd).c_str());

	// Export high res grids
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		PyRun_SimpleString(SMOKE::parseScript(smoke_export_high, smd).c_str());
	}
	PyGILState_Release(gilstate);
}

void* SMOKE::getGridPointer(std::string gridName, std::string solverName)
{
	if ((gridName == "") && (solverName == "")) return NULL;

	PyGILState_STATE gilstate = PyGILState_Ensure();

	// Get pyobject that holds pointer address as string
	PyObject* main = PyImport_AddModule("__main__");
	PyObject* gridObject = PyObject_GetAttrString(main, gridName.c_str());
	PyObject* func = PyObject_GetAttrString(gridObject, (char*) "getDataPointer");
	PyObject* returnedValue = PyObject_CallObject(func, NULL);
	PyObject* encoded = PyUnicode_AsUTF8String(returnedValue);

	// Convert string pointer to void pointer
	std::string pointerString = PyBytes_AsString(encoded);
	std::istringstream in(pointerString);
	void *gridPointer = NULL;
	in >> gridPointer;
	
	Py_DECREF(gridObject);
	Py_DECREF(func);
	Py_DECREF(returnedValue);
	Py_DECREF(encoded);

	PyGILState_Release(gilstate);
	return gridPointer;
}

void SMOKE::updateMeshData(const char* filename)
{
	gzFile gzf;
	float fbuffer[3];
	int ibuffer[3];

	gzf = (gzFile) BLI_gzopen(filename, "rb1"); // do some compression
	if (!gzf)
		std::cout << "readBobj: unable to open file" << std::endl;
	
	// Num vertices
	mNumVertices = 0;
	gzread(gzf, &mNumVertices, sizeof(int));
	
	std::cout << "read mesh , verts " << mNumVertices << std::endl;

	if (mNumVertices)
	{
		mVerticesX.resize(mNumVertices);
		mVerticesY.resize(mNumVertices);
		mVerticesZ.resize(mNumVertices);
		
		// Vertices
		for (int i = 0; i < mNumVertices; i++) {
			gzread(gzf, fbuffer, sizeof(float) * 3);
			
			mVerticesX[i] = fbuffer[0];
			mVerticesY[i] = fbuffer[1];
			mVerticesZ[i] = fbuffer[2];
		}
	}
	
	// Num normals
	mNumNormals = 0;
	gzread(gzf, &mNumNormals, sizeof(float));
	
	if (mNumNormals)
	{
		mNormalsX.resize(mNumNormals);
		mNormalsY.resize(mNumNormals);
		mNormalsZ.resize(mNumNormals);
		
		// Normals
		for (int i = 0; i < mNumNormals; i++) {
			gzread(gzf, fbuffer, sizeof(float) * 3);
			
			mNormalsX[i] = fbuffer[0];
			mNormalsY[i] = fbuffer[1];
			mNormalsZ[i] = fbuffer[2];
		}
	}
	
	// Num triangles
	mNumTriangles = 0;
	gzread(gzf, &mNumTriangles, sizeof(int));
	
	if (mNumTriangles)
	{
		mTrianglesX.resize(mNumTriangles);
		mTrianglesY.resize(mNumTriangles);
		mTrianglesZ.resize(mNumTriangles);
		
		// Triangles
		for (int i = 0; i < mNumTriangles; i++) {
			gzread(gzf, ibuffer, sizeof(int) * 3);
			
			mTrianglesX[i] = ibuffer[0];
			mTrianglesY[i] = ibuffer[1];
			mTrianglesZ[i] = ibuffer[2];
		}
	}

	gzclose( gzf );
}

void SMOKE::updatePointers(SmokeModifierData *smd)
{
	std::cout << "Updating pointers low res" << std::endl;

	mFlags = (int*) getGridPointer("flags", "s");
	
	// Liquid
	if (mUsingLiquid) {
		mPhi        = (float*)         getGridPointer("phi",             "s");
		mPhiInit    = (float*)         getGridPointer("phiInit",         "s");
		mDensity    = (float*)         getGridPointer("density",         "s");
	}
	
	// Smoke
	if (mUsingSmoke) {
		mDensity        = (float*)         getGridPointer("density",     "s");
		mVelocityX      = (float*)         getGridPointer("x_vel",       "s");
		mVelocityY      = (float*)         getGridPointer("y_vel",       "s");
		mVelocityZ      = (float*)         getGridPointer("z_vel",       "s");
		mObVelocityX    = (float*)         getGridPointer("x_obvel",     "s");
		mObVelocityY    = (float*)         getGridPointer("y_obvel",     "s");
		mObVelocityZ    = (float*)         getGridPointer("z_obvel",     "s");
		mForceX         = (float*)         getGridPointer("x_force",     "s");
		mForceY         = (float*)         getGridPointer("y_force",     "s");
		mForceZ         = (float*)         getGridPointer("z_force",     "s");
		mDensityInflow  = (float*)         getGridPointer("inflow_grid", "s");
		mFuelInflow     = (float*)         getGridPointer("fuel_inflow", "s");
		mObstacles      = (unsigned char*) getGridPointer("flags",       "s");
		
		if (mUsingHeat) {
			mHeat       = (float*) getGridPointer("heat",    "s");
		}
		if (mUsingFire) {
			mFlame      = (float*) getGridPointer("flame",   "s");
			mFuel       = (float*) getGridPointer("fuel",    "s");
			mReact      = (float*) getGridPointer("react",   "s");
		}
		if (mUsingColors) {
			mColorR     = (float*) getGridPointer("color_r", "s");
			mColorG     = (float*) getGridPointer("color_g", "s");
			mColorB     = (float*) getGridPointer("color_b", "s");
		}
	}
}

void SMOKE::updatePointersHigh(SmokeModifierData *smd)
{
	std::cout << "Updating pointers high res" << std::endl;

	// Liquid
	if (mUsingLiquid) {
		// TODO (sebbas) phiInitHigh does not exist yet
		// mPhiHigh    = (float*) getGridPointer("phiInitHigh", "xl");
	}
	
	if (mUsingSmoke) {
		mDensityHigh    = (float*) getGridPointer("xl_density", "xl");
		mTextureU       = (float*) getGridPointer("texture_u",  "s");
		mTextureV       = (float*) getGridPointer("texture_v",  "s");
		mTextureW       = (float*) getGridPointer("texture_w",  "s");
		mTextureU2      = (float*) getGridPointer("texture_u2", "s");
		mTextureV2      = (float*) getGridPointer("texture_v2", "s");
		mTextureW2      = (float*) getGridPointer("texture_w2", "s");
		
		if (mUsingFire) {
			mFlameHigh  = (float*) getGridPointer("xl_flame",   "xl");
			mFuelHigh   = (float*) getGridPointer("xl_fuel",    "xl");
			mReactHigh  = (float*) getGridPointer("xl_react",   "xl");
		}
		if (mUsingColors) {
			mColorRHigh = (float*) getGridPointer("xl_color_r", "xl");
			mColorGHigh = (float*) getGridPointer("xl_color_g", "xl");
			mColorBHigh = (float*) getGridPointer("xl_color_b", "xl");
		}
	}
}

void SMOKE::saveMesh(char *filename)
{
	std::string path(filename);
	
	mCommands.clear();
	std::ostringstream save_mesh;
	save_mesh <<  "save_mesh('" << path << "')";
	mCommands.push_back(save_mesh.str());
	
	runPythonString(mCommands);
}

void SMOKE::saveLiquidData(char *pathname)
{
	std::string path(pathname);
	
	mCommands.clear();
	std::ostringstream save_liquid_data;
	save_liquid_data <<  "save_liquid_data('" << path << "')";
	mCommands.push_back(save_liquid_data.str());
	
	runPythonString(mCommands);
}

void SMOKE::loadLiquidData(char *pathname)
{
	std::string path(pathname);
	
	mCommands.clear();
	std::ostringstream load_liquid_data;
	load_liquid_data <<  "load_liquid_data('" <<  path << "')";
	mCommands.push_back(load_liquid_data.str());
	
	runPythonString(mCommands);
}


