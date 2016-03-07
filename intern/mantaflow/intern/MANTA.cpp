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

#include "MANTA.h"

MANTA::MANTA(int *res, float dx, float dtdef, int init_heat, int init_fire, int init_colors, SmokeModifierData *smd) :
_xRes(res[0]), _yRes(res[1]), _zRes(res[2]), _res(0.0f)
{
	std::cout << "MANTA" << std::endl;
	// set simulation consts
	_dt = dtdef;	// just in case. set in step from a RNA factor

	_iterations = 100;
	_tempAmb = 0; 
	_heatDiffusion = 1e-3;
	_totalTime = 0.0f;
	_totalSteps = 0;
	_maxRes = MAX3(_xRes, _yRes, _zRes);
	
	if (!dx)
		_dx = 1.0f / (float)_maxRes;
	else
		_dx = dx;
	_constantScaling = 64.0f / _maxRes;
	_constantScaling = (_constantScaling < 1.0f) ? 1.0f : _constantScaling;
	_vorticityEps = 2.0f / _constantScaling; // Just in case set a default value

	_totalCells   = _xRes * _yRes * _zRes;

	// allocate arrays
	_manta_inflow = NULL;
	_fuel_inflow = NULL;
	_manta_flags = NULL;
	_xVelocity = NULL;
	_yVelocity = NULL;
	_zVelocity = NULL;
	_xVelocityOb = NULL;
	_yVelocityOb = NULL;
	_zVelocityOb = NULL;
	_xForce = NULL;
	_yForce = NULL;
	_zForce = NULL;
	_density = NULL;
	_obstacles = new unsigned char[_totalCells]; // set 0 at end of step

	_using_heat = false;
	_using_fire = false;
	_using_colors = false;
	
	MANTA::start_mantaflow();
	
	// Base setup low res
	std::string setup_script =
		manta_import +
		solver_setup_low +
		uv_setup +
		alloc_base_grids_low +
		prep_domain_low +
		flags +
		manta_step +
		smoke_step_low;
	std::string final_script = MANTA::parse_script(setup_script, smd);
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(final_script.c_str());
	PyGILState_Release(gilstate);
	
	smd->domain->manta = this;
	
	/* heat */
	_heat = NULL;
	if (init_heat) {
		initHeat();
	}
	// Fire simulation
	_flame = NULL;
	_fuel = NULL;
	_react = NULL;
	if (init_fire) {
		initFire();
	}
	// Smoke color
	_color_r = NULL;
	_color_g = NULL;
	_color_b = NULL;
	if (init_colors) {
		initColors(0.0f, 0.0f, 0.0f);
	}
	
	MANTA::update_pointers();
	
	// High resolution
	_xResBig = _amplify * _xRes;
	_yResBig = _amplify * _yRes;
	_zResBig = _amplify * _zRes;
	_slabSizeBig = _xResBig * _yResBig;
	_totalCellsBig = _slabSizeBig * _zResBig;
	
	// allocate high resolution density field
	_densityBig = NULL;
	_flameBig = NULL;
	_fuelBig = NULL;
	_reactBig = NULL;
	
	_color_rBig = NULL;
	_color_gBig = NULL;
	_color_bBig = NULL;
	
	// Base setup high res
	std::string setup_script =
		solver_setup_high +
		alloc_base_grids_high +
		prep_domain_high +
		flags +
		wavelet_turbulence_noise +
		smoke_step_high;
	std::string final_script = Manta_API::parse_script(setup_script);
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(final_script.c_str());
	PyGILState_Release(gilstate);
	
	// Fire grids
	if (init_fire) {
		initFireHigh();
	}
	
	// Color grids
	if (init_colors) {
		initColorsHigh(0.0f, 0.0f, 0.0f);
	}

	// allocate & init texture coordinates
	_tcU = NULL;
	_tcV = NULL;
	_tcW = NULL;
	
	MANTA::update_high_res_pointers();
}

void MANTA::initHeat()
{
	if (!_heat) {
		_using_heat = true;
		PyGILState_STATE gilstate = PyGILState_Ensure();
		PyRun_SimpleString(alloc_heat_low.c_str());
		PyRun_SimpleString(with_heat.c_str());
		PyGILState_Release(gilstate);
		MANTA::update_pointers();
	}
}

void MANTA::initFire()
{
	if (!_flame) {
		_using_fire = true;
		PyGILState_STATE gilstate = PyGILState_Ensure();
		PyRun_SimpleString(alloc_fire_low.c_str());
		PyRun_SimpleString(with_fire.c_str());
		PyGILState_Release(gilstate);
		MANTA::update_pointers();
	}
}

void MANTA::initFireHigh()
{
	if (!_fuelBig) {
		using_fire = true;
		PyGILState_STATE gilstate = PyGILState_Ensure();
		PyRun_SimpleString(alloc_fire_high.c_str());
		PyGILState_Release(gilstate);
		MANTA::update_high_res_pointers();
	}
}

void MANTA::initColors(float init_r, float init_g, float init_b)
{
	if (!_color_r) {
		_using_colors = true;
		PyGILState_STATE gilstate = PyGILState_Ensure();
		std::stringstream ss;
		ss << "manta_color_r = " << init_r << endl;
		ss << "manta_color_g = " << init_g << endl;
		ss << "manta_color_b = " << init_b << endl;
		PyRun_SimpleString(ss.str().c_str());
		PyRun_SimpleString(alloc_colors_low.c_str());
		PyRun_SimpleString(init_colors_low.c_str());
		PyRun_SimpleString(with_colors.c_str());
		PyGILState_Release(gilstate);
		MANTA::update_pointers();
	}
}

void MANTA::initColorsHigh(float init_r, float init_g, float init_b)
{
	if (!_color_rBig) {
		using_colors = true;
		PyGILState_STATE gilstate = PyGILState_Ensure();
		stringstream ss;
		ss << "manta_color_r = " << init_r << endl;
		ss << "manta_color_g = " << init_g << endl;
		ss << "manta_color_b = " << init_b << endl;
		PyRun_SimpleString(ss.str().c_str());
		PyRun_SimpleString(alloc_colors_high.c_str());
		PyRun_SimpleString(init_colors_high.c_str());
		PyRun_SimpleString(with_fire.c_str());
		PyGILState_Release(gilstate);
		MANTA::update_high_res_pointers();
	}
}

void MANTA::initBlenderRNA(SmokeDomainSettings *sds)
{
//	SmokeModifierData *smd = sds->smd;
//	ModifierData *md = (ModifierData*) smd;

	_alpha = &(sds->alpha);
	_beta = &(sds->beta);
	_dtFactor = &(sds->time_scale);
	_vorticity = &(sds->vorticity);
	_border_collisions = &(sds->border_collisions);
	_burning_rate = &(sds->burning_rate);
	_flame_smoke = &(sds->flame_smoke);
	_flame_smoke_color = sds->flame_smoke_color;
	_flame_vorticity = &(sds->flame_vorticity);
	_flame_ignition = &(sds->flame_ignition);
	_flame_max_temp = &(sds->flame_max_temp);
	_flags = &(sds->flags);
	_manta_solver_res = &(sds->manta_solver_res);
	_fps = 24;//md->scene->r.frs_sec / md->scene->r.frs_sec_base;
	_amplify = &(sds->amplify);
	_strength = &(sds->strength);
	_noise_pos_scale = &(sds->noise_pos_scale);
	_noise_time_anim = &(sds->noise_time_anim);	
}

void MANTA::step()
{
	// Blender computes heat buoyancy, not yet impl. in Manta
	//manta_write_effectors(this);

	PyGILState_STATE gilstate = PyGILState_Ensure();
	std::string py_string_0 = string("manta_step()");
	PyRun_SimpleString(py_string_0.c_str());
	PyGILState_Release(gilstate);
	MANTA::update_pointers();
}

void MANTA::processBurn()
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	std::string py_string_0 = string("process_burn_low()");
	PyRun_SimpleString(py_string_0.c_str());
	PyGILState_Release(gilstate);
	MANTA::update_pointers();
}

void MANTA::updateFlame()
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	std::string py_string_0 = string("update_flame_low()");
	PyRun_SimpleString(py_string_0.c_str());
	PyGILState_Release(gilstate);
	MANTA::update_pointers();
}


MANTA::~MANTA()
{
	cout << "~MANTA()" << endl;
	
	PyGILState_STATE gilstate = PyGILState_Ensure();
	if (_using_heat)
		PyRun_SimpleString(del_heat_low.c_str());
	if (_using_fire)
		PyRun_SimpleString(del_fire_low.c_str());
	if (_using_colors)
		PyRun_SimpleString(del_colors_low.c_str());
	PyRun_SimpleString(del_base_grids_low.c_str());
	PyGILState_Release(gilstate);
}

void MANTA::run_python_string(std::vector<std::string> commands)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	for (std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it) {
		std::string command = *it;
		PyRun_SimpleString(command.c_str());
	}
	PyGILState_Release(gilstate);
}

void MANTA::start_mantaflow()
{
	// THIS STILL CAUSES PROBLEMS IN LINUX
//	string filename = "manta_scene.py";
//	std::vector<std::string> fill = std::vector<std::string>();
//	
//	// Initialize extension classes and wrappers
//	srand(0);
//	PyGILState_STATE gilstate = PyGILState_Ensure();
//	
//	if (!_manta_initialized)
//	{	
//		Pb::setup(filename, fill);
//		_manta_initialized = true;
//	}
//	PyGILState_Release(gilstate);
	
	// Using old initialization setup
	vector<string> args;
	args.push_back("manta_scene.py");
	initializeMantaflow(args);
}

std::string MANTA::get_real_value(const std::string& varName)
{
	ostringstream ss;
	bool is2D = (_resolution == 2);
	
	if (varName == "USING_COLORS")
		ss << (_using_colors ? "True" : "False");
	else if (varName == "USING_HEAT")
		ss << (_using_heat ? "True" : "False");
	else if (varName == "USING_FIRE")
		ss << (_using_fire ? "True" : "False");
	else if (varName == "USE_WAVELETS")
		ss << (*_flags & MOD_SMOKE_HIGHRES ? "True" : "False");
	else if (varName == "SOLVER_DIM")
		ss <<  _manta_solver_res;
	else if (varName == "DO_OPEN")
		ss << (*_border_collisions == SM_BORDER_CLOSED ? "False" : "True");
	else if (varName == "BOUNDCONDITIONS")
		if (*_manta_solver_res == 2) {
			if (*_border_collisions == SM_BORDER_OPEN) ss << "xXyY";
			else if (*_border_collisions == SM_BORDER_VERTICAL) ss << "yY";
			else if (*_border_collisions == SM_BORDER_CLOSED) ss << "";
		}
		if (*_manta_solver_res == 3) {
			if (*_border_collisions == SM_BORDER_OPEN) ss << "xXyYzZ";
			else if (*_border_collisions == SM_BORDER_VERTICAL) ss << "zZ";
			else if (*_border_collisions == SM_BORDER_CLOSED) ss << "";
		}
	else if (varName == "RESX")
		ss << _xRes;
	else if (varName == "RESY")
		if (is2D) {	ss << _zRes;}
		else { 		ss << _yRes;}
	else if (varName == "RESZ")
		if (is2D){	ss << 1;}
		else { 		ss << _zRes;}
	else if (varName == "DT_FACTOR")
		ss << _dtFactor;
	else if (varName == "FPS")
		ss << _fps;
	else if (varName == "VORTICITY")
		ss << *_vorticity / _constantScaling;
	else if (varName == "UPRES")
		ss << *_amplify;
	else if (varName == "HRESX")
		ss << _xResBig;
	else if (varName == "HRESY")
		if (is2D) {	ss << _zResBig;}
		else { 		ss << _yResBig;}
	else if (varName == "HRESZ")
		if (is2D) {	ss << 1;}
		else { 		ss << _zResBig;}
	else if (varName == "WLT_STR")
		ss << *_strength;
	else if (varName == "NOISE_POSSCALE")
		ss << *_noise_pos_scale;
	else if (varName == "NOISE_TIMEANIM")
		ss << *_noise_time_anim;
	else if (varName == "ADVECT_ORDER")
		ss << 2;
	else if (varName == "ALPHA")
		ss << *_alpha;
	else if (varName == "BETA")
		ss << *_beta;
	else if (varName == "BURNING_RATE")
		ss << *_burning_rate;
	else if (varName == "FLAME_SMOKE")
		ss << *_flame_smoke;
	else if (varName == "IGNITION_TEMP")
		ss << *_flame_ignition;
	else if (varName == "MAX_TEMP")
		ss << *_flame_max_temp;
	else if (varName == "FLAME_SMOKE_COLOR_X")
		ss << _flame_smoke_color[0];
	else if (varName == "FLAME_SMOKE_COLOR_Y")
		ss << _flame_smoke_color[1];
	else if (varName == "FLAME_SMOKE_COLOR_Z")
		ss << _flame_smoke_color[2];
	else if (varName == "MANTA_EXPORT_PATH") {
		char parent_dir[1024];
		BLI_split_dir_part(_manta_filepath, parent_dir, sizeof(parent_dir));
		ss << parent_dir;
	} else
		cout << "ERROR: Unknown option:" << varName <<endl;
	return ss.str();
}

std::string MANTA::parse_line(const string& line)
{
	if (line.size() == 0) return "";
	string res = "";
	int currPos = 0, start_del = 0, end_del = -1;
	bool readingVar = false;
	const char delimiter = '$';
	while (currPos < line.size()){
		if(line[currPos] == delimiter && ! readingVar) {
			readingVar 	= true;
			start_del	= currPos + 1;
			res 		+= line.substr(end_del + 1, currPos - end_del -1);
		}
		else if(line[currPos] == delimiter && readingVar){
			readingVar 	= false;
			end_del 	= currPos;
			res 		+= get_real_value(line.substr(start_del, currPos - start_del));
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

std::string MANTA::parse_script(const string& setup_string)
{
	std::istringstream f(setup_string);
	ostringstream res;
	string line = "";
	while(getline(f, line)) {
		res << parse_line(line) << "\n";
	}
	return res.str();
}

void MANTA::manta_export_script(SmokeModifierData *smd)
{
	// Setup low
	std::string manta_script =
		manta_import +
		solver_setup_low +
		uv_setup +
		alloc_base_grids_low;
	
	// Add heat grid low if needed
	if (_using_heat) {
		manta_script += alloc_heat_low;
	}
	
	// Add color grids low if needed
	if (_using_colors) {
		manta_script += alloc_colors_low;
	}
	
	// Add fire grids low if needed
	if (_using_fire) {
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
	if (*_flags & MOD_SMOKE_HIGHRES && _using_colors) {
		manta_script += alloc_colors_high;
	}
	
	// Add fire grids high if needed
	if (*_flags & MOD_SMOKE_HIGHRES && _using_fire) {
		manta_script += alloc_fire_high;
	}

	// Rest of high res setup
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += prep_domain_high + wavelet_turbulence_noise;
	}
	
	// Noise low
	// TODO. Maybe drop this grid, because it can only be used for inflow
	
	// Noise high
	// TODO, Same as noise low
	
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
	
	// Step low
	manta_script += smoke_step_low;
	
	// Step high
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		manta_script += smoke_step_high;
	}
	
	// Fill in missing variables in script
	std::string final_script = MANTA::parse_script(manta_script, smd);
	
	// Add standalone mode (for-loop, gui, ...)
	final_script += standalone;
	
	// Write script
	ofstream myfile;
	myfile.open(_manta_filepath);
	myfile << final_script;
	myfile.close();
}

void MANTA::manta_export_grids(SmokeModifierData *smd)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	
	// Export low res grids
	PyRun_SimpleString(MANTA::parse_script(smoke_export_low, smd).c_str());

	// Export high res grids
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		PyRun_SimpleString(MANTA::parse_script(smoke_export_high, smd).c_str());
	}
	PyGILState_Release(gilstate);
}


string MANTA::get_grid_pointer(std::string gridName, std::string solverName)
{
	if ((gridName == "") && (solverName == "")) {
		return "";
	}
	cout << "getting grid pointer " << gridName<< " , " << solverName <<endl;
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyObject *main = PyImport_AddModule("__main__");
	if (main == NULL){cout << "null" << 1 << endl;return "";}
	PyObject *globals = PyModule_GetDict(main);
	if (globals == NULL){cout << "null" << 12 << endl;return "";}
	PyObject *grid_object = PyDict_GetItemString(globals, gridName.c_str());
	if (grid_object == NULL){cout << "null" << 13 << endl;return "";}
	PyObject* func = PyObject_GetAttrString(grid_object,(char*)"getDataPointer");
	if (func == NULL){cout << "null" << 14 << endl;return "";}
	PyObject* retured_value = PyObject_CallObject(func, NULL);
	PyObject* encoded = PyUnicode_AsUTF8String(retured_value);
	if (retured_value == NULL){cout << "null" << 15 << endl;return "";}
	std::string res = PyBytes_AsString(encoded);
	cout << "Pointer on "<< gridName << " " << res << endl;
	PyGILState_Release(gilstate);		
	return res;
}

void* MANTA::pointer_from_string(const std::string& s){
	stringstream ss(s);
	void *gridPointer = NULL;
	ss >> gridPointer;
	return gridPointer;
}

void MANTA::update_pointers()
{
	cout << "Updating pointers" << endl;
	if (_resolution == 2)
	{
//	TODO
//		cout << "2D" << endl;
//		float* manta_fluid_density = (float* )pointer_from_string(get_grid_pointer("density", "s")); 
//		int* manta_fluid_flags = (int* )pointer_from_string(get_grid_pointer("flags", "s"));
//		if (fluid->_density != NULL){
//			for (int cnt(0); cnt < fluid->xRes() * fluid->yRes() * fluid->zRes(); ++cnt) {
//				fluid->_density[cnt] = 0.;
//				fluid->_manta_flags[cnt] = 2;
//			}
//		}
//		int step = 0;
//		for (int cnty(0);cnty<fluid->yRes(); ++cnty)
//			for(int cntz(0);cntz<fluid->zRes(); ++cntz)
//			{
//				step = fluid->xRes() + cnty * fluid->xRes() + cntz * fluid->xRes()*fluid->yRes(); 
//				if ((step < 0) || (step > fluid->_totalCells)){
//					cout << "UpdatePointers: step is larger than cell dim" << step << endl;
//				}
//				fluid->_density[step] = manta_fluid_density[cnty + cntz*fluid->xRes()];
//				fluid->_manta_flags[step] = manta_fluid_flags[cnty + cntz*fluid->xRes()];
//			}		
	}
	else {
		cout << "3D" << endl;
		_density = (float* )pointer_from_string(get_grid_pointer("density", "s"));
		_manta_flags = (int* )pointer_from_string(get_grid_pointer("flags", "s"));
	}
	
	_manta_inflow = (float* )pointer_from_string(get_grid_pointer("inflow_grid", "s"));
	_fuel_inflow = (float* )pointer_from_string(get_grid_pointer("fuel_inflow", "s"));
	
	if (_using_colors) {
		_color_r = (float* )pointer_from_string(get_grid_pointer("color_r", "s"));
		_color_g = (float* )pointer_from_string(get_grid_pointer("color_g", "s"));
		_color_b = (float* )pointer_from_string(get_grid_pointer("color_b", "s"));
	}
	if (_using_heat) {
		_heat = (float* )pointer_from_string(get_grid_pointer("heat", "s"));
	}
	if (_using_fire) {
		_flame = (float* )pointer_from_string(get_grid_pointer("flame", "s"));
		_fuel = (float* )pointer_from_string(get_grid_pointer("fuel", "s"));
		_react = (float* )pointer_from_string(get_grid_pointer("react", "s"));
	}
	
	_xVelocity = (float* )pointer_from_string(get_grid_pointer("x_vel", "s"));
	_yVelocity = (float* )pointer_from_string(get_grid_pointer("y_vel", "s"));
	_zVelocity = (float* )pointer_from_string(get_grid_pointer("z_vel", "s"));
}

void MANTA::update_high_res_pointers()
{
	cout << "Updating pointers high res" << endl;
	_densityBig = (float* )pointer_from_string(get_grid_pointer("xl_density", "xl"));;
	if (_using_colors) {
		_color_rBig = (float* )pointer_from_string(get_grid_pointer("xl_color_r", "xl"));
		_color_gBig = (float* )pointer_from_string(get_grid_pointer("xl_color_g", "xl"));
		_color_bBig = (float* )pointer_from_string(get_grid_pointer("xl_color_b", "xl"));
	}
	if (_using_fire) {
		_flameBig = (float* )pointer_from_string(get_grid_pointer("xl_flame", "xl"));
		_fuelBig = (float* )pointer_from_string(get_grid_pointer("xl_fuel", "xl"));
		_reactBig = (float* )pointer_from_string(get_grid_pointer("xl_react", "xl"));
	}
}


