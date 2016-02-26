#include "MANTA.h"
#include "WTURBULENCE.h"
#include "scenarios/smoke.h"

void Manta_API::start_mantaflow()
{
//	string filename = "manta_scene.py";
//	std::vector<std::string> fill = std::vector<std::string>();
//	
//	// Initialize extension classes and wrappers
//	srand(0);
//	PyGILState_STATE gilstate = PyGILState_Ensure();
//	
//	if (!manta_initialized)
//	{	
//		Pb::setup(filename, fill);
//		manta_initialized = true;
//	}
//	PyGILState_Release(gilstate);
	
	// Using old initialization setup
	vector<string> args;
	args.push_back("manta_scene.py");
	initializeMantaflow(args);
}

std::string Manta_API::get_real_value( const std::string& varName, SmokeModifierData *smd)
{
	ostringstream ss;
	bool is2D = smd->domain->fluid->manta_resoution == 2;
	ModifierData *md = ((ModifierData*) smd);
	if (varName == "UVS_CNT")
		ss << smd->domain->manta_uvs_num ;
	else if (varName == "UPRES")
		ss << smd->domain->amplify;
	else if (varName == "WLT_STR")
		ss << smd->domain->strength ;
	else if (varName == "RES")
		ss << smd->domain->maxres;
	else if (varName == "RESX")
		ss << smd->domain->fluid->_xRes;
	else if (varName == "RESY")
		if (is2D) {	ss <<  smd->domain->fluid->_zRes;}
		else { 		ss <<  smd->domain->fluid->_yRes;}
	else if (varName == "RESZ")
		if (is2D){	ss << 1;}
		else { 		ss << smd->domain->fluid->_zRes;}
	else if (varName == "SOLVER_DIM")
		ss <<  smd->domain->manta_solver_res;
	else if (varName == "NOISE_POSSCALE")
		ss << smd->domain->noise_pos_scale;
	else if (varName == "NOISE_TIMEANIM")
		ss << smd->domain->noise_time_anim;
	else if (varName == "HRESX")
		ss << smd->domain->wt->getResBig()[0];
	else if (varName == "HRESY")
		if (is2D) {	ss << smd->domain->wt->getResBig()[2];}
		else { 		ss << smd->domain->wt->getResBig()[1];}
	else if (varName == "HRESZ")
		if (is2D) {	ss << 1;}
		else { 		ss << smd->domain->wt->getResBig()[2];}
	else if (varName == "TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "XL_TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "USE_WAVELETS")
		ss << ((smd->domain->flags & MOD_SMOKE_HIGHRES) ? "True" : "False");
	else if (varName == "BUYO_X")
		ss << 0.;
	else if (varName == "BUYO_Y")
		ss << 0.;
	else if (varName == "BUYO_Z")
		ss << (-smd->domain->beta);
	else if (varName == "ALPHA")
		ss << (-smd->domain->alpha);
	else if (varName == "BETA")
		ss << (-smd->domain->beta);
	else if (varName == "ADVECT_ORDER")
		ss << 2;
	else if (varName == "MANTA_EXPORT_PATH") {
		char parent_dir[1024];
		BLI_split_dir_part(smd->domain->_manta_filepath, parent_dir, sizeof(parent_dir));
		ss << parent_dir;
	} else if (varName == "VORTICITY") {
		ss << smd->domain->vorticity / smd->domain->fluid->_constantScaling;
	} else if (varName == "BOUNDCONDITIONS") {
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
	}
	else if (varName == "DO_OPEN")
		ss << ((smd->domain->border_collisions == SM_BORDER_CLOSED) ? "False" : "True");
	else if (varName == "GRAVITY")
		ss << "vec3(0,0,-0.981)";
	else if (varName == "ABS_FLOW")
		ss << ((smd->flow->flags & MOD_SMOKE_FLOW_ABSOLUTE)? "True" : "False");
	else if (varName == "DENSITY_MEM")
		ss << smd->domain->fluid->_density;
	else if (varName == "DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->total_cells;
	else if (varName == "XL_DENSITY_MEM")
		ss << smd->domain->wt->_densityBig;
	else if (varName == "XL_DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->wt->_xResBig * smd->domain->wt->_yResBig * smd->domain->wt->_zResBig;
	else if (varName == "BURNING_RATE")
		ss << (smd->domain->burning_rate);
	else if (varName == "FLAME_SMOKE")
		ss << (smd->domain->flame_smoke);
	else if (varName == "IGNITION_TEMP")
		ss << (smd->domain->flame_ignition);
	else if (varName == "MAX_TEMP")
		ss << (smd->domain->flame_max_temp);
	else if (varName == "DT")
		ss << smd->domain->fluid->_dt;
	else if (varName == "FPS")
		ss << md->scene->r.frs_sec / md->scene->r.frs_sec_base;
	else if (varName == "DT_FACTOR")
		ss << smd->domain->time_scale;
	else if (varName == "FLAME_SMOKE_COLOR_X")
		ss << smd->domain->flame_smoke_color[0];
	else if (varName == "FLAME_SMOKE_COLOR_Y")
		ss << smd->domain->flame_smoke_color[1];
	else if (varName == "FLAME_SMOKE_COLOR_Z")
		ss << smd->domain->flame_smoke_color[2];
	else if (varName == "USING_COLORS")
		ss << ((smd->domain->fluid->using_colors) ? "True" : "False");
	else if (varName == "USING_HEAT")
		ss << ((smd->domain->fluid->using_heat) ? "True" : "False");
	else if (varName == "USING_FIRE")
		ss << ((smd->domain->fluid->using_fire) ? "True" : "False");
	else 
		cout << "ERROR: Unknown option:" << varName <<endl;
	return ss.str();
}

std::string Manta_API::parse_line(const string& line, SmokeModifierData *smd)
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
			res 		+= get_real_value(line.substr(start_del, currPos - start_del), smd);
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

std::string Manta_API::parse_script(const string& setup_string, SmokeModifierData *smd)
{
	std::istringstream f(setup_string);
	ostringstream res;
	string line = "";
	while(getline(f, line)) {
		res << parse_line(line, smd) << "\n";
	}
	return res.str();
}

void Manta_API::manta_export_script(SmokeModifierData *smd)
{
	// Setup low
	std::string manta_script =
		manta_import +
		solver_setup_low +
		uv_setup +
		alloc_base_grids_low;
	
	// Add heat grid low if needed
	if (smd->domain->fluid->using_heat) {
		manta_script += alloc_heat_low;
	}
	
	// Add color grids low if needed
	if (smd->domain->fluid->using_colors) {
		manta_script += alloc_colors_low;
	}
	
	// Add fire grids low if needed
	if (smd->domain->fluid->using_fire) {
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
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->fluid->using_colors) {
		manta_script += alloc_colors_high;
	}
	
	// Add fire grids high if needed
	if (smd->domain->flags & MOD_SMOKE_HIGHRES && smd->domain->fluid->using_fire) {
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
	std::string final_script = Manta_API::parse_script(manta_script, smd);
	
	// Add standalone mode (for-loop, gui, ...)
	final_script += standalone;
	
	// Write script
	ofstream myfile;
	myfile.open(smd->domain->_manta_filepath);
	myfile << final_script;
	myfile.close();
}

void Manta_API::manta_export_grids(SmokeModifierData *smd)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	
	// Export low res grids
	PyRun_SimpleString(Manta_API::parse_script(smoke_export_low, smd).c_str());

	// Export high res grids
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		PyRun_SimpleString(Manta_API::parse_script(smoke_export_high, smd).c_str());
	}
	PyGILState_Release(gilstate);
}


string Manta_API::get_grid_pointer(std::string gridName, std::string solverName)
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

void* Manta_API::pointer_from_string(const std::string& s){
	stringstream ss(s);
	void *gridPointer = NULL;
	ss >> gridPointer;
	return gridPointer;
}

void Manta_API::update_pointers(FLUID_3D *fluid)
{
	cout << "Updating pointers" << endl;
	if (fluid->manta_resoution == 2)
	{
		cout << "2D" << endl;
		float* manta_fluid_density = (float* )pointer_from_string(get_grid_pointer("density", "s")); 
		int* manta_fluid_flags = (int* )pointer_from_string(get_grid_pointer("flags", "s"));
		if (fluid->_density != NULL){
			for (int cnt(0); cnt < fluid->xRes() * fluid->yRes() * fluid->zRes(); ++cnt) {
				fluid->_density[cnt] = 0.;
				fluid->_manta_flags[cnt] = 2;
			}
		}
		int step = 0;
		for (int cnty(0);cnty<fluid->yRes(); ++cnty)
			for(int cntz(0);cntz<fluid->zRes(); ++cntz)
			{
				step = fluid->xRes() + cnty * fluid->xRes() + cntz * fluid->xRes()*fluid->yRes(); 
				if ((step < 0) || (step > fluid->_totalCells)){
					cout << "UpdatePointers: step is larger than cell dim" << step << endl;
				}
				fluid->_density[step] = manta_fluid_density[cnty + cntz*fluid->xRes()];
				fluid->_manta_flags[step] = manta_fluid_flags[cnty + cntz*fluid->xRes()];
			}		
	}
	else {
		cout << "3D" << endl;
		fluid->_density = (float* )pointer_from_string(get_grid_pointer("density", "s"));	
		fluid->_manta_flags = (int* )pointer_from_string(get_grid_pointer("flags", "s"));
	}
	
	fluid->_manta_inflow = (float* )pointer_from_string(get_grid_pointer("inflow_grid", "s"));
	fluid->_fuel_inflow = (float* )pointer_from_string(get_grid_pointer("fuel_inflow", "s"));
	if (fluid-> manta_resoution == 2) {return;}
	if (fluid->using_colors) {
		fluid->_color_r = (float* )pointer_from_string(get_grid_pointer("color_r", "s"));
		fluid->_color_g = (float* )pointer_from_string(get_grid_pointer("color_g", "s"));
		fluid->_color_b = (float* )pointer_from_string(get_grid_pointer("color_b", "s"));
	}
	if (fluid->using_heat) {
		fluid->_heat = (float* )pointer_from_string(get_grid_pointer("heat", "s"));
	}
	if (fluid->using_fire) {
		fluid->_flame = (float* )pointer_from_string(get_grid_pointer("flame", "s"));
		fluid->_fuel = (float* )pointer_from_string(get_grid_pointer("fuel", "s"));
		fluid->_react = (float* )pointer_from_string(get_grid_pointer("react", "s"));
	}
	fluid->_xVelocity = (float* )pointer_from_string(get_grid_pointer("x_vel", "s"));
	fluid->_yVelocity = (float* )pointer_from_string(get_grid_pointer("y_vel", "s"));
	fluid->_zVelocity = (float* )pointer_from_string(get_grid_pointer("z_vel", "s"));
}

void Manta_API::update_high_res_pointers(WTURBULENCE *wt)
{
	cout << "Updating pointers high res" << endl;
	wt->_densityBig = (float* )pointer_from_string(get_grid_pointer("xl_density", "xl"));;
	if (wt->using_colors) {
		wt->_color_rBig = (float* )pointer_from_string(get_grid_pointer("xl_color_r", "xl"));
		wt->_color_gBig = (float* )pointer_from_string(get_grid_pointer("xl_color_g", "xl"));
		wt->_color_bBig = (float* )pointer_from_string(get_grid_pointer("xl_color_b", "xl"));
	}
	if (wt->using_fire) {
		wt->_flameBig = (float* )pointer_from_string(get_grid_pointer("xl_flame", "xl"));
		wt->_fuelBig = (float* )pointer_from_string(get_grid_pointer("xl_fuel", "xl"));
		wt->_reactBig = (float* )pointer_from_string(get_grid_pointer("xl_react", "xl"));
	}
}

