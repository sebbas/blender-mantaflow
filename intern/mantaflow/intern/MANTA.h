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

/** \file mantaflow/intern/MANTA.h
 *  \ingroup mantaflow
 */

#ifndef MANTA_H
#define MANTA_H

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include "smoke.h"
#include "Python.h"

#include "DNA_smoke_types.h"
#include "DNA_modifier_types.h"

#include "BLI_path_util.h"
#include "BLI_utildefines.h"

struct MANTA
{
public:
	MANTA(int *res, float dx, float dtdef, int init_heat, int init_fire, int init_colors, struct SmokeModifierData *smd);
	MANTA() {};
	virtual ~MANTA();
	
	// All step functions
	void step();
	void processBurn();
	void updateFlame();
	
	// Initialization functions
	void initHeat();
	void initFire();
	void initColors(float init_r, float init_g, float init_b);
	void initFireHigh();
	void initColorsHigh(float init_r, float init_g, float init_b);
	void initBlenderRNA(SmokeDomainSettings *sds);

	// Getters
	const float* xVelocity() { return _xVelocity; };
	const float* yVelocity() { return _yVelocity; }; 
	const float* zVelocity() { return _zVelocity; }; 

	int xRes() const { return _xRes; };
	int yRes() const { return _yRes; };
	int zRes() const { return _zRes; };
	
	int withHighRes() { return *_flags & MOD_SMOKE_HIGHRES };
	
private:
	void start_mantaflow();
	void update_pointers();
	void update_high_res_pointers();
	void run_python_string(std::vector<std::string> commands);
	std::string get_real_value(const std::string& varName);
	std::string parse_line(const string& line);
	std::string parse_script(const string& setup_string);
	void manta_export_script(struct SmokeModifierData *smd);
	void manta_export_grids(struct SmokeModifierData *smd);
	std::string get_grid_pointer(std::string gridName, std::string solverName);
	void* pointer_from_string(const std::string& s);
	void initializeMantaflow(vector<string>& args);
	void initBlenderRNA(SmokeDomainSettings *sds);
	
	bool _using_heat;
	bool _using_fire;
	bool _using_colors;
	int _resolution;
	bool _manta_initialized;
	
	// dimensions
	int _xRes, _yRes, _zRes, _maxRes;
	size_t _totalCells;
	float _dx;

	// fields
	float* _density;
	float* _heat;
	float* _xVelocity;
	float* _yVelocity;
	float* _zVelocity;
	float* _xVelocityOb;
	float* _yVelocityOb;
	float* _zVelocityOb;
	float* _xForce;
	float* _yForce;
	float* _zForce;
	unsigned char*  _obstacles; /* only used (useful) for static obstacles like domain boundaries */
	unsigned char*  _obstaclesAnim;
	float* _manta_inflow;
	float* _fuel_inflow;

	// fire simulation
	float *_flame;
	float *_fuel;
	float *_react;

	// smoke color
	float *_color_r;
	float *_color_g;
	float *_color_b;

	int *_manta_flags;

	// simulation constants
	float _dt;
	float _vorticityEps;
	float _heatDiffusion;
	float _tempAmb; /* ambient temperature */
	float _constantScaling;

	// RNA pointers
	float *_alpha;
	float *_beta;
	float *_dtFactor;
	float *_vorticity;
	int *_border_collisions;
	float *_burning_rate;
	float *_flame_smoke;
	float *_flame_smoke_color;
	float *_flame_vorticity;
	float *_flame_ignition;
	float *_flame_max_temp;
	int *_flags;
	int *_manta_solver_res;
	int _fps;
	int _fps_base;
	float *_strength;
	float *_noise_pos_scale;
	float *_noise_time_anim;
	char *_manta_filepath;
	int *_amplify;
	
	// High
	int _xResBig;
	int _yResBig;
	int _zResBig;

	float* _densityBig;
	float* _flameBig;
	float* _fuelBig;
	float* _reactBig;

	float* _color_rBig;
	float* _color_gBig;
	float* _color_bBig;

	// texture coordinates for noise
	float* _tcU;
	float* _tcV;
	float* _tcW;
};

#endif