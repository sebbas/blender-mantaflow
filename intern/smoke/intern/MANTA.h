#ifndef MANTA_H
#define MANTA_H
#include "FLUID_3D.h"
#include "zlib.h"
#include "../../../source/blender/makesdna/DNA_scene_types.h"
#include "../../../source/blender/makesdna/DNA_modifier_types.h"
#include "../../../source/blender/makesdna/DNA_smoke_types.h"
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <pthread.h>
#include <Python.h>
#include <vector>

#include "BLI_path_util.h"

#ifdef WIN32
#include "BLI_winstuff.h"
#endif

void export_force_fields(int size_x, int size_y, int size_z, float *f_x, float*f_y, float*f_z);/*defined in pymain.cpp*/
void export_em_fields(float *em_map, float flow_density, int min_x, int min_y, int min_z, int max_x, int max_y, int max_z, int d_x, int d_y, int d_z, float *inf, float *vel);/*defined in pymain.cpp*/
extern "C" void manta_write_effectors(struct FLUID_3D *fluid); /*defined in smoke_api.cpp*/
void initializeMantaflow(vector<string>& args);//defined in manta_pp/pwrapper/pymain.cpp

/*for passing to detached thread*/
struct manta_arg_struct {
	Scene s;
	SmokeModifierData smd;
};

static bool manta_sim_running=true;

extern "C" bool manta_check_grid_size(struct FLUID_3D *fluid, int dimX, int dimY, int dimZ);

extern "C" int read_mantaflow_sim(struct SmokeDomainSettings *sds, char *name, bool read_wavelets);

class Manta_API{
private:	
	Manta_API() {}
	Manta_API(const Manta_API &);	 
	Manta_API & operator=(const Manta_API &);
public:
	~Manta_API();	 
	Manta_API(int *res, float dx, float dtdef, int init_heat, int init_fire, int init_colors, struct SmokeDomainSettings *sds);
	void initBlenderRNA(float *alpha, float *beta, float *dt_factor, float *vorticity, int *border_colli, float *burning_rate,
						float *flame_smoke, float *flame_smoke_color, float *flame_vorticity, float *ignition_temp, float *max_temp);
	int _totalCells;
	int _xRes, _yRes, _zRes;
	float _res;
	int _slabSize;
	float _dt,_dx;
	float* _density;
	float* _xVelocity;
	float* _yVelocity;
	float* _zVelocity;
	float* _xVelocityOb;
	float* _yVelocityOb;
	float* _zVelocityOb;
	float* _xForce;
	float* _yForce;
	float* _zForce;
	float *_alpha; // for the buoyancy density term <-- as pointer to get blender RNA in here
	float *_beta; // was _buoyancy <-- as pointer to get blender RNA in here
	
	float *_dtFactor;
	float *_vorticityRNA;	// RNA-pointer.
	int *_borderColli; // border collision rules <-- as pointer to get blender RNA in here
	float *_burning_rate; // RNA pointer
	float *_flame_smoke; // RNA pointer
	float *_flame_smoke_color; // RNA pointer
	float *_flame_vorticity; // RNA pointer
	float *_ignition_temp; // RNA pointer
	float *_max_temp; // RNA pointer
	
	unsigned char*  _obstacles; /* only used (useful) for static obstacles like domain */
	void step(float dt, float gravity[3]);
//	void runMantaScript(const string&, vector<string>& args);//defined in manta_pp/pwrapper/pymain.cpp
	
	void indent_ss(stringstream& ss, int indent);
	
	void manta_gen_noise(stringstream& ss, char* solver, int indent, char *noise, int seed, bool load, bool clamp, float clampNeg, float clampPos, float valScale, float valOffset, float timeAnim);
	
	void manta_solve_pressure(stringstream& ss, char *flags, char *vel, char *pressure, bool useResNorms, int openBound, int solver_res,float cgMaxIterFac=1.0, float cgAccuracy = 0.01);
	
	void manta_advect_SemiLagr(stringstream& ss, int indent, char *flags, char *vel, char *grid, int order);
	
	/*create solver, handle 2D case*/
	void manta_create_solver(stringstream& ss, char *name, char *nick, char *grid_size_name, int x_res, int y_res, int z_res, int dim);
	
	inline bool file_exists (const std::string& name);
	
	/*blender transforms obj coords to [-1,1]. This method transforms them back*/
	void add_mesh_transform_method(stringstream& ss);
	
	void manta_cache_path(char *filepath);
	
	void create_manta_folder();
	
	void *run_manta_scene_thread(void *threadid);
	
	void run_manta_sim_highRes(WTURBULENCE *wt);
	
	void run_manta_scene(Manta_API * fluid);
	
	void stop_manta_sim();
	
	static void run_manta_sim_file_lowRes(SmokeModifierData *smd);
	
	static void run_manta_sim_file_highRes(SmokeModifierData *smd);
	
	void manta_sim_step(int frame);
	
	static std::string getRealValue(const string& varName, SmokeModifierData *sds);
	
	static std::string parseLine(const string& line, SmokeModifierData *sds);
	
	static std::string parseScript(const string& setup_string, SmokeModifierData *sds);	
	
	pthread_t manta_thread;
	
	static void * pointerFromString(const std::string& s);
	
	static string gridNameFromType(const string& type);
	static void addGrid(void * data,string name, string type, int x, int y, int z, bool is2D);
	static void addAdaptiveGrid(void * data, string gridName, string solverName, string type,int minX, int minY, int minZ, int maxX, int maxY, int maxZ, bool is2D);
	static void export_obstacles(float *data, int x, int y, int z, bool is2D);
	
	static std::string getGridPointer(string gridName, string solverName);
	static void updatePointers(FLUID_3D *fluid);
	static void updateHighResPointers(WTURBULENCE *wt);
	static void manta_export_grids(SmokeModifierData *smd);
	static std::string get_manta_smoke_script(bool highRes, SmokeModifierData *smd);
};



#endif /* MANTA_H */

