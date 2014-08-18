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
#include "../../../source/blender/blenlib/BLI_path_util.h"

void export_force_fields(int size_x, int size_y, int size_z, float *f_x, float*f_y, float*f_z);/*defined in pymain.cpp*/
extern "C" void manta_write_effectors(struct Scene *s, struct SmokeModifierData *smd); /*defined in smoke_api.cpp*/

/*for passing to detached thread*/
struct manta_arg_struct {
	Scene s;
	SmokeModifierData smd;
};

static pthread_t manta_thread;

void runMantaScript(const string&, vector<string>& args);//defined in manta_pp/pwrapper/pymain.cpp

extern "C" bool manta_check_grid_size(struct FLUID_3D *fluid, int dimX, int dimY, int dimZ);

extern "C" int read_mantaflow_sim(struct SmokeDomainSettings *sds, char *name, bool read_wavelets);

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

static bool manta_sim_running=true;

void create_manta_folder();

void *run_manta_scene_thread(void *threadid);

void *run_manta_sim_thread(void *threadid);

void run_manta_scene(Scene *scene, SmokeModifierData *smd);

void stop_manta_sim();

void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd);

void manta_sim_step(int frame);

std::string getRealValue(SmokeModifierData *sds, Scene *s, const string& varName);

std::string parseLine(SmokeModifierData *sds, Scene *s, const string& line);

void parseFile(SmokeModifierData *sds, Scene *s, char *file);

#endif /* MANTA_H */

