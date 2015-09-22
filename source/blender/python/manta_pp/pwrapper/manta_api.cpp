#include "Python.h"
#include "manta_api.h"
#include "manta.h"

using namespace Manta;

PyObject * PyInit_Manta(void)
{
	return Pb::PyInit_Main_Obj();
}

