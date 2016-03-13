#include "Python.h"
#include "manta_python_API.h"
#include "manta.h"

using namespace Manta;

PyObject * PyInit_Manta(void)
{
	return Pb::PyInit_Main_Obj();
}

