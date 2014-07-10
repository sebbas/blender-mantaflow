#include "Python.h"
#include "manta_api.h"
#include "../manta.h"

//#ifdef __cplusplus
//extern "C" {
//#endif
	PyObject * PyInit_Manta(void)
	{
		return Pb::PyInit_Main();
	}
//#ifdef __cplusplus
//}
//#endif