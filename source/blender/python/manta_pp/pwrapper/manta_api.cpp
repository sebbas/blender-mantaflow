#include "Python.h"
//#include "manta_api.h"
#include "manta.h"
#include "../general.h"

using namespace std;
using namespace Manta;

#if PY_MAJOR_VERSION >= 3
typedef wchar_t pyChar;
typedef wstring pyString;
#else
typedef char pyChar;
typedef string pyString;
#endif

PyObject * PyInit_Manta(void)
{
	return Pb::PyInit_Main_Obj();
}
