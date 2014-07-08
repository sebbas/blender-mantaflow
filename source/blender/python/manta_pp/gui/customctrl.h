




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * GUI extension from python
 *
 ******************************************************************************/

#ifndef _CUSTOMCTRL_H__
#define _CUSTOMCTRL_H__

#include <QSlider>
#include <QLabel>
#include <QCheckBox>
#include <QBoxLayout>
#include "manta.h"

namespace Manta {

// fwd decl.
class Mesh;
class GuiThread;
class MainThread;
	
//! Interface for python declared controls
class CustomControl : public PbClass {public:
	CustomControl(); static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "CustomControl::CustomControl" ); { ArgLocker _lock;  obj = new CustomControl(); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"CustomControl::CustomControl" ); return 0; } catch(std::exception& e) { pbSetError("CustomControl::CustomControl",e.what()); return -1; } }
	
	virtual void init(QBoxLayout* layout) {};
 protected: public: PbArgs _args;}
#define _C_CustomControl
;

//! Checkbox with attached text display
class TextCheckbox : public QCheckBox {
	Q_OBJECT
public:
	TextCheckbox(const std::string& name, bool val);
	void attach(QBoxLayout* layout);
	void set(bool v);
	bool get();
	
public slots:
	void update(int v);
		
protected:
	bool mVal;
	QLabel* mLabel;    
	QString mSName;    
};

//! Slider with attached text display
class TextSlider : public QSlider {
	Q_OBJECT
public:
	TextSlider(const std::string& name, float val, float min, float max);
	void attach(QBoxLayout* layout);
	void set(float v);
	float get();
	
public slots:
	void update(int v);
		
protected:
	float mMin, mMax, mScale;
	QLabel* mLabel;    
	QString mSName;    
};
	
//! Links a slider control

class CustomSlider : public CustomControl {public:
	CustomSlider(std::string text, float val, float min, float max); static int _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "CustomSlider::CustomSlider" ); { ArgLocker _lock; std::string text = _args.get<std::string >("text",0,&_lock); float val = _args.get<float >("val",1,&_lock); float min = _args.get<float >("min",2,&_lock); float max = _args.get<float >("max",3,&_lock);  obj = new CustomSlider(text,val,min,max); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"CustomSlider::CustomSlider" ); return 0; } catch(std::exception& e) { pbSetError("CustomSlider::CustomSlider",e.what()); return -1; } }
	virtual void init(QBoxLayout* layout);
	
	float get(); static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CustomSlider* pbo = dynamic_cast<CustomSlider*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "CustomSlider::get"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->get());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CustomSlider::get"); return _retval; } catch(std::exception& e) { pbSetError("CustomSlider::get",e.what()); return 0; } }
	void set(float v); static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CustomSlider* pbo = dynamic_cast<CustomSlider*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "CustomSlider::set"); PyObject *_retval = 0; { ArgLocker _lock; float v = _args.get<float >("v",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->set(v);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CustomSlider::set"); return _retval; } catch(std::exception& e) { pbSetError("CustomSlider::set",e.what()); return 0; } }
	
protected:
	float mMin, mMax, mVal;
	std::string mSName; 	TextSlider* mSlider; public: PbArgs _args;}
#define _C_CustomSlider
;

//! Links a checkbox control

class CustomCheckbox : public CustomControl {public:
	CustomCheckbox(std::string text, bool val); static int _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "CustomCheckbox::CustomCheckbox" ); { ArgLocker _lock; std::string text = _args.get<std::string >("text",0,&_lock); bool val = _args.get<bool >("val",1,&_lock);  obj = new CustomCheckbox(text,val); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"CustomCheckbox::CustomCheckbox" ); return 0; } catch(std::exception& e) { pbSetError("CustomCheckbox::CustomCheckbox",e.what()); return -1; } }
	virtual void init(QBoxLayout* layout);
	
	bool get(); static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CustomCheckbox* pbo = dynamic_cast<CustomCheckbox*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "CustomCheckbox::get"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->get());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CustomCheckbox::get"); return _retval; } catch(std::exception& e) { pbSetError("CustomCheckbox::get",e.what()); return 0; } }
	void set(bool v); static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CustomCheckbox* pbo = dynamic_cast<CustomCheckbox*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "CustomCheckbox::set"); PyObject *_retval = 0; { ArgLocker _lock; bool v = _args.get<bool >("v",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->set(v);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CustomCheckbox::set"); return _retval; } catch(std::exception& e) { pbSetError("CustomCheckbox::set",e.what()); return 0; } }
	
protected:
	bool mVal;
	std::string mSName; 	TextCheckbox* mCheckbox; public: PbArgs _args;}
#define _C_CustomCheckbox
;
	

//! GUI adapter class to call from Python
class Gui : public PbClass {public:
	Gui(); static int _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Gui::Gui" ); { ArgLocker _lock;  obj = new Gui(); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Gui::Gui" ); return 0; } catch(std::exception& e) { pbSetError("Gui::Gui",e.what()); return -1; } }
	
	void setBackgroundMesh(Mesh* m); static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::setBackgroundMesh"); PyObject *_retval = 0; { ArgLocker _lock; Mesh* m = _args.getPtr<Mesh >("m",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setBackgroundMesh(m);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::setBackgroundMesh"); return _retval; } catch(std::exception& e) { pbSetError("Gui::setBackgroundMesh",e.what()); return 0; } }
	void show(); static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::show"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->show();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::show"); return _retval; } catch(std::exception& e) { pbSetError("Gui::show",e.what()); return 0; } }
	void update(); static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::update"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->update();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::update"); return _retval; } catch(std::exception& e) { pbSetError("Gui::update",e.what()); return 0; } }
	void pause(); static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::pause"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->pause();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::pause"); return _retval; } catch(std::exception& e) { pbSetError("Gui::pause",e.what()); return 0; } }
	PbClass* addControl(PbType t); static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::addControl"); PyObject *_retval = 0; { ArgLocker _lock; PbType t = _args.get<PbType >("t",0,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->addControl(t));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::addControl"); return _retval; } catch(std::exception& e) { pbSetError("Gui::addControl",e.what()); return 0; } }
	void screenshot(std::string filename); static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Gui* pbo = dynamic_cast<Gui*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Gui::screenshot"); PyObject *_retval = 0; { ArgLocker _lock; std::string filename = _args.get<std::string >("filename",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->screenshot(filename);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Gui::screenshot"); return _retval; } catch(std::exception& e) { pbSetError("Gui::screenshot",e.what()); return 0; } }
	
protected:
	GuiThread* mGuiPtr; 	MainThread* mMainPtr; public: PbArgs _args;}
#define _C_Gui
;
	
} // namespace

#endif



