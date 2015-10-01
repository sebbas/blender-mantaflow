/****************************************************************************
** Meta object code from reading C++ file 'qtmain.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qtmain.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qtmain.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Manta__GuiThread[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      17,   32,   34,   34, 0x0a,
      35,   34,   34,   34, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__GuiThread[] = {
    "Manta::GuiThread\0sendEvent(int)\0e\0\0"
    "exitApp()\0"
};

void Manta::GuiThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        GuiThread *_t = static_cast<GuiThread *>(_o);
        switch (_id) {
        case 0: _t->sendEvent((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->exitApp(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::GuiThread::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::GuiThread::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Manta__GuiThread,
      qt_meta_data_Manta__GuiThread, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::GuiThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::GuiThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::GuiThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__GuiThread))
        return static_cast<void*>(const_cast< GuiThread*>(this));
    return QObject::qt_metacast(_clname);
}

int Manta::GuiThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_Manta__MainThread[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      18,   33,   39,   39, 0x05,

 // slots: signature, parameters, type, tag, flags
      40,   39,   39,   39, 0x0a,
      49,   39,   39,   39, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__MainThread[] = {
    "Manta::MainThread\0sendToGui(int)\0event\0"
    "\0wakeUp()\0killMe()\0"
};

void Manta::MainThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainThread *_t = static_cast<MainThread *>(_o);
        switch (_id) {
        case 0: _t->sendToGui((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->wakeUp(); break;
        case 2: _t->killMe(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::MainThread::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::MainThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_Manta__MainThread,
      qt_meta_data_Manta__MainThread, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::MainThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::MainThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::MainThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__MainThread))
        return static_cast<void*>(const_cast< MainThread*>(this));
    return QThread::qt_metacast(_clname);
}

int Manta::MainThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void Manta::MainThread::sendToGui(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
