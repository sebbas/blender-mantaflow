/****************************************************************************
** Meta object code from reading C++ file 'qtmain.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qtmain.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qtmain.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_Manta__GuiThread_t {
    QByteArrayData data[5];
    char stringdata[38];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Manta__GuiThread_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Manta__GuiThread_t qt_meta_stringdata_Manta__GuiThread = {
    {
QT_MOC_LITERAL(0, 0, 16), // "Manta::GuiThread"
QT_MOC_LITERAL(1, 17, 9), // "sendEvent"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 1), // "e"
QT_MOC_LITERAL(4, 30, 7) // "exitApp"

    },
    "Manta::GuiThread\0sendEvent\0\0e\0exitApp"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Manta__GuiThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   24,    2, 0x0a /* Public */,
       4,    0,   27,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void,

       0        // eod
};

void Manta::GuiThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GuiThread *_t = static_cast<GuiThread *>(_o);
        switch (_id) {
        case 0: _t->sendEvent((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->exitApp(); break;
        default: ;
        }
    }
}

const QMetaObject Manta::GuiThread::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Manta__GuiThread.data,
      qt_meta_data_Manta__GuiThread,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *Manta::GuiThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Manta::GuiThread::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__GuiThread.stringdata))
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
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 2)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 2;
    }
    return _id;
}
struct qt_meta_stringdata_Manta__MainThread_t {
    QByteArrayData data[6];
    char stringdata[49];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Manta__MainThread_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Manta__MainThread_t qt_meta_stringdata_Manta__MainThread = {
    {
QT_MOC_LITERAL(0, 0, 17), // "Manta::MainThread"
QT_MOC_LITERAL(1, 18, 9), // "sendToGui"
QT_MOC_LITERAL(2, 28, 0), // ""
QT_MOC_LITERAL(3, 29, 5), // "event"
QT_MOC_LITERAL(4, 35, 6), // "wakeUp"
QT_MOC_LITERAL(5, 42, 6) // "killMe"

    },
    "Manta::MainThread\0sendToGui\0\0event\0"
    "wakeUp\0killMe"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Manta__MainThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    0,   32,    2, 0x0a /* Public */,
       5,    0,   33,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    3,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void Manta::MainThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainThread *_t = static_cast<MainThread *>(_o);
        switch (_id) {
        case 0: _t->sendToGui((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->wakeUp(); break;
        case 2: _t->killMe(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (MainThread::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&MainThread::sendToGui)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject Manta::MainThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_Manta__MainThread.data,
      qt_meta_data_Manta__MainThread,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *Manta::MainThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Manta::MainThread::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__MainThread.stringdata))
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
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void Manta::MainThread::sendToGui(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
