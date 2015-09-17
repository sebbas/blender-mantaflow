/****************************************************************************
** Meta object code from reading C++ file 'painter.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.2.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "painter.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'painter.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.2.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_Manta__Painter_t {
    QByteArrayData data[9];
    char stringdata[66];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    offsetof(qt_meta_stringdata_Manta__Painter_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData) \
    )
static const qt_meta_stringdata_Manta__Painter_t qt_meta_stringdata_Manta__Painter = {
    {
QT_MOC_LITERAL(0, 0, 14),
QT_MOC_LITERAL(1, 15, 11),
QT_MOC_LITERAL(2, 27, 0),
QT_MOC_LITERAL(3, 28, 5),
QT_MOC_LITERAL(4, 34, 8),
QT_MOC_LITERAL(5, 43, 5),
QT_MOC_LITERAL(6, 49, 7),
QT_MOC_LITERAL(7, 57, 1),
QT_MOC_LITERAL(8, 59, 5)
    },
    "Manta::Painter\0setViewport\0\0Vec3i\0"
    "gridsize\0paint\0doEvent\0e\0param\0"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Manta__Painter[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   34,    2, 0x06,

 // slots: name, argc, parameters, tag, flags
       5,    0,   37,    2, 0x0a,
       6,    2,   38,    2, 0x0a,
       6,    1,   43,    2, 0x2a,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    7,    8,
    QMetaType::Void, QMetaType::Int,    7,

       0        // eod
};

void Manta::Painter::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Painter *_t = static_cast<Painter *>(_o);
        switch (_id) {
        case 0: _t->setViewport((*reinterpret_cast< const Vec3i(*)>(_a[1]))); break;
        case 1: _t->paint(); break;
        case 2: _t->doEvent((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 3: _t->doEvent((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (Painter::*_t)(const Vec3i & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Painter::setViewport)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject Manta::Painter::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Manta__Painter.data,
      qt_meta_data_Manta__Painter,  qt_static_metacall, 0, 0}
};


const QMetaObject *Manta::Painter::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Manta::Painter::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__Painter.stringdata))
        return static_cast<void*>(const_cast< Painter*>(this));
    return QObject::qt_metacast(_clname);
}

int Manta::Painter::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void Manta::Painter::setViewport(const Vec3i & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
struct qt_meta_stringdata_Manta__LockedObjPainter_t {
    QByteArrayData data[1];
    char stringdata[25];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    offsetof(qt_meta_stringdata_Manta__LockedObjPainter_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData) \
    )
static const qt_meta_stringdata_Manta__LockedObjPainter_t qt_meta_stringdata_Manta__LockedObjPainter = {
    {
QT_MOC_LITERAL(0, 0, 23)
    },
    "Manta::LockedObjPainter\0"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Manta__LockedObjPainter[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

void Manta::LockedObjPainter::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObject Manta::LockedObjPainter::staticMetaObject = {
    { &Painter::staticMetaObject, qt_meta_stringdata_Manta__LockedObjPainter.data,
      qt_meta_data_Manta__LockedObjPainter,  qt_static_metacall, 0, 0}
};


const QMetaObject *Manta::LockedObjPainter::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Manta::LockedObjPainter::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__LockedObjPainter.stringdata))
        return static_cast<void*>(const_cast< LockedObjPainter*>(this));
    return Painter::qt_metacast(_clname);
}

int Manta::LockedObjPainter::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Painter::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
QT_END_MOC_NAMESPACE
