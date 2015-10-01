/****************************************************************************
** Meta object code from reading C++ file 'meshpainter.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "meshpainter.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'meshpainter.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_Manta__MeshPainter_t {
    QByteArrayData data[5];
    char stringdata[48];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Manta__MeshPainter_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Manta__MeshPainter_t qt_meta_stringdata_Manta__MeshPainter = {
    {
QT_MOC_LITERAL(0, 0, 18), // "Manta::MeshPainter"
QT_MOC_LITERAL(1, 19, 17), // "setBackgroundMesh"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 5), // "Mesh*"
QT_MOC_LITERAL(4, 44, 3) // "bgr"

    },
    "Manta::MeshPainter\0setBackgroundMesh\0"
    "\0Mesh*\0bgr"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Manta__MeshPainter[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   19,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,

       0        // eod
};

void Manta::MeshPainter::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MeshPainter *_t = static_cast<MeshPainter *>(_o);
        switch (_id) {
        case 0: _t->setBackgroundMesh((*reinterpret_cast< Mesh*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject Manta::MeshPainter::staticMetaObject = {
    { &LockedObjPainter::staticMetaObject, qt_meta_stringdata_Manta__MeshPainter.data,
      qt_meta_data_Manta__MeshPainter,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *Manta::MeshPainter::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Manta::MeshPainter::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__MeshPainter.stringdata))
        return static_cast<void*>(const_cast< MeshPainter*>(this));
    return LockedObjPainter::qt_metacast(_clname);
}

int Manta::MeshPainter::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LockedObjPainter::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
