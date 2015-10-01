/****************************************************************************
** Meta object code from reading C++ file 'painter.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "painter.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'painter.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Manta__Painter[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      15,   34,   43,   43, 0x05,

 // slots: signature, parameters, type, tag, flags
      44,   43,   43,   43, 0x0a,
      52,   69,   43,   43, 0x0a,
      77,   90,   43,   43, 0x2a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__Painter[] = {
    "Manta::Painter\0setViewport(Vec3i)\0"
    "gridsize\0\0paint()\0doEvent(int,int)\0"
    "e,param\0doEvent(int)\0e\0"
};

void Manta::Painter::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Painter *_t = static_cast<Painter *>(_o);
        switch (_id) {
        case 0: _t->setViewport((*reinterpret_cast< const Vec3i(*)>(_a[1]))); break;
        case 1: _t->paint(); break;
        case 2: _t->doEvent((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 3: _t->doEvent((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::Painter::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::Painter::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Manta__Painter,
      qt_meta_data_Manta__Painter, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::Painter::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::Painter::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::Painter::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__Painter))
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
    }
    return _id;
}

// SIGNAL 0
void Manta::Painter::setViewport(const Vec3i & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_Manta__LockedObjPainter[] = {

 // content:
       6,       // revision
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

static const char qt_meta_stringdata_Manta__LockedObjPainter[] = {
    "Manta::LockedObjPainter\0"
};

void Manta::LockedObjPainter::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObjectExtraData Manta::LockedObjPainter::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::LockedObjPainter::staticMetaObject = {
    { &Painter::staticMetaObject, qt_meta_stringdata_Manta__LockedObjPainter,
      qt_meta_data_Manta__LockedObjPainter, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::LockedObjPainter::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::LockedObjPainter::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::LockedObjPainter::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__LockedObjPainter))
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
