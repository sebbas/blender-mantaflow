/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Manta__MainWnd[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       6,       // signalCount

 // signals: signature, parameters, type, tag, flags
      24,   16,   15,   15, 0x05,
      48,   46,   15,   15, 0x25,
      66,   15,   15,   15, 0x05,
      81,   77,   15,   15, 0x05,
     106,   15,   15,   15, 0x05,
     117,   15,   15,   15, 0x05,

 // slots: signature, parameters, type, tag, flags
     127,   15,   15,   15, 0x0a,
     135,   15,   15,   15, 0x0a,
     142,   15,   15,   15, 0x0a,
     154,  149,   15,   15, 0x0a,
     177,  172,   15,   15, 0x0a,
     219,  197,   15,   15, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__MainWnd[] = {
    "Manta::MainWnd\0\0e,param\0painterEvent(int,int)\0"
    "e\0painterEvent(int)\0wakeMain()\0bgr\0"
    "setBackgroundMesh(Mesh*)\0killMain()\0"
    "exitApp()\0pause()\0play()\0step()\0ctrl\0"
    "addControl(void*)\0file\0screenshot(QString)\0"
    "pos,p0,p1,p2,q0,q1,q2\0"
    "clickLine(QPoint,float,float,float,float,float,float)\0"
};

void Manta::MainWnd::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWnd *_t = static_cast<MainWnd *>(_o);
        switch (_id) {
        case 0: _t->painterEvent((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 1: _t->painterEvent((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->wakeMain(); break;
        case 3: _t->setBackgroundMesh((*reinterpret_cast< Mesh*(*)>(_a[1]))); break;
        case 4: _t->killMain(); break;
        case 5: _t->exitApp(); break;
        case 6: _t->pause(); break;
        case 7: _t->play(); break;
        case 8: _t->step(); break;
        case 9: _t->addControl((*reinterpret_cast< void*(*)>(_a[1]))); break;
        case 10: _t->screenshot((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 11: _t->clickLine((*reinterpret_cast< QPoint(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3])),(*reinterpret_cast< float(*)>(_a[4])),(*reinterpret_cast< float(*)>(_a[5])),(*reinterpret_cast< float(*)>(_a[6])),(*reinterpret_cast< float(*)>(_a[7]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::MainWnd::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::MainWnd::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_Manta__MainWnd,
      qt_meta_data_Manta__MainWnd, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::MainWnd::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::MainWnd::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::MainWnd::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__MainWnd))
        return static_cast<void*>(const_cast< MainWnd*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int Manta::MainWnd::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void Manta::MainWnd::painterEvent(int _t1, int _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 2
void Manta::MainWnd::wakeMain()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void Manta::MainWnd::setBackgroundMesh(Mesh * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void Manta::MainWnd::killMain()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void Manta::MainWnd::exitApp()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}
QT_END_MOC_NAMESPACE
