/****************************************************************************
** Meta object code from reading C++ file 'customctrl.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "customctrl.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'customctrl.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Manta__TextCheckbox[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      20,   32,   34,   34, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__TextCheckbox[] = {
    "Manta::TextCheckbox\0update(int)\0v\0\0"
};

void Manta::TextCheckbox::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        TextCheckbox *_t = static_cast<TextCheckbox *>(_o);
        switch (_id) {
        case 0: _t->update((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::TextCheckbox::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::TextCheckbox::staticMetaObject = {
    { &QCheckBox::staticMetaObject, qt_meta_stringdata_Manta__TextCheckbox,
      qt_meta_data_Manta__TextCheckbox, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::TextCheckbox::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::TextCheckbox::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::TextCheckbox::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__TextCheckbox))
        return static_cast<void*>(const_cast< TextCheckbox*>(this));
    return QCheckBox::qt_metacast(_clname);
}

int Manta::TextCheckbox::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QCheckBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_Manta__TextSlider[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   30,   32,   32, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Manta__TextSlider[] = {
    "Manta::TextSlider\0update(int)\0v\0\0"
};

void Manta::TextSlider::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        TextSlider *_t = static_cast<TextSlider *>(_o);
        switch (_id) {
        case 0: _t->update((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Manta::TextSlider::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Manta::TextSlider::staticMetaObject = {
    { &QSlider::staticMetaObject, qt_meta_stringdata_Manta__TextSlider,
      qt_meta_data_Manta__TextSlider, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Manta::TextSlider::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Manta::TextSlider::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Manta::TextSlider::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Manta__TextSlider))
        return static_cast<void*>(const_cast< TextSlider*>(this));
    return QSlider::qt_metacast(_clname);
}

int Manta::TextSlider::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QSlider::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
