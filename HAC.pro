
TEMPLATE = lib
TARGET = HAC
QMAKE_CXXFLAGS += -std=c++11

!include($$PWD/../../global.pri) {
        WORKDIR = $$PWD/../../Ornament/
        message(File global.pri is not found)
}

BIN_DIR = $$WORKDIR/MinGW-build/bin/
DESTDIR = $$BIN_DIR/
DEPENDPATH += $$BIN_DIR/

message(The project HAC will be installed in $$DESTDIR)
DEFINES +=  HAC_DLL HAC_DLL_EXPORT _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS
#WIN32

HEADERS += headers/bands.h \
           headers/beam.h \
           headers/bellhop.h \
           headers/bellhop_hydrology.h \
           headers/brdy.h \
           headers/clutter.h \
           headers/coord_2d.h \
           headers/field_clutter.h \
           headers/field_gain.h \
           headers/field_reverb.h \
           headers/main.h \
           headers/ray.h \
           headers/reverb.h \
           headers/ssp.h \

SOURCES += src/bands.cpp \
           src/beam.cpp \
           src/bellhop.cpp \
           src/bellhop_hydrology.cpp \
           src/brdy.cpp \
           src/clutter.cpp \
           src/coord_2d.cpp \
           src/field_clutter.cpp \
           src/field_gain.cpp \
           src/field_reverb.cpp \
           src/reverb.cpp \
           src/ssp.cpp

#INCLUDEPATH += .

LIBS += -L$$BIN_DIR/ \
        -lGISDataBase \
