#-------------------------------------------------
#
# Project created by QtCreator 2015-01-21T19:03:11
#
#-------------------------------------------------

QT       += core gui

QT += widgets

TARGET = pgm_som_project
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    somoclu/src/denseCpuKernels.cpp \
    somoclu/src/io.cpp \
    somoclu/src/mapDistanceFunctions.cpp \   
    somoclu/src/somocluWrap.cpp \
    somoclu/src/sparseCpuKernels.cpp \
    somoclu/src/training.cpp \
    somoclu/src/trainOneEpoch.cpp \
    somoclu/src/uMatrix.cpp \
    somoclu/src/Windows/getopt.c \
    qsomthread.cpp \
    cgetdata.cpp \    
    npmpgmthread.cpp \
    realclsnpmpgmthread.cpp

HEADERS  += mainwindow.h \   
    somoclu/src/Windows/getopt.h \
    somoclu/src/somocluWrap.h \   
    ui_mainwindow.h \
    qsomthread.h \  
    cgetdata.h \   
    npmpgmthread.h \
    realclsnpmpgmthread.h

FORMS    += mainwindow.ui

INCLUDEPATH += "./somoclu/src"
INCLUDEPATH += "./armadillo/include"

