#-------------------------------------------------
#
# Project created by QtCreator 2015-01-21T19:03:11
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = pgm_som_project
TEMPLATE = app

DEFINES += HAVE_MPI



#LIBS += "C:/Program Files/Microsoft SDKs/MPI/Lib/x86/msmpi.lib"

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
    qsomthread.cpp

HEADERS  += mainwindow.h \   
    somoclu/src/Windows/getopt.h \
    somoclu/src/Windows/unistd.h \   
    somoclu/src/somocluWrap.h \
    ui_mainwindow.h \
    qsomthread.h

FORMS    += mainwindow.ui

INCLUDEPATH += "C:/Program Files/MPICH2/include"
INCLUDEPATH += "./somoclu/src"
LIBS += "C:/Program Files/MPICH2/lib/mpi.lib"
