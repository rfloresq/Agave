TEMPLATE = app
TARGET = agave
DEPENDPATH += .
INCLUDEPATH += .
LIBS += -lGLU -llapacke -llapack

# Input
QT += opengl 
HEADERS += Agave.h \
           Atom.h \
           Array.h \
           Basis.h \
           Curve.h \
           Element.h \
           Evaluator.h \
           GAPP.h \
           Global.h \
           Grid.h \
	   Integrator.h \
           Math.h \
           Matrix.h \
           Model.h \
           Parameter.h \
           Plotter.h \
           QM.h \
	   SCF.h \
	   Set.h \
	   Shell.h \
	   Species.h \
	   Structure.h \
	   Surface.h \
           System.h \
           Units.h \
           Vector.h \
           Viewer.h \
           WinElement.h \
           WinGeo.h \
           WinPlot.h \
           WinProp.h \
           WinQM.h \
           WinSCF.h \
           WinSpecies.h
SOURCES += Agave.cpp \
           Atom.cpp \
           Array.cpp \
           Basis.cpp \
           Curve.cpp \
           Element.cpp \
           Evaluator.cpp \
           GAPP.cpp \
           Grid.cpp \
	   Integrator.cpp \
           Math.cpp \
           Matrix.cpp \
           Model.cpp \
           main.cpp \
           Plotter.cpp \
           QM.cpp \
	   SCF.cpp \
	   Set.cpp \
	   Shell.cpp \
	   Species.cpp \
	   Structure.cpp \
	   Surface.cpp \
           System.cpp \
           Units.cpp \
           Vector.cpp \
           Viewer.cpp \
           WinElement.cpp \
           WinGeo.cpp \
           WinPlot.cpp \
           WinProp.cpp \
           WinQM.cpp \
           WinSCF.cpp \
           WinSpecies.cpp
