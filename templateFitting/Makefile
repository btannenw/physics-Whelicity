#############################################################################
#									    #	
# Code setup by Philipp Stolte 
#									    #	
#############################################################################

# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

# Programs
CXX          = g++ -pg -ggdb
CXXFLAGS     = -Wall 
LD           = 
LDFLAGS      = -g -O2 -pg
SOFLAGS      = -shared

RM           = rm -f 
MV           = mv 
ECHO         = echo
CINT         = rootcint

# Assign or Add variables
CXXFLAGS    = $(ROOTCFLAGS) 
CXXFLAGS   += -I./include/
CXXFLAGS   += -I. -I./include -I$(BATINSTALLDIR)/include 
LIBS        = $(ROOTLIBS)  -L$(BATINSTALLDIR)/lib -lBAT -lBATmodels 
GLIBS       = $(ROOTGLIBS) -L$(BATINSTALLDIR)/lib -lBAT -lBATmodels 
#CXXFLAGS   += -I. -I./include -I/home/aknue/GoeSelection_11_branch_0410/BAT-0.9.1/include 
#LIBS        = $(ROOTLIBS)  -L/home/aknue/GoeSelection_11_branch_0410/BAT-0.9.1/lib -lBAT -lBATmodels 
#GLIBS       = $(ROOTGLIBS) -L/home/aknue/GoeSelection_11_branch_0410/BAT-0.9.1/lib -lBAT -lBATmodels 

CXSRCS1      = GoeProfiling.c
CXSRCS1b     = ProducePlots.c
CXSRCS1c     = ProduceCalibration.c
CXSRCS1d     = EvaluateExternal.c
CXSRCS1e     = EvaluateExternalUpdate.c
CXSRCS1f     = ProducePlotsCaliExtern.c
CXSRCS2      = src/StatusLogbook.cxx src/AtlasStyle.cxx src/PlotInterpolationCurves.cxx src/ExternalSystematics.cxx src/ValidationClass.cxx src/ProfilingClass.cxx src/interpolation.cxx src/fit_validation.cxx src/helicity.cxx src/plots.cxx src/functions.cxx src/fit.cxx
CXSRCS2b     = src/AtlasStyle.cxx src/helicity.cxx src/plots.cxx src/functions.cxx src/fill_histograms.cxx src/OutputTreeReader.cxx
CXSRCS2c     = src/AtlasStyle.cxx src/StatusLogbook.cxx src/plots.cxx src/CalibrationCurves.cxx
CXSRCS2d     = src/StatusLogbook.cxx src/AtlasStyle.cxx src/ValidationClass.cxx src/interpolation.cxx src/ExternalSystematics.cxx src/ProfilingClass.cxx src/PlotInterpolationCurves.cxx src/helicity.cxx src/functions.cxx src/fit.cxx src/plots.cxx 
CXSRCS2e     = src/StatusLogbook.cxx src/AtlasStyle.cxx src/ValidationClass.cxx src/interpolation.cxx src/ExternalSystematics.cxx src/ProfilingClass.cxx src/PlotInterpolationCurves.cxx src/helicity.cxx src/functions.cxx src/fit.cxx src/plots.cxx
CXSRCS2f     = src/StatusLogbook.cxx src/AtlasStyle.cxx src/helicity.cxx src/functions.cxx src/fit.cxx src/plots.cxx

CXXSRCS      = $(patsubst %.c, obj/%.o, $(CXSRCS1))
CXXSRCS1b    = $(patsubst %.c, obj/%.o, $(CXSRCS1b))
CXXSRCS1c    = $(patsubst %.c, obj/%.o, $(CXSRCS1c))
CXXSRCS1d    = $(patsubst %.c, obj/%.o, $(CXSRCS1d))
CXXSRCS1e    = $(patsubst %.c, obj/%.o, $(CXSRCS1e))
CXXSRCS1f    = $(patsubst %.c, obj/%.o, $(CXSRCS1f))

CXXSRCS     += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2))
CXXSRCS1b   += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2b))
CXXSRCS1c   += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2c))
CXXSRCS1d   += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2d))
CXXSRCS1e   += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2e))
CXXSRCS1f   += $(patsubst src/%.cxx, obj/%.o, $(CXSRCS2f))

CXXOBJS      = obj/AtlasStyle.o obj/StatusLogbook.o obj/PlotInterpolationCurves.o obj/ExternalSystematics.o obj/ProfilingClass.o obj/ValidationClass.o obj/fit_validation.o obj/helicity.o obj/plots.o obj/interpolation.o obj/functions.o obj/fit.o obj/GoeProfiling.o

CXXOBJS1b    = obj/AtlasStyle.o obj/helicity.o obj/plots.o obj/functions.o obj/fill_histograms.o obj/OutputTreeReader.o obj/ProducePlots.o

CXXOBJS1c    = obj/AtlasStyle.o obj/StatusLogbook.o obj/plots.o obj/CalibrationCurves.o obj/ProduceCalibration.o

CXXOBJS1d    = obj/AtlasStyle.o obj/StatusLogbook.o obj/ValidationClass.o obj/interpolation.o obj/PlotInterpolationCurves.o obj/ExternalSystematics.o obj/ProfilingClass.o obj/fit_validation.o obj/helicity.o obj/functions.o obj/fit.o obj/EvaluateExternal.o obj/plots.o

CXXOBJS1e    = obj/AtlasStyle.o obj/StatusLogbook.o obj/ValidationClass.o obj/interpolation.o obj/PlotInterpolationCurves.o obj/ExternalSystematics.o obj/ProfilingClass.o obj/fit_validation.o obj/helicity.o obj/functions.o obj/fit.o obj/EvaluateExternalUpdate.o obj/plots.o

CXXOBJS1f    = obj/AtlasStyle.o obj/StatusLogbook.o obj/ValidationClass.o obj/interpolation.o obj/PlotInterpolationCurves.o obj/ExternalSystematics.o obj/ProfilingClass.o obj/fit_validation.o obj/helicity.o obj/functions.o obj/fit.o obj/EvaluateExternal.o obj/plots.o

CXXOBJS1g    = obj/AtlasStyle.o obj/StatusLogbook.o obj/ValidationClass.o obj/interpolation.o obj/PlotInterpolationCurves.o obj/ExternalSystematics.o obj/ProfilingClass.o obj/fit_validation.o obj/helicity.o obj/functions.o obj/fit.o obj/EvaluateExternalPileup.o obj/plots.o


GARBAGE      = runGoeProfiling runProducePlotsCal runProducePlots runEvaluateExternal runEvaluateExternal runProducePlotsCaliExtern *~ src/*~ obj/* include/*~ *.o scripts/*~

all :main_program plots_cal plots Syst SystU SystP

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

main_program : GoeProfiling.c $(CXXOBJS) 
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS) $(ROOTLIBS) $(LIBS) $(ROOTGLIBS) -o runGoeProfiling

plots_cal : ProduceCalibration.c $(CXXOBJS1c)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1c)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1c) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runProduceCalibration

plots :	ProducePlots.c $(CXXOBJS1b)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1b)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1b) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runProducePlots

Syst  : EvaluateExternal.c $(CXXOBJS1d)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1d)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1d) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runEvaluateExternal

SystU : EvaluateExternalUpdate.c $(CXXOBJS1e)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1e)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1e) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runEvaluateExternalUpdate

SystP : EvaluateExternalPileup.c $(CXXOBJS1g)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1e)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1g) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runEvaluateExternalPileup

plots_cal_ext : ProducePlotsCaliExtern.c $(CXXOBJS1b)
	$(CXX) $(CXXFLAGS) -c $< $(CXSRC1b)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(LDFLAGS) $(CXXOBJS1b) $(ROOTGLIBS) $(LIBS) $(ROOTLIBS) -o runProducePlotsCaliExtern


clean :
	$(RM) $(GARBAGE)   

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ flags : $(CXXFLAGS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)
	echo rootglibs : $(ROOTGLIBS)


