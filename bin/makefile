# MakeFile

CC = g++ 
LIBFLAGS = `root-config --libs` -O3 -lXMLIO -lMLP -lMinuit

LIBDIR = ../lib/
SDIR = ../selection/
TDIR = ../tools/
THREADDIR = ../threadpool/
OBJDIR = ../lib/analyzer_objects/

CFLAGS = -I${OBJDIR} -I${LIBDIR} -I${SDIR} -I${TDIR} -I${THREADDIR} -c -g `root-config --cflags`

# The name of the .cxx file you want to compile into an executable
MAIN = example
#MAIN = categorize

MAINRULES1 = ${LIBDIR}Sample.o ${SDIR}SelectionCuts.o ${SDIR}CategorySelection.o ${SDIR}JetSelectionTools.o ${SDIR}MuonSelectionTools.o 
MAINRULES2 = ${SDIR}ElectronSelectionTools.o ${SDIR}TauSelectionTools.o ${LIBDIR}DiMuPlottingSystem.o ${TDIR}EventTools.o 
MAINRULES3 = ${TDIR}PUTools.o ${TDIR}ParticleTools.o libAnalysisObjects.so ${MAIN}.out 
MAINDEPS   =  ${THREADDIR}ThreadPool.hxx 
DEPS       = ${LIBDIR}Cut.h ${LIBDIR}CutSet.hxx SignificanceMetrics.hxx

OBJH := $(shell find $(OBJDIR) -maxdepth 1 -name '*.h') 
OBJCC := $(shell find $(OBJDIR) -maxdepth 1 -name '*.cc') 

#OBJS = 

# execute all of the rules to make the main program
# $@ is the left of the :, and $^ is everything on the right
${MAIN}: $(MAINRULES1) $(MAINRULES2) $(MAINRULES3)
	$(CC) $(LIBFLAGS) -o $@ $^

# define a rule for each .o rule
# %.o  all files ending in .o, % is the filename without the suffix
# $< is the first item on the right side of the :, and $@ is the left side of the :
%.o: %.cxx %.h $(DEPS)
	$(CC) $(CFLAGS) $< -o $@ 

%.out: %.cxx $(MAINDEPS)
	$(CC) $(CFLAGS) $< -o $@ 

# need this so that ROOT TTrees know how to pull out our analysis objects
MyDict.cxx: ${OBJH} ${OBJDIR}linkdef/Linkdef.h
	rootcint -f $@ -c $(CFLAGS) -p $^

libAnalysisObjects.so: MyDict.cxx ${OBJCC}
	g++ -shared -o$@ `root-config --cflags` -fPIC -I$(ROOTSYS)/include $^

clean:
	rm MyDict*
	rm *.out
	rm *.so

cleanall:
	rm ${LIBDIR}*.o
	rm ${SDIR}*.o
	rm ${TDIR}*.o
	rm MyDict*
	rm *.out
	rm *.so
