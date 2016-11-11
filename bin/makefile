# MakeFile

# Set the Compiler
CC = g++

# Include the root libraries
LIBFLAGS = `root-config --libs` -O3 -lXMLIO -lMLP -lMinuit

# the directory for the base library that enables the analysis
LIBDIR = ../lib/

# the directory for the cut selection objects specific to the analysis
SDIR = ../selection/

# the directory for the cut tool objects specific to the analysis
TDIR = ../tools/

# include the libdir and sdir and the root libraries
CFLAGS = -I${LIBDIR} -I${SDIR} -I${TDIR} -c -g `root-config --cflags`

# The name of the .cxx file you want to compile into an executable
MAIN = categorizeRun1
#MAIN = categorizeRun2
#MAIN = fewzCategorize
#MAIN = fewzOverlay
#MAIN = fewzDYMCGenToRecoSmearing

# link all of the objects to the main program
${MAIN}: ${LIBDIR}Sample.o ${SDIR}SelectionCuts.o ${SDIR}CategorySelection.o ${SDIR}JetSelectionTools.o ${SDIR}MuSelectionTools.o ${SDIR}ElectronSelectionTools.o ${SDIR}TauSelectionTools.o ${LIBDIR}DiMuPlottingSystem.o ${TDIR}EventTools.o ${TDIR}PUTools.o ${TDIR}ParticleTools.o ${MAIN}.o ${LIBDIR}Cut.h ${LIBDIR}CutSet.h
	$(CC) $(LIBFLAGS) ${LIBDIR}Sample.o ${TDIR}EventTools.o ${TDIR}PUTools.o ${TDIR}ParticleTools.o ${SDIR}SelectionCuts.o ${SDIR}CategorySelection.o ${SDIR}JetSelectionTools.o ${SDIR}MuSelectionTools.o ${SDIR}ElectronSelectionTools.o ${SDIR}TauSelectionTools.o ${LIBDIR}DiMuPlottingSystem.o ${MAIN}.o -o ${MAIN}

${MAIN}.o: ${MAIN}.cxx SignificanceMetrics.cxx ${LIBDIR}DataFormats.h 
	$(CC) $(CFLAGS) ${MAIN}.cxx -o ${MAIN}.o

${TDIR}EventTools.o: ${TDIR}EventTools.cxx ${TDIR}EventTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${TDIR}EventTools.cxx -o ${TDIR}EventTools.o

${TDIR}PUTools.o: ${TDIR}PUTools.cxx ${TDIR}PUTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${TDIR}PUTools.cxx -o ${TDIR}PUTools.o

${TDIR}ParticleTools.o: ${TDIR}ParticleTools.cxx ${TDIR}ParticleTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${TDIR}ParticleTools.cxx -o ${TDIR}ParticleTools.o

${SDIR}SelectionCuts.o: ${SDIR}SelectionCuts.cxx ${SDIR}SelectionCuts.h ${LIBDIR}Cut.h ${LIBDIR}CutSet.h ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}SelectionCuts.cxx -o ${SDIR}SelectionCuts.o

${SDIR}CategorySelection.o: ${SDIR}CategorySelection.cxx ${SDIR}CategorySelection.h ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}CategorySelection.cxx -o ${SDIR}CategorySelection.o

${SDIR}JetSelectionTools.o: ${SDIR}JetSelectionTools.cxx ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}JetSelectionTools.cxx -o ${SDIR}JetSelectionTools.o

${SDIR}MuSelectionTools.o: ${SDIR}MuSelectionTools.cxx ${SDIR}MuSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}MuSelectionTools.cxx -o ${SDIR}MuSelectionTools.o

${SDIR}ElectronSelectionTools.o: ${SDIR}ElectronSelectionTools.cxx ${SDIR}ElectronSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}ElectronSelectionTools.cxx -o ${SDIR}ElectronSelectionTools.o

${SDIR}TauSelectionTools.o: ${SDIR}TauSelectionTools.cxx ${SDIR}TauSelectionTools.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${SDIR}TauSelectionTools.cxx -o ${SDIR}TauSelectionTools.o

${LIBDIR}Sample.o: ${LIBDIR}Sample.cxx ${LIBDIR}Sample.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${LIBDIR}Sample.cxx -o ${LIBDIR}Sample.o

${LIBDIR}DiMuPlottingSystem.o: ${LIBDIR}DiMuPlottingSystem.cxx ${LIBDIR}DiMuPlottingSystem.h ${LIBDIR}Sample.h 
	$(CC) $(CFLAGS) ${LIBDIR}DiMuPlottingSystem.cxx -o ${LIBDIR}DiMuPlottingSystem.o


clean:
	rm *.o

cleanall:
	rm *.o
	rm ${LIBDIR}*.o
	rm ${SDIR}*.o
	rm ${TDIR}*.o
