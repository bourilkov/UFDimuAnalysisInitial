# MakeFile

# Set the Compiler
CC = g++

LIBFLAGS = `root-config --libs` -O3 -lXMLIO -lMLP -lMinuit
LIBDIR = ../lib/
SDIR = ../selection/
CFLAGS = -I${LIBDIR} -I${SDIR} -c -g `root-config --cflags`
MAIN = synch

synch: ${LIBDIR}Sample.o ${SDIR}SelectionCuts.o ${SDIR}CategorySelection.o ${SDIR}JetSelectionTools.o ${LIBDIR}DiMuPlottingSystem.o ${MAIN}.o ${LIBDIR}Cut.h ${SDIR}SelectionCuts.h ${SDIR}JetSelectionTools.h
	$(CC) $(LIBFLAGS) ${LIBDIR}Sample.o ${SDIR}SelectionCuts.o ${SDIR}CategorySelection.o ${SDIR}JetSelectionTools.o ${LIBDIR}DiMuPlottingSystem.o ${MAIN}.o -o ${MAIN}

${MAIN}.o: ${MAIN}.cxx 
	$(CC) $(CFLAGS) ${MAIN}.cxx -o ${MAIN}.o

${SDIR}SelectionCuts.o: ${SDIR}SelectionCuts.cxx ${SDIR}SelectionCuts.h ${LIBDIR}Cut.h ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${SDIR}SelectionCuts.cxx -o ${SDIR}SelectionCuts.o

${SDIR}CategorySelection.o: ${SDIR}CategorySelection.cxx ${SDIR}CategorySelection.h ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${SDIR}CategorySelection.cxx -o ${SDIR}CategorySelection.o

${SDIR}JetSelectionTools.o: ${SDIR}JetSelectionTools.cxx ${SDIR}JetSelectionTools.h ${LIBDIR}VarSet.h
	$(CC) $(CFLAGS) ${SDIR}JetSelectionTools.cxx -o ${SDIR}JetSelectionTools.o

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
