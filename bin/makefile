# MakeFile

# Set the Compiler
CC = g++

LIBFLAGS = `root-config --libs` -O3 -lXMLIO -lMLP -lMinuit
LIBDIR = ../lib/
CFLAGS = -I${LIBDIR} -c -g `root-config --cflags`
MAIN = synch

synch: ${LIBDIR}Sample.o ${LIBDIR}DiMuPlottingSystem.o ${MAIN}.o ${LIBDIR}Cut.h ${LIBDIR}CutSet.h
	$(CC) $(LIBFLAGS) ${LIBDIR}Sample.o ${LIBDIR}DiMuPlottingSystem.o ${MAIN}.o -o ${MAIN}

${MAIN}.o: ${MAIN}.cxx 
	$(CC) $(CFLAGS) ${MAIN}.cxx -o ${MAIN}.o

${LIBDIR}Sample.o: ${LIBDIR}Sample.cxx ${LIBDIR}Sample.h ${LIBDIR}VarSet.h ${LIBDIR}DataFormats.h
	$(CC) $(CFLAGS) ${LIBDIR}Sample.cxx -o ${LIBDIR}Sample.o


${LIBDIR}DiMuPlottingSystem.o: ${LIBDIR}DiMuPlottingSystem.cxx ${LIBDIR}DiMuPlottingSystem.h ${LIBDIR}Sample.h
	$(CC) $(CFLAGS) ${LIBDIR}DiMuPlottingSystem.cxx -o ${LIBDIR}DiMuPlottingSystem.o

clean:
	rm *.o

cleanall:
	rm *.o
	rm ${LIBDIR}*.o
