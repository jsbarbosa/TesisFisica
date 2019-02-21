include Makefile.defs

OPT+= -fPIC -DLIBREBOUND

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_ias15.c integrator_sei.c integrator_leapfrog.c integrator_mercurius.c boundary.c input.c binarydiff.c output.c collision.c communication_mpi.c display.c tools.c derivatives.c simulationarchive.c glad.c integrator_janus.c transformations.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so

%.o: %.c $(HEADERS)
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

librebound.so: $(OBJECTS)
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -o $@

clean:
	rm -f *.o
	rm -f *.so
	rm -f *.gcda
	rm -f *.gcno
