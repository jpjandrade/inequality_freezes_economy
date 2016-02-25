CC = gcc
FLAGS = -lgsl -lgslcblas -lm -O2
DFLAGS = -lgsl -lgslcblas -lm -O0 -Wall -g
LFLAGS = -m64

all: gran_econ

gran_econ : gran_econ.c gran_econ.h
	$(CC) -o gran_econ gran_econ.c $(FLAGS) $(LFLAGS)


debug : gran_econ.c gran_econ.h
	$(CC) -o gran_econ gran_econ.c $(DFLAGS) $(LFLAGS)


clean :
	rm $(OBJECTS)
