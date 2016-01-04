# Compiler and options.
COMPILER = gcc
CFLAGS = -std=gnu99
LFLAGS = -lm
EXECUTABLE = 1d_inhomogenous_physical_bcs
OBJECTS = 1d_inhomogenous_physical_bcs.o
SOURCE = 1d_inhomogenous_physical_bcs.c

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(OBJECTS) -o $(EXECUTABLE)

$(OBJECTS): $(SOURCE)
	$(COMPILER) -c $(SOURCE)

.PHONY: run clean

run:
	./$(EXECUTABLE)
clean:
	rm *.o $(EXECUTABLE)
