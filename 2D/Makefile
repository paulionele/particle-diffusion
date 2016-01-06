# Compiler and options.
COMPILER = gcc
CFLAGS = -std=gnu99
LFLAGS = -lm
EXECUTABLE = 2d_inhomogenous
OBJECTS = 2d_inhomogenous.o
SOURCE = 2d_inhomogenous.c

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
