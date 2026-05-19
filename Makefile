CC      = gcc
INCLUDES = -Iinclude
CFLAGS  = -Wall -Wpedantic -std=c11 -O2
LDFLAGS = -lm

SOURCES = \
	main.c \
	src/algorithms/madgwick.c \
	src/algorithms/mahony.c \
	src/algorithms/ekf.c \
	src/algorithms/mathTransform.c

OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = my_prog

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
