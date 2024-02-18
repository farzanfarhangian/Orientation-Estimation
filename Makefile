CC = gcc
INCLUDES = -Iinclude
CFLAGS = -Wall -Wpedantic -std=c11 

SOURCES = \
    main.c \
    src/vn/conv.c \
    src/vn/error.c \
    src/vn/error_detection.c \
    src/vn/sensors.c \
    src/vn/util.c \
    src/vn/math/matrix.c \
    src/vn/math/vector.c \
    src/vn/protocol/spi.c \
    src/vn/protocol/upack.c \
    src/vn/protocol/upackf.c \
    src/vn/sensors/compositedata.c \
    src/vn/sensors/ezasyncdata.c \
    src/vn/sensors/searcher.c \
    src/vn/xplat/criticalsection.c \
    src/vn/xplat/event.c \
    src/vn/xplat/serialport.c \
    src/vn/xplat/thread.c \
    src/vn/xplat/time.c
				
# Set the object file names, with the source directory stripped
# from the path, and the build path prepended in its place			
# Object files
OBJECTS = $(SOURCES:.c=.o)

# Executable name
EXECUTABLE = my_prog

# Rule to compile .c files to .o files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
	
# Main target
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJECTS) -o $@

# Clean rule
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
