INCLUDE_PATH		= -I./include -isystem /usr/include/eigen3 
LIBRARY_PATH		=
BLAS_LIBS			= -lumfpack
SUITESPARSE_LIBS	= -lspqr -lcholmod
OPENGL_LIBS			= -lglut -lGL -lGLU

TARGET = gproshan
CC = g++
LD = g++
CUDA = nvcc
CFLAGS = -O3 -fopenmp $(INCLUDE_PATH) 
CUDAFLAGS = -I./include/cuda -O3 -Xcompiler -fopenmp -D_FORCE_INLINES
LFLAGS = -O3 -fopenmp $(LIBRARY_PATH) -lcublas -lcuda -lcudart -lX11 -lpthread
LIBS = $(OPENGL_LIBS) $(SUITESPARSE_LIBS) $(BLAS_LIBS) -larmadillo -lCGAL

########################################################################################
## !! Do not edit below this line

HEADERS := $(wildcard include/*.h)
SOURCES := $(wildcard src/*.cpp) $(wildcard src/viewer/*.cpp) $(wildcard src/mdict/*.cpp)
CUDA_HEADERS := $(wildcard include/cuda/*.cuh)
CUDA_SOURCES := $(wildcard src/cuda/*.cu)

OBJECTS :=	$(addprefix obj/,$(notdir $(SOURCES:.cpp=.o)))
CUDA_OBJECTS :=	$(addprefix obj/,$(notdir $(CUDA_SOURCES:.cu=_cuda.o)))

all: $(TARGET) | tmp

$(TARGET): $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o
	$(LD) $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o $(TARGET) $(CFLAGS) $(LFLAGS) $(LIBS)

obj/%.o: src/%.cpp | obj
	$(CC) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/viewer/%.cpp | obj
	$(CC) -c $< -o $@ -I./include/viewer $(CFLAGS) 

obj/%.o: src/mdict/%.cpp | obj
	$(CC) -c $< -o $@ -I./include/mdict $(CFLAGS) 

obj/%_cuda.o: src/cuda/%.cu | obj
	$(CUDA) -dc $< -o $@ $(CUDAFLAGS)

obj/link_cuda.o: $(CUDA_OBJECTS) | obj
	$(CUDA) -dlink $(CUDA_OBJECTS) -o obj/link_cuda.o $(CUDAFLAGS)

obj:
	mkdir obj

tmp:
	mkdir tmp

clean:
	rm -f $(OBJECTS) $(CUDA_OBJECTS)
	rm -f $(TARGET)

