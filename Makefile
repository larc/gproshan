INCLUDE_PATH		= -I./include -isystem /usr/include/eigen3 
LIBRARY_PATH		=
BLAS_LIBS			= -lumfpack
SUITESPARSE_LIBS	= -lspqr -lcholmod
OPENGL_LIBS			= -lglut -lGL -lGLU
CGAL_DEFINES		= -DCGAL_EIGEN3_ENABLED -DCGAL_USE_BOOST_PROGRAM_OPTIONS -DCGAL_USE_GMP -DCGAL_USE_MPFR

TARGET = tesiscode
CC = g++
LD = g++
CUDA = nvcc
CFLAGS = -O3 -fopenmp $(INCLUDE_PATH) $(CGAL_DEFINES) 
CUDAFLAGS = --gpu-architecture=sm_50 -I./include/cuda -O3 -Xcompiler -fopenmp -D_FORCE_INLINES
LFLAGS = -O3 -fopenmp $(LIBRARY_PATH) -lcublas -lcuda -lcudart -lX11 -lpthread
LIBS = $(OPENGL_LIBS) $(SUITESPARSE_LIBS) $(BLAS_LIBS) -larmadillo -lCGAL

########################################################################################
## !! Do not edit below this line

HEADERS := $(wildcard include/*.h)
SOURCES := $(wildcard src/*.cpp) $(wildcard src/viewer/*.cpp)
CUDA_HEADERS := $(wildcard include/cuda/*.cuh)
CUDA_SOURCES := $(wildcard src/cuda/*.cu)

OBJECTS :=	$(addprefix obj/,$(notdir $(SOURCES:.cpp=.o)))
CUDA_OBJECTS :=	$(addprefix obj/,$(notdir $(CUDA_SOURCES:.cu=_cuda.o)))

all: $(TARGET)

$(TARGET): $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o
	$(LD) $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o $(TARGET) $(CFLAGS) $(LFLAGS) $(LIBS)

obj/%.o: src/%.cpp
	$(CC) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/viewer/%.cpp
	$(CC) -c $< -o $@ -I./include/viewer $(CFLAGS) 

obj/%_cuda.o: src/cuda/%.cu
	$(CUDA) -dc $< -o $@ $(CUDAFLAGS)

obj/link_cuda.o: $(CUDA_OBJECTS)
	$(CUDA) -dlink $(CUDA_OBJECTS) -o obj/link_cuda.o $(CUDAFLAGS)

clean:
	rm -f $(OBJECTS) $(CUDA_OBJECTS)
	rm -f $(TARGET)

