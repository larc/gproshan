INCLUDE_PATH		= -I./include -isystem /usr/include/eigen3 -I/usr/include/suitesparse 
LIBRARY_PATH		= -L/usr/local/cuda/lib64
BLAS_LIBS			= #-lumfpack
SUITESPARSE_LIBS	= #-lspqr -lcholmod
OPENGL_LIBS			= -lglut -lGL -lGLU

TARGET = gproshan

CC = g++
LD = g++ -no-pie
CUDA = nvcc
CFLAGS = -O3 -fopenmp $(INCLUDE_PATH) 
CUDAFLAGS = -I./include/cuda -O3 -Xcompiler -fopenmp -D_FORCE_INLINES
LFLAGS = -O3 -fopenmp $(LIBRARY_PATH) -lcublas -lcuda -lcudart -lX11 -lpthread
LIBS = $(OPENGL_LIBS) $(SUITESPARSE_LIBS) $(BLAS_LIBS) -larmadillo -lsuperlu -lCGAL

########################################################################################
## !! Do not edit below this line

HEADERS := $(wildcard include/*.h)
SOURCES := $(wildcard src/*.cpp) $(wildcard src/viewer/*.cpp) $(wildcard src/mdict/*.cpp)
CUDA_HEADERS := $(wildcard include/cuda/*.cuh)
CUDA_SOURCES := $(wildcard src/cuda/*.cu)

OBJECTS :=	$(addprefix obj/,$(notdir $(SOURCES:.cpp=.o)))
CUDA_OBJECTS :=	$(addprefix obj/,$(notdir $(CUDA_SOURCES:.cu=_cuda.o)))

all: $(TARGET) | tmp

$(TARGET): obj/$(TARGET).o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o
	$(LD) obj/$(TARGET).o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o $(TARGET) $(CFLAGS) $(LFLAGS) $(LIBS)

test_geodesics: obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o
	$(LD) obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o test_geodesics $(CFLAGS) $(LFLAGS) $(LIBS)

obj/$(TARGET).o: $(TARGET).cpp | obj
	$(CC) -c $< -o $@ $(CFLAGS) 

obj/test_geodesics.o: test_geodesics.cpp | obj
	$(CC) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/%.cpp | obj
	$(CC) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/viewer/%.cpp | obj
	$(CC) -c $< -o $@ -I./include/viewer $(CFLAGS) 

obj/%.o: src/mdict/%.cpp | obj
	$(CC) -c $< -o $@ -I./include/mdict $(CFLAGS) 

obj/%_cuda.o: src/cuda/%.cu | obj
	$(CUDA) -dc $< -o $@ -I./include $(CUDAFLAGS)

obj/link_cuda.o: $(CUDA_OBJECTS) | obj
	$(CUDA) -dlink $(CUDA_OBJECTS) -o obj/link_cuda.o $(CUDAFLAGS)

obj:
	mkdir obj

tmp:
	mkdir tmp

library: $(filter-out obj/main.o, $(OBJECTS) $(CUDA_OBJECTS)) obj/link_cuda.o | lib
	ar -cvq lib/lib$(TARGET).a $(filter-out obj/main.o, $(OBJECTS) $(CUDA_OBJECTS)) obj/link_cuda.o

lib:
	mkdir lib

clean:
	rm -f $(OBJECTS) $(CUDA_OBJECTS)
	rm -f $(TARGET) test_geodesics

