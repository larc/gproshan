INCLUDE_PATH		= -I./include -isystem /usr/include/eigen3 -I/usr/include/suitesparse 
LIBRARY_PATH		= -L/usr/local/cuda/lib64
BLAS_LIBS			= -lumfpack
SUITESPARSE_LIBS	= -lamd -lcholmod -lsuitesparseconfig -lm
OPENGL_LIBS			= -lglut -lGL -lGLU

TARGET = gproshan

# SINGLE_P = -DSINGLE_P to compile with single precision
SINGLE_P = 

CC = g++
LD = g++ -no-pie
CUDA = nvcc -ccbin g++-7
CFLAGS = -O3 -fopenmp $(INCLUDE_PATH) 
CUDAFLAGS = --generate-code arch=compute_50,code=sm_50 \
			--generate-code arch=compute_60,code=sm_60 \
			-I./include/cuda -O3 -Xcompiler -fopenmp -D_FORCE_INLINES
LFLAGS = -O3 -fopenmp $(LIBRARY_PATH) -lcublas -lcusolver -lcusparse -lcuda -lcudart -lX11 -lpthread
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
	$(LD) $(SINGLE_P) obj/$(TARGET).o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o $(TARGET) $(CFLAGS) $(LFLAGS) $(LIBS)

test_geodesics: obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o
	$(LD) $(SINGLE_P) obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) obj/link_cuda.o -o test_geodesics $(CFLAGS) $(LFLAGS) $(LIBS)

obj/$(TARGET).o: $(TARGET).cpp | obj
	$(CC) $(SINGLE_P) -c $< -o $@ $(CFLAGS) 

obj/test_geodesics.o: test_geodesics.cpp | obj
	$(CC) $(SINGLE_P) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/%.cpp | obj
	$(CC) $(SINGLE_P) -c $< -o $@ $(CFLAGS) 

obj/%.o: src/viewer/%.cpp | obj
	$(CC) $(SINGLE_P) -c $< -o $@ -I./include/viewer $(CFLAGS) 

obj/%.o: src/mdict/%.cpp | obj
	$(CC) $(SINGLE_P) -c $< -o $@ -I./include/mdict $(CFLAGS) 

obj/%_cuda.o: src/cuda/%.cu | obj
	$(CUDA) $(SINGLE_P) -dc $< -o $@ -I./include $(CUDAFLAGS)

obj/link_cuda.o: $(CUDA_OBJECTS) | obj
	$(CUDA) $(SINGLE_P) -dlink $(CUDA_OBJECTS) -o obj/link_cuda.o $(CUDAFLAGS)

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

