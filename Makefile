UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	USR = /usr
endif

ifeq ($(UNAME), Darwin)
	USR = /opt/local
endif

INCLUDE_PATH		= -I./include -isystem $(USR)/include/eigen3 -I$(USR)/include/suitesparse 

ifeq ($(UNAME), Darwin)
	INCLUDE_PATH		+= -isystem $(USR)/include
endif

LIBRARY_PATH		= -L$(USR)/lib/
BLAS_LIBS			= -lumfpack
SUITESPARSE_LIBS	= -lamd -lcholmod -lsuitesparseconfig -lm
OPENGL_LIBS			= -lglut -lGL -lGLU

TARGET = gproshan

# SINGLE_P = -DSINGLE_P to compile with single precision
SINGLE_P = 

ifeq ($(UNAME), Linux)
	CC = g++
	LD = $(CC) -no-pie
endif

# macports GNU gcc/g++ compiler
ifeq ($(UNAME), Darwin)
	CC = g++-mp-8
	LD = $(CC)
endif

ifeq (nvcc, $(shell basename $(shell which nvcc)))
	LIBRARY_PATH += -L$(USR)/local/cuda/lib64
	CUDA = nvcc
	CUDA_SUPPORT = -DCUDA_SUPPORT
	CUDA_FLAGS = -arch=sm_50 \
				-gencode=arch=compute_50,code=sm_50 \
				-gencode=arch=compute_52,code=sm_52 \
				-gencode=arch=compute_60,code=sm_60 \
				-gencode=arch=compute_61,code=sm_61 \
				-gencode=arch=compute_70,code=sm_70 \
				-gencode=arch=compute_70,code=compute_70 \
				-I./include/cuda -O3 -D_FORCE_INLINES
	CUDA_LIBS = -lcublas -lcusolver -lcusparse -lcuda -lcudart
endif

CFLAGS = -O3 -fopenmp $(INCLUDE_PATH) $(CUDA_SUPPORT)
LFLAGS = -O3 $(LIBRARY_PATH) $(CUDA_LIBS) -lX11 -lpthread
LIBS = $(OPENGL_LIBS) $(SUITESPARSE_LIBS) $(BLAS_LIBS) -larmadillo -lsuperlu -lCGAL

########################################################################################
## !! Do not edit below this line

HEADERS := $(wildcard include/*.h)
SOURCES := $(wildcard src/*.cpp) $(wildcard src/viewer/*.cpp) $(wildcard src/mdict/*.cpp)
OBJECTS :=	$(addprefix obj/,$(notdir $(SOURCES:.cpp=.o)))

ifeq ($(CUDA), nvcc)
	CUDA_HEADERS := $(wildcard include/cuda/*.cuh)
	CUDA_SOURCES := $(wildcard src/cuda/*.cu)
	CUDA_OBJECTS :=	$(addprefix obj/,$(notdir $(CUDA_SOURCES:.cu=_cuda.o)))
	CUDA_LINK := obj/cuda_link.o
endif

all: $(TARGET) | tmp

$(TARGET): obj/$(TARGET).o $(OBJECTS) $(CUDA_OBJECTS) $(CUDA_LINK)
	$(LD) $(SINGLE_P) obj/$(TARGET).o $(OBJECTS) $(CUDA_OBJECTS) $(CUDA_LINK) -o $(TARGET) $(CFLAGS) $(LFLAGS) $(LIBS)

test_geodesics: obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) $(CUDA_LINK)
	$(LD) $(SINGLE_P) obj/test_geodesics.o $(OBJECTS) $(CUDA_OBJECTS) $(CUDA_LINK) -o test_geodesics $(CFLAGS) $(LFLAGS) $(LIBS)

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
	$(CUDA) $(SINGLE_P) -dc $< -o $@ -I./include/ -I./include/cuda $(CUDAFLAGS)

$(CUDA_LINK): $(CUDA_OBJECTS) | obj
	$(CUDA) $(SINGLE_P) -dlink $(CUDA_OBJECTS) -o $(CUDA_LINK) $(CUDA_FLAGS)

obj:
	mkdir obj

tmp:
	mkdir tmp

library: $(filter-out obj/main.o, $(OBJECTS) $(CUDA_OBJECTS)) $(CUDA_LINK) | lib
	ar -cvq lib/lib$(TARGET).a $(filter-out obj/main.o, $(OBJECTS) $(CUDA_OBJECTS)) $(CUDA_LINK)

lib:
	mkdir lib

clean:
	rm -f $(OBJECTS) $(CUDA_OBJECTS)
	rm -f $(TARGET) test_geodesics

