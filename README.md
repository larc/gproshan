## [gproshan](https://github.com/larc/gproshan): a geometry processing and shape analysis framework

[![Build](https://github.com/larc/gproshan_dev/actions/workflows/build.yml/badge.svg?branch=dev)](https://github.com/larc/gproshan_dev/actions/workflows/build.yml)
[![DOI](https://zenodo.org/badge/88686093.svg)](https://zenodo.org/badge/latestdoi/88686093)


![](https://raw.githubusercontent.com/larc/gproshan/master/gproshan.png)


This framework integrates some algorithms and contributions focus on the areas of computer graphics, geometry processing and computational geometry.


## Build and Run


Install all dependencies and run:

	mkdir build
	cd build
	cmake ..
	make

finally execute:

	./gproshan [mesh_paths.(off,obj,ply)]

### Dependencies (Linux)
g++ >= 9.3, cuda >= 11.0, cmake >= 3.18, armadillo, eigen, cgal, suitesparse, openblas, glew, glfw3, glm, cimg, gnuplot, embree >= 3.13

In Ubuntu you can install them with:

	sudo apt install cmake libarmadillo-dev libcgal-dev libsuitesparse-dev libopenblas-dev libglew-dev libglfw3-dev libglm-dev cimg-dev gnuplot

#### Installing Intel Embree

Intel Embree is a collection of high performance ray tracing kernels that helps graphics application engineers to improve the performance of their photorealistic rendering application ([https://www.embree.org/](https://www.embree.org/)).

##### Ubuntu (Linux)

	# download the key to system keyring
	wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
	| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

	# add signed entry to apt sources and configure the APT client to use Intel repository:
	echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

	sudo apt update
	sudo apt install intel-renderkit

	# configure the enviroment variables
	source /opt/intel/oneapi/setvars.sh

##### MacOS

	brew install embree


## Using gproshan in your project

To use gproshan as library in your project you need to execute the command to install it:

	sudo make install

then you will be able to find it in your cmake project:

	```cmake
	find_package(gproshan REQUIRED)
	include_directories(SYSTEM ${gproshan_INCLUDE_DIR})
	target_link_libraries(your_target gproshan)
	```

An example of a code using gproshan:

	```cpp
	#include <gproshan/app_viewer.h>

	int main(int nargs, const char ** args)
	{
		gproshan::app_viewer app;
		return app.main(nargs, args);
	}
	```


## Contributions

### CHE implementation
We have implemented a [Compact Half-Edge (CHE)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.523.7580) data structure to manipulated triangular meshes, also can be extended for other polygonal meshes. See the paper: [CHE: A scalable topological data structure for triangular meshes](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.523.7580) for more details.

### Geodesics
We proposed a CPU/GPU parallel algorithm to compute geodesics distances on triangular meshes. Our approach is competitive with the current methods and is simple to implement. Please cite our paper:

[A minimalistic approach for fast computation of geodesic distances on triangular meshes](https://doi.org/10.1016/j.cag.2019.08.014)

```bibtex
@Article{RFM19,
  author       = { {Romero Calla}, Luciano A. and {Fuentes Perez}, Lizeth J. and Montenegro, Anselmo A. },
  title        = { A minimalistic approach for fast computation of geodesic distances on triangular meshes },
  issn         = { 0097-8493 },
  year         = { 2019 },
  doi          = { 10.1016/j.cag.2019.08.014 },
  journaltitle = { Computers \& Graphics }
}
```

Also, we have implemented the [Fast Marching algorithm](), and the [Heat method](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/index.html).

### Ray Tracing
We have implemented a ray tracing viewer option for mesh and pointcloud visualization using [Intel Embree](https://www.embree.org/) v3.1x and [Nvidia Optix](https://developer.nvidia.com/optix) v7.1.

### Dictionary Learning
We proposed a Dictionary Learning and Sparse Coding framework, to solve the problems of Denoising, Inpainting, and Multiresolution on triangular meshes. This work is still in process. Please cite our work:

[A Robust Feature-aware Sparse Mesh Representation](https://diglib.eg.org/handle/10.2312/pg20201226)

```bibtex
@InProceedings{FRMMP20,
  booktitle    = { Pacific Graphics Short Papers, Posters, and Work-in-Progress Papers},
  title        = { {A Robust Feature-aware Sparse Mesh Representation} },
  author       = { {Fuentes Perez}, Lizeth J.  and {Romero Calla}, Luciano A. and Montenegro, Anselmo A. and Mura, Claudio and Pajarola, Renato },
  year         = { 2020 },
  publisher    = { The Eurographics Association },
  ISBN         = { 978-3-03868-120-5 },
  DOI          = { 10.2312/pg.20201226 }
}
```

### Hole repairing
We implemented repairing mesh holes in two steps:

1. Generate a mesh to cover the hole. We modified the algorithm presented in the paper: [A robust hole-filling algorithm for triangular mesh](https://doi.org/10.1007/s00371-007-0167-y), in order to
generate a planar triangular mesh using a priority queue.
2. Fit the surface described by the new points in order to minimize the variation of the surface, solving the Poisson equation (see the Chapter 4 of the book [Polygon Mesh Processing](http://www.pmp-book.org/)) or using Biharmonic splines.

Please see and cite our final undergraduate project: [mesh hole repairing report](http://repositorio.unsa.edu.pe/handle/UNSA/2576) (in Spanish).

### Key-Points and Key-Components

We proposed a simple method based on the faces' areas to compute key-points for adaptive meshes.

Please cite our paper (in Spanish):

[Efficient approach for interest points detection in non-rigid shapes](https://doi.org/10.1109/CLEI.2015.7359459)

```bibtex
@InProceedings{LRF15,
  author    = { {Lopez Del Alamo}, Cristian J. and {Romero Calla}, Luciano A. and {Fuentes Perez}, Lizeth J. },
  title     = { Efficient approach for interest points detection in non-rigid shapes },
  booktitle = { Latin American Computing Conference (CLEI) },
  pages     = { 1-8 },
  year      = { 2015 },
  doi       = { 10.1109/CLEI.2015.7359459 }
}
```

Computing key-components depends on the accuracy and definition of the key points. We were inspired by the [work of Ivan Sipiran](https://www.researchgate.net/publication/262350194_Key-component_detection_on_3D_meshes_using_local_features), he defined for the first time the notion of a key-component in meshes. We proposed a method based on geodesics to determine the key components.

Please see and cite our final undergraduate project: [key-components report](http://repositorio.unsa.edu.pe/handle/UNSA/2575) (in Spanish).

### Decimation
We are implementing the algorithm described by the paper [Stellar Mesh Simplification Using Probabilistic Optimization](https://doi.org/10.1111/j.1467-8659.2004.00811.x), to compute a mesh simplification.

### Fairing
We implemented Spectral and Taubin fairing algorithms to smooth a mesh surface. See the Chapter 4 of the book [Polygon Mesh Processing](http://www.pmp-book.org/).

### Laplacian and signatures
Laplace-Beltrami operator and its eigen decomposition, WKS, HKS, GPS signatures.


## Documentation
Execute:

	doxygen Doxyfile

to generate the documentation in *html* and *LaTeX*.


## Viewer
The viewer is done with modern OpenGL and a GUI using Dear ImGui [https://github.com/ocornut/imgui](https://github.com/ocornut/imgui). The viewer was initially based in the viewer of [https://github.com/dgpdec/course](https://github.com/dgpdec/course).


## License

MIT License


## Authors/Contributors
- [Luciano A. Romero Calla](https://github.com/larc)
- [Lizeth J. Fuentes Pérez](https://github.com/lizonly)
