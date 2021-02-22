## [gproshan](https://github.com/larc/gproshan): a geometry processing and shape analysis framework 

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
g++ >= 8.3, cuda >= 10.1, libarmadillo, libeigen, libsuitesparse, libopenblas, opengl, glew, gnuplot, libcgal, libgles2-mesa, cimg

In Ubuntu (>= 18.04) you can install them with:

	sudo apt install libarmadillo-dev libeigen3-dev libopenblas-dev libsuitesparse-dev libglew-dev freeglut3-dev libgles2-mesa-dev cimg-dev libcgal-dev

#### Build Status

##### Ubuntu 18.04 (Bionic)
[![Build Status](https://travis-ci.com/larc/gproshan.svg?branch=main)](https://travis-ci.com/larc/gproshan)

## Contributions

### CHE implementation
We have implemented a [Compact Half-Edge (CHE)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.523.7580) data structure to manipulated triangular meshes, also can be extended for other polygonal meshes.
See the paper: [CHE: A scalable topological data structure for triangular meshes](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.523.7580) for more details.

### Geodesics
We proposed a CPU/GPU parallel algorithm to compute geodesics distances on triangular meshes. Our
approach is competitive with the current methods and is simple to implement. Please cite our paper:

[A minimalistic approach for fast computation of geodesic distances on triangular meshes](https://doi.org/10.1016/j.cag.2019.08.014)

```bibtex
@Article{ROMEROCALLA2019,
  author   = { {Romero Calla}, Luciano A. and {Fuentes Perez}, Lizeth J. and Montenegro, Anselmo A. },
  title    = { A minimalistic approach for fast computation of geodesic distances on triangular meshes },
  journal  = { Computers \& Graphics },
  year     = { 2019 },
  issn     = { 0097-8493 },
  doi      = { https://doi.org/10.1016/j.cag.2019.08.014 },
  keywords = { Geodesic distance, Fast marching, Triangular meshes, Parallel programming, Breadth-first search },
  url      = { http://www.sciencedirect.com/science/article/pii/S0097849319301426 }
}
```

Also, we have implemented the [Fast Marching algorithm](), and the [Heat method](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/index.html).

### Dictionary Learning
We proposed a Dictionary Learning and Sparse Coding framework, to solve the problems of Denoising,
Inpainting, and Multiresolution on triangular meshes. This work is still in process. Please cite
our work:

[A Dictionary Learning-based framework on Triangular Meshes](https://arxiv.org/abs/1810.08266)

```bibtex

@ARTICLE{2018arXiv181008266F,
	author	= { {Fuentes Perez}, L.~J. and {Romero Calla}, L.~A. and {Montenegro}, A.~A. },
	title	= { Dictionary Learning-based Inpainting on Triangular Meshes },
	journal	= { ArXiv e-prints },
	eprint	= { 1810.08266 },
	year	= 2018,
	month	= oct,
	url	= { https://arxiv.org/abs/1810.08266 }
}
```

### Hole repairing
We implemented repairing mesh holes in two steps:

1. Generate a mesh to cover the hole. We modified the algorithm presented in the paper: [A robust hole-filling algorithm for triangular mesh](https://doi.org/10.1007/s00371-007-0167-y), in order to
generate a planar triangular mesh using a priority queue.
2. Fit the surface described by the new points in order to minimize the variation of the surface,
solving the Poisson equation (see the Chapter 4 of the book [Polygon Mesh Processing](http://www.pmp-book.org/)) or using Biharmonic splines.

Please see and cite our final undergraduate project: [mesh hole repairing report](http://repositorio.unsa.edu.pe/handle/UNSA/2576) (in Spanish).

### Key-Points and Key-Components

We proposed a simple method based on the faces' areas to compute key-points for adaptive meshes.

Please cite our paper (in Spanish):

[Efficient approach for interest points detection in non-rigid shapes](https://doi.org/10.1109/CLEI.2015.7359459)

```bibtex
@INPROCEEDINGS{7359459,
	author		= { C. J. Lopez Del Alamo and L. A. Romero Calla and L. J. Fuentes Perez },
	booktitle	= { 2015 Latin American Computing Conference (CLEI) },
	title		= { Efficient approach for interest points detection in non-rigid shapes },
	year		= { 2015 },
	month		= { Oct },
	pages		= { 1-8 },
	doi		= { 10.1109/CLEI.2015.7359459 }
}
```

Computing key-components depends on the accuracy and definition of the key points. We were inspired
by the [work of Ivan Sipiran](https://www.researchgate.net/publication/262350194_Key-component_detection_on_3D_meshes_using_local_features),
he defined for the first time the notion of a key-component in meshes.
We proposed a method based on geodesics to determine the key components.

Please see and cite our final undergraduate project: [key-components report](http://repositorio.unsa.edu.pe/handle/UNSA/2575) (in Spanish).


### Decimation
We are implementing the algorithm described by the paper [Stellar Mesh Simplification Using Probabilistic Optimization](https://doi.org/10.1111/j.1467-8659.2004.00811.x),
to compute a mesh simplification.

### Fairing
We implemented Spectral and Taubin fairing algorithms to smooth a mesh surface.
See the Chapter 4 of the book [Polygon Mesh Processing](http://www.pmp-book.org/).

### Laplacian and signatures
Laplace-Beltrami operator and its eigen decomposition, WKS, HKS, GPS signatures.

## Documentation
Execute:

	doxygen Doxyfile

to generate the documentation in *html* and *LaTeX*.

## Viewer
The viewer was initially based in the viewer of [https://github.com/dgpdec/course](https://github.com/dgpdec/course). The current viewer uses VAO and VBO to render, and the shaders have been modified and upgraded.

## License

MIT License

## Authors
- [Lizeth Joseline Fuentes Pérez](https://github.com/lishh)
- [Luciano Arnaldo Romero Calla](https://github.com/larc)

