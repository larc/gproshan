# gproshan
### geometry processing and shape analysis framework

## Description
This framework include some algorithms of geometry processing and shape analysis as part of our graduate research.

## Build and Run
Install all dependencies and run:

	make

finally execute:

	./gproshan [input mesh paths]

![](gproshan.png) 

## Dependencies (linux)
g++ 7.2, cuda >= 9.1, libarmadillo, libeigen, libsuitesparse, libopenblas, opengl, openmp, gnuplot, libcgal

## Contributions

### CHE implementation
We have implemented the CHE: compact half-edge data structure, to manipulated the meshes.

### Geodesics
Fast Marching, Heat Method.
**We proposed a new parallel algorithm to compute geodesics**.
See our paper...

### Dictionary Learning
**We proposed a Dictionary Learning technique** in order to solve the problems of Denoising, Inpainting and Multiresolution.
See our technical report.

### Hole repairing
We implemented repairing mesh holes in two steps:

1. Generate a mesh to cover the hole (modified algorithm base on ...).
2. Approximate the curvature solving the Poisson equation and using Biharmonic splines.

### Decimation
We are implementing the algorithm ?? to decimate a mesh. 

### Fairing
Spectral and Taubin algorithms.

### Laplacian and signatures
Laplace-Beltrami operator and its eigen decomposition, WKS, HKS, GPS signatures.

## Documentation
Execute:

	doxygen Doxyfile

to generate the documentation in html and latex.

## Viewer
The viewer was initially based in the viewer of [https://github.com/dgpdec/course](https://github.com/dgpdec/course). The current viewer use VAO and VBO to render, and the shaders have been modified and upgraded.

## Authors
- Lizeth Joseline Fuentes Pérez
- Luciano Arnaldo Romero Calla

