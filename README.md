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

## Documentation
Execute:

	doxygen Doxyfile

to generate the documentation in html and latex.

## Viewer
The viewer was initially based in the viewer of [https://github.com/dgpdec/course](https://github.com/dgpdec/course). The current viewer use VAO and VBO to render, and the shaders have been modified and upgraded.

## Authors
- Lizeth Joseline Fuentes Pérez
- Luciano Arnaldo Romero Calla

