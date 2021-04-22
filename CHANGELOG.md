Version History
---------------

### gproshan 3.0.0
- Add module geometry, including a 2D convex hull algorithm implementation.
- Supported file mesh types include: off, obj, ply, ptx, xyz, and any depth image or image file loaded as a mesh.
- Upgraded version of geodesics module: fast marching algorithm, heat method, and parallel topleset propagation algorithm from our published paper.
- Upgraded version of the sparse mesh coding with feature aware sampling module and published paper reference.
- Updated save mesh with options for file types, normals, and point-cloud.
- Added render option using [Embree](https://www.embree.org/) ray tracing mesh, point-cloud disks, and splats.
- Implemented the loading and rendering of point clouds.
- Added heatmap viewer options and loading vertex color from a file.
- New user interface implemented with [ImGui](https://github.com/ocornut/imgui).
- Viewer upgraded using GLEW, GLFW, and GLM.

