Version History
---------------

### gproshan v4.0

- Added KNN module using flann for 3D point clouds.
- Added a simple path tracer implementation using Embree and OptiX.
- Added scene rendering: loader from .obj files and .mtl, handling textures.
- Exporting gproshan as cmake library, use find_package(gproshan) in your project.
- Added Intel Embree as default ray tracing library, for ray casting operations and as rendering option with shadows.
- Adding Scenes module, virtual point cloud scanners and point cloud normals computation for 3D scenes.
- Added render option using [OptiX](https://developer.nvidia.com/optix) ray tracing mesh, shadows.
- Add module geometry, including a 2D convex hull algorithm implementation and connected components detection.
- Supported file mesh types include: off, obj, ply, ptx, xyz, and any depth image or image file loaded as a mesh.
- Upgraded version of geodesics module: fast marching algorithm, heat method, and parallel topleset propagation algorithm from our published paper.
- Upgraded version of the sparse mesh coding with feature aware sampling module and published paper reference.
- Updated save mesh with options for file types, normals, and point-cloud.
- Added render option using [Embree](https://www.embree.org/) ray tracing mesh, shadows, point-cloud disks, and splats.
- Implemented the loading and rendering of point clouds.
- Added heatmap viewer options and loading vertex color from a file.
- New user interface implemented with [ImGui](https://github.com/ocornut/imgui).
- Viewer upgraded using GLEW and GLFW3.

