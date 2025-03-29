Version History
---------------

### gproshan v3.2

- The KNN module was added using flann for 3D point clouds.
- A simple path tracer implementation using Embree and OptiX was added.
- Added scene rendering: loader from .obj files and .mtl, handling textures.
- Exporting gproshan as cmake library, use find_package(gproshan) in your project.
- Added Intel Embree as the default ray tracing library for ray casting operations and as a rendering option with shadows.
- Adding Scenes module, virtual point cloud scanners, and point cloud normals computation for 3D scenes.
- Added render option using [OptiX](https://developer.nvidia.com/optix) ray tracing mesh, shadows.
- Add module geometry, including a 2D convex hull algorithm implementation and connected components detection.
- Supported file mesh types include off, obj, ply, ptx, xyz, and any depth image or image file loaded as a mesh.
- Upgraded version of geodesics module: fast marching algorithm, heat method, and parallel topleset propagation algorithm from our published paper.
- Upgraded the sparse mesh coding version with a feature-aware sampling module and published paper reference.
- Updated save mesh with file types, normals, and point cloud options.
- Implemented the loading and rendering of point clouds.
- Added heatmap viewer options and loaded vertex color from a file.
- A new user interface was implemented with [ImGui](https://github.com/ocornut/imgui).
- Viewer upgraded using GLEW and GLFW3.

