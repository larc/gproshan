#version 410 core

layout (location=0) in vec3 in_position;
layout (location=1) in vec3 in_normal;
layout (location=2) in vec3 in_mesh_color;
layout (location=3) in float in_color;

out vec3 vs_position;
out vec3 vs_normal;
out vec3 vs_mesh_color;
out float vs_color;

uniform mat4 model_view_mat;
uniform mat4 proj_mat;
uniform uint point_size;

void main()
{
	vs_position = in_position;
	vs_normal = in_normal;
	vs_mesh_color = in_mesh_color;
	vs_color = in_color;

	gl_Position = proj_mat * model_view_mat * vec4(in_position, 1);
	gl_PointSize = point_size;
}

