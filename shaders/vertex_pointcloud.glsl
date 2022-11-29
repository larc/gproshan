#version 410 core

layout (location=0) in vec3 in_position;
layout (location=1) in vec3 in_normal;
layout (location=2) in vec3 in_mesh_color;
layout (location=3) in float in_color;
layout (location=4) in vec2 in_texcoord;

out vec3 vs_position;
out vec3 vs_normal;
out vec3 vs_mesh_color;
out float vs_color;
out vec2 vs_texcoord;

uniform mat4 proj_view_mat;
uniform mat4 model_mat;
uniform uint point_size;

void main()
{
	vs_position = vec3(model_mat * vec4(in_position, 1));
	vs_normal = in_normal;
	vs_mesh_color = in_mesh_color;
	vs_color = in_color;
	vs_texcoord = in_texcoord;

	gl_Position = proj_view_mat * vec4(vs_position, 1);
	gl_PointSize = point_size;
}

