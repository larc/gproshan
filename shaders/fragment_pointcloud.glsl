#version 410 core

#include shading.glsl

in vec3 vs_position;
in vec3 vs_normal;
in vec3 vs_mesh_color;
in float vs_color;

layout(location = 0) out vec4 frag_color;

uniform vec3 eye;
uniform vec3 cam_light;
uniform bool point_normals;


void main()
{
	vec3 color = lines_colormap(vs_mesh_color, vs_color);

	if(point_normals)
		color = shading(normalize(vs_normal), normalize(cam_light - vs_position), normalize(eye - vs_position), color);

	frag_color = vec4(color, 1);
}

