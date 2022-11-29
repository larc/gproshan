#version 410 core

#include shading.glsl

in vec3 vs_position;
in vec3 vs_normal;
in vec3 vs_mesh_color;
in float vs_color;
in vec2 vs_texcoord;

layout(location = 0) out vec4 frag_color;

uniform bool point_normals;


void main()
{
	vec3 color = lines_colormap(vs_mesh_color, vs_color);

	if(point_normals)
		color = shading(color, normalize(vs_normal), vs_position, vs_texcoord);

	frag_color = vec4(color, 1);
}

