#version 410 core

#include shading.glsl

in vec3 vs_position;
in vec3 vs_normal;

layout(location = 0) out vec4 frag_color;


void main()
{
	frag_color = vec4(shading(vec3(1, 0, 0), normalize(vs_normal), vs_position, vec2(0)), 1);
}

