#version 410 core

#include shading.glsl

in vec3 vs_position;
in vec3 vs_normal;

layout(location = 0) out vec4 frag_color;

uniform vec3 eye;
uniform vec3 cam_light;


void main()
{
	vec3 color = shading(normalize(vs_normal), normalize(cam_light - vs_position), normalize(eye - vs_position), vec3(1, 0, 0));
	frag_color = vec4(color, 1);
}

