#version 410 core

#include ../shaders/colormap.glsl

in vec3 vs_position;
in vec3 vs_normal;
in float vs_color;

layout(location = 0) out vec4 frag_color;

uniform uint idx_colormap;
uniform vec3 eye;
uniform vec3 light;
uniform bool render_lines;

float diffuse(vec3 N, vec3 L)
{
	return max(0, dot(N, L));
}

float specular(vec3 N, vec3 L, vec3 E)
{
	const float shininess = 4;
	vec3 R = 2 * dot(L, N) * N - L;

	return pow(max(0, dot(R, E)), shininess);
}

float fresnel(vec3 N, vec3 E)
{
	const float sharpness = 10;
	float NE = max(0, dot(N, E));

	return pow(sqrt( 1. - NE * NE ), sharpness);
}

void main()
{
	vec3 color = colormap(idx_colormap, vs_color);

	// lines
	if(render_lines)
	{
		float h = vs_color;
		h = h * 40;
		h = h - floor(h);
		h = (1 / (1 + exp(-100 * (h - .55)))) + (1 / (1 + exp(-100 * (-h + .45))));
		h = 1 - h;
		color = vec3(0) + (1. - h) * color;
	}

 	vec3 N = normalize(vs_normal);
	
	vec3 L = normalize(light - vs_position);
	vec3 E = normalize(eye - vs_position);
	vec3 R = 2 * dot(L, N) * N - L;
	vec3 one = vec3(1);


	color = diffuse(N, L) * color + 0.1 * specular(N, L, E) * one + 0.25 * fresnel(N,E) * one;
	frag_color = vec4(color, 1);
}

