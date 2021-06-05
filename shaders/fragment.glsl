#version 410 core

#include colormap.glsl

in vec3 gs_position;
in vec3 gs_normal;
in vec3 gs_mesh_color;
in float gs_color;

noperspective in vec3 edge_dist;

layout(location = 0) out vec4 frag_color;

uniform uint idx_colormap;
uniform vec3 eye;
uniform vec3 light;
uniform bool render_flat;
uniform bool render_lines;
uniform bool render_wireframe;

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
	vec3 color = idx_colormap > 0 ? colormap(idx_colormap, gs_color) : gs_mesh_color;

	if(render_lines)
	{
		float h = gs_color;
		h = h * 40;
		h = h - floor(h);
		h = (1 / (1 + exp(-100 * (h - .55)))) + (1 / (1 + exp(-100 * (-h + .45))));
		h = 1 - h;
		color = vec3(0) + (1. - h) * color;
	}

 	vec3 N;

 	if(render_flat)
		N = normalize(cross(dFdx(gs_position), dFdy(gs_position)));
	else
		N = normalize(gs_normal);

	vec3 L = normalize(light - gs_position);
	vec3 E = normalize(eye - gs_position);
	vec3 R = 2 * dot(L, N) * N - L;
	vec3 one = vec3(1);


	if(render_wireframe)
	{
		vec3 delta = fwidth(edge_dist);
		vec3 tmp = smoothstep(vec3(0), delta, edge_dist);

		float d = min(min(tmp.x, tmp.y), tmp.z);

		color = mix(vec3(.2), color, d);
	}

	color = diffuse(N, L) * color + 0.1 * specular(N, L, E) * one + 0.25 * fresnel(N,E) * one;
	frag_color = vec4(color, 1);
}

