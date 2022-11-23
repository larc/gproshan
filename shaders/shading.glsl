#include colormap.glsl
#include material.glsl

uniform bool render_lines;
uniform uint idx_colormap;

uniform material mat;


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

	return pow(sqrt(1. - NE * NE ), sharpness);
}

vec3 lines_colormap(vec3 color, float h)
{
	color = idx_colormap > 0 ? colormap(idx_colormap, h) : color;
	color = mat.Kd;	

	if(render_lines)
	{
		h = h * 40;
		h = h - floor(h);
		h = (1 / (1 + exp(-100 * (h - .55)))) + (1 / (1 + exp(-100 * (-h + .45))));
		h = 1 - h;
		color = vec3(0) + (1. - h) * color;
	}

	return color;
}

vec3 shading(vec3 N, vec3 L, vec3 E, vec3 color)
{
	vec3 one = vec3(1);
	return diffuse(N, L) * color + 0.1 * specular(N, L, E) * one + 0.25 * fresnel(N, E) * one;
}

