#version 410 core

in vec3 vs_position;
in vec3 vs_normal;

layout(location = 0) out vec4 frag_color;

uniform vec3 eye;
uniform vec3 light;

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
	vec3 color = vec3(1, 0, 0);

	vec3 N = normalize(vs_normal);
	vec3 L = normalize(light - vs_position);
	vec3 E = normalize(eye - vs_position);
	vec3 R = 2 * dot(L, N) * N - L;
	vec3 one = vec3(1);

	color = diffuse(N, L) * color + 0.1 * specular(N, L, E) * one + 0.25 * fresnel(N,E) * one;
	frag_color = vec4(color, 1);
}

