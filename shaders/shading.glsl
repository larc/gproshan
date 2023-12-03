#include colormap.glsl
#include material.glsl

uniform uint idx_colormap;
uniform bool render_lines;

uniform vec3 eye;

uniform sampler2D tex_Ka;
uniform sampler2D tex_Kd;
uniform sampler2D tex_Ks;

uniform material mat;

struct light
{
	vec3 pos;
	vec3 color;
	float power;
};

uniform light ambient;
uniform light cam_light;

vec3 lines_colormap(vec3 color, float h)
{
	color = idx_colormap > 0 ? colormap(idx_colormap, h) : color;

	if(render_lines)
	{
		h = 50 * h;
		h = h - floor(h);
		h = (1 / (1 + exp(-100 * (h - .55)))) + (1 / (1 + exp(-100 * (-h + .45))));
		h = 1 - h;
		color = vec3(0) + (1. - h) * color;
	}

	return color;
}

vec3 shading(vec3 color, vec3 n, vec3 pos, vec2 texcoord)
{
	vec3 Ka = vec3(1);
	vec3 Kd = color;
	vec3 Ks = vec3(0.2);
	float Ns = 10;

	if(idx_colormap == 5)
	{
		Ka = mat.Ka;
		if(mat.map_Ka != -1)
			Ka = texture(tex_Ka, texcoord).rgb;

		Kd = mat.Kd;
		if(mat.map_Kd != -1)
			Kd = texture(tex_Kd, texcoord).rgb;

		Ks = mat.Ks;
		if(mat.map_Ks != -1)
			Ks = texture(tex_Ks, texcoord).rgb;

		Ns = mat.Ns;
	}

	vec3 l = cam_light.pos - pos;
	float r = length(l);
	l /= r;
	vec3 v = normalize(eye - pos);
	vec3 h = normalize(l + v);
	float lambertian = max(dot(l, n), 0.0);
	float specular = pow(max(dot(h, n), 0.0), Ns);

	return Ka * ambient.color * ambient.power +
			(lambertian * Kd + specular * Ks) * cam_light.color * cam_light.power / (r * r);
}


