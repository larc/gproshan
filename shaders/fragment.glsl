#version 410 core

#include shading.glsl
#include material.glsl

in vec3 gs_position;
in vec3 gs_normal;
in vec3 gs_mesh_color;
in float gs_color;
in vec2 gs_texcoord;

noperspective in vec3 edge_dist;

layout(location = 0) out vec4 frag_color;


uniform vec3 eye;
uniform vec3 cam_light;
uniform bool render_flat;
uniform bool render_wireframe;

uniform sampler2D tex_Ka;
uniform sampler2D tex_Kd;
uniform sampler2D tex_Ks;

uniform material mat;

void main()
{
	vec3 color = lines_colormap(gs_mesh_color, gs_color);

	vec3 Ka = mat.Ka;
	if(mat.map_Ka != -1)
		Ka *= texture(tex_Ka, gs_texcoord.xy).rgb;

	vec3 Kd = mat.Kd;
	if(mat.map_Kd != -1)
		Kd *= texture(tex_Kd, gs_texcoord.xy).rgb;

	vec3 Ks = mat.Ks;
	if(mat.map_Ks != -1)
		Ks *= texture(tex_Ks, gs_texcoord.xy).rgb;


 	vec3 n = render_flat ? normalize(cross(dFdx(gs_position), dFdy(gs_position))) : normalize(gs_normal);

	if(render_wireframe)
	{
		vec3 delta = fwidth(edge_dist);
		vec3 tmp = smoothstep(vec3(0), delta, edge_dist);
		float d = min(min(tmp.x, tmp.y), tmp.z);
		color = mix(vec3(.2), color, d);
	}

	vec3 l = cam_light - gs_position;
	float r = length(l);
	l /= r;
	vec3 v = normalize(eye - gs_position);
	vec3 h = normalize(l + v);
	float lambertian = max(dot(l, n), 0.0);
	float specular = pow(max(dot(h, n), 0.0), mat.Ns);
	float P = 4;
	vec3 La = vec3(0.2, 0.2, 0.2);

	color = Ka * La + (lambertian * Kd + specular * Ks) * P / (r * r);
	//color = shading(n, normalize(cam_light - gs_position), normalize(eye - gs_position), Kd);
	//color = Ka + Kd + Ks;
	frag_color = vec4(color, 1);
}

