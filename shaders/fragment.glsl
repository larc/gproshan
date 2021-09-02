#version 410 core

#include shading.glsl

in vec3 gs_position;
in vec3 gs_normal;
in vec3 gs_mesh_color;
in float gs_color;

noperspective in vec3 edge_dist;

layout(location = 0) out vec4 frag_color;


uniform vec3 eye;
uniform vec3 light;
uniform bool render_flat;
uniform bool render_wireframe;


void main()
{
	vec3 color = lines_colormap(gs_mesh_color, gs_color);
 	vec3 N = render_flat ? normalize(cross(dFdx(gs_position), dFdy(gs_position))) : normalize(gs_normal);

	if(render_wireframe)
	{
		vec3 delta = fwidth(edge_dist);
		vec3 tmp = smoothstep(vec3(0), delta, edge_dist);
		float d = min(min(tmp.x, tmp.y), tmp.z);
		color = mix(vec3(.2), color, d);
	}

	color = shading(N, normalize(light - gs_position), normalize(eye - gs_position), color);
	frag_color = vec4(color, 1);
}

