#version 410 core

layout (triangles) in;
layout (line_strip, max_vertices = 2) out;

in float color[];
in vec3 position[];

uniform float length;
uniform mat4 model_view_mat;
uniform mat4 proj_mat;

void main()
{
	vec3 xi = position[0];
	vec3 xj = position[1];
	vec3 xk = position[2];

	vec3 n = normalize(cross(xj - xi, xk - xi));
	
	vec3 pij = cross(n, xj - xi);
	vec3 pjk = cross(n, xk - xj);
	vec3 pki = cross(n, xi - xk);

	vec3 g = normalize(color[0] * pjk + color[1] * pki + color[2] * pij);

	vec3 a = (xi + xj + xk) / 3.0;
	vec3 b = a + g * 0.3 * length;

	gl_Position = proj_mat * model_view_mat * vec4(a, 1.);
	EmitVertex();
	
	gl_Position = proj_mat * model_view_mat * vec4(b, 1.);
	EmitVertex();
	
	EndPrimitive();
}

