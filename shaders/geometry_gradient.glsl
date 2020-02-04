#version 460 core

layout (triangles) in;
layout (line_strip, max_vertices = 2) out;

in float color[];

uniform float length;

void main()
{
	vec3 xi = vec3(gl_in[0].gl_Position);
	vec3 xj = vec3(gl_in[1].gl_Position);
	vec3 xk = vec3(gl_in[2].gl_Position);

	vec3 n = normalize(cross(xj - xi, xk - xi));
	
	vec3 pij = cross(n, xj - xi);
	vec3 pjk = cross(n, xk - xj);
	vec3 pki = cross(n, xi - xk);

	vec3 g = normalize(color[0] * pjk + color[1] * pki + color[2] * pij);

	vec3 a = (xi + xj + xk) / 3;
	vec3 b = a + g * length;

	gl_Position = vec4(a, 0.);
	EmitVertex();
	
	gl_Position = vec4(b, 0.);
	EmitVertex();
	
	EndPrimitive();
}

