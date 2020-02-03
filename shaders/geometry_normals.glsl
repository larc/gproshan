#version 460 core

layout (triangles) in;
layout (line_strip, max_vertices = 6) out;

in vec3 normal[];

uniform float length;

void line_normal(int i)
{
	gl_Position = gl_in[i].gl_Position;
	EmitVertex();

	gl_Position = gl_in[i].gl_Position + vec4(normal[i], 0.) * length;
	EmitVertex();

	EndPrimitive();
}

void main()
{
	line_normal(0);
	line_normal(1);
	line_normal(2);
}

