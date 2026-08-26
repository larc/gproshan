#version 410 core

layout (points) in;
layout (line_strip, max_vertices = 2) out;

in vec3 normal[];

uniform float length;

void main()
{
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	gl_Position = gl_in[0].gl_Position + vec4(normal[0], 0.) * length;
	EmitVertex();

	EndPrimitive();
}

