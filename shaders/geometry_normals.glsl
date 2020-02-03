#version 460 core

layout (triangles) in;
layout (line_strip, max_vertices = 2) out;

in vec3 normal[];

const float lenght = 0.02;

void main()
{
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	gl_Position = gl_in[0].gl_Position + vec4(normal[0], 0.) * lenght;
	EmitVertex();

	EndPrimitive();
}

