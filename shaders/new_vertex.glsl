#version 450 core

layout (location=0) in vec3 VertexPosition;
layout (location=1) in vec3 VertexNormal;
layout (location=2) in float VertexColor;

out vec3 position;
out vec3 normal;
out float color;

uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;

void main()
{
	position = VertexPosition;
	normal = VertexNormal;
	color = VertexColor;
	gl_Position =  ProjectionMatrix * ModelViewMatrix * vec4(VertexPosition, 1.);
}

