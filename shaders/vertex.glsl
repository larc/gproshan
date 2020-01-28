#version 460 core

layout (location=0) in vec3 in_position;
layout (location=1) in vec3 in_normal;
layout (location=2) in float in_color;

out vec3 position;
out vec3 normal;
out float color;

uniform mat4 model_view_mat;
uniform mat4 proj_mat;

void main()
{
	position = in_position;
	normal = in_normal;
	color = in_color;
	gl_Position =  proj_mat * model_view_mat * vec4(position, 1.);
}

