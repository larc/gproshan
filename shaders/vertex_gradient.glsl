#version 460 core

layout (location=0) in vec3 in_position;
layout (location=2) in float in_color;

out float color;
out vec3 position;

uniform mat4 model_view_mat;
uniform mat4 proj_mat;

void main()
{
	color = in_color;
	position = in_position;
	gl_Position = proj_mat * model_view_mat * vec4(in_position, 1.);
}

