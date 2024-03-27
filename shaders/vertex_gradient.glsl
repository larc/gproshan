#version 410 core

layout (location=0) in vec3 in_position;
layout (location=3) in float in_color;

out float color;
out vec3 position;

uniform mat4 proj_view_mat;
uniform mat4 model_mat;

void main()
{
	color = in_color;
	position = in_position;
	gl_Position = proj_view_mat * model_mat * vec4(in_position, 1.);
}

