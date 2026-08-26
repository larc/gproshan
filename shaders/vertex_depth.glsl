#version 410 core

layout (location=0) in vec3 in_position;

uniform mat4 proj_view_mat;
uniform mat4 model_mat;

void main()
{
	gl_Position = proj_view_mat * model_mat * vec4(in_position, 1);
}

