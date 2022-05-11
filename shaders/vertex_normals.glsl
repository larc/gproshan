#version 410 core

layout (location=0) in vec3 in_position;
layout (location=1) in vec3 in_normal;

out vec3 normal;

uniform mat4 proj_view_mat;
uniform mat4 model_mat;

void main()
{
	normal = normalize(vec3(proj_view_mat * vec4(in_normal, 0.)));
	gl_Position = proj_view_mat * model_mat * vec4(in_position, 1.);
}

