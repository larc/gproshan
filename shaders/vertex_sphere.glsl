#version 410 core

layout (location=0) in vec3 in_position;
layout (location=1) in vec3 in_normal;
layout (location=3) in vec3 in_translation;

out vec3 vs_position;
out vec3 vs_normal;

uniform mat4 proj_view_mat;
uniform mat4 model_mat;
uniform float scale;

void main()
{
	vs_position = scale * in_position + in_translation;
	vs_normal = in_normal;

	gl_Position = proj_view_mat * model_mat * vec4(vs_position, 1);
}

