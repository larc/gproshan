#version 410 core

layout(location = 0) in vec3 in_position;

out vec2 uv;

void main()
{
	gl_Position = vec4(in_position, 1);
	uv = (vec2(in_position.x, in_position.y) + vec2(1)) / 2;
}

