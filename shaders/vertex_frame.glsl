#version 410 core

layout(location = 0) in vec3 in_position;

out vec2 UV;

void main()
{
	gl_Position = vec4(in_position, 1);
	UV = ( vec2(in_position.x, in_position.y) + vec2(1) ) / 2;
}

